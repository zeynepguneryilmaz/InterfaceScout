"""Reproducible downloadable ZIP with map screenshots, PDF, CSV, XLSX, and B-factor PDBs.

Screenshots are captured from the real 3Dmol viewer on the client, NOT invented on
the server. Images and the PDB-coordinate sources are both supplied explicitly.
This is a descriptive research report, not a prediction of adsorption affinity.
"""
from __future__ import annotations

import base64
import csv
from io import BytesIO, StringIO
import json
import hashlib
import math
import re
import zipfile
from typing import Any, Mapping
from xml.sax.saxutils import escape

from PIL import Image as PILImage
from ensemble import build_consensus
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak, KeepTogether

MAP_SPECS = (
    ("anionic", "Anionic surface"),
    ("cationic", "Cationic surface"),
    ("hydrophobic", "Hydrophobic surface"),
    ("pi_carbon", "Aromatic / pi surface"),
    ("hbond_donor", "H-bond donor surface"),
    ("hbond_acceptor", "H-bond acceptor surface"),
    ("oxide", "Metal-oxide surface"),
    ("calcium_phosphate", "Calcium/phosphate sites"),
    ("metal_coord", "Transition-metal coordination"),
    ("soft_metal_sulfur", "Soft-metal / sulfur"),
    ("phosphate", "Phosphate-rich surface"),
)
MAP_NAMES = [name for name, _ in MAP_SPECS]


def safe_name(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip()).strip("._")
    return (out or "structure")[:85]


def _csv_bytes(rows: list[dict[str, Any]], columns: list[str] | None = None) -> bytes:
    cols = list(columns or [])
    for row in rows:
        for key in row:
            if key not in cols:
                cols.append(key)
    out = StringIO(newline="")
    if cols:
        writer = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict))
                             else v for k, v in row.items()})
    return out.getvalue().encode("utf-8-sig")


def _xlsx_bytes(tables: Mapping[str, list[dict[str, Any]]]) -> bytes:
    workbook = Workbook()
    workbook.remove(workbook.active)
    for key, rows in tables.items():
        sheet = workbook.create_sheet(key[:31].replace("/", "_"))
        keys = list(dict.fromkeys(k for row in rows for k in row))
        if not keys:
            sheet.append(["No data"])
        else:
            sheet.append(keys)
            for row in rows:
                sheet.append([json.dumps(row.get(k), ensure_ascii=False)
                              if isinstance(row.get(k), (list, dict))
                              else row.get(k) for k in keys])
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for cell in sheet[1]:
                cell.fill = PatternFill("solid", fgColor="263755")
                cell.font = Font(bold=True, color="FFFFFF")
            for col in sheet.columns:
                letter = col[0].column_letter
                sheet.column_dimensions[letter].width = min(55, max(12,
                    max(len(str(c.value or "")) for c in list(col)[:150]) + 2))
    out = BytesIO()
    workbook.save(out)
    return out.getvalue()


def _image_bytes(data: str) -> bytes:
    if not isinstance(data, str) or not data.startswith("data:image/png;base64,"):
        raise ValueError("Every exported image must be an actual PNG screenshot from the viewer")
    try:
        blob = base64.b64decode(data.split(",", 1)[1], validate=True)
    except Exception as e:
        raise ValueError("Invalid screenshot base64 encoding") from e
    if len(blob) > 5_000_000 or len(blob) < 50 or not blob.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("Invalid or oversized PNG screenshot")
    try:
        im = PILImage.open(BytesIO(blob))
        im.verify()
    except Exception as e:
        raise ValueError("Unable to decode exported screenshot") from e
    return blob


def _bfactor_pdb(original: str, scores: Mapping[str, float], label: str) -> bytes:
    if not original or "\n" not in original or "ATOM" not in original:
        raise ValueError("An original coordinate PDB is required for B-factor export")
    comments = [
        "REMARK 900 INTERFACESCOUT DISPLAY FIELD IN PDB B-FACTOR COLUMN",
        "REMARK 900 SCORE: " + label[:65],
        "REMARK 900 0-100 SCALE; NOT EXPERIMENTAL CRYSTALLOGRAPHIC B FACTORS",
        "REMARK 900 UNMAPPED RESIDUES ARE WRITTEN AS ZERO",
    ]
    out: list[str] = []
    for line in original.replace("\r", "").split("\n"):
        if line.startswith(("ATOM  ", "HETATM")):
            pad = line.ljust(66)
            chain = pad[21:22].strip()
            seq = pad[22:26].strip()
            ic = pad[26:27].strip()
            try:
                key = f"{chain}:{int(seq)}:{ic}"
            except ValueError:
                out.append(line)
                continue
            val = scores.get(key, 0.0)
            val = max(0.0, min(100.0, float(val)))
            out.append(pad[:60] + f"{val:6.2f}" + pad[66:])
        else:
            out.append(line)
    return ("\n".join(comments + out) + "\n").encode("utf-8")


def _map_scores(analysis: Mapping[str, Any], map_key: str) -> dict[str, float]:
    return {
        str(r["key"]): float(r.get("propensity") or 0.0)
        for r in analysis["maps"][map_key].get("residues", [])
    }


def _topk_scores(m: Mapping[str, Any], k: int) -> dict[str, float]:
    # A PDB for top-K region membership, not chemical propensity.
    keys = {str(q) for p in m.get("patches", [])[:k] for q in p.get("members", [])}
    return {key: 100.0 for key in keys}


def structure_tables(name: str, analysis: Mapping[str, Any],
                     combined: Mapping[str, Any] | None = None) -> dict[str, list[dict[str, Any]]]:
    rows = []
    for r in analysis.get("residue_profiles", []):
        rows.append({"Structure": name, **r, "groups": ";".join(r.get("groups", []))})
    maps, favorable_rows, patches, members, interactions = [], [], [], [], []
    all_residues = list(analysis.get("residue_profiles", []))
    for key in analysis["map_order"]:
        m = analysis["maps"][key]
        favorable = {str(r["key"]): r for r in m.get("residues", [])}
        for r in all_residues:
            mapped = favorable.get(str(r["key"]))
            maps.append({
                "Structure": name, "Map": key,
                "ResidueKey": r["key"], "ResidueName": r["res_name"],
                "Chain": r["chain"], "ResidueNumber": r["res_seq"],
                "InsertionCode": r.get("icode", ""),
                "SidechainSASA_A2": r.get("sidechain_sasa"),
                "scRSA_raw": r.get("scrsa_raw"),
                "SurfaceExposed": bool(r.get("surface_exposed")),
                "ChemistrySeed": bool(mapped is not None and float(mapped.get("local_score") or 0) > 0),
                "LocalScore": float(mapped.get("local_score") or 0) if mapped else 0.0,
                "Propensity": float(mapped.get("propensity") or 0) if mapped else 0.0,
                # Undefined for non-favorable residues; do not conflate with a real zero.
                "MultiscalePersistence": mapped.get("multiscale_persistence") if mapped else None,
            })
        favorable_rows.extend({"Structure": name, "Map": key, **r}
                              for r in m.get("residues", []))
        for patch in m.get("patches", []):
            p = dict(patch)
            p.pop("center_residue", None)
            p.pop("member_residues", None)
            p.pop("seed_residues", None)
            patches.append({"Structure": name, "Map": key, **p})
            seeds = set(patch.get("seed_members", []))
            for residue in patch.get("member_residues", []):
                members.append({
                    "Structure": name, "Map": key, "Rank": patch.get("display_rank"),
                    **residue, "chemical_seed": residue.get("key") in seeds,
                })
    if combined:
        for p in combined.get("regions", []):
            interactions.append({
                "Structure": name,
                "Selected_maps": ";".join(combined.get("selected_maps", [])),
                **p,
            })
        intersections = [{
            "Structure": name, "Selected_maps": ";".join(combined.get("selected_maps", [])),
            "Exact_shared_residue": k,
        } for k in combined.get("exact_same_residue_overlap", [])]
    else:
        intersections = []
    metadata = [
        {"Section": section, "Key": k, "Value": v}
        for section in ("input", "method", "structure_preparation", "scope")
        for k, v in analysis.get(section, {}).items()
    ]
    return {
        "Run_Metadata": metadata,
        "Protein_Residues": rows,
        "Map_Residue_Scores": maps,
        "Map_Favorable_Residues": favorable_rows,
        "All_Patches": patches,
        "Patch_Members": members,
        "Map_Colocation": interactions,
        "Exact_Intersections": intersections,
    }


def consensus_tables(consensus: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    rows = [{"Map": key, **r}
            for key in consensus["map_order"]
            for r in consensus["maps"][key]["residues"]]
    summary = []
    for key in consensus["map_order"]:
        m = consensus["maps"][key]
        summary.append({
            "Map": key,
            "Strict_compatible": len(m["strict_compatible_residues"]),
            "Top1_strict": len(m["strict_top1_members"]),
            "Top3_strict": len(m["strict_top3_members"]),
            "Top5_strict": len(m["strict_top5_members"]),
            "Top1_majority": len(m["majority_top1_members"]),
            "Top3_majority": len(m["majority_top3_members"]),
            "Top5_majority": len(m["majority_top5_members"]),
            "Top1_centers": m["top1_centers_by_structure"],
        })
    combined = consensus.get("combined_consensus") or {}
    return {
        "Consensus_Residues": list(consensus.get("residues", [])),
        "Consensus_Map_Scores": rows,
        "Consensus_Summary": summary,
        "Residue_Matching": list(consensus.get("matching", [])),
        "Combined_Consensus": list(combined.get("residues", [])),
    }


def comparison_tables(comparisons: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    info, changes = [], []
    for entry in comparisons:
        a, b = entry["structure_a"], entry["structure_b"]
        for m in entry["map_comparisons"]:
            for k in m["top_k"]:
                info.append({
                    "Structure_A": a, "Structure_B": b,
                    "Map": m["key"],
                    "Center_A": m["top1_center_a"], "Center_B": m["top1_center_b"],
                    "Top_K": k["top_k"], "Jaccard": k["jaccard"],
                    "n_A": k["n_a"], "n_B": k["n_b"], "n_shared": k["n_shared"],
                    "Retained": k["retained"],
                    "Lost_from_B": k["lost_in_b"], "New_in_B": k["new_in_b"],
                    "Coverage_A": k["surface_coverage_a"],
                    "Coverage_B": k["surface_coverage_b"],
                })
        for r in entry.get("residue_changes", []):
            changes.append({"Structure_A": a, "Structure_B": b, **r})
    return {"Comparison_Maps": info, "Comparison_Residues": changes}


STYLES = getSampleStyleSheet()
STYLES.add(ParagraphStyle(name="ISsmall", parent=STYLES["BodyText"], fontSize=8,
                          leading=10, spaceAfter=3))
STYLES.add(ParagraphStyle(name="ISmicro", parent=STYLES["BodyText"], fontSize=6.5,
                          leading=8.2, spaceAfter=2))
STYLES.add(ParagraphStyle(name="ISheading", parent=STYLES["Heading2"], fontSize=13,
                          textColor=colors.HexColor("#283751"), spaceBefore=12, spaceAfter=6))


def _p(v: Any, micro: bool = False):
    s = str(v if v is not None else "")
    s = (s.replace("\u03c0", "pi").replace("\u00c5", "A")
          .replace("\u2013", "-").replace("\u2014", "-")
          .replace("\u2265", ">=").replace("\u2264", "<="))
    return Paragraph(escape(s), STYLES["ISmicro" if micro else "ISsmall"])


def _table(story: list, headers: list[str], rows: list[list[Any]],
           widths: list[float] | None = None):
    records = [[_p(x, True) for x in headers]]
    records += [[_p(x, True) for x in row] for row in rows]
    if len(records) == 1:
        records.append([_p("No data", True)] + [_p("", True)] * (len(headers) - 1))
    table = Table(records, colWidths=widths, hAlign="LEFT", repeatRows=1,
                  splitByRow=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8edf5")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#172033")),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#d4dce8")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(table)


def _screenshot(story: list, blob: bytes):
    im = PILImage.open(BytesIO(blob))
    w, h = im.size
    max_w, max_h = 510, 310
    scale = min(max_w / w, max_h / h)
    story.append(Image(BytesIO(blob), width=w * scale, height=h * scale))
    story.append(Spacer(1, 8))


def _make_pdf(title: str, build) -> bytes:
    out = BytesIO()
    doc = SimpleDocTemplate(out, pagesize=A4,
                            rightMargin=36, leftMargin=36,
                            topMargin=38, bottomMargin=37,
                            title=title, author="InterfaceScout")
    story = [_p(title), Spacer(1, 8)]
    build(story)
    doc.build(story)
    return out.getvalue()


def map_pdf(structure_name: str, analysis: Mapping[str, Any], key: str, png: bytes) -> bytes:
    m = analysis["maps"][key]
    def build(story: list):
        story.append(_p("Structure: " + structure_name))
        story.append(_p("Map: " + m["label"]))
        story.append(_p(m.get("description", "")))
        story.append(_p("Input conditions: " + json.dumps(analysis.get("input", {}))))
        story.append(_p("Image: screenshot of the actual 3D viewer with this channel highlighted."))
        _screenshot(story, png)
        story.append(_p("All ranked regions", False))
        _table(story,
            ["Rank / front", "Center", "Residues / seeds", "Chemistry", "Access.", "Persist.", "Orient."],
            [[f'{p.get("display_rank")}/{p.get("pareto_front")}',
              p.get("center_key"),
              f'{p.get("n_members")}/{p.get("n_chemistry_seeds")}',
              f'{float(p.get("chemistry_support", 0)):.3f}',
              f'{float(p.get("mean_accessibility", 0)):.3f}',
              f'{float(p.get("patch_coherence", 0)):.3f}',
              f'{float(p.get("orientation_coherence", 0)):.3f}']
             for p in m["patches"]], [63, 82, 86, 66, 65, 65, 68])
        story.append(Spacer(1, 8))
        story.append(_p("Complete region membership (rank: region members | chemistry seeds)"))
        for patch in m["patches"]:
            story.append(_p(
                f'{patch["display_rank"]}: {", ".join(patch.get("members", []))} | '
                f'seeds: {", ".join(patch.get("seed_members", []))}', True))
        story.append(Spacer(1, 8))
        story.append(_p("All channel residue measurements"))
        _table(story,
            ["Residue", "Name", "scRSA", "Local score", "Propensity", "Persistence"],
            [[r.get("key"), r.get("res_name"), r.get("scrsa"),
              r.get("local_score"), r.get("propensity"),
              r.get("multiscale_persistence")] for r in m["residues"]],
            [99, 58, 79, 88, 98, 107])
        story.append(Spacer(1, 6))
        story.append(_p("Ranked patches are structural candidates, not validated adsorption sites. "
                        "B-factor exports contain visualization scores, not crystallographic B-factors."))
    return _make_pdf(f"InterfaceScout | {structure_name} | {m['label']}", build)


def report_pdf(structure_name: str, analysis: Mapping[str, Any],
               images: Mapping[str, bytes], combined: Mapping[str, Any] | None,
               combined_image: bytes | None = None) -> bytes:
    def build(story: list):
        story.append(_p("Full 11-channel protein-surface analysis; not an adsorption prediction."))
        for section in ("input", "method", "structure_preparation", "scope"):
            story.append(_p(section.upper()))
            for k, v in analysis.get(section, {}).items():
                story.append(_p(f"{k}: {json.dumps(v, ensure_ascii=True)}", True))
        prop = analysis.get("protein_properties", {})
        story.append(_p("Direct protein properties"))
        story.append(_p(f"All residues: {prop.get('n_residues')}; "
                        f"exposed residues: {prop.get('n_surface_residues')}; "
                        f"total side-chain SASA: {prop.get('total_sidechain_sasa_A2')} A2"))
        for key, value in prop.get("groups", {}).items():
            story.append(_p(f"{key}: {value.get('n_exposed_residues')} exposed, "
                            f"{', '.join(value.get('exposed_residues', []))}", True))
        story.append(_p("Complete direct protein residue properties"))
        _table(story, ["Residue", "Type", "Exposed", "Side-chain SASA", "scRSA", "Charge"],
            [[r.get("key"), r.get("res_name"), r.get("surface_exposed"),
              r.get("sidechain_sasa"), r.get("scrsa_raw"), r.get("charge_descriptor")]
             for r in analysis.get("residue_profiles", [])], [106, 58, 67, 109, 81, 102])
        story.append(PageBreak())
        for index, key in enumerate(analysis["map_order"]):
            m = analysis["maps"][key]
            story.append(Paragraph(escape(m["label"]), STYLES["ISheading"]))
            story.append(_p(m.get("description", "")))
            _screenshot(story, images[key])
            story.append(_p(f"Total ranked regions: {m['n_patches']}; "
                            f"Pareto-front-1 regions: {m['n_pareto_primary_patches']}"))
            _table(story, ["Rank", "Center", "Members", "Seeds", "Chemistry", "Accessibility", "Coherence"],
                [[p.get("display_rank"), p.get("center_key"), p.get("n_members"),
                  p.get("n_chemistry_seeds"), f'{p.get("chemistry_support", 0):.3f}',
                  f'{p.get("mean_accessibility", 0):.3f}', f'{p.get("patch_coherence", 0):.3f}']
                for p in m["patches"]], [48, 85, 63, 54, 87, 99, 81])
            story.append(Spacer(1, 6))
            story.append(_p("All patch-member residue keys, including chemistry seeds"))
            for patch in m["patches"]:
                story.append(_p(
                    f'{patch["display_rank"]} ({patch["center_key"]}): '
                    f'[{", ".join(patch.get("members", []))}] ; '
                    f'seeds [{", ".join(patch.get("seed_members", []))}]', True))
            story.append(Spacer(1, 6))
            story.append(_p("All residue-level map scores"))
            _table(story, ["Residue", "Name", "scRSA", "Local", "Propensity", "Persistence"],
                [[r.get("key"), r.get("res_name"), r.get("scrsa"),
                  r.get("local_score"), r.get("propensity"),
                  r.get("multiscale_persistence")] for r in m["residues"]],
                [99, 58, 79, 88, 98, 107])
            if index < len(analysis["map_order"]) - 1:
                story.append(PageBreak())
        if combined:
            story.append(PageBreak())
            story.append(Paragraph("Selected-map intersections and co-localization", STYLES["ISheading"]))
            story.append(_p("Selected channels: " + ", ".join(combined.get("selected_maps", []))))
            story.append(_p("Exact same-residue overlap: " +
                            ", ".join(combined.get("exact_same_residue_overlap", []))))
            if combined_image is not None:
                story.append(_p("Selected channels: actual 3D viewer composite screenshot"))
                _screenshot(story, combined_image)
            for region in combined.get("regions", []):
                story.append(_p(
                    f'{region.get("source_map")} {region.get("source_center_key")}: '
                    f'maps [{", ".join(region.get("represented_maps", []))}]; '
                    f'members [{", ".join(region.get("member_keys", []))}]; '
                    f'exact all-map residues [{", ".join(region.get("exact_all_map_seed_overlap", []))}]',
                    True))
        story.append(_p("All original inputs and full-precision structured data are retained "
                        "in the ZIP JSON and CSV tables."))
    return _make_pdf(f"InterfaceScout | Complete report | {structure_name}", build)


def consensus_pdf(consensus: Mapping[str, Any], images: Mapping[str, bytes]) -> bytes:
    def build(story: list):
        story.append(_p(consensus["method"]))
        story.append(_p(consensus["limits"]))
        story.append(_p("Input structures: " + ", ".join(consensus["structure_names"])))
        story.append(_p(f"Identity-matched residues: {consensus['n_matched_residues']}"))
        _table(story, ["Residue", "Exposed count", "Fraction", "scRSA min", "scRSA max"],
            [[r["key"], f'{r["n_exposed"]}/{r["n_structures"]}',
              r["exposure_fraction"], r["scRSA_min"], r["scRSA_max"]]
             for r in consensus["residues"]], [130, 100, 90, 90, 99])
        combined = consensus.get("combined_consensus")
        if combined:
            story.append(PageBreak())
            story.append(Paragraph("Selected-map consensus", STYLES["ISheading"]))
            story.append(_p("Selected maps: " + ", ".join(combined["selected_maps"])))
            story.append(_p(combined["note"]))
            story.append(_p("Strict same-residue overlap: " +
                            ", ".join(combined["strict_exact_intersection"])))
            story.append(_p("Strict spatial co-location membership: " +
                            ", ".join(combined["strict_spatial_colocation_members"])))
            _table(story, ["Residue", "Exact overlap N", "Exact fraction", "Co-location N", "Co-location fraction"],
                [[r["key"], r["exact_intersection_n"], r["exact_intersection_fraction"],
                  r["co_location_n"], r["co_location_fraction"]]
                 for r in combined["residues"]], [120, 103, 97, 102, 101])
        for key in consensus["map_order"]:
            m = consensus["maps"][key]
            story.append(PageBreak())
            story.append(Paragraph(escape(m["label"]), STYLES["ISheading"]))
            story.append(_p("Color shows fraction of input PDB structures where a residue belongs "
                            "to the Top-5 union for this channel. Not adsorption likelihood."))
            _screenshot(story, images[key])
            story.append(_p("Strict Top-1: " + ", ".join(m["strict_top1_members"])))
            story.append(_p("Strict Top-3: " + ", ".join(m["strict_top3_members"])))
            story.append(_p("Strict Top-5: " + ", ".join(m["strict_top5_members"])))
            story.append(_p("Majority Top-5: " + ", ".join(m["majority_top5_members"])))
            story.append(_p("Full residue membership frequencies"))
            _table(story, ["Residue", "Compatible", "Top-1", "Top-3", "Top-5"],
                [[r["key"], r["compatible_fraction"], r["top1_membership_fraction"],
                  r["top3_membership_fraction"], r["top5_membership_fraction"]]
                 for r in m["residues"]], [132, 105, 94, 94, 94])
    return _make_pdf("InterfaceScout | Multi-PDB consensus report", build)


def comparison_pdf(rows: list[dict[str, Any]]) -> bytes:
    def build(story: list):
        story.append(_p("Descriptive comparison of matching residues and ranked-region membership. "
                        "Jaccard overlap is not prediction accuracy."))
        for item in rows:
            story.append(Paragraph(
                escape(str(item["structure_a"]) + " versus " + str(item["structure_b"])),
                STYLES["ISheading"]))
            summary = item.get("accessibility", {})
            story.append(_p("Accessibility class changes: " +
                            str(summary.get("n_exposure_status_changes"))))
            _table(story, ["Map", "Center A", "Center B", "Top-1 J", "Top-3 J", "Top-5 J"],
                [[m["key"], m["top1_center_a"], m["top1_center_b"],
                  m["top_k"][0]["jaccard"], m["top_k"][1]["jaccard"],
                  m["top_k"][2]["jaccard"]] for m in item["map_comparisons"]],
                [103, 92, 91, 73, 74, 76])
            story.append(_p("All residue-level changes"))
            _table(story, ["Residue", "Exposed A/B", "scRSA A", "scRSA B", "Delta scRSA"],
                [[r["key"], f'{int(r["exposed_a"])}/{int(r["exposed_b"])}',
                  r["scRSA_a"], r["scRSA_b"], r["delta_scRSA"]]
                 for r in item["residue_changes"]], [120, 115, 95, 95, 89])
            story.append(PageBreak())
    return _make_pdf("InterfaceScout | Multi-PDB comparison report", build)


def _add_csv_and_xlsx(zipfile_: zipfile.ZipFile, prefix: str,
                      tables: Mapping[str, list[dict[str, Any]]]) -> None:
    for table_name, rows in tables.items():
        zipfile_.writestr(f"{prefix}/csv/{table_name}.csv", _csv_bytes(rows))
    zipfile_.writestr(f"{prefix}/complete_results.xlsx", _xlsx_bytes(tables))


def build_export_zip(payload: Mapping[str, Any]) -> bytes:
    structures = payload.get("structures", [])
    if not isinstance(structures, list) or not 1 <= len(structures) <= 20:
        raise ValueError("Export requires one to 20 previously analyzed PDB structures")
    consensus = payload.get("consensus")
    comparisons = payload.get("comparison") or []
    # Rebuild consensus from the exact per-structure analyses and currently
    # exported selected-map intersections. Never export a stale consensus
    # when the user changes the chosen maps after the initial multi-PDB run.
    if consensus is not None:
        if len(structures) < 2:
            raise ValueError("Consensus requires multiple structures")
        consensus = build_consensus([
            {
                "name": str(item.get("name") or f"Structure_{i + 1}"),
                "analysis": item["analysis"],
                "combined_maps": item.get("combined_maps"),
            }
            for i, item in enumerate(structures)
        ])
    master = BytesIO()
    used: set[str] = set()
    with zipfile.ZipFile(master, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for index, item in enumerate(structures):
            name = str(item.get("name") or f"Structure_{index + 1}")
            safe = safe_name(name)
            if safe in used:
                raise ValueError("Duplicate normalized structure names")
            used.add(safe)
            analysis = item["analysis"]
            if list(analysis.get("map_order", [])) != MAP_NAMES:
                raise ValueError("Analysis is missing one or more canonical 11 maps")
            shots = item.get("screenshots", {})
            images = {key: _image_bytes(shots[key]) for key in MAP_NAMES}
            pdb_text = item.get("pdb_text")
            if not isinstance(pdb_text, str) or not pdb_text:
                raise ValueError(f"Missing original PDB for {name}")
            recorded = analysis.get("input", {}).get("source_pdb_sha256")
            if recorded and hashlib.sha256(pdb_text.encode("utf-8")).hexdigest() != recorded:
                raise ValueError(
                    f"{name}: original PDB no longer matches the analyzed input; "
                    "run the analysis again before exporting results"
                )
            base = f"structures/{safe}"
            z.writestr(f"{base}/original.pdb", pdb_text)
            prepared_pdb = analysis.get("prepared_pdb")
            if not isinstance(prepared_pdb, str) or not prepared_pdb.strip():
                raise ValueError(f"{name}: missing analyzed prepared PDB for model-safe export")
            z.writestr(f"{base}/prepared_analyzed.pdb", prepared_pdb)
            z.writestr(f"{base}/analysis_complete.json",
                       json.dumps(analysis, ensure_ascii=False, indent=2))
            for i, key in enumerate(MAP_NAMES, start=1):
                subname = f"{i:02d}_{key}"
                z.writestr(f"{base}/screenshots/{subname}.png", images[key])
                z.writestr(f"{base}/pdf_maps/{subname}.pdf",
                           map_pdf(name, analysis, key, images[key]))
                z.writestr(f"{base}/pdb_bfactor/{subname}_Bfactor.pdb",
                           _bfactor_pdb(prepared_pdb, _map_scores(analysis, key),
                                        f"{key}: canonical propensity"))
                for k in (1, 3, 5):
                    z.writestr(f"{base}/pdb_bfactor/{subname}_Top{k}_Bfactor.pdb",
                               _bfactor_pdb(prepared_pdb,
                                            _topk_scores(analysis["maps"][key], k),
                                            f"{key}: Top-{k} region membership"))
            combined = item.get("combined_maps")
            if combined:
                z.writestr(f"{base}/combined/combined_maps.json",
                           json.dumps(combined, ensure_ascii=False, indent=2))
                z.writestr(f"{base}/combined/selected_exact_overlap_Bfactor.pdb",
                           _bfactor_pdb(prepared_pdb,
                               {key: 100.0 for key in combined.get("exact_same_residue_overlap", [])},
                               "exact intersection of selected channel chemistry seeds"))
                if item.get("combined_screenshot"):
                    blob = _image_bytes(item["combined_screenshot"])
                    z.writestr(f"{base}/combined/selected_combination.png", blob)
            _add_csv_and_xlsx(z, base, structure_tables(name, analysis, combined))
            z.writestr(f"{base}/full_report.pdf",
                       report_pdf(name, analysis, images, combined,
                                  _image_bytes(item['combined_screenshot'])
                                  if item.get('combined_screenshot') else None))

        if consensus:
            base = "consensus"
            if len(structures) < 2:
                raise ValueError("Consensus requires multiple PDB structures")
            consensus_images = payload.get("consensus_screenshots", {})
            images = {key: _image_bytes(consensus_images[key]) for key in MAP_NAMES}
            z.writestr(f"{base}/consensus_complete.json",
                       json.dumps(consensus, ensure_ascii=False, indent=2))
            _add_csv_and_xlsx(z, base, consensus_tables(consensus))
            representative = structures[0]["analysis"].get("prepared_pdb")
            if not representative:
                raise ValueError("Consensus export lacks the prepared reference structure")
            for i, key in enumerate(MAP_NAMES, start=1):
                label = f"{i:02d}_{key}"
                z.writestr(f"{base}/screenshots/{label}_Top5_frequency.png", images[key])
                z.writestr(f"{base}/pdb_bfactor/{label}_Top5_frequency_Bfactor.pdb",
                           _bfactor_pdb(representative, {
                               r["key"]: 100.0 * float(r["top5_membership_fraction"])
                               for r in consensus["maps"][key]["residues"]
                           }, f"{key}: Top-5 membership percentage across input structures"))
                z.writestr(f"{base}/pdb_bfactor/{label}_compatible_frequency_Bfactor.pdb",
                           _bfactor_pdb(representative, {
                               r["key"]: 100.0 * float(r["compatible_fraction"])
                               for r in consensus["maps"][key]["residues"]
                           }, f"{key}: compatible-residue percentage across input structures"))
            combined_consensus = consensus.get("combined_consensus")
            if combined_consensus:
                z.writestr(f"{base}/pdb_bfactor/selected_maps_exact_intersection_frequency_Bfactor.pdb",
                    _bfactor_pdb(representative, {
                        r["key"]: 100.0 * float(r["exact_intersection_fraction"])
                        for r in combined_consensus["residues"]
                    }, "selected maps: exact overlap percentage across structures"))
                z.writestr(f"{base}/pdb_bfactor/selected_maps_colocation_frequency_Bfactor.pdb",
                    _bfactor_pdb(representative, {
                        r["key"]: 100.0 * float(r["co_location_fraction"])
                        for r in combined_consensus["residues"]
                    }, "selected maps: spatial co-location percentage across structures"))
            z.writestr(f"{base}/consensus_report.pdf", consensus_pdf(consensus, images))

        if comparisons:
            base = "comparison"
            _add_csv_and_xlsx(z, base, comparison_tables(comparisons))
            z.writestr(f"{base}/comparison_report.pdf", comparison_pdf(comparisons))
            z.writestr(f"{base}/comparison_complete.json",
                       json.dumps(comparisons, ensure_ascii=False, indent=2))

        warnings = {
            "consensus_unavailable_reason": payload.get("consensus_unavailable_reason"),
            "comparison_skipped": payload.get("comparison_skipped", []),
        }
        z.writestr("analysis_warnings.json",
                   json.dumps(warnings, ensure_ascii=False, indent=2))
        if warnings["comparison_skipped"]:
            z.writestr("comparison/skipped_comparisons.csv",
                       _csv_bytes(list(warnings["comparison_skipped"])))
        z.writestr("README_EXPORT.txt", (
            "InterfaceScout research data package\n"
            "Every original PDB has 11 full map PDFs, viewer screenshots, "
            "11 map B-factor PDBs, Top-1/3/5 membership PDBs, all-residue CSVs, "
            "an Excel workbook, a full report PDF and raw JSON results.\n"
            "Consensus B-factor percentages are across supplied structures only "
            "and must not be read as adsorption probabilities.\n"
            "B-factor fields are display values, not experimental B-factors.\n"
            "No trajectory, GNM, ANM or RINpy analysis is performed.\n"
        ))
    return master.getvalue()
