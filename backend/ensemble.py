"""Multi-PDB analysis and descriptive consensus using the existing 11-map engine.

No trajectory, elastic network, GNM, ANM or residue-network dependency.
Structure frequencies reflect the supplied PDB *set*, not equilibrium populations.
"""
from __future__ import annotations

from itertools import combinations
from typing import Any, Mapping

from interface_engine import analyze_all_maps
from research_workflows import compare_results, combine_maps


def analyze_structure_set(
    structures: list[dict[str, Any]], *,
    pH: float = 7.4,
    ionic_mM: float = 150.0,
    temp_K: float = 298.0,
    combined_maps: list[str] | None = None,
) -> dict[str, Any]:
    if not 1 <= len(structures) <= 20:
        raise ValueError("Provide between 1 and 20 PDB structures")
    used: set[str] = set()
    runs: list[dict[str, Any]] = []
    for index, spec in enumerate(structures):
        name = str(spec.get("name") or f"Structure_{index + 1}").strip()
        if not name or name in used:
            raise ValueError("Each structure must have a unique nonempty name")
        used.add(name)
        try:
            analysis = analyze_all_maps(
                pdb_id=spec.get("pdb_id"),
                pdb_text=spec.get("pdb_text"),
                chain=spec.get("chain"),
                pH=pH, ionic_mM=ionic_mM, temp_K=temp_K,
            )
            combined = (
                combine_maps(analysis, combined_maps, max_regions=None)
                if combined_maps and len(combined_maps) >= 2
                else None
            )
        except Exception as exc:
            raise ValueError(f"{name}: {exc}") from exc
        runs.append({"name": name, "analysis": analysis, "combined_maps": combined})

    out = {
        "n_structures": len(runs),
        "structures": runs,
        "conditions": {"pH": pH, "ionic_mM": ionic_mM, "temp_K": temp_K},
        "map_order": runs[0]["analysis"]["map_order"],
        "comparison": [],
        "consensus": None,
        "disclaimer": (
            "Structural comparison is descriptive. Counts across provided structures "
            "are not equilibrium probabilities or kinetic predictions."
        ),
    }
    # Independent analyses always succeed for valid inputs. Consensus/comparison
    # are additional operations requiring matched residues, not prerequisites.
    out["comparison_skipped"] = []
    out["consensus_unavailable_reason"] = None
    if len(runs) > 1:
        try:
            out["consensus"] = build_consensus(runs)
        except ValueError as exc:
            out["consensus_unavailable_reason"] = str(exc)
        for a, b in combinations(runs, 2):
            try:
                item = compare_results(a["analysis"], b["analysis"])
            except ValueError as exc:
                out["comparison_skipped"].append({
                    "structure_a": a["name"], "structure_b": b["name"],
                    "reason": str(exc),
                })
            else:
                # compare_results() includes original input metadata under
                # structure_a/b. Keep it, but never let it override the names
                # used for pairwise reporting, CSV and PDF identity matching.
                out["comparison"].append({
                    **item,
                    "input_a": item.get("structure_a", {}),
                    "input_b": item.get("structure_b", {}),
                    "structure_a": a["name"],
                    "structure_b": b["name"],
                })
    return out


def _patch_members(map_row: Mapping[str, Any], top_k: int,
                   allowed: set[str]) -> set[str]:
    return {
        str(key)
        for patch in map_row.get("patches", [])[:top_k]
        for key in patch.get("members", [])
        if str(key) in allowed
    }


def _identity_matched(runs: list[dict[str, Any]]) -> tuple[set[str], list[dict[str, Any]]]:
    identities = [
        {str(r["key"]): str(r["res_name"]) for r in run["analysis"]["residue_profiles"]}
        for run in runs
    ]
    matched = set.intersection(*(set(d) for d in identities))
    matched = {key for key in matched if len({d[key] for d in identities}) == 1}
    details = []
    for i, d in enumerate(identities):
        coverage = len(matched) / max(len(d), 1)
        details.append({
            "structure": runs[i]["name"],
            "residues": len(d),
            "matched_positions": len(matched),
            "coverage": round(coverage, 5),
            "excluded_residues": sorted(set(d) - matched),
        })
        if coverage < 0.80:
            raise ValueError(
                "Consensus requires at least 80% identity-matched residue keys for every "
                f"structure; {runs[i]['name']} has {coverage:.1%}. Align chain IDs and "
                "residue numbering or analyze separately."
            )
    return matched, details


def build_consensus(runs: list[dict[str, Any]]) -> dict[str, Any]:
    if len(runs) < 2:
        raise ValueError("Consensus requires at least two structures")
    names = [r["name"] for r in runs]
    if len(set(names)) != len(names):
        raise ValueError("Structure names must be distinct")
    inputs = [r["analysis"].get("input", {}) for r in runs]
    if len({float(i.get("pH", 7.4)) for i in inputs}) != 1:
        raise ValueError("Use one common pH for structure consensus")
    map_order = list(runs[0]["analysis"]["map_order"])
    if any(list(r["analysis"]["map_order"]) != map_order for r in runs[1:]):
        raise ValueError("Inconsistent map definitions")

    allowed, mapping_details = _identity_matched(runs)
    n = len(runs)
    profiles = [
        {str(row["key"]): row for row in r["analysis"]["residue_profiles"]}
        for r in runs
    ]
    residue_results = []
    for key in sorted(allowed):
        rows = [r[key] for r in profiles]
        n_exposed = sum(bool(x["surface_exposed"]) for x in rows)
        scrsas = [float(x["scrsa_raw"]) for x in rows]
        residue_results.append({
            "key": key, "res_name": rows[0]["res_name"],
            "n_structures": n, "n_exposed": n_exposed,
            "exposure_fraction": round(n_exposed / n, 5),
            "scRSA_min": round(min(scrsas), 5),
            "scRSA_max": round(max(scrsas), 5),
            "scRSA_mean": round(sum(scrsas) / n, 5),
            "scRSA_range": round(max(scrsas) - min(scrsas), 5),
            "strictly_exposed": n_exposed == n,
            "exposed_in_any": n_exposed > 0,
        })

    consensus_maps: dict[str, dict[str, Any]] = {}
    for name in map_order:
        map_rows = [r["analysis"]["maps"][name] for r in runs]
        compatible = [
            {str(x["key"]) for x in m["residues"] if float(x.get("local_score") or 0) > 0}
            & allowed
            for m in map_rows
        ]
        at_top = {k: [_patch_members(m, k, allowed) for m in map_rows]
                  for k in (1, 3, 5)}
        rrows = []
        for key in sorted(allowed):
            count = sum(key in cs for cs in compatible)
            record = {
                "key": key, "res_name": profiles[0][key]["res_name"],
                "compatible_in_n": count, "compatible_fraction": round(count / n, 5),
                "compatible_in_all": count == n,
            }
            for k in (1, 3, 5):
                c = sum(key in ps for ps in at_top[k])
                record[f"top{k}_membership_n"] = c
                record[f"top{k}_membership_fraction"] = round(c / n, 5)
            rrows.append(record)

        def strict(k: int) -> list[str]:
            return [row["key"] for row in rrows
                    if row[f"top{k}_membership_n"] == n]

        consensus_maps[name] = {
            "label": map_rows[0]["label"],
            "n_structures": n,
            "residues": rrows,
            "strict_compatible_residues": sorted(set.intersection(*compatible)),
            "any_compatible_residues": sorted(set.union(*compatible)),
            "strict_top1_members": strict(1),
            "strict_top3_members": strict(3),
            "strict_top5_members": strict(5),
            "majority_top1_members": [row["key"] for row in rrows
                if row["top1_membership_n"] >= (n // 2 + 1)],
            "majority_top3_members": [row["key"] for row in rrows
                if row["top3_membership_n"] >= (n // 2 + 1)],
            "majority_top5_members": [row["key"] for row in rrows
                if row["top5_membership_n"] >= (n // 2 + 1)],
            "top1_centers_by_structure": {
                run["name"]: (
                    map_row["patches"][0]["center_key"]
                    if map_row.get("patches") else None
                )
                for run, map_row in zip(runs, map_rows)
            },
            "note": "Strict/majority consensus refers to RESIDUE MEMBERSHIP in ranked "
                    "patch sets. It does not imply one spatially contiguous consensus patch.",
        }

    combined_consensus = None
    combinations = [run.get("combined_maps") for run in runs]
    if all(combinations):
        selections = [tuple(c.get("selected_maps", [])) for c in combinations]
        if len(set(selections)) == 1 and len(selections[0]) >= 2:
            exact_sets = [
                set(c.get("exact_same_residue_overlap", [])) & allowed
                for c in combinations
            ]
            spatial_sets = [
                {
                    str(key)
                    for region in c.get("regions", [])
                    if region.get("all_maps_colocated")
                    for key in region.get("member_keys", [])
                    if str(key) in allowed
                }
                for c in combinations
            ]
            rows = []
            for key in sorted(allowed):
                exact_n = sum(key in entry for entry in exact_sets)
                spatial_n = sum(key in entry for entry in spatial_sets)
                rows.append({
                    "key": key,
                    "res_name": profiles[0][key]["res_name"],
                    "exact_intersection_n": exact_n,
                    "exact_intersection_fraction": round(exact_n / n, 5),
                    "co_location_n": spatial_n,
                    "co_location_fraction": round(spatial_n / n, 5),
                })
            combined_consensus = {
                "selected_maps": list(selections[0]),
                "residues": rows,
                "strict_exact_intersection": [
                    row["key"] for row in rows if row["exact_intersection_n"] == n
                ],
                "majority_exact_intersection": [
                    row["key"] for row in rows if row["exact_intersection_n"] >= (n // 2 + 1)
                ],
                "strict_spatial_colocation_members": [
                    row["key"] for row in rows if row["co_location_n"] == n
                ],
                "majority_spatial_colocation_members": [
                    row["key"] for row in rows if row["co_location_n"] >= (n // 2 + 1)
                ],
                "note": (
                    "Co-localization frequency is descriptive membership in at least "
                    "one same-patch multi-chemistry candidate per structure, "
                    "NOT an experimentally observed interface or one persistent 3D region."
                ),
            }

    return {
        "n_structures": n,
        "structure_names": names,
        "n_matched_residues": len(allowed),
        "matching": mapping_details,
        "residues": residue_results,
        "combined_consensus": combined_consensus,
        "maps": consensus_maps,
        "map_order": map_order,
        "method": (
            "Identity-matched residue positions across user-supplied structures; "
            "scRSA and each original 11-channel map re-evaluated separately. "
            "A residue's consensus frequency is the number of input structures "
            "where that feature is present divided by the total number of structures."
        ),
        "limits": (
            "Input structures are not necessarily independent, equally populated, or "
            "experimental replicates. Frequencies must not be interpreted as "
            "thermodynamic populations, adsorption likelihoods, or binding affinities."
        ),
    }
