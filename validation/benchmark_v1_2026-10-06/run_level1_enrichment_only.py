from __future__ import annotations
import json, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import run_level1_scoring as m

def main():
    proteins=m.parse_csv(m.BENCH/"proteins.csv")
    pcdb,total=m.resolve_pcdb_names(proteins)
    payne=m.resolve_payne(proteins)
    identities=pcdb+payne
    enriched=m.enrich_identities_batch(identities)
    enriched.sort(key=lambda r:(r["source_database"],r["source_identifier"]))
    m.write_csv(m.OUT/"level1_structure_enrichment_only.csv",enriched)
    audit={
        "human_uniprot_entries_downloaded":total,
        "pcdb_source_entities":len(pcdb),
        "pcdb_identity_resolved":sum(bool(r.get("uniprot_accession")) for r in pcdb),
        "pcdb_unique_exact":sum(r.get("identity_resolution_status")=="unique_exact_normalized" for r in pcdb),
        "pcdb_unique_reviewed_exact":sum(r.get("identity_resolution_status")=="unique_reviewed_exact_normalized" for r in pcdb),
        "pcdb_ambiguous_exact":sum(r.get("identity_resolution_status")=="ambiguous_exact_normalized" for r in pcdb),
        "pcdb_no_exact":sum(r.get("identity_resolution_status")=="no_exact_normalized_match" for r in pcdb),
        "payne_source_entities":len(payne),
        "payne_identity_resolved":sum(bool(r.get("uniprot_accession")) for r in payne),
        "payne_identity_unresolved":sum(not bool(r.get("uniprot_accession")) for r in payne),
        "experimental_pdb_selected":sum(r.get("structure_origin")=="experimental_pdb" for r in enriched),
        "alphafold_selected":sum(r.get("structure_origin")=="alphafold_predicted" for r in enriched),
        "no_structure_or_identity":sum(not r.get("structure_origin") for r in enriched),
        "scoreable_total":sum(bool(r.get("structure_origin")) for r in enriched),
    }
    (m.OUT/"level1_enrichment_only_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(json.dumps(audit,indent=2))
if __name__=="__main__":
    main()
