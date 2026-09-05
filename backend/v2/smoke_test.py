"""Minimal non-network smoke tests for InterfaceScout V2-alpha."""

from .patches import assemble_candidate_patches, nonredundant_top_patches
from .scoring import rank_patches
from .surface_profiles import get_surface_profile


def run_smoke_test() -> None:
    profile = get_surface_profile("polystyrene")
    w = profile.normalized_weights()
    assert abs(sum(w.values()) - 1.0) < 1e-9

    fake = {
        "hydrophobic": {
            "patch_centers": [
                {"center_key": "A:10:", "res_name": "LEU", "res_seq": 10, "icode": "", "chain": "A", "multiscale_persistence": 80.0, "multiscale_geomean": 82.0, "compatible_members_8A": ["A:10:", "A:11:"]},
                {"center_key": "A:30:", "res_name": "VAL", "res_seq": 30, "icode": "", "chain": "A", "multiscale_persistence": 50.0, "multiscale_geomean": 55.0, "compatible_members_8A": ["A:30:"]},
            ]
        },
        "pi_carbon": {
            "patch_centers": [
                {"center_key": "A:10:", "res_name": "LEU", "res_seq": 10, "icode": "", "chain": "A", "multiscale_persistence": 20.0, "multiscale_geomean": 25.0, "compatible_members_8A": ["A:12:"]},
                {"center_key": "A:30:", "res_name": "VAL", "res_seq": 30, "icode": "", "chain": "A", "multiscale_persistence": 90.0, "multiscale_geomean": 91.0, "compatible_members_8A": ["A:31:"]},
            ]
        },
    }
    candidates = assemble_candidate_patches(fake, w)
    ranked = rank_patches(candidates)
    assert ranked[0]["center_key"] == "A:10:"
    top = nonredundant_top_patches(ranked, top_n=2)
    assert len(top) == 2


if __name__ == "__main__":
    run_smoke_test()
    print("V2-alpha smoke test: OK")
