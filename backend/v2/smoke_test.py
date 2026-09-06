"""Minimal non-network smoke tests for the current coarse InterfaceScout V2."""

from __future__ import annotations

import numpy as np

from .coarse_patch import build_coarse_patches
from .surface_modes import get_surface_mode


def _fake_inputs():
    keys = [f"A:{i}:" for i in range(1, 7)]
    coords = [
        np.array([0.0, 0.0, 0.0]),
        np.array([3.8, 0.0, 0.0]),
        np.array([7.5, 0.0, 0.0]),
        np.array([22.0, 0.0, 0.0]),
        np.array([25.8, 0.0, 0.0]),
        np.array([29.5, 0.0, 0.0]),
    ]
    nodes = [
        {"key": k, "chain": "A", "res_seq": i + 1, "icode": "", "res_name": "LYS", "coord": coords[i]}
        for i, k in enumerate(keys)
    ]
    corr = np.eye(len(keys), dtype=float)
    for i in range(len(keys) - 1):
        corr[i, i + 1] = corr[i + 1, i] = 0.5
    gnm = {
        "nodes": nodes,
        "index": {k: i for i, k in enumerate(keys)},
        "correlation_matrix": corr,
    }
    surface = [
        {"key": k, "scrsa": 0.5, "scrsa_raw": 0.5}
        for k in keys
    ]
    residues = [
        {"key": k, "local_score": 0.5, "multiscale_persistence": p}
        for k, p in zip(keys, [1.0, 0.8, 0.5, 1.0, 0.7, 0.4])
    ]
    centers = [
        {"center_key": k, "res_name": "LYS", "res_seq": i + 1, "icode": "", "chain": "A",
         "multiscale_persistence": p * 100.0, "multiscale_geomean": p * 100.0}
        for i, (k, p) in enumerate(zip(keys, [1.0, 0.8, 0.5, 0.9, 0.7, 0.4]))
    ]
    v1_result = {
        "surface_residues": surface,
        "chemistries": {
            "anionic": {"residues": residues, "patch_centers": centers}
        },
    }
    return v1_result, gnm


def run_smoke_test() -> None:
    assert get_surface_mode("silica").chemistry == "anionic"
    v1_result, gnm = _fake_inputs()
    patches = build_coarse_patches(v1_result=v1_result, chemistry="anionic", gnm=gnm)
    assert len(patches) >= 2
    assert all(p["n_members"] <= 3 for p in patches), "patch growth must remain local/non-transitive"
    assert all("dynamic_coupling_abs" in p for p in patches)
    assert all("orientation_coherence" in p for p in patches)
    assert all("patch_coherence" in p for p in patches)
    assert all("pareto_front" in p for p in patches)


if __name__ == "__main__":
    run_smoke_test()
    print("InterfaceScout V2 coarse smoke test: OK")
