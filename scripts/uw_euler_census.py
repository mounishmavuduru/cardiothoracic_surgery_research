"""Full Euler-characteristic census over every UW pre-ablation mesh.

The manuscript states that removing tag 164 drops the Euler characteristic "by a
median of 7 (range 2-12)" and, in the same paragraph, names ID040 and ID049 as
having odd topology before removal.  The source comment behind that figure says
those two "go the other way", and the regression test skips them -- so the quoted
range is over 80 meshes, not the 82 the sentence appears to cover.  Its own
docstring says both "measured over all 82" and "they are excluded from this
guard", which cannot both be true.

This recomputes the census over all 82 with nothing excluded, so the manuscript
can state the real distribution and the sign of the two exceptions.  Written
because this exact paragraph is the one the paper already corrects once for
over-generalising from a single mesh; a second over-generalisation in the
correction itself would be worse than the original.

Writes results/uw_euler_census.json.
"""
from __future__ import annotations

import glob
import json
import os
import statistics as st
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from asb.substrate.uw_boyle import read_vtk_unstructured_bin  # noqa: E402

OUT = "results/uw_euler_census.json"
DROP_TAG = 164


def _euler_characteristic(faces: np.ndarray) -> int:
    """chi = V - E + F over the vertices actually referenced by ``faces``.

    Identical to the helper in tests/test_uw_boyle.py, reproduced here so the census
    and the regression guard cannot drift apart silently.
    """
    v = np.unique(faces).size
    e = np.unique(np.sort(np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0).shape[0]
    return int(v - e + faces.shape[0])


def main() -> None:
    meshes = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))
    if not meshes:
        raise SystemExit("no UW pre-ablation meshes on disk")

    rows = []
    for path in meshes:
        parsed = read_vtk_unstructured_bin(path)
        faces = np.asarray(parsed["faces"])
        tags = np.asarray(parsed["cell_scalars"]["elemTag"]).astype(int)
        before = _euler_characteristic(faces)
        after = _euler_characteristic(faces[tags != DROP_TAG])
        rows.append({
            "subject": os.path.basename(path),
            "chi_before": int(before),
            "chi_after": int(after),
            "delta_chi": int(before - after),
        })

    deltas = [r["delta_chi"] for r in rows]
    falls = [r for r in rows if r["delta_chi"] > 0]
    rises = [r for r in rows if r["delta_chi"] < 0]
    flat = [r for r in rows if r["delta_chi"] == 0]
    fall_d = [r["delta_chi"] for r in falls]

    out = {
        "description": (
            "Euler characteristic before and after removing tag 164, over every UW "
            "pre-ablation mesh with nothing excluded. chi_before is the surface as "
            "released; delta_chi = chi_before - chi_after, so a positive value means "
            "removal opened the surface."
        ),
        "drop_tag": DROP_TAG,
        "n_meshes": len(rows),
        "all_meshes": {
            "delta_chi_median": st.median(deltas),
            "delta_chi_min": min(deltas),
            "delta_chi_max": max(deltas),
        },
        "n_falls": len(falls),
        "n_rises": len(rises),
        "n_unchanged": len(flat),
        "falling_subset": {
            "n": len(falls),
            "delta_chi_median": st.median(fall_d) if fall_d else None,
            "delta_chi_min": min(fall_d) if fall_d else None,
            "delta_chi_max": max(fall_d) if fall_d else None,
        },
        "exceptions": [r for r in rows if r["delta_chi"] <= 0],
        "chi_before_min": min(r["chi_before"] for r in rows),
        "chi_before_max": max(r["chi_before"] for r in rows),
        "rows": rows,
    }

    os.makedirs("results", exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    print(f"UW Euler census over {len(rows)} meshes (nothing excluded)\n")
    print(f"  delta chi over ALL      : median {out['all_meshes']['delta_chi_median']}, "
          f"range {out['all_meshes']['delta_chi_min']} to {out['all_meshes']['delta_chi_max']}")
    print(f"  opened  (delta > 0)     : {len(falls)}")
    print(f"  closed  (delta < 0)     : {len(rises)}")
    print(f"  unchanged (delta == 0)  : {len(flat)}")
    if fall_d:
        print(f"  delta chi among openers : median {out['falling_subset']['delta_chi_median']}, "
              f"range {out['falling_subset']['delta_chi_min']} to {out['falling_subset']['delta_chi_max']}")
    print(f"  chi before removal      : {out['chi_before_min']} to {out['chi_before_max']}")
    if out["exceptions"]:
        print("\n  exceptions (removal does not open the surface):")
        for r in out["exceptions"]:
            print(f"    {r['subject']}: chi {r['chi_before']} -> {r['chi_after']} "
                  f"(delta {r['delta_chi']})")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
