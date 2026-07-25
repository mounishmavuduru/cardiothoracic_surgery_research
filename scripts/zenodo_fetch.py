"""Download a Zenodo record's files (no proof-of-work wall, unlike Dryad).

Used for the Roney LA virtual cohort, Zenodo 5801337 (CC-BY-4.0), which ships the
continuous LGE ``IIR`` field and real ``fiber_endo``/``fiber_epi`` vectors -- the
fibre-carrying control against which the UW/Boyle substrate defects are measured.

Usage:  python scripts/zenodo_fetch.py 5801337 data/roney [--limit N]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import requests


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("record")
    ap.add_argument("dest")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--pattern", default="")
    a = ap.parse_args()

    meta = requests.get(f"https://zenodo.org/api/records/{a.record}", timeout=120).json()
    files = meta.get("files", [])
    if a.pattern:
        files = [f for f in files if a.pattern in f["key"]]
    # Match the loader's deterministic ordering: smallest first, then by name.
    files.sort(key=lambda f: (f["size"], f["key"]))
    if a.limit:
        files = files[: a.limit]

    os.makedirs(a.dest, exist_ok=True)
    total = sum(f["size"] for f in files)
    print(f"[zenodo {a.record}] {len(files)} files, {total/1e9:.2f} GB -> {a.dest}", flush=True)

    t0 = time.time()
    got = 0
    for i, f in enumerate(files, 1):
        out = os.path.join(a.dest, f["key"])
        if os.path.exists(out) and os.path.getsize(out) == f["size"]:
            got += f["size"]
            print(f"  [skip {i}/{len(files)}] {f['key']}", flush=True)
            continue
        with requests.get(f["links"]["self"], stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(out, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
        got += f["size"]
        el = time.time() - t0
        print(f"  [{i}/{len(files)}] {f['key']} {f['size']/1e6:.1f} MB   "
              f"{got/1e9:.2f}/{total/1e9:.2f} GB  {el/60:.1f} min", flush=True)

    print(f"ZENODO_DONE {len(files)} files in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    sys.exit(main())
