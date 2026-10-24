#!/usr/bin/env python3
"""Download pinned public inputs; validate all files before replacing any of them."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def download(destination, verify_only=False):
    sources = json.loads((ROOT / "data/sources.json").read_text())
    pending = {}
    for name, source in sources.items():
        target = destination / name
        if target.exists():
            content = target.read_bytes()
        elif verify_only:
            raise FileNotFoundError(target)
        else:
            with urlopen(source["url"], timeout=60) as response:
                content = response.read()
        if hashlib.sha256(content).hexdigest() != source["sha256"]:
            raise ValueError(f"SHA-256 mismatch: {name}; no files written")
        pending[name] = content
    if not verify_only:
        destination.mkdir(parents=True, exist_ok=True)
        for name, content in pending.items():
            if not (destination / name).exists():
                temporary = destination / (name + ".tmp")
                temporary.write_bytes(content)
                temporary.replace(destination / name)
    print(f"Verified {len(pending)}/{len(sources)} inputs against pinned SHA-256 hashes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=ROOT / "data/raw")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    download(args.destination, args.verify_only)
