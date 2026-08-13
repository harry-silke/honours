"""compare_sets.py — diff frame coverage between the external drive and OneDrive."""

import re
import done
from pathlib import Path


LOCAL = Path("/Volumes/Expansion/Heat_stress_experiment_photos/07-05-2026")
CLOUD = Path.home() / "Library/CloudStorage/OneDrive-AustralianNationalUniversity/Honours/kelp_photos_070525"

FRAME_RE = re.compile(r"DSC(\d+)")




odd = [p.name for p in CLOUD.rglob("*.ARW")
       if not p.name.startswith("._") and not FRAME_RE.search(p.stem)]
print(len(odd), odd[:20])

def frames(root):
    found = {}
    for p in root.rglob("*.ARW"):
        if p.name.startswith("._"):
            continue
        m = FRAME_RE.search(p.stem)
        if m:
            found.setdefault(int(m.group(1)), []).append(p)
    return found


def runs(nums):
    """Collapse a sorted list of ints into contiguous ranges."""
    out = []
    for n in sorted(nums):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return out


def main():
    local, cloud = frames(LOCAL), frames(CLOUD)
    print(f"local {len(local)} frames, cloud {len(cloud)} frames")

    only_local = set(local) - set(cloud)
    only_cloud = set(cloud) - set(local)
    both = set(local) & set(cloud)
    print(f"{len(both)} in both, {len(only_local)} local only, {len(only_cloud)} cloud only")

    for label, s in (("local only", only_local), ("cloud only", only_cloud)):
        if not s:
            continue
        print(f"\n{label}:")
        for lo, hi in runs(s)[:20]:
            span = f"DSC{lo:05d}" if lo == hi else f"DSC{lo:05d}-DSC{hi:05d}"
            print(f"  {span} ({hi - lo + 1})")

    dupes = {n: v for n, v in local.items() if len(v) > 1}
    if dupes:
        print(f"\n{len(dupes)} frames appear more than once locally")


if __name__ == "__main__":
    main()