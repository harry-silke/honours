import csv, os
from datetime import datetime, timezone

root = os.path.expanduser("~/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Honours/kelp_photos_070525")

def iso(ts):
    return datetime.fromtimestamp(ts, timezone.utc).astimezone().isoformat()

with open("photos.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "created", "modified"])
    for dirpath, _, _ in os.walk(root):
        with os.scandir(dirpath) as it:
            for e in it:
                if e.is_file(follow_symlinks=False) and e.name.lower().endswith(".arw"):
                    st = e.stat(follow_symlinks=False)
                    w.writerow([e.name, iso(st.st_birthtime), iso(st.st_mtime)])