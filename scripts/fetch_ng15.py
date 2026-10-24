"""Download the NANOGrav 15-year HD free-spectrum KDEs (PTArcade data, Zenodo 10495907).

Writes data/ng15/ng15_30f_fs{hd}_ceffyl/{density,log10rhogrid,freqs,bandwidths}.npy and the text
files, after checking the archive's MD5 against the Zenodo record. Only plain .npy/.txt members
are extracted (no pickles; the arrays are loaded with allow_pickle=False).

Usage: python scripts/fetch_ng15.py [--dest data/ng15]
"""

import argparse
import urllib.request
import zipfile
from pathlib import Path

from ngpbh.pta import DATA_MD5, DATA_SUBDIR, DATA_URL, DEFAULT_DIR, md5sum

ALLOWED = {"density.npy", "log10rhogrid.npy", "freqs.npy", "bandwidths.npy", "log10rholabels.txt",
           "log.txt", "pulsar_list.txt"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=str(DEFAULT_DIR))
    args = ap.parse_args()
    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)
    archive = dest / "ng15_30f_fs_hd_ceffyl.zip"
    if not archive.exists() or md5sum(archive) != DATA_MD5:
        print(f"downloading {DATA_URL}")
        urllib.request.urlretrieve(DATA_URL, archive)
    digest = md5sum(archive)
    if digest != DATA_MD5:
        raise SystemExit(f"MD5 mismatch: {digest} != {DATA_MD5}")
    out = dest / DATA_SUBDIR
    out.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as zf:
        for name in zf.namelist():
            parts = Path(name).parts
            if len(parts) == 2 and parts[0] == DATA_SUBDIR and parts[1] in ALLOWED:
                (out / parts[1]).write_bytes(zf.read(name))
    print(f"wrote {sorted(p.name for p in out.iterdir())} to {out}")


if __name__ == "__main__":
    main()
