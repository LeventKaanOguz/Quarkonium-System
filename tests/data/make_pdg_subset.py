#!/usr/bin/env python3
"""Cut the small PDG fixture used by the tests (and so by the CI) out of the
full PDG SQLite file.

    python3 tests/data/make_pdg_subset.py data/pdg-2026.0.sqlite tests/data/pdg-subset.sqlite

The fixture keeps, with the original schema, every row of pdgparticle, pdgid
and pdgdata that belongs to a particle listed in config/*.states (the c cbar,
b bbar and B_c states the fits use). The reader (src/data/pdg.cpp) touches
nothing else, so every number it returns from the fixture equals the number
from the full file. The data are the PDG's, CC BY 4.0 (see README.md here).
"""
import re
import sqlite3
import sys
from pathlib import Path

TABLES = ("pdgparticle", "pdgid", "pdgdata")


def mcids(config_dir):
    out = set()
    for f in sorted(Path(config_dir).glob("*.states")):
        for line in f.read_text().splitlines():
            line = line.split("#")[0].strip()
            if line:
                out.add(int(line.split()[0]))
    return out


def main(src, dst, config_dir):
    ids = sorted(mcids(config_dir))
    if Path(dst).exists():
        Path(dst).unlink()
    s, d = sqlite3.connect(src), sqlite3.connect(dst)
    for t in TABLES:
        ddl = s.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name=?",
                        (t,)).fetchone()[0]
        d.execute(ddl)
    q = ",".join("?" * len(ids))
    parts = s.execute(f"SELECT * FROM pdgparticle WHERE mcid IN ({q})", ids).fetchall()
    cols = [r[1] for r in s.execute("PRAGMA table_info(pdgparticle)")]
    pid = cols.index("pdgid")
    d.executemany(f"INSERT INTO pdgparticle VALUES ({','.join('?' * len(cols))})", parts)
    parents = [p[pid] for p in parts]
    qp = ",".join("?" * len(parents))
    rows = s.execute(f"SELECT * FROM pdgid WHERE parent_pdgid IN ({qp})", parents).fetchall()
    cols_i = [r[1] for r in s.execute("PRAGMA table_info(pdgid)")]
    d.executemany(f"INSERT INTO pdgid VALUES ({','.join('?' * len(cols_i))})", rows)
    idc = cols_i.index("id")
    keep = [r[idc] for r in rows]
    n_data = 0
    for i in range(0, len(keep), 500):
        chunk = keep[i:i + 500]
        qq = ",".join("?" * len(chunk))
        data = s.execute(f"SELECT * FROM pdgdata WHERE pdgid_id IN ({qq})", chunk).fetchall()
        cols_d = [r[1] for r in s.execute("PRAGMA table_info(pdgdata)")]
        d.executemany(f"INSERT INTO pdgdata VALUES ({','.join('?' * len(cols_d))})", data)
        n_data += len(data)
    d.commit()
    d.execute("VACUUM")
    d.close()
    print(f"{len(parts)} particles, {len(rows)} pdgid rows, {n_data} pdgdata rows "
          f"({Path(dst).stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    src, dst = sys.argv[1:3]
    main(src, dst, sys.argv[3] if len(sys.argv) > 3 else str(Path(__file__).parents[2] / "config"))
