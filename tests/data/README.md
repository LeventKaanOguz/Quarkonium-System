# PDG fixture for the tests

`pdg-subset.sqlite` is a cut of the PDG 2026 SQLite file (`pdg-2026.0.sqlite`, pdg.lbl.gov/api): the
`pdgparticle`, `pdgid` and `pdgdata` rows of the c c̄, b b̄ and B_c states listed in `config/*.states`, with the
original schema. The reader in `src/data/pdg.cpp` touches nothing else, so it returns the same numbers from the
fixture as from the full file; that is what lets the fit tests run in the CI.

Data: Particle Data Group, *Review of Particle Physics*, 2026 edition, licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The subset is a modified extract of that data.

Remake it with `python3 tests/data/make_pdg_subset.py data/pdg-2026.0.sqlite tests/data/pdg-subset.sqlite`.
The analysis (`analysis/validity.py`, the sweeps) needs the full file in `data/`.
