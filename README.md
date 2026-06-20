# plume-rmt-deep-earth

**Level Repulsion in the Deep Earth: A Random Matrix Theory Analysis of Mantle-Plume Eruption Rhythms**

Source isolation recovers GOE statistics that spectral superposition conceals.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

## Core Finding

The **thermal-shadow hypothesis** holds — a single deep plume depletes the
local thermal boundary layer at the core-mantle boundary and must re-accumulate
heat before its next eruption, imposing level repulsion — **but only after
source isolation.**

| Configuration | n | ⟨r⟩ | Result |
|---|---|---|---|
| **Single-source hotspot tracks (pooled)** | 73 | **0.630** | **GOE repulsion** — Poisson rejected, CI [0.422, 0.563] |
| Pacific LLSVP (Jason) | 19 | 0.461 | intermediate |
| African LLSVP (Tuzo) | 23 | 0.331 | clustered (CV=1.39, p=0.017) |
| Global LIP mixed | 42 | 0.397 | Poisson |

**A single deep plume shows GOE level repulsion. Mixing many plumes collapses
to Poisson** by the spectral superposition theorem. The recovery of repulsion
upon source isolation is itself the evidence that the signal is real.

## Four Single-Source Hotspot Tracks
Hawaii–Emperor, Louisville, Tristan–Gough, Rurutu–Arago — each fed by one
deep-rooted mantle plume.

> Note: the Aleutian Islands are **excluded** — they are a subduction-related
> volcanic arc, not a plume source. (The Emperor chain's oldest end subducts
> beneath the Aleutian Trench and is already covered by the Hawaii–Emperor track.)

## The Two Superplumes Differ
- **African (Tuzo):** significantly clustered (CV=1.39, p=0.017, γ=0.61, Fano=4.08)
  — pulse-like plume-head eruptions (e.g., 5 LIPs within 16 Myr at 118–134 Ma)
- **Pacific (Jason):** intermediate / mild regularity

## Repository Structure
```
plume-rmt-deep-earth/
├── README.md · LICENSE · requirements.txt · CITATION.cff · .zenodo.json
├── paper/
│   ├── paper.tex · paper.pdf      # 7 pp.
│   └── figs/                      # figures embedded by LaTeX
├── code/
│   ├── lip_rmt_pipeline.py        # Experiments A (global) + B (LLSVP domains)
│   └── hotspot_analysis.py        # 4 single-source tracks + pooled GOE test
├── data/
│   └── lip_database.csv           # 44 LIPs, 17–1880 Ma, LLSVP-assigned
├── figures/                       # 4 standalone PDFs
└── results/                       # JSON outputs
```

## Method
1. Compile LIP absolute ages (Ernst 2014) + hotspot 40Ar/39Ar ages
2. Assign LIPs to African/Pacific LLSVP via plate reconstruction (Torsvik 2010)
3. Unfold inter-eruption intervals, compute ⟨r⟩, KS tests, CV, bootstrap CI
4. Contrast mixed-source (→ Poisson) vs single-source (→ GOE)

## Reproduce
```bash
pip install -r requirements.txt
cd code
python lip_rmt_pipeline.py     # Experiments A + B
python hotspot_analysis.py     # Four single-source tracks + pooled GOE
```

## Three-Racetrack RMT Program
This is the third study in a unified program. A theme runs through all three:
**GOE repulsion appears wherever a single long-memory process is isolated;
Poisson appears wherever many independent sources superpose.**

1. Geological boundaries (Myr) → GOE repulsion — [zenodo 20766310](https://zenodo.org/records/20766310)
2. Seismotectonic rhythms → scale-dependent — [zenodo 20768130](https://zenodo.org/records/20768130)
3. Mantle plumes (this work) → single-source GOE, superposition-concealed

## Citation
```bibtex
@misc{chen2026plume,
  author = {Chen, Ruqing},
  title  = {Level Repulsion in the Deep Earth: A Random Matrix Theory
            Analysis of Mantle-Plume Eruption Rhythms},
  year   = {2026},
  publisher = {GitHub},
  url    = {https://github.com/Ruqing1963/plume-rmt-deep-earth}
}
```

## License
[MIT](LICENSE). LIP ages from Ernst (2014); hotspot ages from published 40Ar/39Ar studies.
