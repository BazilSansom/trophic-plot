# Trophic Layout for Directed Network Visualisation

This repository contains the research code for the **Trophic Layout for Directed Network Visualisation** paper (forthecoming/insert).

The trophic layout method provides a principled and interpretable way to visualise **global direction, hierarchy, and feedback** in weighted directed networks. It is built on *generalised trophic levels* (MacKay et al., [2020](https://doi.org/10.1098/rsos.201138)), a quantitative measure of how “upstream” or “downstream” each node is in the flow structure of the network.

In this approach:

- Nodes are placed vertically according to their trophic level, revealing the network’s **directional backbone**.
- Horizontal positioning is optimised to reduce edge crossings while preserving flow structure.
- Feedback loops and incoherent flow patterns naturally appear as **horizontal** or **backwards** edges.
- The resulting plots offer a consistent way to compare directionality, identify hierarchy, and expose cycles across different networks.

This repository provides MATLAB (and forthcoming Python) implementations of:

- Trophic-level computation and trophic incoherence measures  
- Single-component and multi-component layout algorithms
- Optional layout refinements (e.g. automatic barycentre sweeps with x-smoothing for coherent networks).  
- Publication-quality visualisation via `plotTFL`  


### Main user-facing MATLAB functions

- **`tflPlot`** — Runs the full layout pipeline and produces a plot
- **`trophicLayoutMulti`** — General layout funciton -> returns layout for networks with multiple components
- **`trophicLayoutSoft`** — Core layout function -> returns layout for a single strongly connected component  
- **`plotTFL`** - Main network plotting function
- **`trophic_levels`** — Trophic analysis function -> returns trophic levels + trophic incoherence

All other `.m` files in `matlab/` are internal helpers used by these core functions.

Python equivalents will follow the same API structure.

## Repository Structure

```
tfl-layout-paper/
├─ README.md
├─ LICENSE
│
├─ matlab/
│   ├─ core/                 % core MATLAB implementation
│   ├─ plotting/
│   ├─ utils/
│   ├─ examples/
│   ├─ tests/
│   └─ data/                % MATLAB-friendly example data (MAT, CSV)
│
├─ python/
│   ├─ tfl/                 % python package (future)
│   │    ├─ core/
│   │    ├─ plotting/
│   │    ├─ utils/
│   │    └─ __init__.py
│   ├─ examples/
│   ├─ tests/
│   └─ data/                % Python-friendly data (CSV, JSON, NumPy arrays)
│
└─ shared_data/             % Optional: if large datasets used by both

```

## Requirements

### MATLAB
- MATLAB R20xx (tested with R20xxb and later)
- Uses built-in `graph`, `digraph`, `conncomp`, `subgraph`, `eigs`, etc.
- No additional toolboxes required.

Add the code to your MATLAB path:

```matlab
addpath(genpath('matlab'));
```

## Citing

If you use this code, please cite:

- The Trophic Layout paper (details forthcoming)
- MacKay, R.S., Johnson, S., Sansom, B. (2020). How directed is a directed network? (Royal Society Open Science) [https://doi.org/10.1098/rsos.201138](https://doi.org/10.1098/rsos.201138)

A CITATION.cff file will be added once publication details are final.


## MIT License

Copyright (c) <YEAR>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software...

