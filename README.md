# Trophic Layout for Directed Network Visualisation

This repository contains research code for **Trophic Layout for Directed Network Visualisation**.

The trophic layout method provides a principled and interpretable way to visualise **global direction, hierarchy, and feedback** in weighted directed networks. It is built on the notion of *generalised trophic levels* introduced by MacKay, Johnson, and Sansom ([2020](https://doi.org/10.1098/rsos.201138)), which quantify how “upstream” or “downstream” each node is in the directed flow structure of a network. The repository also includes visualisation tools related to recent work on directional and circular structure in weighted directed networks, including Circular Directional Flow Decomposition (CDFD) by Homs-Dones, MacKay, Sansom, and Zhou ([2025, arXiv:2506.12546](https://arxiv.org/abs/2506.12546)).

In this approach:

- nodes are placed vertically according to trophic level, revealing the network’s **directional backbone**
- horizontal positions are optimised for readability while preserving directional structure
- feedback loops and incoherent flow patterns appear naturally as **horizontal** or **backward** edges
- circular flow structure can be visually highlighted using the balanced-flow-forwarding (BFF) CDFD of Homs-Dones et al. (2025)
- the resulting layouts provide a consistent way to compare hierarchy, directionality, and cyclic structure across networks

This repository currently provides a MATLAB implementation of:

- trophic-level computation and related directionality diagnostics, including trophic incoherence and CDFD-based directionality/circularity measures
- single-component and multi-component trophic layout algorithms
- publication-quality directed-network visualisation
- circular–directional flow decomposition (CDFD) plotting utilities used in the associated paper workflow

Python support is planned, but is not yet included in this repository.

## Main MATLAB functions

The main user-facing MATLAB functions are:

- **`tflPlot`** — run the full trophic layout pipeline and produce a plot
- **`trophicLayoutMulti`** — compute a layout for a directed network with multiple weakly connected components
- **`trophicLayout`** — compute a trophic layout for a single connected network or component
- **`plotTFL`** — render a network from supplied coordinates using the TFL plotting layer
- **`trophic_levels`** — compute generalised trophic levels and, optionally, incoherence/coherence diagnostics
- **`plotCDFD`** — visualise a circular–directional flow decomposition
- **`cdfd_bff`** — compute a balanced-flow-forwarding CDFD decomposition

Additional analysis and plotting helpers are included under `toolbox/src/` and are used internally by these main functions.

## Repository structure

```text
trophic-plot/
├─ toolbox/
│  └─ src/
│     ├─ analysis/      % trophic analysis, incoherence, Levine levels, CDFD
│     ├─ layout/        % trophic layout algorithms
│     ├─ plotting/      % plotting functions and plotting utilities
│     └─ tflPlot.m      % high-level plotting entry point
│
├─ paper/
│  ├─ build/            % reproducible figure-build pipeline
│  ├─ figures/          % one callable function per paper figure
│  ├─ utils/            % paper-specific data loaders and helpers
│  ├─ data/             % source data used by the current figure pipeline
│  └─ outputs/          % generated figure outputs (typically not versioned)
│
└─ .gitignore

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
If you also want to reproduce the paper figures, add the paper code as well.

## Citing

If you use this code, please cite:

- Sansom, B. (2026) *Trophic Flow Layout for Directed Network Visualisation.* Working paper / preprint, 2026. Preprint link to be added.
- MacKay, R.S., Johnson, S., Sansom, B. (2020). How directed is a directed network? (Royal Society Open Science) [https://doi.org/10.1098/rsos.201138](https://doi.org/10.1098/rsos.201138)
- Homs-Dones, M., MacKay, R. S., Sansom, B., and Zhou, Y. (2025). *Circular Directional Flow Decomposition of Networks.* arXiv:2506.12546. [https://doi.org/10.48550/arXiv.2506.12546] Forthcoming in *Royal Society Open Science*.

A CITATION.cff file will be added once publication details are final.


## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

