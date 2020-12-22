# `Koopmans code`
CP-based code for Γ-only calculations with Koopmans-compliant functionals

## Directories
`cp_koopmans` contains the main CP-based (v4.1) Koopmans code \
`qe_koopmans` a fork of QE that contains also the wannier2odd interface (see below for installation details) \
`unfold_and_interpolate` submodule containing the Python code for band structure unfolding and interpolation

## Installation
You will need to install the git repositories `qe_koopmans` and `unfold_and_interpolate`. You can do that by running:

```
git submodule update --init --recursive
```

the `--recursive` allows to install also the submodules within `qe_koopmans` like, in particular, `wannier90`.
