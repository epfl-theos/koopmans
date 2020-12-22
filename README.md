# `Koopmans code`
CP-based code for Γ-only calculations with Koopmans-compliant functionals

## Directories
`cp_koopmans` contains the main CP-based (v4.1) Koopmans code \
`qe_koopmans` a fork of QE that contains also the wannier2odd interface (see below for installation details) \
`unfold_and_interpolate` submodule containing the Python code for band structure unfolding and interpolation (you might want to check also its README)

## Installation
You will need to install the git repositories `qe_koopmans` and `unfold_and_interpolate`. You can do that by running:

```
git submodule init
git submodule update
```

In order to have also the `Wannier90` code within `qe_koopmans` you need to enter `qe_koopmans` and then type:

```
git submodule update --init external/wannier90
```
