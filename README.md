# `Koopmans code`
CP-based code for Γ-only calculations with Koopmans-compliant functionals

## Directories

## Installation
You will need to install the git repositories `qe_koopmans` and `unfold_and_interpolate`. You can do that by running:

```
git submodule init
git submodule update --remote --merge
```

In order to have also the `Wannier90` code within `qe_koopmans` you need to enter `qe_koopmans` and then type:

```
git submodule update --init external/wannier90
```
