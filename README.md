# dispuτauble

[![](https://img.shields.io/badge/arXiv-2504.16932%20-red.svg)](https://arxiv.org/abs/2504.16932)

This repository contains scripts used to perform the analysis cited above, including a modified version of `CAMB` that extrapolates to negative neutrino mass. If you make use of this code please cite [Sailer, Farren, Ferraro, White (2025)](https://inspirehep.net/literature/2915153). 

## Extrapolation to negative neutrino mass

The modified version of `CAMB` can be installed with

```
python3 -m pip install -v git+https://github.com/NoahSailer/disputauble --user
```

To use it in `Cobaya` include `disputauble.CobayaCAMB_mnuEff` in your `theory` block. See `defaults/camb_mnu-eff_nuDegenerate.yaml` for an example.

## Public chains and scripts

Chains are publicly available on [zenodo](https://zenodo.org/records/15298950). Download the chains and reproduce the figures with:

```
wget https://raw.github.com/NoahSailer/disputauble/main/do_analysis.py
wget https://zenodo.org/records/15298950/files/chains.zip
mkdir figures
unzip chains.zip
python do_analysis.py
```

Should you wish to rerun the chains, do so with:

```
git clone https://github.com/NoahSailer/disputauble
cd disputauble
bash do_analysis.sh --submit=true --njobs=3
```