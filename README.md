# dispuτauble

[![](https://img.shields.io/badge/arXiv-2504.16932%20-red.svg)](https://arxiv.org/abs/2504.16932)

This repository contains scripts used to perform the analysis cited above, including a modified version of `CAMB` that extrapolates to negative neutrino masses. If you make use of this code please cite [Sailer, Farren, Ferraro, White (2025)](https://ui.adsabs.harvard.edu/abs/2025arXiv250416932S/exportcitation). 

The modified version of `CAMB` can be installed with

```
python3 -m pip install -v git+https://github.com/NoahSailer/disputauble --user
```

To use it in `Cobaya` include `disputauble.CobayaCAMB_mnuEff` in your `theory` block. See `defaults/camb_mnu-eff_nuDegenerate.yaml` for an example.

-------

Chains are publicly available on [zenodo](https://www.zenodo.com). Download the chains and remake the figures with:

```
mkdir figures
wget ... (will upload chains to zenodo shortly)
wget https://github.com/NoahSailer/disputauble/do_analysis.py
python do_analysis.py
```

Should you wish to rerun these chains, do so with:

```
git clone https://github.com/NoahSailer/disputauble
cd disputauble
bash do_analysis.sh --submit=true
```
