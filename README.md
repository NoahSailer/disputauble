# dispuτauble

Reproduce the analysis with:
```
git clone https://github.com/NoahSailer/disputauble
cd disputauble
bash do_analysis.sh
```
**WARNING:** will submit loads of jobs for running chains and minimizers to the queue by default if they haven't already been submitted or downloaded from elsewhere. If you want to avoid this add the flag:
```
bash do_analysis.sh --submit false
```