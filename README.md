# dispuτauble

```
git clone https://github.com/NoahSailer/disputauble
cd disputauble
mkdir chains
mkdir log
bash create_cobaya_env.sh
sbatch minimize.sh pact_tau\=0.06
sbatch run_chains.sh pact_tau\=0.06
sbatch add_data.sh add_bao pact_tau\=0.06
sbatch add_data.sh add_bao_sne pact_tau\=0.06
```
