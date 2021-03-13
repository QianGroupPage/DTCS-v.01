#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=60
#SBATCH --nodes=30
#SBATCH --constraint=haswell

python3 runner.py | srun --no-kill --ntasks=30 --wait=0 payload.sh
