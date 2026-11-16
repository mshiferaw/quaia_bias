#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=pool_ndim3_zmin3.0zmax4.6_Lmax-27.0
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo -e "Starting batch job\n"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

# python -u sample.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 2000 --log_thetamin -2.8
# python -u sample.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 500 --log_thetamin -2.0
# python -u sample.py 3 3.0 4.6 -27.0 --plot --dlogz 0.01 --nlive 1500 --log_thetamin -2.8
# python -u sample.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8
# python -u sample.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8
# python -u sample.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --sigma 0.4
python -u sampler.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 10
python -u sampler.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 15
python -u sampler.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 20
python -u sampler.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 25
python -u sampler.py 3 3.0 4.6 -27.0 --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 30
