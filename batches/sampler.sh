#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=sampler
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
conda activate pyccl-env
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches

for corrfunc in 'xi' 'wp'
do 
    for penalty in '_soft' '_hard' 
    do 
        for ndim in '2' '3'
        do
            cat > sampler_${corrfunc}${penalty}_ndim${ndim}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sampler_${corrfunc}${penalty}_ndim${ndim}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

python -u sampler.py ${corrfunc} ${penalty} ${ndim}
EOF

                sbatch sampler_${corrfunc}${penalty}_ndim${ndim}.sh
            done
        done
    done
done