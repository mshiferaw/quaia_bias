#!/bin/bash
#SBATCH --output=logs/%x.out
L_bins=(-32.0 -31.0 -30.0 -29.0 -28.0 -27.0 -26.0 -25.0 -24.0 -23.0 -22.0 -21.0 -20.0)
G_max="20.0"
save_tag=""

for ((bb=0; bb<$((${#L_bins[@]}-1)); bb++)); do
    fn_gcat_Lbin="../data/quaia_G${G_max}_Lmin${L_bins[bb]}Lmax${L_bins[bb+1]}${save_tag}.fits"
    fn_selfunc="../data/maps/selection_function_NSIDE64_G${G_max}_Lmin${L_bins[bb]}Lmax${L_bins[bb+1]}${save_tag}.fits"
    
    cat > sel_func_Lmin${L_bins[bb]}Lmax${L_bins[bb+1]}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac
#SBATCH --job-name=sel_func_G${G_max}_Lmin${L_bins[bb]}Lmax${L_bins[bb+1]}${save_tag}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
##SBATCH --mem=175GB
#SBATCH --mem=360GB
#SBATCH --time=72:00:00

echo "Starting batch job"
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env
python selection_function_map.py $fn_gcat_Lbin $fn_selfunc -p ../data/quaia_G${G_max}.fits
EOF
    
    sbatch sel_func_Lmin${L_bins[bb]}Lmax${L_bins[bb+1]}.sh
done
