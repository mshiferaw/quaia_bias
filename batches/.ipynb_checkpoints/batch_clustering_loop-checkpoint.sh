#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=clustering_loop
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

#conda init 
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
conda activate quaia-env
G=20.5
 
# # for b in '0' '10' '20' '30' '40' '50' '60'
# for mask in '20' '30' '40' '50' '60' '70' '80'
# do 
#     # fname="_b${b}"
#     fname="_mask${mask}"
#     cat > G${G}${fname}_wp.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac
# #SBATCH --job-name=G${G}${fname}_wtheta_clustering
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=72:00:00

# echo "Starting batch job"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate quaia-env

# # python -u clustering.py wp 20.5 1 200 14 log --L max --mask ${mask} _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6
# # python -u clustering.py wp 20.5 0.15 240 20 log --L max --mask ${mask} _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# # python -u clustering.py wp 20.5 0.15 240 20 log --L max --mask ${mask} _MC_ Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# # python -u clustering.py wp 20.5 30 180 30 linear --L max --mask ${mask} _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# # python -u clustering.py wp 20.5 30 180 30 linear --L minmax --mask ${mask} _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# python -u clustering.py wtheta 20.5 0.1 25 10 log --L max --mask ${mask} _MC_ Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# EOF

#     sbatch G${G}${fname}_wp.sh
# done

# for b in '0' '10' '20' '30' '40' '50' '60'
# do 
#     fname="_b${b}"
#     cat > G${G}${fname}_wp.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=G${G}${fname}_wtheta_clustering
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo "Starting batch job"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate quaia-env

# # python -u clustering.py wp 20.5 1 200 14 log --L max _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6
# # python -u clustering.py wp 20.5 0.15 240 20 log --L max _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# # python -u clustering.py wp 20.5 0.15 240 20 log --L max _MC_ --b ${b} Quaia --z_bins 0.0 1.0 2.0 3.0 4.6 #1 200 14 log
# # python -u clustering.py wp 20.5 30 180 30 linear --L max _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# # python -u clustering.py wtheta 20.5 0.1 25 10 log --L max _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# python -u clustering.py wtheta 20.5 0.001 10 20 log --L max _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 
# EOF

#     sbatch G${G}${fname}_wp.sh
# done

# for pimax in '100' '200' '300' '400' '500' '600' '700' '800' '900'
# do 
#     fname="_pimax${pimax}"
#     cat > G${G}${fname}_wp_log_linear.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=G${G}${fname}_log_linear
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo "Starting batch job"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate quaia-env

# # # python -u clustering.py wp 20.5 30 180 30 linear --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --pimax ${pimax}
# python -u clustering.py wp 20.5 1 200 14 log --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --pimax ${pimax}
# python -u clustering.py wp 20.0 30 180 30 linear --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --pimax ${pimax}
# python -u clustering.py wp 20.0 1 200 14 log --L max _MC_ --b 30 Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --pimax ${pimax}
# EOF

#     sbatch G${G}${fname}_wp_log_linear.sh
# done

n_zbins=(4 4 4 3 3 2)
n_Lbins=(4 3 2 3 2 2)
G=20.5

for b in '0' '10' '20' '30' '40' '50' '60'
do 
    for j in "${!n_zbins[@]}" 
    do
        for L_method in 'split' 'maxsplit'
        do
            fname="_b${b}_zsplit${n_zbins[j]}_L${L_method}${n_Lbins[j]}"
            cat > G${G}${fname}_wtheta.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G${G}${fname}_wtheta
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# python -u clustering.py wtheta ${G} 0.001 10 30 log --L_method maxsplit _MC_ --b ${b} Planck --n_zbins ${zsplit} --n_Lbins ${Lsplit}
# python -u clustering.py wtheta ${G} 0.001 10 30 log --L_method split _MC_ --b ${b} Planck --n_zbins ${zsplit} --n_Lbins ${Lsplit}
python -u clustering.py wtheta ${G} 0.001 10 30 log --L_method ${L_method} _MC_ --b ${b} Planck --n_zbins ${n_zbins[j]} --n_Lbins ${n_Lbins[j]}
EOF

            sbatch G${G}${fname}_wtheta.sh
        done
    done
done

for b in '0' '10' '20' '30' '40' '50' '60'
do 
    for n_Lbins in '4' '3' '2'
    do
        for L_method in 'split' 'maxsplit'
        do
            fname="_b${b}_L${L_method}${n_Lbins}"
            cat > G${G}${fname}_wtheta.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=G${G}${fname}_wtheta
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate quaia-env

# python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method maxsplit _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --n_Lbins ${Lsplit}
# python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method split _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --n_Lbins ${Lsplit}
python -u clustering.py wtheta 20.5 0.001 10 30 log --L_method ${L_method} _MC_ --b ${b} Planck --z_bins 0.0 1.0 2.0 3.0 4.6 --n_Lbins ${n_Lbins}
EOF

            sbatch G${G}${fname}_wtheta.sh
        done
    done
done