#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=sample_loop
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo "Starting batch job"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/batches
conda activate pyccl-env

zmin=(0.0 1.0 2.0 3.0)
zmax=(1.0 2.0 3.0 4.6)
Lmax=(-20.0 -25.0 -26.0 -27.0)

# for ndim in '2' '3'
# do
#     for i in "${!Lmax[@]}"
#     do 
#         for j in "${!zmin[@]}" 
#         do
#             for dlogz in '0.01' '0.001'
#             do 
#                 fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_dlogz${dlogz}"
#                 cat > sample${fname}.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=pool${fname}
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo -e "Starting batch job\n"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate pyccl-env

# python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz ${dlogz}
# EOF

#                 sbatch sample${fname}.sh
#             done
#         done
#     done
# done

# for ndim in '2' '3'
# do
#     for i in "${!Lmax[@]}"
#     do 
#         for j in "${!zmin[@]}" 
#         do
#             for sample in 'unif' 'rwalk' 'slice' 'rslice' #'auto'
#             do 
#                 fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_${sample}"
#                 cat > sample${fname}.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=pool${fname}
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo -e "Starting batch job\n"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate pyccl-env

# python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --sample ${sample}
# EOF

#                 sbatch sample${fname}.sh
#             done
#         done
#     done
# done

# for ndim in '2' '3'
# do
#     for i in "${!Lmax[@]}"
#     do 
#         for j in "${!zmin[@]}" 
#         do
#             for thetamin in '-2.8' '-2.4' '-2.0' '-1.6' '-1.2' #'0.001' '0.002' '0.003' '0.004' '0.005' '0.006' '0.007' '0.008' '0.009' #'0.01'#'auto'
#             # for thetamin in '-3.0' '-2.6' '-2.2' '-1.8' '-1.4' '-1.0' 
#             do 
#                 fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_thetamin${thetamin}"
#                 cat > sample${fname}.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=pool${fname}
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo -e "Starting batch job\n"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate pyccl-env

# python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --log_thetamin ${thetamin}
# EOF

#                 sbatch sample${fname}.sh
#             done
#         done
#     done
# done

# for ndim in '2' '3'
# do
#     for i in "${!Lmax[@]}"
#     do 
#         for j in "${!zmin[@]}" 
#         do
#             for dlogz in '0.01' '0.001'
#             do
#                 for nlive in '1000' '1500' '2000' 
#                 do 
#                     fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_dlogz${dlogz}_nlive${nlive}"
#                     cat > sample${fname}.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=pool${fname}
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo -e "Starting batch job\n"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate pyccl-env

# python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz ${dlogz} --nlive ${nlive}
# EOF

#                     sbatch sample${fname}.sh
#                 done
#             done
#         done
#     done
# done

# for ndim in '2' '3'
# do
#     for i in "${!Lmax[@]}"
#     do 
#         for j in "${!zmin[@]}" 
#         do
#             fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_fsat"
#             cat > sample${fname}.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=pool${fname}
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo -e "Starting batch job\n"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate pyccl-env

# # python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 2000 --log_thetamin -2.8
# # python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 500 --log_thetamin -2.0
# # python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.01 --nlive 1500 --log_thetamin -2.8
# # python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8
# # python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8
# # python -u sample.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --sigma 0.4
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 10
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 15
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 20
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 25
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --fsat 30
# EOF

#             sbatch sample${fname}.sh
#         done
#     done
# done

# for ndim in '2' '3'
# do
#     for i in "${!Lmax[@]}"
#     do 
#         for j in "${!zmin[@]}" 
#         do
#             fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_frac"
#             cat > sample${fname}.sh << EOF
# #!/bin/bash
# #SBATCH --partition=kipac,hns,normal
# #SBATCH --job-name=pool${fname}
# #SBATCH --output=logs/%x.out
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=48
# #SBATCH --mem=360GB
# #SBATCH --time=48:00:00

# echo -e "Starting batch job\n"
# source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
# cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
# conda activate pyccl-env

# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 1
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 1.4
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 1.8
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 2.1
# python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --frac 2.5
# EOF

#             sbatch sample${fname}.sh
#         done
#     done
# done

for ndim in '2' '3'
do
    for i in "${!Lmax[@]}"
    do 
        for j in "${!zmin[@]}" 
        do
            fname="_ndim${ndim}_zmin${zmin[j]}zmax${zmax[j]}_Lmax${Lmax[i]}_diag"
            cat > sample${fname}.sh << EOF
#!/bin/bash
#SBATCH --partition=kipac,hns,normal
#SBATCH --job-name=pool${fname}
#SBATCH --output=logs/%x.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360GB
#SBATCH --time=48:00:00

echo -e "Starting batch job\n"
source /home/users/mahlet/miniconda3/etc/profile.d/conda.sh
cd /oak/stanford/orgs/kipac/users/mahlet/quaia_bias/scripts
conda activate pyccl-env

python -u sampler.py ${ndim} ${zmin[j]} ${zmax[j]} ${Lmax[i]} --plot --dlogz 0.001 --nlive 1500 --log_thetamin -2.8 --diag
EOF

            sbatch sample${fname}.sh
        done
    done
done