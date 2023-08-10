#PBS -S /bin/bash
#PBS -N merge_fastq
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -o ~/HRJ_monocytes/CHiC/hicup_combinations/OU
#PBS -e ~/HRJ_monocytes/CHiC/hicup_combinations/ER


BASE=~/HRJ_monocytes/CHiC/hicup_combinations

cd ${BASE}/S025TA-05_tech2
cat S025TA-05_tech2A_R1.fastq.gz S025TA-05_tech2B_R1.fastq.gz > S025TA-05_tech2_R1.fastq.gz
cat S025TA-05_tech2A_R2.fastq.gz S025TA-05_tech2B_R2.fastq.gz > S025TA-05_tech2_R2.fastq.gz

cd ${BASE}/S025V6-07
cat S025V6-07_A_R1.fastq.gz S025V6-07_B_R1.fastq.gz > S025V6-07_R1.fastq.gz
cat S025V6-07_A_R2.fastq.gz S025V6-07_B_R2.fastq.gz > S025V6-07_R2.fastq.gz

cd ${BASE}/S0262N-12
cat S0262N-12_A_R1.fastq.gz S0262N-12_B_R1.fastq.gz > S0262N-12_R1.fastq.gz
cat S0262N-12_A_R2.fastq.gz S0262N-12_B_R2.fastq.gz > S0262N-12_R2.fastq.gz

cd ${BASE}/S0263L-13
cat S0263L-13_A_R1.fastq.gz S0263L-13_B_R1.fastq.gz > S0263L-13_R1.fastq.gz
cat S0263L-13_A_R2.fastq.gz S0263L-13_B_R2.fastq.gz > S0263L-13_R2.fastq.gz

cd ${BASE}/S0264J-14
cat S0264J-14_A_R1.fastq.gz S0264J-14_B_R1.fastq.gz > S0264J-14_R1.fastq.gz
cat S0264J-14_A_R2.fastq.gz S0264J-14_B_R2.fastq.gz > S0264J-14_R2.fastq.gz

cd ${BASE}/S0267D-16
cat S0267D-16_A_R1.fastq.gz S0267D-16_B_R1.fastq.gz > S0267D-16_R1.fastq.gz
cat S0267D-16_A_R2.fastq.gz S0267D-16_B_R2.fastq.gz > S0267D-16_R2.fastq.gz

cd ${BASE}/S026C3-20
cat S026C3-20_A_R1.fastq.gz S026C3-20_B_R1.fastq.gz S026C3-20_C_R1.fastq.gz > S026C3-20_R1.fastq.gz
cat S026C3-20_A_R2.fastq.gz S026C3-20_B_R2.fastq.gz S026C3-20_C_R2.fastq.gz > S026C3-20_R2.fastq.gz

cd ${BASE}/S026T6-30_tech2
cat S026T6-30_tech2A_R1.fastq.gz S026T6-30_tech2B_R1.fastq.gz > S026T6-30_tech2_R1.fastq.gz
cat S026T6-30_tech2A_R2.fastq.gz S026T6-30_tech2B_R2.fastq.gz > S026T6-30_tech2_R2.fastq.gz

cd ${BASE}/S026D1-21
cat S026D1-21_A_R1.fastq.gz S026D1-21_B_R1.fastq.gz S026D1-21_C_R1.fastq.gz > S026D1-21_R1.fastq.gz
cat S026D1-21_A_R2.fastq.gz S026D1-21_B_R2.fastq.gz S026D1-21_C_R2.fastq.gz > S026D1-21_R2.fastq.gz
