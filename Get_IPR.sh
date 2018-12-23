#!/bin/bash
for dis in 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.8 0.9 1.0 1.2 1.4 2.0 3.0 4.0 5.0
do printf "${dis} "
std_dev=0.0
mean=0.0

for conf in {1..8}
do val=$(awk 'NR == 30 {print $7}' Runs/J_1.15/Dis_${dis}/unbiased_OP_seed_101/Disorder_seed_${conf}/out_run.txt)
printf "${val} "
mean=$(echo "$mean + $val" | bc)
done
mean=$(echo "${mean}*(0.125)" | bc)
printf "${mean} "


for conf in {1..8}
do val=$(awk 'NR == 30 {print $7}' Runs/J_1.15/Dis_${dis}/unbiased_OP_seed_101/Disorder_seed_${conf}/out_run.txt)
std_dev=$(echo "${std_dev} + ((${val} - ${mean})*(${val} - ${mean}))" | bc)
done
std_dev=$(echo "${std_dev}*(0.125)" | bc)
std_dev=$(echo "sqrt(${std_dev})" | bc)


printf "${std_dev} "
echo ""
done
