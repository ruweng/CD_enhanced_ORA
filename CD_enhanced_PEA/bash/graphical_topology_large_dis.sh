#PBS -l walltime=24:00:00
#PBS -l select=4:ncpus=32:mem=60gb
#PBS -N topology
#PBS -q med-bio

cd /rds/general/user/rjc4717/home/CD_enhanced_PEA/R

module load anaconda3/personal

Rscript graphical_topology_large_dis.R