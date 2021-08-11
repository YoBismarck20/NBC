#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=rosenMRIPrj
#SBATCH --time=48:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
start=`date +%s`
bash genome_download.sh species true
end=`date +%s`
runtime=$((end-start))
echo "Download time is $runtime"
