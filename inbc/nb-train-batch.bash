#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=rosenMRIPrj
#SBATCH --time=12:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
NTHREADS=1
EXEC_ROOT=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/c++/NB.run
GEN_SRC=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/lustre/scratch/zz374/nbc_paper_data/genomes/fold1
OUTPUT=./NB_save
KMERSIZE=15

MEMLIM=$(($NTHREADS * 1500000)) # 3 GB/core - 1.5 GB/class

$EXEC_ROOT train $GEN_SRC/ -s $OUTPUT -t $NTHREADS -k $KMERSIZE -m $MEMLIM
