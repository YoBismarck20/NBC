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

GEN_SRC=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/lustre/scratch/zz374/nbc_paper_data/ncbi_reads_noerr/fold1/
MODEL=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/c++/inbc/NB_save/
#MODEL=1
KMERSIZE=15

MEMLIM=$((NTHREADS * 1400000)) # 48GB for 16 cores: 3 GB/core - 1.6 GB/class = 1400000

echo 'TESTING_READS', $GEN_SRC
echo 'MODEL_INDEX', $MODEL

$EXEC_ROOT classify $GEN_SRC -s $MODEL -t $NTHREADS -k $KMERSIZE -m $MEMLIM
