NTHREADS=16
EXEC_ROOT=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/c++/NB.run

#GEN_SRC=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/lustre/scratch/zz374/nbc_paper_data/ncbi_reads_noerr/fold1/
GEN_SRC=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/c++/read2

MODEL=/ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/c++/inbc/NB_save/
#MODEL=1
KMERSIZE=15

MEMLIM=$((NTHREADS * 1400000)) # 48GB for 16 cores: 3 GB/core - 1.6 GB/class = 1400000

echo 'TESTING_READS', $GEN_SRC
echo 'MODEL_INDEX', $MODEL

gdb --args $EXEC_ROOT classify $GEN_SRC -s $MODEL -t $NTHREADS -k $KMERSIZE -m $MEMLIM


