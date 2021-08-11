#!/bin/bash
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


. /etc/profile.d/modules.sh

### These four modules must ALWAYS be loaded

start=`date +%s`
#bash jellyfish_gen_new.bash /ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/lustre/scratch/zz374/nbc_paper_data/ncbi_reads_noerr/fold1/ 15 true
bash jellyfish_gen_new.bash /ifs/groups/rosenMRIGrp/wc492/timingStuff/Lastne/lustre/scratch/zz374/nbc_paper_data/ncbi_reads_noerr/fold1_computed/fold1/ 15 true

end=`date +%s`
runtime=$((end-start))
echo "jellyfish runtime is $runtime"
