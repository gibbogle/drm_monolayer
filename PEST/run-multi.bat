rjm_batch_submit.exe -c "module load intel/2017a" -c "chmod a+x *.sh" -c "bash script.sh" -m 100M -j mpi:4 -w 24:00:0 -f pest-dirs-multi.txt -ll info
rjm_batch_wait.exe -f pest-dirs-multi.txt -z 60 -ll info
