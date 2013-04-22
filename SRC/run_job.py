
import sys
import math
import os

if len(sys.argv) < 2:
    sys.exit('Usage: run_job.py <nb_procs> ' )





nb_tasks = int(sys.argv[1])
node = math.ceil(nb_tasks/16.)



job= """
#!/bin/bash
## submit from ./bin directory with "llsubmit run.ll"

#@ class            = clallmds
#@ job_name         = ZEPHYR
#@ total_tasks      = %d
#@ node             = %d
#@ node_usage       = not_shared
#@ wall_clock_limit = 0:10:00
#@ output           = $(job_name).$(jobid)
#@ error            = $(job_name).$(jobid)
#@ environment      = COPY_ALL
#@ job_type         = mpich
#@ queue

module load intel-env/13.0.1 intelmpi/4.0.3

export OMP_NUM_THREADS=4
RUN=`pwd`


\rm -rf ../RES/${LOADL_TOTAL_TASKS} 
mkdir -p ../RES/${LOADL_TOTAL_TASKS} 
cd ../RES/${LOADL_TOTAL_TASKS} 
cp ${RUN}/input .

mpirun -np ${LOADL_TOTAL_TASKS} ${RUN}/zephyr > output
""" % (nb_tasks, node)

print "creating and submitting job for %d tasks running on %d nodes" %(nb_tasks,node)

f = open("job2submit","w")
f.write(job)
f.close()

os.system("llsubmit job2submit")



