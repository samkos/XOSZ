
import sys
import math
import os
import getopt

import exceptions, traceback


LOCAL = False
TEST = False
nb_tasks = False
nb_nodes = False

class MyError(Exception):
    """
    class reporting own exception message
    """
    def __init__(self, value):
        self.value  =  value
    def __str__(self):
        return repr(self.value)

def except_print():

    exc_type, exc_value, exc_traceback = sys.exc_info()
    
    print "!!!!!!!!!!!!!!!!!!!!!!!!!"
    print exc_type
    print exc_value

    print "!!!!!!!!!!!!!!!!!!!!!!!!!"
    traceback.print_exception(exc_type, exc_value, exc_traceback,
                              file=sys.stderr)


def print_debug(s,r="nulll"):
  if debug:
    if r=="nulll":
      print s
    else:
      print s,":",r
      

#########################################################################
# welcome message
#########################################################################

def welcome_message():
    """ welcome message"""
    
    print """
                   #################################################
                   #                                               #
                   #      Welcome to run_job v 0.1 !               #
                   #                                               #
                   #################################################

                   
     """

#########################################################################
# usage ...
#########################################################################

def usage(message = None):
    """ helping message"""
    if message:
        print "  Error : %s \n" % message
    print "  usage: \n \t python run_job.py \
             \n\t\t--tasks=<nb_tasks>      # total number of tasks \
             \n\t\t[ --local]              # launch on the frontend\
             \n\t\t[ --test]  \
             \n\t\t[ --debug]  \
             \n\t\t[ --help]               # print this message \
           \n"  

    sys.exit(1)

#########################################################################
# command line parsing...
#########################################################################

def parse(args=sys.argv[1:]):
    """ parse the command line and set global _flags according to it """

    global LOCAL, TEST, DEBUG, nb_tasks, nb_nodes
    
    #print "parsing vishnu parameter : ",args
    try:
        opts, args = getopt.getopt(args, "h", 
                          ["help","debug", "local", "test", "tasks="])
    except getopt.GetoptError, err:
        # print help information and exit:
        usage(err)

    # first scan opf option to get prioritary one first
    # those who sets the state of the process
    # especially those only setting flags are expected here
    for option, argument in opts:
        if option in ("-h", "--help"):
            usage("")
        elif option in ("--local"):
            LOCAL = True
        elif option in ("--test"):
            TEST = True
        elif option in ("--debug"):
            DEBUG = True
        elif option in ("--tasks"):
            nb_tasks = int(argument)
            nb_nodes = math.ceil(nb_tasks/16.)

    if not(nb_tasks):
        usage()


scorep_filter = """
SCOREP_REGION_NAMES_BEGIN EXCLUDE
binvcrhs*
matmul_sub*
matvec_sub*
exact_solution*
binvrhs*
lhs*init*
timer_*
"""


job_scalasca =     job= """
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


module use /gpfslocal/pub/vihps/UNITE/local
module load intel-env/13.0.1 intelmpi/4.0.3
module load UNITE scalasca




#export OMP_NUM_THREADS=4
RUN=`pwd`

export DEST=../RES/${LOADL_TOTAL_TASKS} 

#\\rm -rf ${DEST}
mkdir -p ${DEST}
cd ${DEST}
cp ${RUN}/input .

#export SCOREP_EXPERIMENT_DIRECTORY=scorep_trace
#export SCOREP_METRIC_PAPI=PAPI_L2_TCM,PAPI_FP_OPS
#export SCOREP_ENABLE_TRACING=true
#export SCOREP_ENABLE_PROFILING=false
#export SCOREP_TOTAL_MEMORY=500m
#export SCOREP_FILTERING_FILE=./config/scorep.filt



scan mpirun -np ${LOADL_TOTAL_TASKS} ${RUN}/zephyr > output

#export SCOREP_EXPERIMENT_DIRECTORY=scorep_trace
scan -t  mpirun -np ${LOADL_TOTAL_TASKS} ${RUN}/zephyr > output_scan
""" 

job_tau =  """
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


module use /gpfslocal/pub/vihps/UNITE/local
module load UNITE
module load VI-HPS-TW
module load intel-env/13.0.1 intelmpi/4.0.3
module load tau
module load papi



module list

papi_avail

RUN=`pwd`

export DEST=./RES/${LOADL_TOTAL_TASKS} 

\\rm -rf ${DEST}
mkdir -p ${DEST}
cd ${DEST}
cp ${RUN}/input .


export TAU_MAKEFILE=$TAU/Makefile.tau-mpi-pdt
export TAU_MAKEFILE=$TAU/Makefile.tau-icpc-papi-mpi-pdt

wc $TAU_MAKEFILE 

#export TAU_TRACK_SIGNALS=1
export TAU_SAMPLING=1
export TAU_COMM_MATRIX=1
export TAU_METRICS=TIME:PAPI_FP_OPS:PAPI_L2_TCM

mpirun -np ${LOADL_TOTAL_TASKS} tau_exec ${RUN}/zephyr > output

""" 

def create_job(job_template):
    global nb_tasks, nb_nodes
    job = job_template % (nb_tasks, nb_nodes)
    return job



welcome_message()
parse()
job=create_job(job_tau)

print "creating and submitting job for %d tasks running on %d nodes" %(nb_tasks,nb_nodes)

f = open("job2submit","w")
f.write(job)
f.close()

if TEST:
    print job
elif LOCAL:
    os.system("sh job2submit")
else:
    os.system("llsubmit job2submit")



