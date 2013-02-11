mpirun.openmpi  -np $1  xterm -hold -e gdb -ex "run" -q ./zephyr 
