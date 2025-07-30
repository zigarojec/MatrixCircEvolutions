This is a tutorial on the usage of the MPI interface. 
Study it in the order in which the files are numbered. 

Before you run the examples, edit the host files (hosts.openmpi for OpenMPI. 

MPI
---
To run the first example with 4 task slots (using openmpi) type: 

  mpirun -n 4 python3 01-config.py
  
The hostlist can be specified in the following way
  
  mpirun -n 12 -hostfile hosts.openmpi python3 01-config.py

Under Windows you should use Microsoft MPI. mpiexec and python should be in 
the system path. 

  mpiexec /machinefile hosts.msmpi /np <number of processes> python <example>.py