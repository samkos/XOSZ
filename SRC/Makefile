# 5765-D93 (C) COPYRIGHT IBM CORP 1993, 1999                
# All Rights Reserved
# Licensed Materials - Property of IBM
#  
# US Government Users Restricted Rights - Use, duplication or disclosure
# restricted by GSA ADP Schedule Contract with IBM Corp.
#
# makefile for poetest.bcast, test of point to point bandwidth
   
TARGET = zephyr
OBJECTS =  zparam.o zuncol_mpi_timed.o    zpara.o zgarbage.o \
	   zdac.o   zcompact.o zinit.o zblas.o  zio.o  zmg.o \
	   zcalcul.o  znavier.o ztest.o zmainone.o  second.o \
	   zplot.o

COMPILER  = mpiifort
LINKER = mpiifort
FFLAGS    =  
#FFLAGS    =    -g -fdump-core -fbacktrace 
LIBS = -pthread -I/home/samy/OMPI/lib -L/home/samy/OMPI/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl 
all: $(TARGET)

$(TARGET) : $(OBJECTS)
	$(LINKER) $(FFLAGS) $(OBJECTS) -o zephyr

zplot_control : zplot_control.c
	mpicc -o zplot_control zplot_control.c

%.o : %.c
	mpiicc -c $<

zparam.o : zparam.f90 
	$(COMPILER) $(MPI_INC) $(FFLAGS) -DMPI -cpp -c $<

%.o : %.f90 zparam.o
	$(COMPILER) $(MPI_INC) $(FFLAGS) -DMPI -cpp -c $<

clean: 
	rm -f z*.o second.o zephyr *.mod  *.lst
# end
