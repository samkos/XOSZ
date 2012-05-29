# 5765-D93 (C) COPYRIGHT IBM CORP 1993, 1999                
# All Rights Reserved
# Licensed Materials - Property of IBM
#  
# US Government Users Restricted Rights - Use, duplication or disclosure
# restricted by GSA ADP Schedule Contract with IBM Corp.
#
# makefile for poetest.bcast, test of point to point bandwidth
   
TARGET = zephyr zplot_control
OBJECTS =  zuncol_fake.o zparam.o   zpara.o zgarbage.o \
	   zdac.o   zcompact.o zinit.o zblas.o  zio.o  zmg.o \
	   zcalcul.o  znavier.o ztest.o zmainone.o  second.o zplot.o

COMPILER  = gfortran
FFLAGS    =  -O3
#FFLAGS    =  -g

all: $(TARGET) 


zplot_control : zplot_control.c
	cc -o zplot_control zplot_control.c

zephyr : $(OBJECTS)
	$(COMPILER) -o $@ $(FFLAGS) $(OBJECTS) $(MPI_LIB) 

%.o : %.f90
	$(COMPILER) $(MPI_INC) $(FFLAGS) -c $<

clean: 
	rm -f z*.o second.o zephyr *.mod  *.lst
# end
