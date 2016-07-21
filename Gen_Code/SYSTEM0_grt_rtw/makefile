program_C_SRCS := $(wildcard *.c)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_OBJS := $(program_C_OBJS)

MATLAB_ROOT := /usr/local/MATLAB/R2016a

program_INCLUDE_DIRS := ./ $(MATLAB_ROOT)/simulink/include \
	$(MATLAB_ROOT)/simulink/include \
	$(MATLAB_ROOT)/extern/include \
	$(MATLAB_ROOT)/rtw/c/src \
	$(MATLAB_ROOT)/rtw/c/src/ext_mode/common

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))

CFLAGS += -g

all: $(program_C_SRCS)
	$(CC) $(CPPFLAGS) $(program_C_SRCS)

SYSTEM0.o: SYSTEM0.c
	$(CC) -I $(MATLAB_ROOT)/simulink/include -I $(MATLAB_ROOT)/extern/include -I ./ -c SYSTEM0.c -lm
	
Dynamics_data.o: Dynamics_data.c
	$(CC) -I $(MATLAB_ROOT)/simulink/include/ -I $(MATLAB_ROOT)/extern/include/ -I ./ -c Dynamics_data.c

 
