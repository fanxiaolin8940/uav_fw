CPPFLAGS += -I. -I mavlink/include/mavlink/v1.0 -I ptask/src -I Gen_Code/DynModel_grt_rtw/
DBFLAG += -g
LIBS += -lpthread -lrt -lptask -lm
MAIN_SOURCE = main_routing.cpp
OBJECTS = time_utils.o serial_port.o udp_port.o autopilot_interface.o \
		gs_interface.o sim_interface.o DynModel.o DynModel_data.o

MATLAB_ROOT := /usr/local/MATLAB/R2016a
MATLABPATH := -I $(MATLAB_ROOT)/simulink/include -I $(MATLAB_ROOT)/extern/include

SUBDIR := Gen_Code/DynModel_grt_rtw

all: main_routing.cpp $(OBJECTS) 
	$(CXX) -o main_routing $(MAIN_SOURCE) $(CPPFLAGS)  $(DBFLAG) $(MATLABPATH)  \
	-L ptask/src  $(OBJECTS) $(LIBS) 
	$(info )
	$(info ************    COMPILING     ***********)
	$(info *****************************************)
	rm -rf *o *~ .*.swn .*.swo .*.swp

DynModel.o: $(SUBDIR)/DynModel.c
	$(MAKE) -C $(SUBDIR)
	
DynModel_data.o: $(SUBDIR)/DynModel_data.c
	$(MAKE) -C $(SUBDIR)

time_utils.o: time_utils.c time_utils.h
	$(CXX) -c $(DBFLAG) time_utils.c

serial_port.o: serial_port.cpp serial_port.h time_utils.h
	$(CXX) -c $(DBFLAG) $(LIBS) serial_port.cpp 

udp_port.o: udp_port.cpp udp_port.h
	$(CXX) -c $(CPPFLAGS) $(DBFLAG) $(LIBS) udp_port.cpp

autopilot_interface.o: autopilot_interface.cpp autopilot_interface.h serial_port.h
	$(CXX) -c $(CPPFLAGS) $(DBFLAG) $(LIBS) autopilot_interface.cpp

gs_interface.o: gs_interface.cpp gs_interface.h udp_port.h
	$(CXX) -c $(CPPFLAGS) $(DBFLAG) $(LIBS) gs_interface.cpp

sim_interface.o: sim_interface.cpp sim_interface.h udp_port.h
	$(CXX) -c $(CPPFLAGS) $(DBFLAG) $(LIBS) sim_interface.cpp


clean:
	 rm -rf *o *~ mavlink_control .*.swn .*.swo .*.swp

clean_txt:
	rm -rf *.txt
