CC = g++
CFLAGS = -I. -I mavlink/include/mavlink/v1.0 -I ptask/src
DBFLAG = -g
LIBS = -lpthread -lrt -lptask
MAIN_SOURCE = main_routing.cpp
OBJECTS = time_utils.o serial_port.o udp_port.o autopilot_interface.o gs_interface.o sim_interface.o# queue.o

all: main_routing.cpp $(OBJECTS) 
	$(CC) -o main_routing $(MAIN_SOURCE) $(CFLAGS) $(DBFLAG) -L ptask/src  $(OBJECTS) $(LIBS) 
	$(info )
	$(info ************    COMPILING     ***********)
	$(info *****************************************)
	rm -rf *o *~ .*.swn .*.swo .*.swp


time_utils.o: time_utils.c time_utils.h
	$(CC) -c $(DBFLAG) time_utils.c

#queue.o: queue.c queue.h
#	$(CC) -c $(DBFLAG) $(CFLAGS) queue.c

serial_port.o: serial_port.cpp serial_port.h time_utils.h
	$(CC) -c $(DBFLAG) $(LIBS) serial_port.cpp 

udp_port.o: udp_port.cpp udp_port.h
	$(CC) -c $(CFLAGS) $(DBFLAG) $(LIBS) udp_port.cpp

autopilot_interface.o: autopilot_interface.cpp autopilot_interface.h serial_port.h
	$(CC) -c $(CFLAGS) $(DBFLAG) $(LIBS) autopilot_interface.cpp

gs_interface.o: gs_interface.cpp gs_interface.h udp_port.h
	$(CC) -c $(CFLAGS) $(DBFLAG) $(LIBS) gs_interface.cpp

sim_interface.o: sim_interface.cpp sim_interface.h udp_port.h
	$(CC) -c $(CFLAGS) $(DBFLAG) $(LIBS) sim_interface.cpp


#main_routing: main_routing.cpp autopilot_interface.cpp serial_port.cpp udp_port.cpp time_utils.cpp
#	g++ -o main_routing -I mavlink/include/mavlink/v1.0 -I. serial_port.cpp ptime.c pmutex.c ptask.c \
#	autopilot_interface.cpp udp_port.cpp time_utils.cpp gs_interface.cpp sim_interface.cpp main_routing.cpp \
#	-lpthread -lrt -lptask \
	

clean:
	 rm -rf *o *~ mavlink_control .*.swn .*.swo .*.swp
