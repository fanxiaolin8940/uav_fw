/*
 *
 *  file: ue_interface.h
 *
 *  This code constitutes the interface for the communication with a Sim
 *  instance.
 *
 */

#include "udp_port.h"
#include <time.h>
#include <poll.h>

extern "C" {
#include "ptask.h"
}


struct UE_SendData
{
	// Drone Index
	int Id;

	// Position
	float X;
	float Y;
	float Z;

	// Orientation
	float r;
	float p;
	float y;
};

struct UE_RecData
{
	// Drone Index
	int Id;

	// Position
	float X;
	float Y;
	float Z;

	// Orientation
	float r;
	float p;
	float y;
};




class UE_Interface {

    public:

        // Constructor
        UE_Interface();
        UE_Interface(char* ip, uint32_t r_port, uint32_t w_port);

        // Destructor
        ~UE_Interface();

        // Settings
        int setReadPort(unsigned int port);
        int setWritePort(unsigned int port);

        // Communication 
        int sendData();
		int setData(struct UE_SendData);

		int receiveData();
        int getData(struct UE_RecData* data);

        int started;

        Udp_Port udp_port;
        

    private:

        uint32_t r_port, w_port;

        pthread_mutex_t mut_sendData;
        pthread_mutex_t mut_recData;

        struct pollfd fdsR[1];
        struct pollfd fdsW[1];
		
		char rbuff[512];
};

