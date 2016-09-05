/*
 *
 *  file: gs_interface.h
 *
 *  This code constitutes the interface for the communication with a Sim
 *  instance.
 *
 */

#include "udp_port.h"
#include <time.h>
#include "common/mavlink.h"
#include <poll.h>
#include <queue>
//#include "queue.h"


extern "C" {
#include "ptask.h"
}

class GS_Interface {

    public:

        // Constructor
        GS_Interface();
        GS_Interface(char* ip, uint32_t r_port, uint32_t w_port);

        // Destructor
        ~GS_Interface();

        // Settings
        int setReadPort(unsigned int port);
        int setWritePort(unsigned int port);

        // Communication 
        int sendMessage();
		int sendMessage(mavlink_message_t* message);
        int receiveMessage();

        int pushMessage(mavlink_message_t* message);
        int getMessage(mavlink_message_t* message);

        int started;

        Udp_Port udp_port;
        

    private:

        uint32_t r_port, w_port;

        //struct mess_queue sendQueue;
        //struct mess_queue recQueue;

        std::queue<mavlink_message_t> sendQueue;
        std::queue<mavlink_message_t> recQueue;
        
        pthread_mutex_t mut_sendQueue;
        pthread_mutex_t mut_recQueue;

        struct pollfd fdsR[1];
        struct pollfd fdsW[1];
		
		char rbuff[512];

};

