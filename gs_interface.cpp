/*
 * file: gs_interface.cpp
 *
 */

#include "gs_interface.h"


// ------------------------------------------------------------------
// Ground Station Interface
// ------------------------------------------------------------------
/*
 *                     +----------------+
 *  UDP (14551)  <-----|                |
 *                     |  GS (x.x.x.x)  |
 *  UDP (14550)  ----->|                |
 *                     +----------------+
 */
GS_Interface::GS_Interface():
    udp_port((const char *)"127.0.0.1", (uint32_t)14551, (uint32_t)14550) 
{
    // Inizialize UDP
    setReadPort(14551);
    setWritePort(14550);

    // Define the structure for the polling
    fdsR[0].fd = udp_port.sock;
    fdsR[0].events = POLLIN;

    fdsW[0].fd = udp_port.sock;
    fdsW[0].events = POLLOUT;

    // Initialize Queues
    //init_mess_queue(&sendQueue);
    //init_mess_queue(&recQueue);

    // Initialize Mutexes
    pthread_mutex_init(&mut_sendQueue, 0);
    pthread_mutex_init(&mut_recQueue, 0);
    started = 1;
}

GS_Interface::GS_Interface(char *ip, uint32_t r_port, uint32_t w_port):
    udp_port(ip, r_port, w_port)
{
    // Initialize UDP
    setReadPort(r_port);
    setWritePort(w_port);

    // Define the structure for the polling
    fdsR[0].fd = udp_port.sock;
    fdsR[0].events = POLLIN;

    fdsW[0].fd = udp_port.sock;
    fdsW[0].events = POLLOUT;

    // Initialize Queues
    //init_mess_queue(&sendQueue);
    //init_mess_queue(&recQueue);

    // Initilize Mutexes
    pthread_mutex_init(&mut_sendQueue, 0);
    pthread_mutex_init(&mut_recQueue, 0);
    started = 1;
}

GS_Interface::~GS_Interface()
{
    printf("GS Destructor\n");
}

// ----------------------------------------------------------------
//  FUNCTIONS
// ----------------------------------------------------------------

//
// setReadPort
//
int GS_Interface::setReadPort(unsigned int port)
{
    r_port = port;
}

//
// setWritePort
//
int GS_Interface::setWritePort(unsigned int port)
{
    w_port = port;
}

//
// sendMessage
//
int GS_Interface::sendMessage()
{
    int bytes_sent;
    mavlink_message_t sendMessage;
    char buf[254];
    int len;
    
    // Take the message from the queue
    
    //int ret = extract(&sendQueue, &sendMessage);
    while (!sendQueue.empty())
    {
        pthread_mutex_lock(&mut_sendQueue);
        sendMessage = sendQueue.front();
        sendQueue.pop();
        pthread_mutex_unlock(&mut_sendQueue);

        //printf("sendQueue # = %d\n", sendQueue.size());
        bytes_sent = udp_port.send_mav_mess(&sendMessage);
    }
    return bytes_sent;
}

//
// receiveMessage
//
int GS_Interface::receiveMessage()
{
    uint16_t timeout = 0;  // ms
    int8_t read_bytes = 0;
    mavlink_message_t recMessage;

    int ret = poll(fdsR, 1, timeout);

    if (ret < 0) 
    { 
        printf("GS_Interface::receiveMessage : ERROR IN READING UDP PORT\n");
        return -1;
    }
    
    if (ret == 0)
    {
        //printf("GS_Interface::receiveMessage : NO DATA RETURNED\n");
        return 0;
    }
    else
    {
        // Receive data num bytes over UDP and put them at the sensors address
        read_bytes = udp_port.receive_bytes((char *)&recMessage, 
                                            sizeof(mavlink_message_t));
        if ( read_bytes =! sizeof(mavlink_message_t))
        {
            printf("GS_Interface::receiveMessage : WRONG DATA RETURNED\n");
            return -1;
        }
            // Retrieving complete, now copy the data into the variable
        pthread_mutex_lock(&mut_recQueue);
        recQueue.push(recMessage);
        //printf("recQueue # = %d\n", recQueue.size());
        //printf("Message id %d\n", recMessage.msgid);
        pthread_mutex_unlock(&mut_recQueue);
    }
    return 1;
}


//
// pushMessage
//
int GS_Interface::pushMessage(mavlink_message_t* msg)
{
    pthread_mutex_lock(&mut_sendQueue);
    sendQueue.push(*msg);
    pthread_mutex_unlock(&mut_sendQueue);
    return 1;
}


// 
// getMessage
//
int GS_Interface::getMessage(mavlink_message_t* msg)
{
    pthread_mutex_lock(&mut_recQueue);
    if (!recQueue.empty())
    {
        *msg = recQueue.front();
        recQueue.pop();
    }
    
    //printf("recQueue # = %d\n", recQueue.size());
    pthread_mutex_unlock(&mut_recQueue);

    return 1;
}
