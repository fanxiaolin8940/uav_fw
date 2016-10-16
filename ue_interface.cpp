/*
 * file: ue_interface.cpp
 *
 */

#include "ue_interface.h"

struct UE_SendData UEDataOut;
struct UE_RecData UEDataIn;


// ------------------------------------------------------------------
// Unreal Engine Interface
// ------------------------------------------------------------------
/*
 *                    +----------------+
 *  UDP (8001)  <-----|                |
 *                    |  UE (x.x.x.x)  |
 *  UDP (8000)  ----->|                |
 *                    +----------------+
 */
UE_Interface::UE_Interface():
	udp_port((const char *)"127.0.0.1", (uint32_t)8001, (uint32_t)8000) 
{
	int i;

	// Inizialize UDP
	setReadPort(8000);
	setWritePort(8001);

	// Define the structure for the polling
	fdsR[0].fd = udp_port.sock;
	fdsR[0].events = POLLIN;

	fdsW[0].fd = udp_port.sock;
	fdsW[0].events = POLLOUT;

	// Initialize Mutexes
	pthread_mutex_init(&mut_sendData, 0);
	pthread_mutex_init(&mut_recData, 0);
	started = 1;

	for (i = 0; i < 512; i++)
		rbuff[i] = 0;
}

UE_Interface::UE_Interface(char *ip, uint32_t r_port, uint32_t w_port):
	udp_port(ip, r_port, w_port)
{
	int i;

	// Initialize UDP
	setReadPort(r_port);
	setWritePort(w_port);

	// Define the structure for the polling
	fdsR[0].fd = udp_port.sock;
	fdsR[0].events = POLLIN;

	fdsW[0].fd = udp_port.sock;
	fdsW[0].events = POLLOUT;

	// Initilize Mutexes
	pthread_mutex_init(&mut_sendData, 0);
	pthread_mutex_init(&mut_recData, 0);
	started = 1;

	for (i = 0; i < 512; i++)
		rbuff[i] = 0;
}

UE_Interface::~UE_Interface()
{
	printf("UE Destructor\n");
}

// ----------------------------------------------------------------
//  FUNCTIONS
// ----------------------------------------------------------------

//
// setReadPort
//
int UE_Interface::setReadPort(unsigned int port)
{
	r_port = port;
}

//
// setWritePort
//
int UE_Interface::setWritePort(unsigned int port)
{
	w_port = port;
}

//
// sendData
// Send Data to the UE through the UDP port 
//
int UE_Interface::sendData()
{
	int bytes_sent;

	struct UE_SendData tempData;

	// Comping the data in a temp structure
	pthread_mutex_lock(&mut_sendData);
	tempData = UEDataOut;
	pthread_mutex_unlock(&mut_sendData);

	// Sending the data
	bytes_sent = udp_port.send_bytes((char* )&tempData, sizeof(struct UE_SendData));

	return bytes_sent;
}

// 
// setData
// Set the data to send 
//
int UE_Interface::setData(struct UE_SendData Data)
{
	pthread_mutex_lock(&mut_sendData);
	UEDataOut = Data;
	pthread_mutex_unlock(&mut_sendData);
		
	return 1;
}


//
// receiveMessage
// Wait for data on the UDP and receive one message from the GroundStation 
//
int UE_Interface::receiveData()
{
	int i;
	uint16_t timeout = 0;  // ms
	int8_t read_bytes = 0;

	int ret = poll(fdsR, 1, timeout);

	if (ret < 0) 
	{ 
		printf("UE_Interface::receiveMessage : ERROR IN READING UDP PORT\n");
		return -1;
	}

	if (ret == 0)
	{
		//printf("UE_Interface::receiveMessage : NO DATA RETURNED\n");
		return 0;
	}
	else
	{
		// Receive data num bytes over UDP and put them at the sensors address
		read_bytes = udp_port.receive_bytes(rbuff, sizeof(struct UE_RecData));
	}
	return 1;
}

// 
// getData
//
int UE_Interface::getData(struct UE_RecData* data)
{
	pthread_mutex_lock(&mut_recData);

	//printf("recQueue # = %d\n", recQueue.size());
	pthread_mutex_unlock(&mut_recData);

	return 1;
}

