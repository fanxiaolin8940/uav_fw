/**
 * upd_port.h
 * 
 * @brief UDP interface definition
 * 
 * Functions for opening, closing , reading and writing via serial ports
 *
 */

#ifndef UDP_PORT_H_
#define UDP_PORT_H_

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <errno.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>

#include <common/mavlink.h>


#define BUFFER_LENGTH 2041



// ---------------------------------------------------------------------------
// UDP Port Manager Class
// ---------------------------------------------------------------------------
/*
 * UDP Port Class
 *
 *
 */
class Udp_Port
{

public:

  Udp_Port();
  Udp_Port(const char* target_ip, uint16_t r_port, uint16_t w_port);
  ~Udp_Port();
  
  int send_mav_mess(mavlink_message_t* message);
  int send_bytes(char* data, unsigned int len);
  int receive_bytes(char* data, unsigned int len);

  int sock;
private:
  
  char target_ip[100];

  struct sockaddr_in remAddr;
  struct sockaddr_in locAddr;

  unsigned char buf[BUFFER_LENGTH];

  socklen_t fromlen;
  int bytes_sent;
  mavlink_message_t msg;
  
};


#endif // UDP_PORT_H_
