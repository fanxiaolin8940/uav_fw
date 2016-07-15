// ---------------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------------
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>


#include "udp_port.h"

// ---------------------------------------------------------------------------
// Costructor / Destructor
// ---------------------------------------------------------------------------

Udp_Port::Udp_Port()
{
    memset(&locAddr, 0, sizeof(locAddr));
    memset(&remAddr, 0, sizeof(remAddr));
}

Udp_Port::Udp_Port(const char* ip_addr, uint16_t r_port, uint16_t w_port)
{
    // Endpoint for communication file descriptor
    sock = socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP);


    // Initialize local sockaddr_in structure
    memset(&locAddr, 0, sizeof(locAddr));

    locAddr.sin_family = AF_INET;
    locAddr.sin_port = htons(r_port);
    locAddr.sin_addr.s_addr = INADDR_ANY;

    // Initialize remote sockaddr_in structure
    memset(&remAddr, 0, sizeof(remAddr));

    remAddr.sin_family = AF_INET;
    remAddr.sin_port = htons(w_port);
    remAddr.sin_addr.s_addr = inet_addr(ip_addr);

    // Binding the socket to read
    if( -1 == bind(sock,(struct sockaddr *)&locAddr,sizeof(struct sockaddr)))
    {
        fprintf(stderr,"Udp_Port::Udp_Port error: bind failed!\n");
        close(sock);
        exit(EXIT_FAILURE);
    }

    
    // Perfoming a non blocking access
    if (fcntl(sock, F_SETFL, O_NONBLOCK | FASYNC) < 0)
    {
    	fprintf(stderr, "Udp_Port::Udp_Port error: unable to set nonblocking %s\n",
    	strerror(errno));
    	close(sock);
    	exit(EXIT_FAILURE);
    } 
}

Udp_Port::~Udp_Port()
{
    close(sock);
}

int Udp_Port::send_mav_mess(mavlink_message_t* message)
{
    int bytes = 0;
    unsigned int len = mavlink_msg_to_send_buffer(buf, message);
    bytes = sendto(sock, buf, len, 0,
            (struct sockaddr*)&remAddr, sizeof(struct sockaddr_in));

    return bytes;
}

int Udp_Port::send_bytes(char * data, unsigned int len)
{
    int bytes = 0;
    bytes = sendto(sock,data,len,0,(struct sockaddr*)&remAddr,sizeof(struct sockaddr_in));
    return bytes;
}

int Udp_Port::receive_bytes(char* data, unsigned int len)
{
    int i;
    char rec_buf[len];
    memset(rec_buf,0,len);

    ssize_t recsize;

    //recsize = recvfrom(sock,(void* )buf,BUFFER_LENGTH, 0, (struct sockaddr *)&remAddr, &fromlen)  
    recsize = recvfrom(sock,(void* )rec_buf, len, 0, 0, 0);

    /*
    while (nrec < len)
    {
        recsize = recvfrom(sock,(void* )buf,len,MSG_WAITALL,0,0);
        if( recsize == 0 )
        {
            return 0;
        }

        for(i = 0 ; i<recsize ; i++)
            received_data[i+nrec] = buf[i];

        nrec += recsize;
        
           if (recsize != len)
           {
           printf("Udp_Port::receive_bytes : Number of bytes mismatch!\n %d != %d \n",recsize,len);
           return -1;
           }
    }
    */
   
    if ( recsize > 0 )
    {
        for(i = 0 ; i < recsize ; i++)
        {
            data[i] = rec_buf[i];
            //printf(" %x ",data[i]);
        }
        //printf("\n");
    }
    
    if (recsize == -1)
        printf("Error reading UDP\n");

    return recsize;

}




