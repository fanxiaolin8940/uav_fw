/*
 *
 * file: sim_interface.cpp
 * @author : Luigi Pannocchi <l.pannocchi@gmail.com>
 *
 */

#include "sim_interface.h"

// --------------------------------------
// Simulation Interface Class
// --------------------------------------

Sim_Interface::Sim_Interface(): 
    udp_port((char*)"127.0.0.1", (uint32_t)49000, (uint32_t)49001 )
{
    // Set the ports
    setReadPort(49000);
    setWritePort(49001);

    // Define the structure for the polling
    fdsR[0].fd = udp_port.sock;
    fdsR[0].events = POLLIN;

    fdsW[0].fd = udp_port.sock;
    fdsW[0].events = POLLOUT;

    // Initialize Mutexes
    pthread_mutex_init(&mut_act_controls,0);
    pthread_mutex_init(&mut_sensors,0);

    pthread_cond_init(&sim_contr_cond,0);
    ActUpdated = 0;

    memset(&sensors,0,NUM_FLOAT_SENSORS * sizeof(float));
    memset(&copy_sensors,0,NUM_FLOAT_SENSORS*sizeof(float));

    memset(&act_controls,0,NUM_FLOAT_ACT_CONTROLS*sizeof(float));
    memset(&copy_act_controls,0,NUM_FLOAT_ACT_CONTROLS*sizeof(float));

    F_simulator_period_time_R = fopen("F_simulator_period_time_R.txt","w");
    F_simulator_period_time_S = fopen("F_simulator_period_time_S.txt","w");

    simulator_period_time_R = 0;
    simulator_period_time_S = 0;

}

Sim_Interface::Sim_Interface(char* &ip, unsigned int r_port, unsigned int w_port):
    udp_port(ip,r_port,w_port)
{
    setReadPort(r_port);
    setWritePort(w_port);

    // Define the structure for the polling
    fdsR[0].fd = udp_port.sock;
    fdsR[0].events = POLLIN;

    fdsW[0].fd = udp_port.sock;
    fdsW[0].events = POLLOUT;

    // Initialize Mutexes
    pthread_mutex_init(&mut_act_controls,0);
    pthread_mutex_init(&mut_sensors,0);

    ActUpdated = 0;

    memset(&sensors,0,NUM_FLOAT_SENSORS * sizeof(float)); 
    memset(&copy_sensors,0,NUM_FLOAT_SENSORS*sizeof(float));

    memset(&act_controls,0,NUM_FLOAT_ACT_CONTROLS*sizeof(float));
    memset(&copy_act_controls,0,NUM_FLOAT_ACT_CONTROLS*sizeof(float));

    F_simulator_period_time_R = fopen("F_simulator_period_time_R.txt","w");
    F_simulator_period_time_S = fopen("F_simulator_period_time_S.txt","w");

    simulator_period_time_R = 0;
    simulator_period_time_S = 0;

}

Sim_Interface::~Sim_Interface()
{

    fclose(F_simulator_period_time_S);
    fclose(F_simulator_period_time_R);

    printf("Destructor of Sim_Interface/n");
}

// ------------------------------------------------------
//  FUNCTIONS
// ------------------------------------------------------ 


//
// setReadPort

int Sim_Interface::setReadPort(unsigned int port)
{
    r_port = port;
}

//
// setWritePort
//
int Sim_Interface::setWritePort(unsigned int port)
{
    w_port = port;
}


//
// sendActuatorCommand
// 
int Sim_Interface::sendActuatorCommand(float* values, const unsigned int N)
{
    int i = 0;
    uint16_t timeout = 0; // ms
    int8_t bytes_written  = 0;

    float send_buf[NUM_FLOAT_ACT_CONTROLS] = {0,0,0,0};

    for (i=0 ; i<NUM_FLOAT_ACT_CONTROLS ; i++)
        send_buf[i] = values[i];

    // Check if the port is free
    int ret = poll(fdsW,1,timeout);
    if ( ret < 0 )
    {
        printf("Sim_Interface::sendActuatorOutput : ERROR IN WRITING UDP PORT\n");
        return -1;
    }

    // Write byte over UDP

    //printf("Sim_Interface::sendActuatorOutput : Writing on the UDP\n");

    bytes_written = udp_port.send_bytes((char*)send_buf, NUM_FLOAT_ACT_CONTROLS*sizeof(float));
    if ( bytes_written != NUM_FLOAT_ACT_CONTROLS*sizeof(float))
    {
        printf("Sim_Interface::sendActuatorOutput : ERROR  WHILE SENDING DATA\n");
        printf("Written %d bytes\n",bytes_written);
        return 0;
    }
    
    simulator_period_time_S = ptask_gettime(MICRO);
    fprintf(F_simulator_period_time_S,"%llu \n",simulator_period_time_S);

    return bytes_written;
}




//
// set_act_controls
//
/*
   int Sim_Interface::set_act_controls(float* values, unsigned int N)
   {
   if(N != NUM_FLOAT_ACT_CONTROLS)
   {
   printf("Sim_Interface::set_act_controls : NUMBER OF FLOAT TO SET DOES NOT CORRESPOND!\n %d != %d\n",N,NUM_FLOAT_ACT_CONTROLS);
   return -1;
   }

   int i = 0;
   pthread_mutex_lock(&mut_act_controls);

// Set the updated flag
ActUpdated = 1;
if ( SIM_INTERFACE_DBG ) printf("Act_controls Set to: [ "); 

for (i=0 ; i<NUM_FLOAT_ACT_CONTROLS ; i++) 
{
act_controls[i] = values[i];
if ( SIM_INTERFACE_DBG ) printf("%1.2f ",act_controls[i]);    
}

if ( SIM_INTERFACE_DBG ) printf("]\n");


pthread_mutex_unlock(&mut_act_controls); 
return 1;
}
*/


//
// fetchSensors
//
int Sim_Interface::fetchSensors(float* dest_sensors, const unsigned int N)
{
    uint16_t timeout = 0;  // ms
    int8_t read_bytes = 0;

    if (N != NUM_FLOAT_SENSORS)
    {
        printf("Sim_Interface::get_sensors : NUMBER OF FLOAT VALUES REQUESTED FOR THE SENSOR DOES NOT MATCH!\n %d != %d\n", N, NUM_FLOAT_SENSORS);
        return -1;
    }   

    int ret = poll(fdsR, 1, timeout);
    int i= 0;	

    // Something wrong happened
    if ( ret < 0 ) 
    { 
        printf("Sim_Interface::fetchSensors : ERROR IN READING UDP PORT\n");
        return -1;
    }

    // No data returned
    if (ret == 0)
    {
    //    printf("Sim_Interface::fetchSensors : NO DATA RETURNED IN THE TIMEOUT\n");
        return 0;
    }

    // Fetch data num bytes over UDP 
    if ( ret > 0 )
    {
        read_bytes = udp_port.receive_bytes((char *)&copy_sensors, NUM_FLOAT_SENSORS * sizeof(float));

        if ( read_bytes =! NUM_FLOAT_SENSORS * sizeof(float))
        {
            printf("Sim_Interface::fetchSensors : WRONG DATA RETURNED : read %d bytes\n", read_bytes);
            return -1;
        }

        // Retrieving complete, now copy the data into the variable
        for (i = 0 ; i < NUM_FLOAT_SENSORS ; i++)
            dest_sensors[i] = copy_sensors[i]; 

        // Debug printing
        
        //if ( SIM_INTERFACE_DBG && ( read_bytes == NUM_FLOAT_SENSORS*sizeof(float))) 
        {
            printf("READ BYTES:  %d\n",read_bytes);
            printf("Time     =   %1.2f \n",sensors[0]); 
            printf("Acc      = [ %1.2f %1.2f %1.2f ] \n", sensors[1], sensors[2], sensors[3]); 
            printf("Gyro     = [ %1.2f %1.2f %1.2f ] \n", sensors[4], sensors[5], sensors[6]); 
            printf("Magn     = [ %1.2f %1.2f %1.2f ] \n", sensors[7], sensors[8], sensors[9]); 
            printf("Abs Pres =   %1.2f \n",sensors[10]); 
            printf("Dif Pres =   %1.2f \n",sensors[11]); 
            printf("Pres Alt =   %1.2f \n",sensors[12]); 
            printf("Temp     =   %1.2f \n",sensors[13]); 
            printf("LLA      = [ %1.2f %1.2f %1.2f ] \n",sensors[14],sensors[15],sensors[16]); 

            printf("V_gps_mod=   %1.2f \n",sensors[17]); 
            printf("V_gps    = [ %1.2f %1.2f %1.2f ] \n", sensors[18], sensors[19], sensors[20]); 
            printf("COG      =   %1.2f \n",sensors[21]);
        } 
        

        simulator_period_time_R = ptask_gettime(MICRO);

        fprintf(F_simulator_period_time_R,"%llu \n",simulator_period_time_R);
    } 
    return ret;
}
