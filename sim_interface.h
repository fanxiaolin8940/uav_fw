/*
 *
 *  file: sim_interface.h
 *
 *  This code constitutes the interface for the communication with a Sim
 *  instance.
 *
 *
 */

#include "udp_port.h"
#include <time.h>
#include <poll.h>

extern "C" {
#include "ptask.h"
}

#define SIM_INTERFACE_DBG 1
#define NUM_FLOAT_ACT_CONTROLS 4
#define NUM_FLOAT_SENSORS 22
 


class Sim_Interface {

    public:

        // Constructor
        Sim_Interface();
        Sim_Interface(char*& ip, unsigned int r_port, unsigned int w_port);

        // Destructor
        ~Sim_Interface();

        // Settings
        int setReadPort(unsigned int port);
        int setWritePort(unsigned int port);

        /*
         * Send actuator values to the simulator over UDP
         */
        //int sendActuatorControls();

        /*
         * Set Actuator Control variables
         */
        //int set_act_controls(float* values, const unsigned int N);

        /*
         * Send Actuator Command 
         */
        int sendActuatorCommand(float* values, const unsigned int N);

        /*
         * Fetch sensor data from Simulator
         */
        int fetchSensors(float* dest_sens, const unsigned int N);

        /*
         * Get Sensor variables
         */
        //int get_sensors(float* dest_sens, const unsigned int N);


        ptime simulator_period_time_S;
        ptime simulator_period_time_R;

    private:

        FILE* F_simulator_period_time_S;
        FILE* F_simulator_period_time_R;

        /*
         * Communication Structures
         */

        // Ports from/to the Simulator 
        unsigned int r_port, w_port; 
        // UDP port class
        Udp_Port udp_port;
        //Structure for the reading/writing poll()
        struct pollfd fdsR[1];
        struct pollfd fdsW[1];


        // Time structs
        tspec act_timestamp, sens_timestamp;
        
        // Condition Variable 
        int ActUpdated;
        pthread_cond_t sim_contr_cond;

        // Data structures mutex and condition variables
        pthread_mutex_t mut_act_controls;
        pthread_mutex_t mut_sensors;

        float copy_act_controls[NUM_FLOAT_ACT_CONTROLS];
        float act_controls[NUM_FLOAT_ACT_CONTROLS];

        float sensors[NUM_FLOAT_SENSORS];
        float copy_sensors[NUM_FLOAT_SENSORS];

};

