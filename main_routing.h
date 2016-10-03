/**
 * @file mavlink_control.h
 *
 * This application perform the routing of messages between Ground Station
 * Simulator and Autopilot Board. 
 *
 *
 *
 * @author Luigi Pannocchi, <l.pannocchi@gmail.com>
 */


// -----------------------------------------------------------------------
//   Includes
// -----------------------------------------------------------------------

#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <string.h>
#include <inttypes.h>
#include <fstream>
#include <signal.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>

using namespace std;

#include <common/mavlink.h>

#include "sim_interface.h"
#include "gs_interface.h"
#include "ue_interface.h"
#include "autopilot_interface.h"

extern "C" {
#include <ptask.h>
#include <tstat.h>
#include <pbarrier.h>
}


// File descriptor
FILE* file_TFtcSns;
FILE* file_TSndSns;

FILE* file_TFtcComm;
FILE* file_TSndComm;

FILE* file_TGS;

// -----------------------------------------------------------------------
//   Prototypes
// -----------------------------------------------------------------------
int main(int argc, char **argv);

void parse_commandline(int argc, char **argv, char *&gs_ip, 
						unsigned int &gs_r_port, unsigned int &gs_w_port);

void routing_messages(mavlink_message_t *msg, struct Interfaces* p);
void unpack_ctr_mess(mavlink_message_t *msg, struct Interfaces* p);

// Threads Bodies
//
void inflow_thread();
void simulator_thread();
void gs_thread();
void ue_thread();

// Threads Indexes
int inflowT_id;
int simT_id;
int gsT_id;
int ueT_id;

// Struct with the pointes to the interfaces
struct Interfaces 
{
    GS_Interface* gs;
    Autopilot_Interface* aut;
	UE_Interface* ue;
};  


// Global Variables
uint8_t UAV_base_mode; // Mode of the UAV

pthread_mutex_t mut_first_heartbeat;
pthread_cond_t cond_first_heartbeat;

float hil_ctr[4];


// Flags
bool autopilot_connected = false;
bool inflow_thread_active = false;
bool simulator_thread_active = false;
bool gs_thread_active = false;
bool ue_thread_active = false;


// Time variables
ptime hil_ctr_time;
ptime hil_ctr_old_time;

tspec sim_period;
tspec aut_period;
tspec gs_period;
tspec ue_period;

ptime inflow_drift_time;
ptime inflow_end_time;
ptime simulator_drift_time;
ptime gs_time;

pbarrier barrier;


// quit handler
GS_Interface* gs_interface_quit;
Autopilot_Interface *autopilot_interface_quit;
void quit_handler( int sig );



