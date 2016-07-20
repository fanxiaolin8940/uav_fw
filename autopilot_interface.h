/**
 * @file autopilot_interface.h
 *
 * @brief Autopilot interface definition
 *
 * Functions for sending and recieving commands to an autopilot via MAVlink
 *
 * @author Luigi Pannocchi, <l.pannocchi@gmail.com>
 *
 */


#ifndef AUTOPILOT_INTERFACE_H_
#define AUTOPILOT_INTERFACE_H_


// -----------------------------------------------------------------------
//   Includes
// -----------------------------------------------------------------------

#include "serial_port.h"

#include <signal.h>
#include <sys/time.h>
#include <stdint.h>
#include <common/mavlink.h>
#include <queue>
extern "C" {
#include <ptask.h>
}
// ------------------------------------------------------------------------
//   Defines
// ------------------------------------------------------------------------


// ------------------------------------------------------------------------
//   Data Structures
// ------------------------------------------------------------------------
// Struct containing information on the MAV we are currently connected to
struct Mavlink_Messages {

	int sysid;
	int compid;

	std::queue<mavlink_message_t> messages;

	// Time Stamps
	long unsigned int time_stamps[256];
};


// ---------------------------------------------------------------------
//   Autopilot Interface Class
// ---------------------------------------------------------------------
/*
 * Autopilot Interface Class
 *
 */
class Autopilot_Interface
{

	public:

		Autopilot_Interface();
		Autopilot_Interface(char *&uart_name_, int &baudrate);
		~Autopilot_Interface();

		ptime write_count;
		ptime tot_byte;

		// Counters
		uint64_t highres_imu_count;
		uint64_t heartbeat_count;
		uint64_t hil_controls_count;

		// Variable to compute the frequency of reception of mavlink_messages
		ptime read_highres_imu_old;
		ptime read_heartbeat_old;
		ptime read_hil_controls_old;

		// State of the vehicle
		int system_id;
		int autopilot_id;
		int companion_id;
		uint8_t mav_type;
		uint8_t system_status;
		uint32_t custom_mode;
		uint8_t base_mode;

		// Structure to stock the messages
		Mavlink_Messages current_messages;
		std::queue<int> queueIndexFetched;

		// COMMUNICATION
		// Function to read messages from the Serial Port, manage them and
		// retrieve them.
		int fetch_data();
		int get_message(mavlink_message_t* req_mess);

		// Write mavlink message on the serial interface
		int send_message(mavlink_message_t* msg);

		//        int write_messages();
		void write_hilsensors();

		// Interface Handles
		void start_hil();
		void stop_hil();
		bool is_hil();

		void enable_offboard_control();
		void disable_offboard_control();

		Serial_Port uart_port; 

	private:
		FILE* f_aut_THilCtr;
		FILE* f_aut_TSens;
		FILE* f_aut_TSens_before;
		int control_status;

		int handle_message(mavlink_message_t* message);

		int toggle_offboard_control( bool flag );
		void write_ping();
		void write_raw();

		// Synchronization Part


		pthread_mutex_t mut_Messages;
		pthread_mutex_t mut_queueIndex;

		// Mutex and condition variable for the sending queues
		pthread_mutex_t mut_sendMessage;
		pthread_cond_t cond_empty;

		//mutex for the access to the base_mode variable 
		bool new_heartbeat;
		pthread_mutex_t mut_heartbeat;
		pthread_cond_t cond_heartbeat;

		bool hil_mode;

		std::queue<mavlink_message_t> LPsendQueue;
		std::queue<mavlink_message_t> HPsendQueue;

		struct pollfd fdsR[1];
		struct pollfd fdsW[1];

};



#endif // AUTOPILOT_INTERFACE_H_


