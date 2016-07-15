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

struct Time_Stamps
{
    ptime heartbeat;
    ptime sys_status;
    ptime battery_status;
    ptime radio_status;
    ptime local_position_ned;
    ptime global_position_int;
    ptime position_target_local_ned;
    ptime position_target_global_int;
    ptime highres_imu;
    ptime attitude;
    ptime hil_controls;
    ptime hud;
    ptime encapsulated_data;
    ptime nav_controller_output;
    ptime attitude_target;
    ptime log_entry;
    ptime log_data;
    ptime param_value;
    ptime set_position_target_local_ned;
    ptime data_transmission_handshake;
    ptime gps_status;
    ptime gps_raw_int;
    ptime global_vision_position_estimate;
    ptime attitude_quaternion;
};

/*
void reset_timestamps(struct Time_Stamps* ts)
{
	ts->heartbeat = 0;
	ts->sys_status = 0;
	ts->battery_status = 0;
	ts->radio_status = 0;
	ts->local_position_ned = 0;
	ts->global_position_int = 0;
	ts->position_target_local_ned = 0;
	ts->position_target_global_int = 0;
	ts->highres_imu = 0;
	ts->attitude = 0;
	ts->hil_controls = 0;
}
*/


// Struct containing information on the MAV we are currently connected to
struct Mavlink_Messages {

	int sysid;
	int compid;

	// Heartbeat
	std::queue<mavlink_heartbeat_t> heartbeat;
	std::queue<mavlink_message_t> heartbeat_m;

	// System Status
	std::queue<mavlink_sys_status_t> sys_status;
	std::queue<mavlink_message_t> sys_status_m;

	// Hud
	std::queue<mavlink_vfr_hud_t> hud;
	std::queue<mavlink_message_t> hud_m;    
    
	// Battery Status
	std::queue<mavlink_battery_status_t> battery_status;
	std::queue<mavlink_message_t> battery_status_m;

    // Log Entry
	std::queue<mavlink_log_entry_t> log_entry;
	std::queue<mavlink_message_t> log_entry_m;
    
    // Log Data
	std::queue<mavlink_log_data_t> log_data;
	std::queue<mavlink_message_t> log_data_m;
    
    // Data Trasmission Handshake
	std::queue<mavlink_data_transmission_handshake_t> data_transmission_handshake;
	std::queue<mavlink_message_t> data_transmission_handshake_m;
    
    // Encapsulated_Data
	std::queue<mavlink_encapsulated_data_t> encapsulated_data;
	std::queue<mavlink_message_t> encapsulated_data_m;
    
    // Set Position target local ned
	std::queue<mavlink_set_position_target_local_ned_t> set_position_target_local_ned;
	std::queue<mavlink_message_t> set_position_target_local_ned_m;
    
    // Navigation Controller Output
	std::queue<mavlink_nav_controller_output_t> nav_controller_output;
	std::queue<mavlink_message_t> nav_controller_output_m;
    
	// Radio Status
	std::queue<mavlink_radio_status_t> radio_status;
	std::queue<mavlink_message_t> radio_status_m;

    // Parameters Values
    std::queue<mavlink_param_value_t> param_value;
	std::queue<mavlink_message_t> param_value_m;
    
	// Local Position
	std::queue<mavlink_local_position_ned_t> local_position_ned;
	std::queue<mavlink_message_t> local_position_ned_m;

	// Global Position
	std::queue<mavlink_global_position_int_t> global_position_int;
	std::queue<mavlink_message_t> global_position_int_m;
    
    // GPS Raw Int
	std::queue<mavlink_gps_raw_int_t> gps_raw_int;
	std::queue<mavlink_message_t> gps_raw_int_m;
    
    // GPS Raw Int
	std::queue<mavlink_gps_status_t> gps_status;
	std::queue<mavlink_message_t> gps_status_m;
    
    // Attitude Quaternion
    std::queue<mavlink_attitude_quaternion_t> attitude_quaternion;
	std::queue<mavlink_message_t> attitude_quaternion_m;
    
	// Local Position Target
	std::queue<mavlink_position_target_local_ned_t> position_target_local_ned;
	std::queue<mavlink_message_t> position_target_local_ned_m;
    
    // Attitude Target
	std::queue<mavlink_attitude_target_t> attitude_target;
	std::queue<mavlink_message_t> attitude_target_m;

	// Global Position Target
	std::queue<mavlink_position_target_global_int_t> position_target_global_int;
	std::queue<mavlink_message_t> position_target_global_int_m;

    // Global Vision Position Estimate
	std::queue<mavlink_global_vision_position_estimate_t> global_vision_position_estimate;
	std::queue<mavlink_message_t> global_vision_position_estimate_m;
    
	// HiRes IMU
	std::queue<mavlink_highres_imu_t> highres_imu;
	std::queue<mavlink_message_t> highres_imu_m;

	// Attitude
	std::queue<mavlink_attitude_t> attitude;
	std::queue<mavlink_message_t> attitude_m;

	// Hil Controls
	std::queue<mavlink_hil_controls_t> hil_controls;
	std::queue<mavlink_message_t> hil_controls_m;

	// System Parameters?


	// Time Stamps
	Time_Stamps time_stamps;
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
		int fetch_message();
		void get_message(mavlink_message_t* req_mess);
		int retrieve_message_id();
		int retrieve_nread_messages();

		// Write mavlink message on the serial interface
		int send_message(mavlink_message_t* msg);
		int send_sensorData(mavlink_message_t* msg);

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


