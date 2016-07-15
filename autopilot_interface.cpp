/**
 * @file autopilot_interface.cpp
 *
 * @brief Autopilot interface functions
 *
 * Functions for sending and recieving commands to an autopilot 
 * via MAVlink
 *
 * @author Luigi Pannocchi, <l.pannocchi@gmail.com>
 */


// ---------------------------------------------------------------------
//   Includes
// ---------------------------------------------------------------------

#include "autopilot_interface.h"

#define AUT_INTERFACE_DBG 1 



// ---------------------------------------------------------------------
//   Autopilot Interface Class
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
//   Con/De structors
// ---------------------------------------------------------------------
Autopilot_Interface::Autopilot_Interface(char *&uart_name_, int &baudrate):
   uart_port(uart_name_,baudrate) 
{

    // Initialize attributes (State of the Vehicle)
    write_count = 0;
    tot_byte = 0;

    control_status = 0;
    
    system_id    = 1; // system id
    autopilot_id = 1; // autopilot component id
    companion_id = 255; // companion computer component id
    system_status = 0;
    mav_type = 2;

    base_mode = 0;
    custom_mode = 0;

    current_messages.sysid  = system_id;
    current_messages.compid = autopilot_id;

    hil_mode = false;

    // Initialization of mutexes
    
    pthread_mutex_init(&mut_sendMessage,0);
    pthread_cond_init(&cond_empty,0);

    pthread_mutex_init(&mut_Messages,0);
    pthread_mutex_init(&mut_queueIndex,0);
    
    // Conditional Mutex
    pthread_mutex_init(&mut_heartbeat,0);
    pthread_cond_init(&cond_heartbeat,0);

    f_aut_THilCtr = fopen("./T_aut_HilCtr.txt","w");
    f_aut_TSens = fopen("./T_aut_Sens.txt","w");

    // Initialize the serial port giving device name and 
    // baudrate.
    // Initialize the structure for the poll() 
    fdsR[0].fd = uart_port.fd;
    fdsR[0].events = POLLIN;

    fdsW[0].fd = uart_port.fd;
    fdsW[0].fd = POLLOUT;

    read_heartbeat_old = 0;

}

Autopilot_Interface::~Autopilot_Interface() 
{
    fclose(f_aut_THilCtr);
    fclose(f_aut_TSens);

    printf("DESTRUCTOR\n");
}

// -------------------------------------------------------- 
//  FUNCTIONS
// --------------------------------------------------------

//
// read_message 
//
int Autopilot_Interface::fetch_message()
{
    //printf("Autopilot_Interface::fetch_message() \n");
    int i = 0;

    // Allocate Space for the Reading Buffer
    uint8_t NBytes = 128;
    char rbuff[NBytes];
    
    // Allocate Mavlink message variables
    mavlink_status_t status;
    mavlink_message_t recMessage;

    // Flag for the Message reception
    int msgReceived = 0; 

    // Number of messages parsed in read operation
    int NMessages = 0;

    // Id of the received message
    int message_Id = -1;

     
    // Until I don't receive a message...
    while (!msgReceived)
    {
        // Lock the device and try to read at most NBytes 
        int nread = uart_port.readBytes(rbuff, NBytes); 

        //printf("Autopilot_Interface::fetch_message() [Inside the while_loop()] 
		//		line %d \n read %d bytes from serial\n", __LINE__, nread);
    
        // Check for errors
        if (nread < 0)
        {
            printf("Autopilot_Interface::read_message : ERROR IN READING UART\n");
            return -1;
        }

        // No data on the device
        if (nread == 0)
        {
            return 0;
        }

        // Seek for mavlink messages in the received data stream
        for (i = 0; i < nread; i++)
        {
            // Parse 1 byte at time
            if (mavlink_parse_char(MAVLINK_COMM_1, rbuff[i], &recMessage, &status))
            {
                pthread_mutex_lock(&mut_Messages);

                // Handle the message and save in the Stock Structure 
                message_Id = handle_message(&recMessage);
                // Take trace of the received messages
                queueIndexFetched.push(message_Id);

                pthread_mutex_unlock(&mut_Messages);
                
                NMessages++;
                // Set the flag to 1 to signal that a full message has been retrieved 
                msgReceived = 1;
            }
        }
    }

    /* Return the number of read messages */
    return NMessages;
}


//
//  handle_message 
//
// In the Autopilot interface the handling consists in stocking the data into 
// queues for a future retrieval. 
//
int Autopilot_Interface::handle_message(mavlink_message_t* message)
{
    Time_Stamps this_timestamps;
    // Handle Message ID
    switch (message->msgid)
    {
        case MAVLINK_MSG_ID_HEARTBEAT:
            {
                // Copy the message to the structure both in the form of mavlink message
                // and of data structure
                current_messages.heartbeat_m.push(*message); 

                //printf("MAVLINK_MSG_ID_HEARTBEAT\n");
                mavlink_heartbeat_t heartbeat;
                mavlink_msg_heartbeat_decode(message, &heartbeat);
                current_messages.heartbeat.push(heartbeat);

                // Look the mutex for accessing the vehicle state information.
                // This information could be shared by other threads.
                // Wakeup eventually blocked thread waiting for new state information
                pthread_mutex_lock(&mut_heartbeat);
                base_mode = heartbeat.base_mode;
                custom_mode = heartbeat.custom_mode;
                mav_type = heartbeat.type;
                system_status = heartbeat.system_status;
                custom_mode = heartbeat.custom_mode;

                //printf("base_mode = %u\n",base_mode);
                pthread_cond_signal(&cond_heartbeat);
                pthread_mutex_unlock(&mut_heartbeat);
                
                current_messages.time_stamps.heartbeat = ptask_gettime(MICRO);

                this_timestamps.heartbeat = current_messages.time_stamps.heartbeat;
                /*
                heartbeat_count++;
                if ( AUT_INTERFACE_DBG )
                {
                    if ( (current_messages.time_stamps.heartbeat - read_heartbeat_old) > 10000000 )
                    {
                        printf("HEARTBEAT frequency :   %lu \n",heartbeat_count/10);
                        heartbeat_count = 0;
                        read_heartbeat_old = current_messages.time_stamps.heartbeat;
                    }
                }
                */
                break;
            }

        case MAVLINK_MSG_ID_SYS_STATUS:
            {
                current_messages.sys_status_m.push(*message);
                //printf("MAVLINK_MSG_ID_SYS_STATUS\n");
                
                /* Decode the message and push it in the respective queue */
                mavlink_sys_status_t sys_status; 
                mavlink_msg_sys_status_decode(message, &sys_status);
                current_messages.sys_status.push(sys_status);

                current_messages.time_stamps.sys_status = ptask_gettime(MICRO);
                this_timestamps.sys_status = current_messages.time_stamps.sys_status;
            }
            break;

        case MAVLINK_MSG_ID_VFR_HUD:
            {
                current_messages.hud_m.push(*message);
                //printf("MAVLINK_MSG_ID_VFR_HUD\n");
                
                /* Decode the message and push it in the respective queue */
                mavlink_vfr_hud_t hud;
                mavlink_msg_vfr_hud_decode(message, &hud);
                current_messages.hud.push(hud);
                
                current_messages.time_stamps.hud = ptask_gettime(MICRO);
                this_timestamps.hud = current_messages.time_stamps.hud;
            }
            break;
            
        case MAVLINK_MSG_ID_BATTERY_STATUS:
            {
                current_messages.battery_status_m.push(*message);
                //printf("MAVLINK_MSG_ID_BATTERY_STATUS\n");

                /* Decode the message and push it in the respective queue */
                mavlink_battery_status_t battery_status; 
                mavlink_msg_battery_status_decode(message, &battery_status);
                current_messages.battery_status.push(battery_status);
                
                current_messages.time_stamps.battery_status = ptask_gettime(MICRO);
                this_timestamps.battery_status = current_messages.time_stamps.battery_status;
                break;
            }

        case MAVLINK_MSG_ID_RADIO_STATUS:
            {
                current_messages.radio_status_m.push(*message);
                //printf("MAVLINK_MSG_ID_RADIO_STATUS\n");
                
                /* Decode the message and push it in the respective queue */
                mavlink_radio_status_t radio_status; 
                mavlink_msg_radio_status_decode(message, &radio_status);
                current_messages.radio_status.push(radio_status);
                
                current_messages.time_stamps.radio_status = ptask_gettime(MICRO);
                this_timestamps.radio_status = current_messages.time_stamps.radio_status;
            }
            break;
            
        case MAVLINK_MSG_ID_ENCAPSULATED_DATA:
        {
            current_messages.encapsulated_data_m.push(*message);
            //printf("MAVLINK_MSG_ID_ENCAPSULATED_DATA\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_encapsulated_data_t img;
            mavlink_msg_encapsulated_data_decode(message, &img);
            current_messages.encapsulated_data.push(img);
                
            current_messages.time_stamps.encapsulated_data = ptask_gettime(MICRO);
            this_timestamps.encapsulated_data = current_messages.time_stamps.encapsulated_data;
        }
        break;
        
        case MAVLINK_MSG_ID_NAV_CONTROLLER_OUTPUT:
        {
            current_messages.nav_controller_output_m.push(*message);
            //printf("MAVLINK_MSG_ID_NAV_CONTROLLER_OUTPUT\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_nav_controller_output_t p;
            mavlink_msg_nav_controller_output_decode(message,&p);
            current_messages.nav_controller_output.push(p);
                
            current_messages.time_stamps.nav_controller_output = ptask_gettime(MICRO);
            this_timestamps.nav_controller_output = current_messages.time_stamps.nav_controller_output;
        }
        break;
            
        case MAVLINK_MSG_ID_LOG_ENTRY:
        {
            current_messages.log_entry_m.push(*message);
            //printf("MAVLINK_MSG_ID_LOG_ENTRY\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_log_entry_t log;
            mavlink_msg_log_entry_decode(message, &log);
            current_messages.log_entry.push(log);
                
            current_messages.time_stamps.log_entry = ptask_gettime(MICRO);
            this_timestamps.log_entry = current_messages.time_stamps.log_entry;
        }
        break;
        
        case MAVLINK_MSG_ID_PARAM_VALUE:
        {
            current_messages.param_value_m.push(*message);
            //printf("MAVLINK_MSG_ID_PARAM_VALUE\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_param_value_t rawValue;
            mavlink_msg_param_value_decode(message, &rawValue);
            current_messages.param_value.push(rawValue);
                
            current_messages.time_stamps.param_value = ptask_gettime(MICRO);
            this_timestamps.param_value = current_messages.time_stamps.param_value;
        }
        break;
        
        case MAVLINK_MSG_ID_LOG_DATA:
        {
            current_messages.log_data_m.push(*message);
            //printf("MAVLINK_MSG_ID_LOG_DATA\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_log_data_t log;
            mavlink_msg_log_data_decode(message, &log);
            current_messages.log_data.push(log);
                
            current_messages.time_stamps.log_data = ptask_gettime(MICRO);
            this_timestamps.log_data = current_messages.time_stamps.log_data;
        }
        break;
        
        case MAVLINK_MSG_ID_DATA_TRANSMISSION_HANDSHAKE:
        {
            current_messages.data_transmission_handshake_m.push(*message);
            //printf("MAVLINK_MSG_ID_DATA_TRANSMISSION_HANDSHAKE\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_data_transmission_handshake_t p;
            mavlink_msg_data_transmission_handshake_decode(message, &p);
            current_messages.data_transmission_handshake.push(p);
                
            current_messages.time_stamps.data_transmission_handshake = ptask_gettime(MICRO);
            this_timestamps.data_transmission_handshake = current_messages.time_stamps.data_transmission_handshake;
        }
        break;
        
        case MAVLINK_MSG_ID_SET_POSITION_TARGET_LOCAL_NED:
        {
            current_messages.set_position_target_local_ned_m.push(*message);
            //printf("MAVLINK_MSG_ID_SET_POSITION_TARGET_LOCAL_NED\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_set_position_target_local_ned_t p;
            mavlink_msg_set_position_target_local_ned_decode(message, &p);
            current_messages.set_position_target_local_ned.push(p);
                
            current_messages.time_stamps.set_position_target_local_ned = ptask_gettime(MICRO);
            this_timestamps.set_position_target_local_ned = current_messages.time_stamps.set_position_target_local_ned;
        }
        break;
        
        case MAVLINK_MSG_ID_ATTITUDE_TARGET:
        {
            current_messages.attitude_target_m.push(*message);
            //printf("MAVLINK_MSG_ID_ATTITUDE_TARGET\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_attitude_target_t out;
            mavlink_msg_attitude_target_decode(message, &out);
            current_messages.attitude_target.push(out);
                
            current_messages.time_stamps.attitude_target = ptask_gettime(MICRO);
            this_timestamps.attitude_target = current_messages.time_stamps.attitude_target;
        }
        break;

        case MAVLINK_MSG_ID_GPS_RAW_INT:
        {
            current_messages.gps_raw_int_m.push(*message);
            //printf("MAVLINK_MSG_ID_GPS_RAW_INT\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_gps_raw_int_t pos;
            mavlink_msg_gps_raw_int_decode(message, &pos);
            current_messages.gps_raw_int.push(pos);
                
            current_messages.time_stamps.gps_raw_int = ptask_gettime(MICRO);
            this_timestamps.gps_raw_int = current_messages.time_stamps.gps_raw_int;
        }
        break;

        case MAVLINK_MSG_ID_GPS_STATUS:
        {
            current_messages.gps_status_m.push(*message);
            //printf("MAVLINK_MSG_ID_GPS_STATUS\n");
                
            /* Decode the message and push it in the respective queue */
            mavlink_gps_status_t pos;
            mavlink_msg_gps_status_decode(message, &pos);
            current_messages.gps_status.push(pos);
                
            current_messages.time_stamps.gps_status = ptask_gettime(MICRO);
            this_timestamps.gps_status = current_messages.time_stamps.gps_status;
        }
        break;
            
        case MAVLINK_MSG_ID_LOCAL_POSITION_NED:
            {
                current_messages.local_position_ned_m.push(*message);
                
                //printf("MAVLINK_MSG_ID_LOCAL_POSITION_NED\n");

                /* Decode the message and push it in the respective queue */
                mavlink_local_position_ned_t local_position_ned;
                mavlink_msg_local_position_ned_decode(message, &local_position_ned);
                current_messages.local_position_ned.push(local_position_ned);

                current_messages.time_stamps.local_position_ned = ptask_gettime(MICRO);
                this_timestamps.local_position_ned = current_messages.time_stamps.local_position_ned;
                break;
            }

        case MAVLINK_MSG_ID_GLOBAL_POSITION_INT:
            {
                current_messages.global_position_int_m.push(*message);
                //printf("MAVLINK_MSG_ID_GLOBAL_POSITION_INT\n");

                /* Decode the message and push it in the respective queue */
                mavlink_global_position_int_t global_position_int;
                mavlink_msg_global_position_int_decode(message, &global_position_int);
                current_messages.global_position_int.push(global_position_int);
                
                current_messages.time_stamps.global_position_int = ptask_gettime(MICRO);
                this_timestamps.global_position_int = current_messages.time_stamps.global_position_int;
                break;
            }

        case MAVLINK_MSG_ID_POSITION_TARGET_LOCAL_NED:
            {
                current_messages.position_target_local_ned_m.push(*message);
                //printf("MAVLINK_MSG_ID_POSITION_TARGET_LOCAL_NED\n");

                /* Decode the message and push it in the respective queue */
                mavlink_position_target_local_ned_t position_target_local_ned;
                mavlink_msg_position_target_local_ned_decode(message, &position_target_local_ned);
                current_messages.position_target_local_ned.push(position_target_local_ned);
                
                current_messages.time_stamps.position_target_local_ned = ptask_gettime(MICRO);
                this_timestamps.position_target_local_ned = current_messages.time_stamps.position_target_local_ned;
                break;
            }

        case MAVLINK_MSG_ID_POSITION_TARGET_GLOBAL_INT:
            {
                current_messages.position_target_global_int_m.push(*message);
                //printf("MAVLINK_MSG_ID_POSITION_TARGET_GLOBAL_INT\n");

                /* Decode the message and push it in the respective queue */
                mavlink_position_target_global_int_t position_target_global_int;
                mavlink_msg_position_target_global_int_decode(message, &position_target_global_int);
                current_messages.position_target_global_int.push(position_target_global_int);

                current_messages.time_stamps.position_target_global_int = ptask_gettime(MICRO);
                this_timestamps.position_target_global_int = current_messages.time_stamps.position_target_global_int;
                break;
            }
            
        case MAVLINK_MSG_ID_GLOBAL_VISION_POSITION_ESTIMATE:
        {
            current_messages.global_vision_position_estimate_m.push(*message);
            
            mavlink_global_vision_position_estimate_t pos;
            mavlink_msg_global_vision_position_estimate_decode(message, &pos);
            current_messages.global_vision_position_estimate.push(pos);
            
            current_messages.time_stamps.global_vision_position_estimate = ptask_gettime(MICRO);
                this_timestamps.global_vision_position_estimate = current_messages.time_stamps.global_vision_position_estimate;
        }
        break;

        case MAVLINK_MSG_ID_HIGHRES_IMU:
            {
                current_messages.highres_imu_m.push(*message);
                //printf("MAVLINK_MSG_ID_HIGHRES_IMU\n");

                /* Decode the message and push it in the respective queue */
                mavlink_highres_imu_t highres_imu;
                mavlink_msg_highres_imu_decode(message, &highres_imu);
                current_messages.highres_imu.push(highres_imu);

                current_messages.time_stamps.highres_imu = ptask_gettime(MICRO);
                this_timestamps.highres_imu = current_messages.time_stamps.highres_imu;
               
                break;
            }

        case MAVLINK_MSG_ID_ATTITUDE:
            {
                current_messages.attitude_m.push(*message);
                //printf("MAVLINK_MSG_ID_ATTITUDE\n");

                /* Decode the message and push it into the respective queue */
                mavlink_attitude_t attitude;
                mavlink_msg_attitude_decode(message, &attitude);
                current_messages.attitude.push(attitude);

                current_messages.time_stamps.attitude = ptask_gettime(MICRO);
                this_timestamps.attitude = current_messages.time_stamps.attitude;
            }
            break;
            
        case MAVLINK_MSG_ID_ATTITUDE_QUATERNION:
        {
            current_messages.attitude_quaternion_m.push(*message);
            
            mavlink_attitude_quaternion_t attitude_q;
            mavlink_msg_attitude_quaternion_decode(message, &attitude_q);
            current_messages.attitude_quaternion.push(attitude_q);
            
            current_messages.time_stamps.attitude_quaternion = ptask_gettime(MICRO);
            this_timestamps.attitude_quaternion = current_messages.time_stamps.attitude_quaternion;
        }
        break;

        case MAVLINK_MSG_ID_HIL_CONTROLS:
            {

                current_messages.time_stamps.hil_controls = ptask_gettime(MICRO);
                fprintf(f_aut_THilCtr,"%lu \n",current_messages.time_stamps.hil_controls);

                current_messages.hil_controls_m.push(*message);

                //printf("MAVLINK_MSG_ID_HIL_CONTROLS\n");

                mavlink_hil_controls_t hil_controls;
                mavlink_msg_hil_controls_decode(message, &hil_controls);
                current_messages.hil_controls.push(hil_controls);

                this_timestamps.hil_controls = current_messages.time_stamps.hil_controls;
                hil_controls_count++;
                if ( AUT_INTERFACE_DBG ) 
                {
                    uint64_t curr = current_messages.time_stamps.hil_controls;
                    //printf("IMU Sensors Timestamp %u          \n",current_messages.time_stamps.highres_imu/1000);
                    if ( (curr - read_hil_controls_old) > 10000000 )
                    {
                        printf("HIL_CONTROLS frequency :   %lu Hz\n",hil_controls_count/10);
                        hil_controls_count = 0;
                        read_hil_controls_old = curr;
                        printf("%1.2f | ",hil_controls.roll_ailerons); 
                        printf("%1.2f | ",hil_controls.pitch_elevator);
                        printf("%1.2f | ",hil_controls.yaw_rudder);
                        printf("%1.2f | \n",hil_controls.throttle);
                    }
                } 
            }
            break;

        default:
            {
                // printf("Warning, did not handle message id %i\n",message.msgid);
                break;
            }


    } // end: switch msgid

    return message->msgid;
}

// 
// get_message
//
void Autopilot_Interface::get_message(mavlink_message_t* rqmsg)
{
    //pthread_mutex_lock(&mut_Messages);
    uint8_t msgid = retrieve_message_id();
    // Handle Message ID
    switch (msgid)
    {

        case MAVLINK_MSG_ID_HEARTBEAT:
            {
                //printf("MAVLINK_MSG_ID_HEARTBEAT\n");

                *rqmsg = current_messages.heartbeat_m.front();
                current_messages.heartbeat_m.pop();
                current_messages.heartbeat.pop();
                break;
            }

        case MAVLINK_MSG_ID_SYS_STATUS:
            {
                *rqmsg = current_messages.sys_status_m.front();
                current_messages.sys_status_m.pop();
                current_messages.sys_status.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_VFR_HUD:
            {
                *rqmsg = current_messages.hud_m.front();
                current_messages.hud_m.pop();
                current_messages.hud.pop();
                break;
            }

        case MAVLINK_MSG_ID_BATTERY_STATUS:
            {
                *rqmsg = current_messages.battery_status_m.front();
                current_messages.battery_status_m.pop();
                current_messages.battery_status.pop();
                break;
            }

        case MAVLINK_MSG_ID_RADIO_STATUS:
            {
                *rqmsg = current_messages.radio_status_m.front();
                current_messages.radio_status_m.pop();
                current_messages.radio_status.pop();
                break;
            }

        case MAVLINK_MSG_ID_ENCAPSULATED_DATA:
            {
                *rqmsg = current_messages.encapsulated_data_m.front();
                current_messages.encapsulated_data_m.pop();
                current_messages.encapsulated_data.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_LOCAL_POSITION_NED:
            {
                *rqmsg = current_messages.local_position_ned_m.front();
                current_messages.local_position_ned_m.pop();
                current_messages.local_position_ned.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_GLOBAL_VISION_POSITION_ESTIMATE:
            {
                *rqmsg = current_messages.global_vision_position_estimate_m.front();
                current_messages.global_vision_position_estimate_m.pop();
                current_messages.global_vision_position_estimate.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_NAV_CONTROLLER_OUTPUT:
            {
                *rqmsg = current_messages.nav_controller_output_m.front();
                current_messages.nav_controller_output_m.pop();
                current_messages.nav_controller_output.pop();
                break;
            }

        case MAVLINK_MSG_ID_LOG_ENTRY:
            {
                *rqmsg = current_messages.log_entry_m.front();
                current_messages.log_entry_m.pop();
                current_messages.log_entry.pop();
                break;
            }    
        
        case MAVLINK_MSG_ID_PARAM_VALUE:
            {
                *rqmsg = current_messages.param_value_m.front();
                current_messages.param_value_m.pop();
                current_messages.param_value.pop();
                break;
            }
        
        case MAVLINK_MSG_ID_LOG_DATA:
            {
                *rqmsg = current_messages.log_data_m.front();
                current_messages.log_data_m.pop();
                current_messages.log_data.pop();
                break;
            }
            
            case MAVLINK_MSG_ID_ATTITUDE_QUATERNION:
            {
                *rqmsg = current_messages.attitude_quaternion_m.front();
                current_messages.attitude_quaternion_m.pop();
                current_messages.attitude_quaternion.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_DATA_TRANSMISSION_HANDSHAKE:
            {
                *rqmsg = current_messages.data_transmission_handshake_m.front();
                current_messages.data_transmission_handshake_m.pop();
                current_messages.data_transmission_handshake.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_SET_POSITION_TARGET_LOCAL_NED:
            {
                *rqmsg = current_messages.set_position_target_local_ned_m.front();
                current_messages.set_position_target_local_ned_m.pop();
                current_messages.set_position_target_local_ned.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_ATTITUDE_TARGET:
            {
                *rqmsg = current_messages.attitude_target_m.front();
                current_messages.attitude_target_m.pop();
                current_messages.attitude_target.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_GPS_RAW_INT:
            {
                *rqmsg = current_messages.gps_raw_int_m.front();
                current_messages.gps_raw_int_m.pop();
                current_messages.gps_raw_int.pop();
                break;
            }
            
        case MAVLINK_MSG_ID_GLOBAL_POSITION_INT:
            {
                *rqmsg = current_messages.global_position_int_m.front();
                current_messages.global_position_int_m.pop();
                current_messages.global_position_int.pop();
                break;
            }

        case MAVLINK_MSG_ID_POSITION_TARGET_LOCAL_NED:
            {
                *rqmsg = current_messages.position_target_local_ned_m.front();
                current_messages.position_target_local_ned_m.pop();
                current_messages.position_target_local_ned.pop();
                break;
            }

        case MAVLINK_MSG_ID_POSITION_TARGET_GLOBAL_INT:
            {
                *rqmsg = current_messages.position_target_global_int_m.front();
                current_messages.position_target_global_int_m.pop();
                current_messages.position_target_global_int.pop();
                break;
            }

            
            case MAVLINK_MSG_ID_GPS_STATUS:
            {
                *rqmsg = current_messages.gps_status_m.front();
                current_messages.gps_status_m.pop();
                current_messages.gps_status.pop();

                //printf("MAVLINK_MSG_ID_HIGHRES_IMU\n");
                break;
            }
            
        case MAVLINK_MSG_ID_HIGHRES_IMU:
            {
                *rqmsg = current_messages.highres_imu_m.front();
                current_messages.highres_imu_m.pop();
                current_messages.highres_imu.pop();

                //printf("MAVLINK_MSG_ID_HIGHRES_IMU\n");
                break;
            }

        case MAVLINK_MSG_ID_ATTITUDE:
            {
                *rqmsg = current_messages.attitude_m.front();
                current_messages.attitude_m.pop();
                current_messages.attitude.pop();

                break;
            }

        case MAVLINK_MSG_ID_HIL_CONTROLS:
            {
                *rqmsg = current_messages.hil_controls_m.front();
                current_messages.hil_controls_m.pop();
                current_messages.hil_controls.pop();

                //printf("MAVLINK_MSG_ID_HIL_CONTROLS\n");
                break;
            }

        default:
            {
                // printf("Warning, did not handle message id %i\n",message.msgid);
                break;
            }


    } // end: switch msgid
    //pthread_mutex_unlock(&mut_Messages);
}

//
// retrieve_message_id
//
int Autopilot_Interface::retrieve_message_id()
{
    //pthread_mutex_lock(&mut_Messages);
    int index = queueIndexFetched.front();
    if (index < 0)
        printf("Error retrieving message id\n");
    queueIndexFetched.pop();
    //pthread_mutex_unlock(&mut_Messages);

    return index;
}

//
// retrieve_n_read_messages
//
int Autopilot_Interface::retrieve_nread_messages()
{
    pthread_mutex_lock(&mut_Messages);
    int dimension = queueIndexFetched.size();
    pthread_mutex_unlock(&mut_Messages);
    
    return dimension;
}


//
// send_message
//
int Autopilot_Interface::send_message(mavlink_message_t* message)
{
    uint16_t len = 0;   
    
    char buf[300];
    len = mavlink_msg_to_send_buffer((uint8_t *)buf, message);

    // Write buffer to serial port, locks port while writing
    int writtenB = uart_port.write_bytes(buf,len);
    
    if (len != writtenB)
    {
        printf("Autopilot_Interface::write_message: ERROR WHILE WRITING"); 
        return -1;
    }
    // Done!
    return writtenB;
}

//
// send_sensorData
//
int Autopilot_Interface::send_sensorData(mavlink_message_t* message)
{
    uint16_t len = 0;

    char buf[300];
    len = mavlink_msg_to_send_buffer((uint8_t *)buf, message);

    // Write buffer to serial port, locks port while writing
    int writtenB = uart_port.write_bytes(buf,len);
    
    if (len != writtenB)
    {
        printf("Autopilot_Interface::write_message: ERROR WHILE WRITING"); 
        return -1;
    }

    // Done!
    return writtenB;
}



// -----------------------------------------------------------------------
//   Write SensorHil message
// -----------------------------------------------------------------------
void Autopilot_Interface::write_hilsensors()
{
    // --------------------------------------------------------------------
    //   PACK PAYLOAD
    // --------------------------------------------------------------------

    // Allocate Space 
    mavlink_message_t msg;
    mavlink_hil_sensor_t sensors;


    ptime timestamp = ptask_gettime(MICRO);

    // Pack the message
    /*
       mavlink_msg_hil_sensor_pack(system_id, autopilot_id, &msg, timestamp,
       0.0f, 0.0f, -9.8f,
       0.0f, 0.0f, 0.0f,
       0.2f, 0.0f, 0.5f,
       1023.0f, 0.0f,0.0f,
       20.0f, 0xFFF);
       */
    sensors.time_usec = timestamp; ///< Timestamp (microseconds)
    sensors.xacc = 0.0f; ///< X acceleration (m/s^2)
    sensors.yacc = 0.0f; ///< Y acceleration (m/s^2)
    sensors.zacc = -9.8f; ///< Z acceleration (m/s^2)
    sensors.xgyro = 0.0f; ///< Angular speed around X axis in body frame (rad / sec)
    sensors.ygyro = 0.0f; ///< Angular speed around Y axis in body frame (rad / sec)
    sensors.zgyro = 0.0f; ///< Angular speed around Z axis in body frame (rad / sec)
    sensors.xmag = 0.2f; ///< X Magnetic field (Gauss)
    sensors.ymag = 0.0f; ///< Y Magnetic field (Gauss)
    sensors.zmag = 0.5f; ///< Z Magnetic field (Gauss)
    sensors.abs_pressure = 1023.0f; ///< Absolute pressure in millibar
    sensors.diff_pressure = 0.0f; ///< Differential pressure (airspeed) in millibar
    sensors.pressure_alt = 0.0f; ///< Altitude calculated from pressure
    sensors.temperature = 20.0f; ///< Temperature in degrees celsius
    sensors.fields_updated = 0xFF; ///< Bitmask for fields that have updated 

    mavlink_msg_hil_sensor_encode(system_id,companion_id,&msg,&sensors);

    // -------------------------------------------------------------------
    //   WRITE
    // -------------------------------------------------------------------

    char buf[300]; 
    uint16_t len = mavlink_msg_to_send_buffer((uint8_t*)buf, &msg);

    int writtenB = uart_port.write_bytes(buf,len);

    // check the write
    if ( writtenB < 0 )
        fprintf(stderr,"WARNING: could not send HIL_SENSORS \n");

    return;
}




// -----------------------------------------------------------------------
//   Write Ping message
// -----------------------------------------------------------------------
void
Autopilot_Interface::
write_ping()
{
    // -------------------------------------------------------------------
    //   PACK PAYLOAD
    // -------------------------------------------------------------------


    return;
}


void
Autopilot_Interface::
write_raw()
{
    int i;
    char send_buff[5];
    for (i = 0 ; i<sizeof(send_buff);i++)
    {
        send_buff[i] = i+1;
    }

    int len = uart_port.write_bytes(send_buff,sizeof(send_buff));
    printf("Written  %d bytes \n",len);
    printf("%d %d %d %d %d \n", send_buff[0], send_buff[1], send_buff[2], send_buff[3], send_buff[4]);
    return ;
}


// -----------------------------------------------------------------------
//   Start Hil Mode
// -----------------------------------------------------------------------
void Autopilot_Interface::start_hil()
{
    char buf[300];
    uint8_t newBaseMode;
    uint32_t newCustomMode;

    newBaseMode = ( base_mode | MAV_MODE_FLAG_HIL_ENABLED );
    newBaseMode &= ~MAV_MODE_FLAG_SAFETY_ARMED;

    newBaseMode |= ( base_mode & MAV_MODE_FLAG_SAFETY_ARMED );
    
    newCustomMode = custom_mode;

    printf("BaseMode    : %u \n", base_mode);
    printf("newBaseMode : %u \n", newBaseMode);
    newBaseMode = 113;
    mavlink_message_t msg;
    mavlink_msg_set_mode_pack(system_id, 0, &msg, (uint8_t)1, newBaseMode, newCustomMode);
    //mavlink_msg_set_mode_pack(255, 0, &msg, (uint8_t)1, newBaseMode, 65536); 
    uint16_t len = mavlink_msg_to_send_buffer((uint8_t*)buf, &msg);
    printf("Setting HIL mode \n");
    int attempts = 0;
    bool condition = 1;
    while( condition )
    {
        attempts++;
        
        // Write the request message to serial
        printf("Sending Request #%d \n",attempts);
        int writtenB = uart_port.write_bytes(buf,len);
      

        // We must wait for an acknowledgment for the new state, that is
        // a new heartbeat message where the state of the AUV is 
        // contained
        pthread_mutex_lock(&mut_heartbeat);
        new_heartbeat = 0;
        if(!new_heartbeat)
        {
            printf("Waiting for update...\n");
            pthread_cond_wait(&cond_heartbeat,&mut_heartbeat);
        }

        // Check the condition for the re-try
        condition = (base_mode & MAV_MODE_FLAG_HIL_ENABLED) == 0;
        // Release the resource 
        pthread_mutex_unlock(&mut_heartbeat);
        
    }

    printf("After Check BaseMode    : %u \n", base_mode);

    hil_mode = true;
    printf("Message = \n");
    printf("len     : %u \n", msg.len);
    printf("seq     : %u \n", msg.seq);
    printf("sysid   : %u \n", msg.sysid);
    printf("compid  : %u \n", msg.compid);
    printf("msgid   : %u \n", msg.msgid);
    char* pchar = (char*)(&msg.payload64[0]);
    printf("payload : %u %u %u %u %u %u \n", *(pchar+5), *(pchar+4), *(pchar+3), *(pchar+2), *(pchar+1), *(pchar+0));
    printf("HIL mode set!  [%d attempts]\n",attempts);
}


// ----------------------------------------------------------------------
//   Stop Hil Mode
// ----------------------------------------------------------------------
void Autopilot_Interface::stop_hil()
{
    
    char buf[300];

    uint8_t newBaseMode;
    uint32_t newCustomMode = 0;
    
    newBaseMode = ( base_mode & ~MAV_MODE_FLAG_HIL_ENABLED );

    newBaseMode &= ~MAV_MODE_FLAG_SAFETY_ARMED;
    newBaseMode |= ( base_mode & MAV_MODE_FLAG_SAFETY_ARMED );

    printf("BaseMode    : %u \n", base_mode);
    printf("newBaseMode : %u \n", newBaseMode);

    // Create the message and put into the buffer
    mavlink_message_t msg;
    mavlink_msg_set_mode_pack(system_id, 0, &msg, (uint8_t)1, newBaseMode, newCustomMode);
    uint16_t len = mavlink_msg_to_send_buffer((uint8_t*)buf, &msg);

    printf("Unsetting HIL mode \n");
    int attempts = 0;
    bool condition = 1;
    while( condition )
    { 
        attempts++;

       // Write the request message to serial
        printf("Sending Request #%d \n",attempts);
        int writtenB = uart_port.write_bytes(buf,len);
      
        // We must wait for an acknowledgment for the new state, that is
        // a new heartbeat message where the state of the AUV is 
        // contained
        pthread_mutex_lock(&mut_heartbeat);
        new_heartbeat = 0;
        if(!new_heartbeat)
        {
            printf("Waiting for update...\n");
            pthread_cond_wait(&cond_heartbeat,&mut_heartbeat);
        }

       
        // Check if it is necessary to re-try
        condition = (base_mode & MAV_MODE_FLAG_HIL_ENABLED) != 0;
        printf("After Check BaseMode    : %u \n", base_mode);
      // Release the resource 
        pthread_mutex_unlock(&mut_heartbeat);
    }

    hil_mode = false;

    printf("Message = \n");
    printf("len     : %u \n", msg.len);
    printf("seq     : %u \n", msg.seq);
    printf("sysid   : %u \n", msg.sysid);
    printf("compid  : %u \n", msg.compid);
    printf("msgid   : %u \n", msg.msgid);
    char* pchar = (char*)(&msg.payload64[0]);
    printf("payload : %u %u %u %u %u %u \n", *(pchar+5), *(pchar+4), *(pchar+3), *(pchar+2), *(pchar+1), *(pchar+0));
    printf("HIL mode unset!  [%d attempts]\n",attempts);
}

// ----------------------------------------------------------------------
//   Check if it is in HIL mode
// ----------------------------------------------------------------------
bool Autopilot_Interface::is_hil()
{
    return hil_mode;    
}


// -----------------------------------------------------------------------
//   Start Off-Board Mode
// -----------------------------------------------------------------------
void
Autopilot_Interface::
enable_offboard_control()
{
    // Should only send this command once
    if ( control_status == false )
    {
        printf("ENABLE OFFBOARD MODE\n");

        // --------------------------------------------------------------
        //   TOGGLE OFF-BOARD MODE
        // --------------------------------------------------------------

        // Sends the command to go off-board
        int success = toggle_offboard_control( true );

        // Check the command was written
        if ( success )
            control_status = true;
        else
        {
            fprintf(stderr,"Error: off-board mode not set, could not write message\n");
            //throw EXIT_FAILURE;
        }

        printf("\n");

    } // end: if not offboard_status

}


// ------------------------------------------------------------------------------
//   Stop Off-Board Mode
// ------------------------------------------------------------------------------
void
Autopilot_Interface::
disable_offboard_control()
{

    // Should only send this command once
    if ( control_status == true )
    {
        printf("DISABLE OFFBOARD MODE\n");

        // ----------------------------------------------------------------------
        //   TOGGLE OFF-BOARD MODE
        // ----------------------------------------------------------------------

        // Sends the command to stop off-board
        int success = toggle_offboard_control( false );

        // Check the command was written
        if ( success )
            control_status = false;
        else
        {
            fprintf(stderr,"Error: off-board mode not set, could not write message\n");
            //throw EXIT_FAILURE;
        }

        printf("\n");

    } // end: if offboard_status

}


// ------------------------------------------------------------------------------
//   Toggle Off-Board Mode
// ------------------------------------------------------------------------------
int Autopilot_Interface::toggle_offboard_control( bool flag )
{
    // Prepare command for off-board mode
    mavlink_command_long_t com;
    com.target_system    = system_id;
    com.target_component = autopilot_id;
    com.command          = MAV_CMD_NAV_GUIDED_ENABLE;
    com.confirmation     = true;
    com.param1           = (float) flag; // flag >0.5 => start, <0.5 => stop

    // Encode
    mavlink_message_t message;
    mavlink_msg_command_long_encode(system_id, companion_id, &message, &com);

    // Send the message
    char buf[300];
    
    uint16_t len = mavlink_msg_to_send_buffer((uint8_t*)buf, &message);
    // Write buffer to serial port, locks port while writing
    int writtenB = uart_port.write_bytes(buf,len);

    // Done!
    return writtenB;
}


// End Autopilot_Interface
