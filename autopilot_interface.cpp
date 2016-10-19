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

#define AUT_INTERFACE_DBG 0 



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

    //f_aut_THilCtr = fopen("./T_aut_HilCtr.txt","w");
    //f_aut_TSens = fopen("./T_aut_Sens.txt","w");

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
    //fclose(f_aut_THilCtr);
    //fclose(f_aut_TSens);

    printf("DESTRUCTOR\n");
}

// -------------------------------------------------------- 
//  FUNCTIONS
// --------------------------------------------------------

//
// Fetch Data
//
int Autopilot_Interface::fetch_data()
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
                // queueIndexFetched.push(message_Id);

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
    int message_id = message->msgid;
    // Push the message in the queue
    current_messages.messages.push(*message);
    // Record the time
    current_messages.time_stamps[message_id] = ptask_gettime(MICRO);
    
    // Handle Message ID
    switch (message_id)
    {
        case MAVLINK_MSG_ID_HEARTBEAT:
            {
                mavlink_heartbeat_t heartbeat;
                mavlink_msg_heartbeat_decode(message, &heartbeat);
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
            }
            break;

        case MAVLINK_MSG_ID_HIL_CONTROLS:
            {
                //fprintf(f_aut_THilCtr,"%lu \n",current_messages.time_stamps[message_id]);
                //printf("MAVLINK_MSG_ID_HIL_CONTROLS\n");

                mavlink_hil_controls_t hil_controls;
                mavlink_msg_hil_controls_decode(message, &hil_controls);

                hil_controls_count++;
                if ( AUT_INTERFACE_DBG ) 
                {
                    uint64_t curr = current_messages.time_stamps[message_id];
                    //printf("IMU Sensors Timestamp %u          \n",current_messages.time_stamps.highres_imu/1000);
                    if ( (curr - read_hil_controls_old) > 10000000 )
                    {
                        printf("HIL_CONTROLS frequency :   %lu Hz\n",hil_controls_count/10);
                        hil_controls_count = 0;
                        read_hil_controls_old = curr;
                    }
                } 
            }
            break;

        default:
            break;
    } 
	
//	printf("Message ID from AUV : %d\n", message_id);
    return message_id;
}

// 
// get_message
//
int Autopilot_Interface::get_message(mavlink_message_t* rqmsg)
{
    
    // Extract the message from the front of the queue
    mavlink_message_t message = current_messages.messages.front();
    current_messages.messages.pop();
    
    *rqmsg = message;
    return message.msgid;
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
    //newBaseMode = 32;
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
