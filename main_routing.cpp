/*
 * file: main_routing.cpp
 *
 */


#include "main_routing.h"
extern "C"
{
    #include "DynModel.h"
    #include "DynModel_private.h"
}

int time_to_exit = 0;

int main(int argc, char *argv[])
{
	int i; 
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Ports for the communication with GS
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	unsigned int gs_r_port = 14551;
	unsigned int gs_w_port = 14550;
	char* gs_ip = (char *)"127.0.0.1";
	/*
	 *                  +--------+
	 * UDP (14550) ---> |   GS   | ----> (14551)  UDP
	 *                  +--------+
	 */

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Port for the communication with the Board
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	char* uart_name = (char *)"/dev/ttyUSB0";
	int baudrate = 921600;
	/*
	 *                           +---------+
	 *                           |         |
	 *                      ---->|         |
	 * UART /dev/ttyUSB0         | Pixhawk |
	 *                      <----|         |
	 *                           |         |
	 *                           +---------+
	 */

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Ports for the communication with Unreal Engine
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	unsigned int ue_r_port = 8001;
	unsigned int ue_w_port = 8000;
	char* ue_ip = (char *)"10.30.3.136";
	/*
	 *                  +--------+
	 * UDP (8000)  ---> |   UE   | ----> (8001)  UDP
	 *                  +--------+
	 */


    
    //======================================================================
    // Initialization
    //======================================================================
	// Open Files descriptors to write inn
	/*
	file_TFtcSns = fopen("Times_FtcSens.txt", "w");
	file_TSndSns = fopen("Times_SndSens.txt", "w");
	file_TFtcComm = fopen("Times_FtcCom.txt", "w");
	file_TSndComm = fopen("Times_SndCom.txt", "w");
	file_TGS = fopen("Times_GS.txt", "w");
    */
	UAV_base_mode = 0;

    // Initialize the global Mutex and Condition Variable
    // for the first heartbeat check
	pthread_mutex_init(&mut_first_heartbeat, 0);
	pthread_cond_init(&cond_first_heartbeat, 0);
    
	autopilot_connected = false; // Variable associated with the first heartbeat

	pbarrier_init(&barrier, 2); // Barrier for the synch of simulator/inflow tasks

    // Initialize the autogenerated code
    DynModel_initialize();

	// --------------------------------------------------------------------
	//   START PTASK ENVIRONMENT AND RUN THREADS 
	// --------------------------------------------------------------------
	ptask_init(SCHED_OTHER, GLOBAL, NO_PROTOCOL);


	// --------------------------------------------------------------------
	//   PARSE THE COMMANDS
	// --------------------------------------------------------------------
	// Parse the Command Line 
	parse_commandline(argc, argv, gs_ip, gs_r_port, gs_w_port);

	// --------------------------------------
	//    INSTANTIATE CLASSES
	// --------------------------------------	
	/*
	 * Instantiate an autopilot interface object
	 * This object handles the communication with the Autopilot board.
	 * The communication is accomplished over the Serial. The Serial
	 * port management is done through Serial port class instantiated
	 * inside the Autopilot_Interface object.
	 */
	Autopilot_Interface autopilot_interface(uart_name, baudrate);

	/*
	 * Instantiate an ground station interface object
	 * This object handles the communication with the ground station.
	 * The communication is accomplished over the UDP. The UDP 
	 * port management is done through Serial port class instantiated
	 * inside the GS_Interface object.
	 */
	GS_Interface gs_interface(gs_ip, gs_r_port, gs_w_port);


	/*
 	 * Instatiate the interface for the management of the 
 	 * Unreal Engine 
 	 */ 
	UE_Interface ue_interface(ue_ip, ue_r_port, ue_w_port);


	// --------------------------------------
	//    POINTERS TO CLASSES AND HANDLERS
	// --------------------------------------
	/*
	 * Setup interrupt signal handler
	 *
	 * Responds to early exits signaled with Ctrl-C.  The handler will 
	 * command the UAV to return the standar operative mode, and close 
	 * threads and the port.
	 *
	 */
	autopilot_interface_quit    = &autopilot_interface;
	signal(SIGINT,quit_handler);

	struct Interfaces point_to_interfaces;
	point_to_interfaces.gs = &gs_interface;
	point_to_interfaces.aut = &autopilot_interface;
	point_to_interfaces.ue = &ue_interface;


	//======================================================================
	//    THREADS  
	//======================================================================
	// Setting periods
	aut_period = tspec_from(4, MILLI);
	sim_period = tspec_from(4, MILLI);
	gs_period = tspec_from(4, MILLI);
	ue_period = tspec_from(20, MILLI);

	tspec_init();

	/************ INFLOW THREAD *************/

	// Defining the Structure containing the thread Properties 
	tpars params_inflow = TASK_SPEC_DFL;
	params_inflow.period = aut_period;
	params_inflow.rdline = aut_period;
	params_inflow.priority = 90;
	params_inflow.act_flag = NOW;
	params_inflow.measure_flag = 0;
	params_inflow.processor = 0;
	params_inflow.arg = &point_to_interfaces;

	// Start to retrieve data from the Autopilot Interface
	//
	// Lock the resource autopilot_connected, then wait for the first heartbeat
	// before going on...
	printf("Creating the inflow_thread....\n"); 
	pthread_mutex_lock(&mut_first_heartbeat);
	inflowT_id = ptask_create_param(inflow_thread, &params_inflow);
	if (inflowT_id == -1)
	{
		pthread_mutex_unlock(&mut_first_heartbeat);
		printf("Inflow Thread not created!");
		return -1;
	}
	if(!autopilot_connected)
	{
		printf("Waiting for autopilot connection...\n");
		pthread_cond_wait(&cond_first_heartbeat, &mut_first_heartbeat);
	}
	printf("Connected!\n");
	pthread_mutex_unlock(&mut_first_heartbeat);

	// Now we are receiving messages from the autopilot board

	/************* SIMULATOR THREAD ***************/
	// Defining the Structure containing the thread Properties 
	tpars params_sim = TASK_SPEC_DFL;
	params_sim.period = sim_period;
	params_sim.rdline = sim_period;
	params_sim.priority = 90;
	params_sim.act_flag = NOW;
	params_sim.measure_flag = 0;
	params_sim.processor = 0;
	params_sim.arg = &point_to_interfaces;


	simT_id = ptask_create_param(simulator_thread, &params_sim);
	if(simT_id == -1)
	{ 
		printf("Error creating the simulator thread\n");
	}

	/************ SWITCHING TO HIL MODE *************/
	autopilot_interface.start_hil();

	// Change the priority of the task related to the Ground Station
	tpars params_gs = TASK_SPEC_DFL;
	params_gs.period = gs_period;
	params_gs.rdline = gs_period;
	params_gs.priority = 80;
	params_gs.act_flag = NOW;
	params_gs.measure_flag = 0;
	params_gs.processor = 0;
	params_gs.arg = &point_to_interfaces;
	gsT_id = ptask_create_param(gs_thread, &params_gs);
	if(gsT_id == -1)
	{
		printf("Error creating the Ground Station thread\n");
	}

	// Change the priority of the task related to the Unreal Engine 
	tpars params_ue = TASK_SPEC_DFL;
	params_ue.period = ue_period;
	params_ue.rdline = ue_period;
	params_ue.priority = 70;
	params_ue.act_flag = NOW;
	params_ue.measure_flag = 0;
	params_ue.processor = 0;
	params_ue.arg = &point_to_interfaces;
	ueT_id = ptask_create_param(ue_thread, &params_ue);
	if(ueT_id == -1)
	{
		printf("Error creating the Ground Station thread\n");
	}


	for(;;)
		usleep(500000);

	return 0;

}

// -------------------------------------------------------
//  Function to manage the retrieved message
//
//  Unpack the control message
//
// -------------------------------------------------------
void unpack_ctr_mess(mavlink_message_t *msg)
{
    // Extract the timestamp of the control generated by the board
    uint64_t rec_time = mavlink_msg_hil_controls_get_time_usec(msg);
	hil_ctr[0] = mavlink_msg_hil_controls_get_roll_ailerons(msg);
    hil_ctr[1] = mavlink_msg_hil_controls_get_pitch_elevator(msg);
	hil_ctr[2] = mavlink_msg_hil_controls_get_yaw_rudder(msg);
	hil_ctr[3] = mavlink_msg_hil_controls_get_throttle(msg);
    
    // Update the Input Structure 
    DynModel_U.PWM1 = hil_ctr[2];
    DynModel_U.PWM2 = hil_ctr[0];
    DynModel_U.PWM3 = hil_ctr[3];
    DynModel_U.PWM4 = hil_ctr[1];
}


// ----------------------------------------------------------------------
//    INFLOW THREAD
// ----------------------------------------------------------------------
/*
 * This thread takes the information from the PX4 and directs them 
 * towards the simulator and the ground station. 
 *
 *
 *   PX4 >-----+----> Simulator
 *             |
 *             +----> Ground Station
 *
 *
 */
void inflow_thread()
{
	inflow_thread_active = true;
	struct Interfaces* p = (struct Interfaces*)ptask_get_argument();

	int counter = 0;
	int i;
	int NMessRead = 0;  // Number of messages read
	int rec_message_id; // Id of the current message
	mavlink_message_t rec_message; // Current received mavlink message

	// Calculate the number of bytes to represent the controls
	// originally represented as floats
	const int NFloatCont = 4;
	const int NBytesCont = NFloatCont * sizeof(float);
	const unsigned int DIM_BUFF = 256;

	int tid = ptask_get_index();

	bool first = true;  // Flag for the first synchronization

	// Times for the control sending loop
	ptime hil_ctr_rec_time = ptask_gettime(MICRO);
	ptime hil_ctr_rec_time_old = ptask_gettime(MICRO); 

	// Initialize the control vector
	for (i = 0; i < 4; i++)
		hil_ctr[i] = 0;

	inflow_end_time = 0;

	printf("***  Starting Inflow Thread  ***\n");
	// The autopilot interface has been configured
	while (!time_to_exit)
	{
        //ptask_set_deadline(tid, ptask_get_period(tid, MICRO), MICRO);
		// Access the Autopilot interface and obtain the 
		// number of read messages 
		NMessRead = p->aut->fetch_data();

		// If something has been found, then manage it:
		if (NMessRead > 0)
		{
			// Retrieve the messages
			for (i = 0; i < NMessRead; i++)
			{
				p->aut->get_message(&rec_message);
				// Handle the message selecting the target 
				routing_messages(&rec_message, p);
			}
		}

		// If not synchonized and we have already received info from the
		// the autopilot board
		if (first && autopilot_connected)
		{
			pbarrier_wait(&barrier, 0);
			first = false;
		}

		// Send controls
		// (Don't send messages if you don't have passed the synchronization step)
		if (!first && p->aut->is_hil() && 
				simulator_thread_active && 
				autopilot_connected)
		{
            // Update the Input Structure 
            /*
			if (p->aut->is_hil() && counter > 500) 
			{
            	printf("HIL_CTR : %1.2f | %1.2f | %1.2f | %1.2f\n", 
						hil_ctr[0], hil_ctr[1], hil_ctr[2], hil_ctr[3]);
				counter = 0;
            }
			counter++;
			
			//p->sim->sendActuatorCommand(hil_ctr, NFloatCont);
			hil_ctr_time = ptask_gettime(MICRO);
			//fprintf(file_TSndComm,"%lu \n", hil_ctr_time);
			*/
		}

		/*
        if (ptask_deadline_miss())
        {
			printf("Inflow Thread Deadline Miss!  Abs deadline = %lu || Time = %lu\n", 
                   task_absdl(tid), ptask_gettime(MILLI));
        }
        */
		ptask_wait_for_period();


	}
	p->aut->stop_hil();
}

// -------------------------------------------------------
//  Function to manage the retrieved messages
//
//  Send messages to the GS or to the SIM
//
// -------------------------------------------------------
void routing_messages(mavlink_message_t *msg, struct Interfaces* p)
{
	switch (msg->msgid) 
	{
		// TO SIMULATOR
		case MAVLINK_MSG_ID_HIL_CONTROLS:
			if (p->aut->is_hil() && simulator_thread_active)
			{
				// Extract the timestamp of the control generated by the board
				uint64_t rec_time = mavlink_msg_hil_controls_get_time_usec(msg);
				unpack_ctr_mess(msg);
				
			}
			break;

			// TO GROUND STATION
		case MAVLINK_MSG_ID_HEARTBEAT:
			UAV_base_mode = mavlink_msg_heartbeat_get_base_mode(msg); 
			pthread_mutex_lock(&mut_first_heartbeat);
			autopilot_connected = true;
			pthread_cond_signal(&cond_first_heartbeat);
			pthread_mutex_unlock(&mut_first_heartbeat);
			if (gs_thread_active)
				p->gs->pushMessage(msg);
			break;

		default:
		if (gs_thread_active)
				p->gs->pushMessage(msg);
            break;
	}
}


// ----------------------------------------------------------------------
//    SIMULATION THREAD
// ----------------------------------------------------------------------
/*
 *
 *
 *    Simulator >--------------> PX4  
 *                         
 *
 */
void simulator_thread()
{
    int i = 0;
	int nread = 0;

	int tid = ptask_get_index();
	struct Interfaces* p = (struct Interfaces*)ptask_get_argument();

	int ret_sens;
	int ret_gps;

	uint8_t system_id = p->aut->system_id;
	uint8_t component_id = p->aut->autopilot_id;

	mavlink_message_t sensor_msg;
	mavlink_message_t gps_msg;
	mavlink_status_t status;

	bool first = true;
	const unsigned int NFloatSens = 22;

	float      xacc;
	float      yacc;
	float      zacc;
	float      xgyro;
	float      ygyro;
	float      zgyro;
	float      xmag;
	float      ymag;
	float      zmag;
	float      abs_pressure;
	float      diff_pressure;
	float      pressure_alt;
	float      temperature;
	uint32_t   fields_updated;

	uint8_t    fix_type;
	int32_t    lat;
	int32_t    lon;
	int32_t    alt;
	uint16_t   eph;
	uint16_t   epv;
	uint16_t   vel;
	int16_t    vn;
	int16_t    ve;
	int16_t    vd;
	int16_t    cog;
	uint8_t    satellites_visible;
    
    uint64_t   time_usec;

	int counter = 0;

	printf("***  Starting Simulator Thread  ***\n");
	simulator_thread_active = true;
    
    ptime old_sent_time = 0;

	// Check the initialization of the necessary classes
	while (!time_to_exit)
	{
        //ptask_set_deadline(tid, ptask_get_period(tid, MICRO), MICRO);
		if (first)
		{
			first = false;
			pbarrier_wait(&barrier, 0);
			printf("Simulation Thread STARTED! \n");
		}

		DynModel_step();

		/*
		if (counter > 200)
		{
		printf("R = %3.2f\n", DynModel_Y.RPY[0] * 180.0 / 3.14);
		printf("P = %3.2f\n", DynModel_Y.RPY[1] * 180.0 / 3.14);
		printf("Y = %3.2f\n\n", DynModel_Y.RPY[2] * 180.0 / 3.14);
		printf("Thrusts = [%3.2f, %3.2f, %3.2f, %3.2f]\n\n", 
				DynModel_Y.Thrusts[0],
				DynModel_Y.Thrusts[1],
				DynModel_Y.Thrusts[2],
				DynModel_Y.Thrusts[3]);
		printf("Torques = [%3.2f, %3.2f, %3.2f]\n\n", 
				DynModel_Y.Torques[0],
				DynModel_Y.Torques[1],
				DynModel_Y.Torques[2]);

		
		counter = 0;
		}
		counter++;
		*/
        time_usec = ptask_gettime(MICRO);

        xacc = (float)DynModel_Y.Accelerometer[0];
        yacc = (float)DynModel_Y.Accelerometer[1]; 
        zacc = (float)DynModel_Y.Accelerometer[2];
        xgyro = (float)DynModel_Y.Gyro[0];
        ygyro = (float)DynModel_Y.Gyro[1];
        zgyro = (float)DynModel_Y.Gyro[2];
        xmag = (float)DynModel_Y.Magn[0];
        ymag = (float)DynModel_Y.Magn[1];
        zmag = (float)DynModel_Y.Magn[2];
        abs_pressure = (float)DynModel_Y.Press;
        diff_pressure = (float)DynModel_Y.diff_Pres;
        pressure_alt = (float)DynModel_Y.Baro_Alt;
        temperature = (float)DynModel_Y.Temp;
        fields_updated = (uint32_t)0xFF;
        fix_type = 3;
        lat = (int32_t)(DynModel_Y.Gps_Lat * 1e7);
        lon = (int32_t)(DynModel_Y.Gps_Lon * 1e7);
        alt = (int32_t)(DynModel_Y.Gps_Alt * 1e3);
        eph = 1;
        epv = 1;
        vel = (uint16_t)(DynModel_Y.Gps_V_Mod * 100); // cm/s
        vn  = (int16_t)(DynModel_Y.Gps_V[0] * 100); 
        ve  = (int16_t)(DynModel_Y.Gps_V[1] * 100);
        vd  = (int16_t)(DynModel_Y.Gps_V[2] * 100);
        cog = (int16_t)(DynModel_Y.COG * 100);  
        satellites_visible = 8;


		printf("Accelerometers : %3.2f | %3.2F | %3.2f\n", xacc, yacc, zacc);
        // Composition of the mavlink messages  
        //  Sensors Message 
        mavlink_msg_hil_sensor_pack(system_id, component_id, &sensor_msg, time_usec, 
                xacc, yacc, zacc, xgyro, ygyro, zgyro, xmag, ymag, zmag, abs_pressure, 
                diff_pressure, pressure_alt, temperature, fields_updated);
        //  GPS Message
        mavlink_msg_hil_gps_pack(system_id, component_id, &gps_msg, 
                time_usec, fix_type, lat, lon, alt, eph, epv, 
                vel, vn, ve, vd, cog, satellites_visible);
		if (p->aut->is_hil())
		{
            // Send Sensor Data to Board
			ret_sens = p->aut->send_message(&sensor_msg);
            
            // Record Sending Time
			ptime sendTime = ptask_gettime(MICRO);
			//fprintf(file_TSndSns,"%lu \n", sendTime);

            // Send GPS data to Board
			if ( (time_usec - old_sent_time) > 400000)
			{
				ret_gps = p->aut->send_message(&gps_msg);
				old_sent_time = time_usec;
			}
		}
		/*
        if (ptask_deadline_miss())
        {
			printf("Simulator Thread Deadline Miss!  Abs deadline = %lu || Time = %lu\n", 
                   task_absdl(tid), ptask_gettime(MILLI));
        }
        */
        ptask_wait_for_period();
		
	}

}



// ----------------------------------------------------------------------
//    GROUND STATION THREAD
// ----------------------------------------------------------------------
/*
 * This thread takes the information from the simulator and ground 
 * station and directs them towards the PX4. 
 *
 *
 *    GroundStation >--------------> PX4  
 *                         
 * 
 *
 *
 */
void gs_thread()
{
	printf("***  Starting Ground Station Thread  ***\n");
	const unsigned int DIM_BUFF = 256;

	int tid = ptask_get_index();
	struct Interfaces* p = (struct Interfaces*)ptask_get_argument();

	mavlink_message_t msg_message;
	mavlink_status_t status;

	int chan_sel = -1;

	int i;
	int nread = 0;
    
    int first = 1;
    
	gs_thread_active = true;

	// Check the initialization of the necessary classes
	while (!time_to_exit)
	{
        if (first)
		{
			first = false;
			printf("GS Thread STARTED! \n");
		}
		// Send all the pending data to the Ground Station
		p->gs->sendMessages();

		//Wait for data from the Ground Station 
		p->gs->receiveMessages();
		
		// Send data to the Autopilot
		// We send data until the recQueue is empty
		while (p->gs->getDimrecQueue() > 0)
		{
        	// Retrieve message from the Ground Station
			p->gs->getMessage(&msg_message);

			p->aut->send_message(&msg_message);
        }
        // Record Sending Time
        gs_time = ptask_gettime(MICRO); 
		//fprintf(file_TGS,"%lu \n",gs_time);
        
        ptask_wait_for_period();
	}

}


// ----------------------------------------------------------------------
//    UNREAL ENGINE THREAD
// ----------------------------------------------------------------------
/*
 * This thread sends data to the Unreal Engine for visualization purposes 
 *
 *
 *    Pose of the Drone >--------------> UE 
 *                         
 * 
 *
 *
 */
void ue_thread()
{
	printf("***  Starting UE Communicator Thread  ***\n");

	struct UE_SendData dataOut;
	struct UE_RecData dataIn;
	
	int tid = ptask_get_index();
	struct Interfaces* p = (struct Interfaces*)ptask_get_argument();
		
	dataOut.Id = 1;
    
    int first = 1;
    
	ue_thread_active = true;

	// Check the initialization of the necessary classes
	while (!time_to_exit)
	{
        if (first)
		{
			first = false;
			printf("UE Thread STARTED! \n");
		}

		dataOut.X = (float)DynModel_Y.Xe[0];
		dataOut.Y = (float)DynModel_Y.Xe[1];
		dataOut.Z = (float)DynModel_Y.Xe[2];
	
		dataOut.r = (float)DynModel_Y.RPY[0];
		dataOut.p = (float)DynModel_Y.RPY[1];
 		dataOut.y = (float)DynModel_Y.RPY[2];

		p->ue->setData(dataOut);
		// Send all the pending data to Unreal Engine
		p->ue->sendData();
	
		
		p->ue->receiveData();
		p->ue->getData(&dataIn);
/*
			
		if (dataIn.PenDepth != 0.0)	
		{
			printf("Nx = %1.2f | Ny = %1.2f | Nz = %1.2f\n", dataIn.Nx, dataIn.Ny, dataIn.Nz);
			printf("Penetration = %2.2f\n\n", dataIn.PenDepth);
			printf("Reaction = [%2.2f | %2.2f | %2.2f]", DynModel_Y.Freact[0], DynModel_Y.Freact[1], DynModel_Y.Freact[2]);
		}
*/
		
		//DynModel_U.n_collision[0] = dataIn.Nx;
		//DynModel_U.n_collision[1] = -dataIn.Ny;
		//DynModel_U.n_collision[2] = -dataIn.Nz;
		//DynModel_U.pen_collision = dataIn.PenDepth;	

		
/*
		DynModel_U.n_collision[0] = 0.0;
		DynModel_U.n_collision[1] = 0.0;	
		DynModel_U.n_collision[2] = 0.0;
		DynModel_U.pen_collision = 0.0;
*/
        ptask_wait_for_period();
	}

}




// ----------------------------------------------------------------------
//   Parse Command Line
// ----------------------------------------------------------------------
// This function parse the option used to launch the application
// throws EXIT_FAILURE if could not open the port
//
// Inputs
// - gs_ip:	Ground Station IP <"x.x.x.x">
// - gs_r_port: Port number for reading incoming data from GS <"x">
// - gs_w_port: Port number for writing outgoing data to GS <"x">
void parse_commandline(int argc, char **argv, char *&gs_ip, unsigned int &gs_r_port, unsigned int &gs_w_port)
{

	// string for command line usage
	const char *commandline_usage = "usage: routing -gs_ip <GSip> -gs_rd <GSreadPort> - gs_wr <GSwritePort>";

	// Read input arguments
	for (int i = 1; i < argc; i++) { // argv[0] is "mavlink"

		// Help
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("%s\n",commandline_usage);
			throw EXIT_FAILURE;
		}

		// Ground Station Instance IP
		if (strcmp(argv[i], "-gs_ip") == 0) {
			if (argc > i + 1) {
				gs_ip = argv[i + 1];
			}
			else {
				printf("%s\n",commandline_usage);
				throw EXIT_FAILURE;
			}
		}

		// From Ground Station Read Port
		if (strcmp(argv[i], "-gs_rp") == 0) {
			if (argc > i + 1) {
				gs_r_port = (unsigned int) atoi(argv[i + 1]);
			}
			else {
				printf("%s\n",commandline_usage);
				throw EXIT_FAILURE;
			}
		}

		// To Ground Station Write Port
		if (strcmp(argv[i], "-gs_wp") == 0) {
			if (argc > i + 1) {
				gs_w_port = (unsigned int) atoi(argv[i + 1]);
			}
			else {
				printf("%s\n",commandline_usage);
				throw EXIT_FAILURE;
			}
		}

	}
	// end: for each input argument

	// Done!
	return;
}


// ----------------------------------------------------------------------
//   Quit Signal Handler
// ----------------------------------------------------------------------
// this function is called when you press Ctrl-C
	void
quit_handler( int sig )
{
	printf("\n");
	printf("TERMINATING AT USER REQUEST\n");
	printf("\n");

	time_to_exit = 1;

	try {
		// Check if it is in HIL mode

		if(autopilot_interface_quit->is_hil())
		{
			printf("The AUV is in hil mode: Try disabling...\n");
			// Disable HIL mode
			//autopilot_interface_quit->stop_hil();
		}

		//Close the Serial Port
		autopilot_interface_quit->uart_port.handle_quit(sig);

		/*
		   tspec wcet, acet;
		   int instances;
		   printf("Now printing the stats\n");
		   wcet = ptask_get_wcet(inflowT_id);
		   acet = ptask_get_avg(inflowT_id);
		   instances = ptask_get_numinstances(inflowT_id); 
		   printf("TASK INFLOW: NINST = %d  WCET = %ld\t ACET = %ld\t \n",instances,tspec_to(&wcet, MICRO), tspec_to(&acet, MICRO));

		   wcet = ptask_get_wcet(simT_id);
		   acet = ptask_get_avg(simT_id);
		   instances = ptask_get_numinstances(simT_id); 
		   printf("TASK SIMULATOR: NINST = %d WCET = %ld\t ACET = %ld\t \n",instances,tspec_to(&wcet, MICRO),tspec_to(&acet, MICRO));
		   */ 


		/*
		printf("Closing Files...\n\n");

		printf("close(file_TSndSns\n");
		fclose(file_TSndSns);
		printf("file_TSndComm\n");
		fclose(file_TSndComm);
		printf("file_TFtcSns\n");
		fclose(file_TFtcSns);
		printf("file_TFtcComm\n");
		fclose(file_TFtcComm);

		printf("file_TGS\n");
		fclose(file_TGS);
		*/
	} 
	catch (int error){}


	// end program here
	exit(0);

}
