#!/bin/bash
#
gdb --args ./main_routing -d /dev/ttyUSB0 -b 921600 -sim_ip 10.30.3.136 -sim_rp 49000 -sim_wr 49001 -gs_ip 127.0.0.1 -gs_rd 14551 -gs_wr 14550  
