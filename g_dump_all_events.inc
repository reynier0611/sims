*
* variables that help dump all GENERATED events in the HUT ntuple 
*
      real*8 stop_where
      real*8 x_stop,y_stop
      real*8 good_evt
      logical dump_all_in_ntuple, dump_all_this_ev, temp_success

      common /g_dump_all/ stop_where,x_stop,y_stop,good_evt,
     > 	     dump_all_in_ntuple,dump_all_this_ev,temp_success


