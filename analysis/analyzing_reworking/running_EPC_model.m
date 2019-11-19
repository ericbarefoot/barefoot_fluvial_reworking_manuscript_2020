in_str.nits = 1000; 
in_str.mw = 100; 
in_str.avul = 2; 
in_str.ch_z = 1; 
in_str.ch_w_i = 11; 
in_str.IR = 0.75; 
in_str.randint = 0; 
in_str.vertagg_rate = 0.1; 
in_str.fp_lat = 5; 
in_str.walk_mag = 5; 
in_str.wwide = 10; 
in_str.wthick = 10; 
in_str.plot_bool = true; 
in_str.latmob = 1;
in_str.fp_filling = 'everywhere';

out_str = StickModel_BarPres(in_str);

close
plot_strat_str(out_str)
