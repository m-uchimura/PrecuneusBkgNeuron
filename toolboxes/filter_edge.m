%filtering
%edge cut
%filterd/sti 81*61
%raw spike map 61*41
%filtered spike map 61*41
%filterd/dot 61*41
function [Spike_raw_L_filter,Spike_raw,Spike_raw_filter,Spike_dot_filter]=filter_edge(Spike_raw_L_bin,dot_L_fil,global_params)
F_real=global_params.F_real;
stim_center=global_params.stim_center;
Spike_raw_filter =filter2(F_real,Spike_raw_L_bin);
Spike_dot_filter=Spike_raw_filter./dot_L_fil;
%filterd/sti 81*61
Spike_raw_L_filter=Spike_dot_filter;
%delete the edge
Spike_raw=Spike_raw_L_bin(stim_center(1)-30:stim_center(1)+30,stim_center(2)-20:stim_center(2)+20);
Spike_raw_filter=Spike_raw_filter(stim_center(1)-30:stim_center(1)+30,stim_center(2)-20:stim_center(2)+20);
Spike_dot_filter=Spike_dot_filter(stim_center(1)-30:stim_center(1)+30,stim_center(2)-20:stim_center(2)+20);





    
    

