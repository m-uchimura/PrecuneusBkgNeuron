
function [spike_peridot,spike_peribkg, bkgpos_map]=spike_map_bkg(bkg_pos_x,bkg_pos_y,spike,trial_spike, photo,bin_frame,bin_dot,params)
spike_peridot=zeros(12,bin_dot,51,41);
dot_start=params.dot_start;
dot_end=params.dot_end;
calc_time_start=params.calc_time_start;
calc_time_end=params.calc_time_end;
spike_peribkg = zeros(length(trial_spike),bin_frame,51,41); %trial*bin_frame*51*41
bkgpos_map=zeros(51,41);

for i=1:length(trial_spike)
    spike_trial=spike{i};  
    stim_on_time_dot=photo{i};
    bkgpos_map((26+bkg_pos_x(i)),21+bkg_pos_y(i))=bkgpos_map((26+bkg_pos_x(i)),21+bkg_pos_y(i))+1;
    bkg_on=0.8+stim_on_time_dot(1)-1;
    for bin_fix=1:bin_frame
        bin_spikes=sum(((bkg_on-0.1)+params.bin_step*(bin_fix-1))<spike_trial&spike_trial<=(bkg_on+params.bin_step*(bin_fix-1)));
        spike_peribkg(i,bin_fix,round(26+bkg_pos_x(i)),21+bkg_pos_y(i))=bin_spikes;
    end
    
    spike_peridot_temp=zeros(12,36,51,41);
    for k=dot_start:dot_end
        for fra_bin=1:bin_dot
            bin_spikes_dot=sum((stim_on_time_dot(k)+params.bin_step*(fra_bin-1)+calc_time_start)<spike_trial&spike_trial<=(stim_on_time_dot(k)+params.bin_step*(fra_bin-1)+calc_time_end));
            spike_peridot_temp(k,fra_bin,round(26+bkg_pos_x(i)),21+bkg_pos_y(i))=bin_spikes_dot;
        end
    end   
    spike_peridot=spike_peridot+spike_peridot_temp;
end