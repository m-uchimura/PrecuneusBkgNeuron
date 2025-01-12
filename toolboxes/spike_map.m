function [Spike_raw_full]=spike_map(dot_x,dot_y,Spike,trial_spike,Photo,bin_number,params)
dot_start=params.dot_start;
dot_end=params.dot_end;
calc_time_start=params.calc_time_start;
calc_time_end=params.calc_time_end;
bin_step=params.bin_step;
stim_center=params.stim_center;

stim_dist=[81,61];
Spike_map_temp= zeros(stim_dist(1),stim_dist(2),'uint8');
trial_spike_sim=1:length(trial_spike);
Spike_raw_full=zeros(bin_number,max(trial_spike_sim),stim_dist(1),stim_dist(2),'uint8');
for bin=1:bin_number
    Spike_raw_full_temp= zeros(max(trial_spike_sim),stim_dist(1),stim_dist(2),'uint8');
    for i=trial_spike_sim
        spike_trial_bin_shifta=Spike{i};
        Spike_map=Spike_map_temp;
        stim_on_time_bin_shift=squeeze(Photo{i});
        for k=dot_start:dot_end            
            dot_pos_x=stim_center(1)+dot_x(i,k);
            dot_pos_y=stim_center(2)+dot_y(i,k);
            spike_count_stim(i,k)=numel(find((spike_trial_bin_shifta>stim_on_time_bin_shift(k)+calc_time_start+(bin-1)*bin_step)&(spike_trial_bin_shifta<stim_on_time_bin_shift(k)+calc_time_end+(bin-1)*bin_step)));
            if ((dot_pos_x>0)&&(dot_pos_x<=stim_dist(1)))
                if ((dot_pos_y>0)&&(dot_pos_y<=stim_dist(2)))
                    Spike_map(dot_pos_x,dot_pos_y)= spike_count_stim(i,k);
                end
            end
        end
        Spike_raw_full_temp(i,:,:)=Spike_map;
    end
    %raw spike map 81*61
    Spike_raw_full(bin,:,:,:)=Spike_raw_full_temp;   
end


