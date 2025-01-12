function trial_spike=get_trial_spike(dat)
tk=0;
for i=1:length(dat.photo_ts)
    if sum((dat.spike_ts>dat.photo_ts{i}(1)-1)&(dat.spike_ts<dat.photo_ts{i}(1)+2.5))
        if ~isnan(dat.fixation(i,1))
            tk=tk+1;
            trial_spike(tk)=i;
        end
    end
end