function [I_bkgpos,ChiV_bkgpos]=calc_entropy_bkgposition(bkgpos_map,spike_peribkg,temp_for_bkgentropy,bin_n)

spikes_bkg_row=temp_for_bkgentropy;
spikes_bkg_row(:,1)= reshape(bkgpos_map,51*41,1);% %1:51,1;1:51;51
spikes_bkg_row( spikes_bkg_row(:,1)==0,:)=[];
spikes_bkg_row(:,2)=spikes_bkg_row(:,2)-26;
spikes_bkg_row(:,3)=spikes_bkg_row(:,3)-21;

spike_bkgpos_temp=NaN(bin_n,size(temp_for_bkgentropy,1),size(temp_for_bkgentropy,2));
I_bkgpos_rot=NaN(bin_n,9);
ChiV_bkgpos_rot=NaN(bin_n,9);

for cal_bin=1:bin_n 
    resp_frafix_bin= squeeze(spike_peribkg(cal_bin,:,:));
    resp_frafix_row=temp_for_bkgentropy;
    resp_frafix_row(:,1)= reshape(resp_frafix_bin,51*41,1);% %1:51,1;1:51;51
    spike_bkgpos_temp(cal_bin,:,:)=resp_frafix_row;
end

for t=0:8
    th=(t*10)*pi/180;
    R=[cos(th) -sin(th);sin(th) cos(th)];
    spike_bkgret_bkg_rot=spikes_bkg_row;
    spike_bkgret_bkg_rot(:,2:3)=spike_bkgret_bkg_rot(:,2:3)*R;
    center_x_bkgret=sum(spike_bkgret_bkg_rot(:,1).*spike_bkgret_bkg_rot(:,2))/sum(spike_bkgret_bkg_rot(:,1));
    center_y_bkgret=sum(spike_bkgret_bkg_rot(:,1).*spike_bkgret_bkg_rot(:,3))/sum(spike_bkgret_bkg_rot(:,1));

    for calc_bin=1:bin_n
        spike_bkg_tem=squeeze(spike_bkgpos_temp(calc_bin,:,:));
        spike_bkg_tem(spike_bkg_tem(:,1)==0,:)=[];
        spike_bkg_tem(:,2)=spike_bkg_tem(:,2)-26;
        spike_bkg_tem(:,3)=spike_bkg_tem(:,3)-21;
        spike_bkg_tem(:,2:3)=spike_bkg_tem(:,2:3)*R;
        [I_bkgpos_rot(calc_bin,t+1),ChiV_bkgpos_rot(calc_bin,t+1),~]=calc_chi2_2_rotate(spike_bkg_tem,spike_bkgret_bkg_rot,center_x_bkgret,center_y_bkgret);
    end
end
I_bkgpos=nanmean(I_bkgpos_rot,2);
ChiV_bkgpos=nanmean(ChiV_bkgpos_rot,2);
