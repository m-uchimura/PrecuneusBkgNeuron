function [Spike_dot_filter_ALL,I,ChiP,ChiV,GchiV_mn,Gp_mn,I_bkg_gain]=calc_gainI(Spike_raw_L_full,dot_L_full,g,F_real,stim_center,params,dot_bkg_orig,Spike_raw_bkg_orig)
for gain_i=1:4
A=Spike_raw_L_full(:,g==gain_i,:,:);
N_gain_bkg=dot_L_full(g==gain_i,:,:); 
%%
dot_L=squeeze(sum(N_gain_bkg));
dot_L_fil=double(filter2(F_real,dot_L));
%delete the edge
dot_bkg=dot_L(stim_center(1)-30:stim_center(1)+30,stim_center(2)-20:stim_center(2)+20);
Spike_raw_L=squeeze(sum(A,2));
%delete edge
% bin=14;
% gain_i=2;
for bin=1:36%(bin_number)
    Spike_raw_L_bin=squeeze(Spike_raw_L(bin,:,:));
    [Spike_raw_L_filter{bin},Spike_raw{bin},Spike_raw_filter{bin},Spike_dot_filter{bin}]=...
        filter_edge(Spike_raw_L_bin,dot_L_fil,params) ;
    %%%%%%%
    dot_bkg_for_entropy=repmat(dot_bkg,3,3);
    dot_bkg_orig_for_entropy=repmat(dot_bkg_orig,3,3);
    spike_bkg_for_entropy=repmat(Spike_raw{bin},3,3);
    spike_bkg_orig_for_entropy=repmat(Spike_raw_bkg_orig{bin},3,3);
    stim_center_L=floor(size(spike_bkg_for_entropy)/2)+1;
    %%%%%%%
   
    for m=1:5
        for n=1:5
            spike_bkg_for_entropy_mn=spike_bkg_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3));
            spike_bkg_orig_for_entropy_mn=spike_bkg_orig_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3));
            dot_bkg_for_entropy_mn=dot_bkg_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3));
            dot_bkg_orig_for_entropy_mn=dot_bkg_orig_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3));
            [GchiV_mn(gain_i,bin,m,n),Gp_mn(gain_i,bin,m,n)]=calc_chi2_gain(spike_bkg_orig_for_entropy_mn,dot_bkg_orig_for_entropy_mn,spike_bkg_for_entropy_mn,dot_bkg_for_entropy_mn);
        
           [I_bkg_gain(gain_i,bin,m,n),ChiV_bkg_gain(gain_i,bin,m,n),ChiP_bkg_gain(gain_i,bin,m,n)]=calc_chi2(spike_bkg_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3)),dot_bkg_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3)));
         end
    end
end

%calculate entropy
dot_bkg_for_entropy=repmat(dot_bkg,3,3);
for bin=1:length(Spike_raw)
    spike_bkg_for_entropy=repmat(Spike_raw{bin},3,3);
    stim_center_L=floor(size(spike_bkg_for_entropy)/2)+1;
    for m=1:5
        for n=1:5
       [I_bkg_full(bin,m,n),ChiV_bkg_full(bin,m,n),ChiP_bkg_full(bin,m,n)]=calc_chi2(spike_bkg_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3)),dot_bkg_for_entropy(stim_center_L(1)-30+5*(m-3):stim_center_L(1)+30+5*(m-3),stim_center_L(2)-20+5*(n-3):stim_center_L(2)+20+5*(n-3)));
        end
    end
end
for bin=1:length(Spike_raw)
    I(gain_i,bin)=[nanmean(nanmean(squeeze(I_bkg_full(bin,:,:)))) ];
    ChiP(gain_i,bin)=[nanmean(nanmean(squeeze(ChiP_bkg_full(bin,:,:)))) ];
    ChiV(gain_i,bin)=[nanmean(nanmean(squeeze(ChiV_bkg_full(bin,:,:))))];
end
Spike_dot_filter_ALL{gain_i}=Spike_dot_filter;
end