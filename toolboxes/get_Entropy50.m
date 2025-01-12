function [Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100,N_dotbkg_50,N_dotret_50]=get_Entropy50
load ../Data_open/Entropy_data.mat

p_thre=0.004;
Entropy_thre_tmp_50(:,:,1:2)=(P_CHIV_50(:,:,1:2)<=p_thre).*(P_CHIV_dotsh_50(:,:,1:2)<=p_thre).*Entropy_50(:,:,1:2);
Entropy_thre_tmp_50(:,:,3)=(P_CHIV_50(:,:,3)<=p_thre).*Entropy_50(:,:,3);
for i=1:size(Entropy_thre_tmp_50,1)
Entropy_thre_50(i,:,1)=Suc5(Entropy_thre_tmp_50(i,:,1));
Entropy_thre_50(i,:,2)=Suc5(Entropy_thre_tmp_50(i,:,2));
Entropy_thre_50(i,:,3)=Suc5(Entropy_thre_tmp_50(i,:,3));
end
N_info_50=squeeze(sum(Entropy_thre_50,2)>0);
N_dotbkg_50=find(N_info_50(:,1));
N_dotret_50=find(N_info_50(:,2));


Entropy_thre_tmp_100(:,:,1:2)=(P_CHIV_100(:,:,1:2)<=p_thre).*(P_CHIV_dotsh_100(:,:,1:2)<=p_thre).*Entropy_100(:,:,1:2);
Entropy_thre_tmp_100(:,:,3)=(P_CHIV_100(:,:,3)<=p_thre).*Entropy_100(:,:,3);
for i=1:size(Entropy_thre_tmp_100,1)
Entropy_thre_100(i,:,1)=Suc5(Entropy_thre_tmp_100(i,:,1));
Entropy_thre_100(i,:,2)=Suc5(Entropy_thre_tmp_100(i,:,2));
Entropy_thre_100(i,:,3)=Suc5(Entropy_thre_tmp_100(i,:,3));
end
N_info_100=squeeze(sum(Entropy_thre_100,2)>0);
N_dotbkg_100=find((N_info_100(:,1)-N_info_50(:,1))==1);
N_dotret_100=find((N_info_100(:,2)-N_info_50(:,2))==1);
N_bkgret_100=find((N_info_100(:,3)-N_info_50(:,3))==1);
% N_dotbkg_100=false(size(N_dotbkg_100));
% N_dotret_100=false(size(N_dotbkg_100));
% N_bkgret_100=false(size(N_dotbkg_100));

Entropy=Entropy_50;
P_CHIV=P_CHIV_50;
P_CHIV_dotsh=P_CHIV_dotsh_50;
I_SHUFF=I_SHUFF_50;
Entropy(N_dotbkg_100,:,1)=Entropy_100(N_dotbkg_100,:,1);
P_CHIV(N_dotbkg_100,:,1)=P_CHIV_100(N_dotbkg_100,:,1);
P_CHIV_dotsh(N_dotbkg_100,:,1)=P_CHIV_dotsh_100(N_dotbkg_100,:,1);
I_SHUFF(N_dotbkg_100,:,:,1)=I_SHUFF_100(N_dotbkg_100,:,:,1);


Entropy(N_dotret_100,:,2)=Entropy_100(N_dotret_100,:,2);
P_CHIV(N_dotret_100,:,2)=P_CHIV_100(N_dotret_100,:,2);
P_CHIV_dotsh(N_dotret_100,:,2)=P_CHIV_dotsh_100(N_dotret_100,:,2);
I_SHUFF(N_dotret_100,:,:,2)=I_SHUFF_100(N_dotret_100,:,:,2);


Entropy(N_bkgret_100,:,3)=Entropy_100(N_bkgret_100,:,3);
P_CHIV(N_bkgret_100,:,3)=P_CHIV_100(N_bkgret_100,:,3);
%%
Spike_dot_bkg_filter=Spike_dot_bkg_filter_50;
Spike_dot_bkg_filter(N_dotbkg_100)=Spike_dot_bkg_filter_100(N_dotbkg_100);
Spike_dot_ret_filter=Spike_dot_ret_filter_50;
Spike_dot_ret_filter(N_dotret_100)=Spike_dot_ret_filter_100(N_dotret_100);
spike_peridot_full=spike_peridot_full_50;
spike_peridot_full(N_dotret_100)=spike_peridot_full_100(N_dotret_100);
bkgpos_map=bkgpos_map_50;
bkgpos_map(N_bkgret_100)=bkgpos_map_100(N_bkgret_100);


N_dotret_50=find((N_info_50(:,2)-N_info_100(:,2))==1);
N_dotret_50=find((N_info_50(:,2)))';

N_dotbkg_50=find(N_info_50(:,1))';
find((N_info_50(:,1)-N_info_100(:,1))==1);

