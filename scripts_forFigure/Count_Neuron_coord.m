%%
%this code calculate the numbers of neurons of each information types in hist_coordinate
%this is used to visualze Venn Diagram in Figure 13
%%
clear 
close all
addpath('../toolboxes/')
%%
load ../Data_open/data.mat
p_thre=0.004;
mon=[data_bkg1(:).monkey];
monkey_no=0;%2,4,5 5
if monkey_no==0
monkey_flt=true(1,length(mon));
else
monkey_flt=mon==monkey_no;
end
%%
%load neual information, Pvalue etc
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100]=get_Entropy50;
%%
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);
%%
%Set info to 0 when the p-value was below the threshold 
Entropy_thre_tmp(:,:,1:2)=(P_CHIV(:,:,1:2)<=p_thre).*(P_CHIV_dotsh(:,:,1:2)<=p_thre).*Entropy(:,:,1:2);
Entropy_thre_tmp(:,:,3)=(P_CHIV(:,:,3)<=p_thre).*Entropy(:,:,3);
for i=1:size(Entropy_thre_tmp,1)
Entropy_thre(i,:,1)=Suc5(Entropy_thre_tmp(i,:,1));
Entropy_thre(i,:,2)=Suc5(Entropy_thre_tmp(i,:,2));
Entropy_thre(i,:,3)=Suc5(Entropy_thre_tmp(i,:,3));
end
%subtract spread of dot/ret info from dot/bkg info
ent_shuff_mean=squeeze(nanmean(I_SHUFF,2));
entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0).*(Entropy_thre(:,:,2)>0);
entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
Entropy_thre(:,:,1)=entropy_thre_minus_spread;

%%
cell_coordinate=squeeze(sum(Entropy_thre(monkey_flt,:,:),2)>0);
sum(cell_coordinate)
coord_2bit=cell_coordinate(:,1)+cell_coordinate(:,2)*2+cell_coordinate(:,3)*4;
hist_coordinate=histcounts(coord_2bit);
%0 no 1:dot/bkg 2:dot/ret 3:dot/bkg+dot/ret 4:bkg/ret 5:bkg/ret+dot/bkg
%6:dot/ret+bkg/ret 7:triple
totalN=sum(cell_coordinate);%1:total of dot/bkg %2 dot/ret %3 bkg/ret


