
%this code draw figure 13 C
%%
clear 
close all
addpath(genpath('../toolboxes'))
%%
save1=0;
lr=3;
p_thre=0.004;%p_thre=0.049;
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;

%%
load ../Data_open/data.mat
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100]=get_Entropy50;
%%
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);

%%
Entropy_thre_tmp(:,:,1:2)=(P_CHIV(:,:,1:2)<=p_thre).*(P_CHIV_dotsh(:,:,1:2)<=p_thre).*Entropy(:,:,1:2);
Entropy_thre_tmp(:,:,3)=(P_CHIV(:,:,3)<=p_thre).*Entropy(:,:,3);
for i=1:size(Entropy_thre_tmp,1)
Entropy_thre(i,:,1)=Suc5(Entropy_thre_tmp(i,:,1));
Entropy_thre(i,:,2)=Suc5(Entropy_thre_tmp(i,:,2));
Entropy_thre(i,:,3)=Suc5(Entropy_thre_tmp(i,:,3));
end
ent_shuff_mean=squeeze(nanmean(I_SHUFF,2));
entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0).*(Entropy_thre(:,:,2)>0);
entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
Entropy_thre(:,:,1)=entropy_thre_minus_spread;
Neuron_information=squeeze(sum(Entropy_thre,2)>0);
%%
load ../Data_open/cell_POS.mat
cell_xy=cell_xy_all;
cell_center=cell_POS_all;
load ../Data_open/anatomy Anat
lateral=Anat.X+1; 
dep=Anat.dep;

%%
xy=cell_xy;
pos=cell_center;
edge=[-Inf,-1000:250:250,Inf];
y=cell_xy(:,1)-cell_center(:,1);
for x=1:5
flt=(lateral==x);
hist_all(x,:)=histcounts(y(flt),edge);
end

%count and draw the histograme of anatomical positons of dot/ret neurons (Figure 8C left)
figure('position',[100 100 floor(1.25*16000*asp_rat_x/res_x) floor(1.25*16000*asp_rat_y/res_y)])
for x=1:5
flt=(lateral'==x)&Neuron_information(:,2)==1;
hist_dotRetino(x,:)=histcounts(y(flt),edge);
end
imagesc(((hist_dotRetino(:,:)./hist_all(:,:))),[0 0.5])
set(gca,'xtick',[],'ytick',[],'box','of')

if save1
print('hist_dot_ret.ai', '-dpdf', '-painters')
end

%count and draw the histograme of anatomical positons of blg/ret neurons (Figure 8C left)
figure('position',[100 100 floor(1.25*16000*asp_rat_x/res_x) floor(1.25*16000*asp_rat_y/res_y)])
for x=1:5
flt=(lateral'==x)&Neuron_information(:,3)==1;
hist_bkgRetino(x,:)=histcounts(y(flt),edge);
end
imagesc(((hist_bkgRetino(:,:)./hist_all(:,:))),[0 1])
set(gca,'xtick',[],'ytick',[],'box','of')

if save1
print('hist_bkg_ret.ai', '-dpdf', '-painters')
end

%count and draw the histograme of anatomical positons of dot/bkg neurons (Figure 8C left)
figure('position',[100 100 floor(1.25*16000*asp_rat_x/res_x) floor(1.25*16000*asp_rat_y/res_y)])
for x=1:5
flt=(lateral'==x)&Neuron_information(:,1)==1;
hist_dotBkg(x,:)=histcounts(y(flt),edge);
end
imagesc(((hist_dotBkg(:,:)./hist_all(:,:))),[0 0.25])
set(gca,'xtick',[],'ytick',[],'box','of')

if save1
print('hist_dot_bkg.ai', '-dpdf', '-painters')
end

