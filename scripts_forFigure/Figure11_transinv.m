%%
%this file draw the resuls of translation invariance neurons for dot/bkg
%(Figure 11)
%%
clear 
colormap default
close all
addpath(genpath('../Scripts_forTransInv'))
addpath(genpath('../Scripts_forGain'))
addpath(genpath('../toolboxes'))

res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
time_end=30;
p_thre=0.004;
p_thre_gain=0.00125;
fs=12;
%%
%load neual information, Pvalue etc
load ../Data_open/data.mat
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100]=get_Entropy50;
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);

%%
% replace entropy if time window of 100ms make significant information, 
N_bkg100=false(size(Entropy,1),1);
N_bkg100(N_dotbkg_100)=true;
Entropy(~N_bkg100,:,1)=20*Entropy(~N_bkg100,:,1);
Entropy(N_bkg100,:,1)=10*Entropy(N_bkg100,:,1);
N_ret100=false(size(Entropy,1),1);
N_ret100(N_dotret_100)=true;
Entropy(~N_ret100,:,2)=20*Entropy(~N_ret100,:,2);
Entropy(N_ret100,:,2)=10*Entropy(N_ret100,:,2);

N_bkgret100=false(size(Entropy,1),1);
N_bkgret100(N_bkgret_100)=true;
Entropy(~N_bkgret100,:,3)=20*Entropy(~N_bkgret100,:,3);
Entropy(N_bkgret100,:,3)=10*Entropy(N_bkgret100,:,3);

I_SHUFF(~N_bkg100,:,:,1)=20*I_SHUFF(~N_bkg100,:,:,1);
I_SHUFF(N_bkg100,:,:,1)=10*I_SHUFF(N_bkg100,:,:,1);

%set entropy zerop if p-value was abobe the threshold 
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

%%
Ni='m4c186r1s1n1'; Bin1=15; flip=0;
for i=1:length(data_bkg1)
    if strcmp(data_bkg1(i).file,Ni)
        iCell=i;
    end
end
dat=data_bkg1(iCell);
Entropy_cell=squeeze(Entropy(iCell,:,:));
P_CHIV_cell=squeeze(P_CHIV(iCell,:,:));
P_CHIV_dotsh_cell=squeeze(P_CHIV_dotsh(iCell,:,:));
I_SHUFF_cell=squeeze(I_SHUFF(iCell,:,:,:));
RF_dot_bkg=Spike_dot_bkg_filter{iCell};
RF_dot_retino=Spike_dot_ret_filter{iCell};
spikes_bkg_ret=spike_peridot_full{iCell};
stim_bkg_ret=bkgpos_map{iCell};

load ../Data_open/trans_bkg.mat GchiV Spike_dot_filter_g
use_neurons=sum(Entropy_thre(:,:,1),2)>0;%dot/bkg

use_N=find(find(use_neurons)==iCell);
GchiV_cell=squeeze(GchiV(use_N,:,:));
Gp_cell=chi2cdf(GchiV_cell,4*3-1,'upper');
Spike_dot_filter_ALL=Spike_dot_filter_g{use_N};
for gain_i=1:4
max_sp1(gain_i)=nanmax(nanmax(Spike_dot_filter_ALL{gain_i,Bin1}));
min_sp1(gain_i)=nanmin(nanmin(Spike_dot_filter_ALL{gain_i,Bin1}));
end
color_min_bin_cand= min(min_sp1);
color_max_bin_cand=max(max_sp1);
%draw receptive fields by using all trials (Figure 11 A)
subplot('position',[0.16 0.75+0.1 0.12 0.09]) 
show_RF(RF_dot_bkg{Bin1}',color_min_bin_cand,color_max_bin_cand,[],Entropy_cell(Bin1,1),'bkg',fs,flip)
set(gca,'XTickLabel',[])

subplot('position',[0.30 0.75+0.1 0.015 0.09])
axis off
colbar_leg=[0 10 20];%i
col_diff=10*(color_max_bin_cand-color_min_bin_cand);
col_min=color_min_bin_cand*10;
col_max=color_max_bin_cand*10;
colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],...
    'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',12,'FontName','Arial narrow')

%draw receptive fields when background grame was presented in each quadrant
for gain_i=1:4
 switch gain_i%1:3quadrant  2:4quadrant  3:2quadrant   4:1quadrant
     case 3, subplot('position',[0.1 0.46 0.12 0.09])   %2quadrant
     case 4, subplot('position',[0.235 0.46 0.12 0.09]) %1quadrant
     case 1, subplot('position',[0.1 0.33 0.12 0.09])   %3quadrant
     case 2, subplot('position',[0.235 0.33 0.12 0.09]) %4quadrant
 end
 max_sp1(isinf(max_sp1))=10;
 show_RF(Spike_dot_filter_ALL{gain_i,Bin1}',color_min_bin_cand,color_max_bin_cand,Gp_cell(gain_i,Bin1),GchiV_cell(gain_i,Bin1),'bkg',fs,flip)
 set(gca,'XTickLabel',[])
end
%%
%draw the position of background presentation (Figure 11B)
load spine
for quad_i=1:4
switch quad_i
     case 2, ax1=subplot('position',[0.1 0.73-0.005 0.12 0.09]);   
     case 1, ax1=subplot('position',[0.235 0.73-0.005 0.12 0.09]); 
     case 3, ax1=subplot('position',[0.1 0.6-0.005 0.12 0.09]); 
     case 4, ax1=subplot('position',[0.235 0.6-0.005 0.12 0.09]); 
 end
conv_resp(quad_i,:,:)=draw_BkgRet_pos(dat,quad_i);
colormap(ax1,map)
end

%%
%draw the time course of neural information dot/bkg neurons with translation variant(Up) or invariant(Bottom) Figure 11D
use_neurons=sum(Entropy_thre(:,:,1),2)>0;
Entropy_thre_tmp=Entropy_thre(:,:,1);
subplot('position',[0.1 0.18 0.25 0.08])
Entropy_thre_use=Entropy_thre_tmp(use_neurons,:,1);
for i=1:4
GchiV(:,i,:)=squeeze(GchiV(:,i,:)).*(Entropy_thre_use>0);
end
Gp=chi2cdf(GchiV,4*3-1,'upper');
[~,maxim]=max(Entropy_thre_use,[],2);
N_bkg=Entropy_thre_use>0;

Gp_under_thre=(Gp<p_thre_gain);
Gain_Neuron_bkg=sum(sum(Gp_under_thre,2),3)>0;

for i=1:size(N_bkg,1)
Gp_under_thre_sig(i,:,:)=squeeze(Gp_under_thre(i,:,:)).*repmat(N_bkg(i,:),4,1);
end
Gp_thre_period=squeeze(sum(Gp_under_thre_sig,2))>0;

T_gain=Entropy_thre_use(Gain_Neuron_bkg,:);
N_gain_bkg=sum(Gain_Neuron_bkg);
T_gain=T_gain/N_gain_bkg;
amp=1;
show_SumEntropy(T_gain,amp,time_end,fs,'m')
T_gain_period=amp*Entropy_thre_use.*Gp_thre_period;

plot([1,1],[0 2],'k:','Linewidth',1)
plot([6,6],[0 2],'k--','Linewidth',1)
set(gca,'ytick',[0 .5 1],'yticklabel',{'0','0.5','1'},'TickDir','out','box','off')
ylim([0 1.5])
ylim([0 1.2])

subplot('position',[0.1 0.06 0.25 0.08])
T_nogain=Entropy_thre_use(~Gain_Neuron_bkg,:);
N_nogain_bkg=sum(~Gain_Neuron_bkg);
T_nogain=T_nogain./N_nogain_bkg;
show_SumEntropy(T_nogain,amp,time_end,fs,'m')
plot([1,1],[0 1],'k:','Linewidth',1)
plot([6,6],[0 1],'k--','Linewidth',1)
set(gca,'ytick',[0 0.2 0.4],'yticklabel',{'0','0.2','0.4'},'TickDir','out','box','off')
ylim([0 0.5])

%%
if 0
 cd figs
 print('fig11_trans.ai', '-dpdf', '-painters','-bestfit')
 cd ../
end
