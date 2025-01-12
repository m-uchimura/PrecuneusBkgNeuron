%%
%this file draw the resuls of gain field neurons for dot/ret or dot/bkg
%(Figure 10)
%%
clear 
close all
addpath(genpath('../toolboxes'))
%%
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
figure('position',[0 0 floor(100000*asp_rat_x/res_x) floor(100000*asp_rat_y/res_y)])
time_end=30;
p_thre=0.004;
p_thre_gain=0.00125;
fs=12;
%%
%load neual information, Pvalue etc
load ../Data_open/data.mat
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100]=get_Entropy50;
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
%%
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);

%%
I_SHUFF(~N_bkg100,:,:,1)=20*I_SHUFF(~N_bkg100,:,:,1);
I_SHUFF(N_bkg100,:,:,1)=10*I_SHUFF(N_bkg100,:,:,1);
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
%set parameters
iter_r=0;
iter=0;
fig=1;
rng(iter_r);
params.calc_time_start=-50/1000;
params.calc_time_end=0;
params.bin_step=10/1000;
params.dot_start=1;
params.dot_end=12;
bin_number=36;
F_real=zeros(10,10);
for f_i=1:10
    for f_j=1:10
        if norm([f_i,f_j]-[5,5])<5
            F_real(f_i,f_j)=1;
        end
    end
end
params.F_real=F_real;
map_x=81; map_y=61;
stim_dist=[map_x,map_y];
stim_center=(stim_dist+1)/2;
params.stim_center=stim_center;

%%
%draw receptive fields of gain field neurons for dot/bkg(Figure 10A) and dot/ret (Figure 10B) 
use_co=1;
for N_i=1:2
 if N_i==1,Ni='m4c109r1s1n1'; Bin1=9; coordinate=1;flip=0;
 else Ni='m4c137r1s1n6';Bin1=11;coordinate=2;flip=0;
 end

for i=1:length(data_bkg1)
    if strcmp(data_bkg1(i).file,Ni)
        iCell=i;
    end
end
dat=data_bkg1(iCell);
sig_period=Entropy_thre(iCell,:,use_co)>0;
sig_b=Entropy_thre(iCell,:,1)>0;
sig_r=Entropy_thre(iCell,:,2)>0;

Entropy_cell=squeeze(Entropy(iCell,:,:));
P_CHIV_cell=squeeze(P_CHIV(iCell,:,:));
P_CHIV_dotsh_cell=squeeze(P_CHIV_dotsh(iCell,:,:));
I_SHUFF_cell=squeeze(I_SHUFF(iCell,:,:,:));
RF_dot_bkg=Spike_dot_bkg_filter{iCell};
RF_dot_retino=Spike_dot_ret_filter{iCell};
spikes_bkg_ret=spike_peridot_full{iCell};
stim_bkg_ret=bkgpos_map{iCell};


dot_x=round(squeeze(dat.stim(:,:,1)));
dot_y=round(squeeze(dat.stim(:,:,2)));
bkg_x=dat.frame(:,1);
bkg_y=dat.frame(:,2);
fix_x=dat.fixation(:,1);
fix_y=dat.fixation(:,2);
trial_number=length(bkg_x);
%
clear trial_spike ts_cell
t=0;
for i=1:trial_number
   if sum(( dat.spike_ts>dat.photo_ts{i}(1)-1)&( dat.spike_ts<dat.photo_ts{i}(1)+2.5))
         if ~isnan(fix_x(i))
           t=t+1;
           trial_spike(t)=i;
           ts_cell{i}=dat.spike_ts((dat.spike_ts>dat.photo_ts{i}(1)-1)&( dat.spike_ts<dat.photo_ts{i}(1)+2));
         end
   end
end
Spikes= ts_cell(trial_spike);
Photo=dat.photo_ts(trial_spike);
%
gain_x=fix_x>=0;
gain_y=fix_y>=0;
gain=1+gain_x+gain_y*2;%1:3quadrant  2:4quadrant  3:2quadrant   4:1quadrant 
g=gain(trial_spike);

%count the number of dot appearance
[dot_bkg_L_full,dot_ret_L_full,~]=...
    stim_map(dot_x(trial_spike,:),dot_y(trial_spike,:),bkg_x(trial_spike),bkg_y(trial_spike),fix_x(trial_spike),fix_y(trial_spike),1:length(trial_spike),params);
if coordinate==1 %dot/bkg
dot_L_full=dot_bkg_L_full(:,:,:);
dot_x_frame=dot_x(trial_spike,:)-repmat(bkg_x(trial_spike),1,12);
dot_y_frame=dot_y(trial_spike,:)-repmat(bkg_y(trial_spike),1,12);
[ Spike_raw_L_full]=...
    spike_map(dot_x_frame,dot_y_frame,Spikes,trial_spike,Photo,bin_number,params);
elseif coordinate==2  %dot/ret
dot_L_full=dot_ret_L_full(:,:,:);
dot_x_fix=dot_x(trial_spike,:)-repmat(fix_x(trial_spike),1,12);
dot_y_fix=dot_y(trial_spike,:)-repmat(fix_y(trial_spike),1,12);
[ Spike_raw_L_full]=...
    spike_map(dot_x_fix,dot_y_fix,Spikes,trial_spike,Photo,bin_number,params);
end
%
dot_L=squeeze(sum(dot_L_full));
dot_L_fil=double(filter2(F_real,dot_L));
dot=dot_L(stim_center(1)-30:stim_center(1)+30,stim_center(2)-20:stim_center(2)+20);
Spike_raw_L=squeeze(sum(Spike_raw_L_full,2));
for bin=1:(bin_number)
    Spike_raw_L_bin=squeeze(Spike_raw_L(bin,:,:));
    [Spike_raw_L_filter{bin},Spike_raw{bin},Spike_raw_filter{bin},Spike_dot_filter{bin}]=...
        filter_edge(Spike_raw_L_bin,dot_L_fil,params) ;
end
Spike_raw_orig=Spike_raw;
dot_orig=dot;
[Spike_dot_filter_ALL,I,ChiP,ChiV,GchiV_mn,Gp_mn,I_bkg_gain_mn]=calc_gainI(Spike_raw_L_full,dot_L_full,g,F_real,stim_center,params,dot_orig,Spike_raw_orig);
for gain_i=1:4
max_sp1(gain_i)=nanmax(nanmax(Spike_dot_filter_ALL{gain_i}{Bin1}));
min_sp1(gain_i)=nanmin(nanmin(Spike_dot_filter_ALL{gain_i}{Bin1}));
end
Spike_dot_filter_ALL2{N_i}=Spike_dot_filter_ALL;

color_max_bin_cand=nanmax(nanmax(RF_dot_bkg{Bin1}));
color_min_bin_cand=nanmin(nanmin(RF_dot_bkg{Bin1}));
color_min_bin_cand= min(min_sp1);
color_max_bin_cand=max(max_sp1);

if N_i==1
%draw receptive fields by using all trials (Figure 10 A,B up)
subplot('position',[0.16 0.75+0.04 0.12 0.09])
RF_dot_bkg=Spike_dot_bkg_filter{iCell};
show_RF(RF_dot_bkg{Bin1}',color_min_bin_cand,color_max_bin_cand,[],Entropy_cell(Bin1,1),'bkg',fs,flip)
set(gca,'XTickLabel',[])
subplot('position',[0.30 0.75+0.04 0.015 0.09])
axis off
colbar_leg=[0 10 20];
col_diff=10*(color_max_bin_cand-color_min_bin_cand);
col_min=color_min_bin_cand*10;
col_max=color_max_bin_cand*10;
colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],...
    'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',12,'FontName','Arial narrow')

GchiV=nanmean(nanmean(GchiV_mn,4),3);
Gp=chi2cdf(GchiV,4*3-1,'upper');
I_bkg_gain=nanmean(nanmean(I_bkg_gain_mn,4),3);
%draw receptive fields when fixation point was presented in each quadrant
for gain_i=1:4
 switch gain_i%1:3quadrant  2:4quadrant  3:2quadrant   4:1quadrant 
     case 3, subplot('position',[0.1 0.66 0.12 0.09])   %2quadrant
     case 4, subplot('position',[0.235 0.66 0.12 0.09]) %1quadrant
     case 1, subplot('position',[0.1 0.53 0.12 0.09])   %3quadrant
     case 2, subplot('position',[0.235 0.53 0.12 0.09]) %4quadrant
 end
 max_sp1(isinf(max_sp1))=10;
 show_RF(Spike_dot_filter_ALL2{1}{gain_i}{Bin1}',color_min_bin_cand,color_max_bin_cand,Gp(gain_i,Bin1),GchiV(gain_i,Bin1),'bkg',fs,flip)
 set(gca,'XTickLabel',[])
end
else
subplot('position',[0.16+0.3 0.75+0.04 0.12 0.09]) 
RF_dot_ret=Spike_dot_ret_filter{iCell};
show_RF(RF_dot_ret{Bin1}',color_min_bin_cand,color_max_bin_cand,[],Entropy_cell(Bin1,2),'retino',fs,flip)
set(gca,'XTickLabel',[])
subplot('position',[0.305+0.3 0.75+0.04 0.015 0.09])
axis off

colbar_leg=[5 10 15];
col_diff=10*(color_max_bin_cand-color_min_bin_cand);
col_min=color_min_bin_cand*10;
col_max=color_max_bin_cand*10;
colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],...
    'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',12,'FontName','Arial narrow')
GchiV=nanmean(nanmean(GchiV_mn,4),3);
Gp=chi2cdf(GchiV,4*3-1,'upper');
for gain_i=1:4
 switch gain_i
     case 3, subplot('position',[0.1+0.3 0.66 0.12 0.09])
     case 4, subplot('position',[0.235+0.3 0.66 0.12 0.09])
     case 1, subplot('position',[0.1+0.3 0.53 0.12 0.09])
     case 2, subplot('position',[0.235+0.3 0.53 0.12 0.09])
 end
show_RF(Spike_dot_filter_ALL{gain_i}{Bin1}',color_min_bin_cand,color_max_bin_cand,Gp(gain_i,Bin1),GchiV(gain_i,Bin1),'retino',fs,flip)
set(gca,'XTickLabel',[])
end
end
GchiV(:,Bin1)
Gp(:,Bin1)
end

%%
%draw the time course of neural information dot/bkg neurons with gain
%modulation(Up) or without gain modulation (Bottom) Figure 10C
use_neurons=sum(Entropy_thre(:,:,1),2)>0;
Entropy_thre_tmp=Entropy_thre(:,:,1);
subplot('position',[0.1 0.35 0.25 0.08])

load ../Data_open/gain_bkg.mat GchiV

Entropy_thre_use=Entropy_thre_tmp(use_neurons,:,use_co);
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
set(gca,'ytick',[0 1 2],'yticklabel',{'0','1','2'},'TickDir','out','box','off')
ylim([0 2])
subplot('position',[0.1 0.23 0.25 0.08])

T_nogain=Entropy_thre_use(~Gain_Neuron_bkg,:);
N_nogain_bkg=sum(~Gain_Neuron_bkg);
T_nogain=T_nogain./N_nogain_bkg;
show_SumEntropy(T_nogain,amp,time_end,fs,'m')
plot([1,1],[0 1],'k:','Linewidth',1)
plot([6,6],[0 1],'k--','Linewidth',1)
set(gca,'ytick',[0 0.2 0.4],'yticklabel',{'0','0.2','0.4'},'TickDir','out','box','off')
ylim([0 0.5])

%%
%draw the time course of neural information dot/ret neurons with gain
%modulation(Up) or without gain modulation (Bottom) Figure 10D
use_neurons=sum(Entropy_thre(:,:,2),2)>0;%dot/ret
subplot('position',[0.1+0.3 0.35 0.25 0.08])
load ../Data_open/gain_retino.mat GchiV

Entropy_thre_use=Entropy_thre(use_neurons,:,2);
Entropy_thre_use_nan=Entropy_thre_use>0;
for i=1:4
GchiV(:,i,:)=squeeze(GchiV(:,i,:)).*(Entropy_thre_use>0);
end

Gp=chi2cdf(GchiV,4*3-1,'upper');
[~,maxim]=max(Entropy_thre_use,[],2);
Gp_under_thre=(Gp<p_thre_gain);
Gain_Neuron_ret=sum(sum(Gp_under_thre,2),3)>0;
G_Neuron_r=repmat(Gain_Neuron_ret,1,36);
T_gain=Entropy_thre_use.*G_Neuron_r;
N_gain_ret=sum(sum(G_Neuron_r,2)>0);
T_gain=T_gain/N_gain_ret;
show_SumEntropy(T_gain,amp,time_end,fs,'b')
plot([1,1],[0 8],'k:','Linewidth',1)
plot([6,6],[0 8],'k--','Linewidth',1)
set(gca,'ytick',[0 3 6],'yticklabel',{'0','3','6'},'TickDir','out','box','off')
ylim([0 7.5])

subplot('position',[0.1+0.3 0.23 0.25 0.08])
T_nogain=Entropy_thre_use.*~G_Neuron_r;
T_all_nogain=Entropy_thre(use_neurons,:,2).*~G_Neuron_r;
N_nogain_ret=sum(sum(T_all_nogain,2)>0);
T_nogain=T_nogain./N_nogain_ret;
show_SumEntropy(T_nogain,amp,time_end,fs,'b')
plot([1,1],[0 5],'k:','Linewidth',1)
plot([6,6],[0 5],'k--','Linewidth',1)
set(gca,'ytick',[0 2 4],'yticklabel',{'0','2','4'},'TickDir','out','box','off')
ylim([0 4])


%%
if 0
cd figs%/material
print('fig10_gain.ai', '-dpdf', '-painters','-bestfit')
cd ../
end