%%
%this file draw the resuls of dot/ret or dot/bkg neurons for Figure 4 or 5
%%
clear
close all
addpath(genpath('../toolboxes'))
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
figure('position',[0 0 floor(100000*asp_rat_x/res_x) floor(100000*asp_rat_y/res_y)])
%%
%choose coordinate to draw
coordinate='dot/retino'; %for figure 3
%coordinate='dot/bkg';   %for figure 4
%%
switch coordinate
    case 'dot/retino',coord=1;
    case 'dot/bkg',coord=2;
end
fs=8;
time_end=30;
p_thre=0.004;
t=6;
ti=0:20;
myu=10;
sd=5;
con=exp(-(ti-myu).^2/(2*sd^2))/(sqrt(2*pi)*sd);
con=con/sum(con);
frame_con=zeros(30,20);
frame_con(1,:)=1;
frame_con(30,:)=1;
frame_con(:,1)=1;
frame_con(:,20)=1;
frame_con=frame_con/sum(sum(frame_con));
switch coord
    case 1 
    iCell=788;data_bin=10;%file=m5c128r1s1n8' 
    case 2
    iCell=404;data_bin=15;%file='m4c186r1s1n1';
end
load ../Data_open/data.mat
dat=data_bkg1(iCell);
file=dat.file;
mon=str2double(file(2));
cite=str2double(file(4:6));
tk=0;
dat.spike_correct=[];
for i=1:length(dat.photo_ts)
    if sum((dat.spike_ts>dat.photo_ts{i}(1)-1)&(dat.spike_ts<dat.photo_ts{i}(1)+2.5))
        if ~isnan(dat.fixation(i,1))
            tk=tk+1;
            trial_spike(tk)=i;
            dat.spike_correct{i}=dat.spike_ts((dat.spike_ts>dat.photo_ts{i}(1)-1)&(dat.spike_ts<dat.photo_ts{i}(1)+2.5));
        end
    end
end

%%
%load neual information, Pvalue etc
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100]=get_Entropy50;
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
%%
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);
%%
Entropy_cell=squeeze(Entropy(iCell,:,:));
P_CHIV_cell=squeeze(P_CHIV(iCell,:,:));
I_SHUFF_cell=squeeze(I_SHUFF(iCell,:,:,:));
P_CHIV_dotsh_cell=squeeze(P_CHIV_dotsh(iCell,:,:));
RF_dot_bkg=Spike_dot_bkg_filter{iCell};
RF_dot_retino=Spike_dot_ret_filter{iCell};
spikes_bkg_ret=spike_peridot_full{iCell};
stim_bkg_ret=bkgpos_map{iCell};
%%
%draw raster plot Figure4/5 A
subplot('position',[0.3 0.75+0.04 0.296774193548387 0.02]) %raster and psth
trial_number=length(trial_spike);
raster_mat=zeros(trial_number*12,1000);
i=0;
for i_trial=trial_spike
    i=i+1;
    for j=1:12
        spike_peri_ras=[dat.spike_correct{i_trial}-dat.photo_ts{i_trial}(j)]*1000;
        spike_peri_ras(spike_peri_ras<-70)=[];
        spike_peri_ras(spike_peri_ras>710)=[];
        i_fig=(i-1)*11+j;
        spike_peri_ras(spike_peri_ras<-20)=[];
        raster_mat(i_fig,floor(spike_peri_ras)+20+1)=1;
        hold on
    end
end
switch coord
    case 1 ,fil_dot=ones(40,1);
    case 2, fil_dot=ones(20,1);
end

raster_mat_fil = conv2(raster_mat,fil_dot);
imagesc(raster_mat_fil==0,[0 1])
colormap(gca,gray)
xlim([1 (time_end-7)*10])
ylim([0 trial_number*12])
switch coord
case 1, set(gca,'ytick',[1 200*12 ],'yticklabel',{'1','200'})
case 2, set(gca,'ytick',[1,400*12],'yticklabel',{'1','400'})
end
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow','tickdir','out','box','off')
ylh=ylabel('Trials ','Fontsize',fs,'FontName','Arial narrow');

%% 
%calculate psth
edges_kthdot=-0.02:0.001:0.34;
spike_kthdot=cell(1,12);
N_use_spike=0;
for i=1:length(dat.spike_correct)
N_use_spike=N_use_spike+1;
for k=1:12
    spike_kthdot{k}=[spike_kthdot{k} ;dat.spike_correct{i}-dat.photo_ts{i}(k)];
end
end
for k=1:12
   psth_kthdot(k,:)=histcounts( spike_kthdot{k},edges_kthdot);
end
psth_conv=conv(nanmean(psth_kthdot),con);
psth_conv=psth_conv(11:310);
psth_conv=psth_conv*1000/length(trial_spike);
%%
%draw psth Figure4/5 A
subplot('position',[0.3 0.75 0.296774193548387 0.04]) %raster and psth
trial_number=1;
psth_con_sample=2*trial_number*psth_conv;
psth_con_dot_fig=-2*trial_number+2*trial_number*psth_conv/max(psth_conv)/1.2;
plot(psth_con_dot_fig,'-','Color','k','Linewidth',2)
hold on
plot([20,20],[-2*trial_number 0],'k:','Linewidth',1)
plot([70,70],[-2*trial_number 0],'k--','Linewidth',0.5)
plot([170,170],[-2*trial_number 0],'k:','Linewidth',0.5)
plot([220,220],[-2*trial_number 0],'k--','Linewidth',0.25)
plot([320,320],[-2*trial_number 0],'k:','Linewidth',0.5)
ylim([-2*trial_number 0])
xlim([1 (time_end-7)*10])
hold on
switch coord
    case 1
        legend_psth_1=5;
        legend_psth_2=8;
    case 2
        legend_psth_1=5;
        legend_psth_2=10;
end
legend_psth_1_norm=-2*trial_number+2*trial_number*(legend_psth_1)/max(psth_conv)/1.2;
legend_psth_2_norm=-2*trial_number+2*trial_number*(legend_psth_2)/max(psth_conv)/1.2;
legend_psth_0_norm=-2*trial_number+2*trial_number*(0)/max(psth_conv)/1.2;
set(gca,'xtick',[20 70 170 220],'xticklabel',{''},'TickDir','out','box','off')
set(gca,'ytick',[legend_psth_0_norm,legend_psth_2_norm],'yticklabel',{0,legend_psth_1})
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
ylh=ylabel('sp/s ','Fontsize',fs,'FontName','Arial narrow');
%%
%Preparation to draw Recptive Fields
for bin=1:length(RF_dot_bkg)
color_min_bin_cand(bin)=nanmin(nanmin(nanmin([RF_dot_bkg{bin} RF_dot_retino{bin}])));
color_max_bin_cand(bin)=nanmax(nanmax(nanmax([RF_dot_bkg{bin} RF_dot_retino{bin}])));
end
%%
%draw receptrive field for dot/ret Figure4/5 C
subplot('position',[0.3 0.45 0.12 0.09])
if mon==4
    flip=false;
else
    flip=true;
end
show_RF(RF_dot_retino{data_bin}',color_min_bin_cand(data_bin),color_max_bin_cand(data_bin),P_CHIV_cell(data_bin,2),Entropy_cell(data_bin,2),'retino',fs,flip)
%%
%draw receptrive field for dot/bkg Figure4/5 D
subplot('position',[0.57 0.45 0.12 0.09])
show_RF(RF_dot_bkg{data_bin}',color_min_bin_cand(data_bin),color_max_bin_cand(data_bin),P_CHIV_cell(data_bin,1),Entropy_cell(data_bin,1),'bkg',fs,flip)
%%
%draw receptrive field for bkg.ret Figure4/5 E
subplot('position',[0.435 0.45 0.12 0.09])
%convolution with bkg shape
conv_resp_frafix=conv2(squeeze(nansum(spikes_bkg_ret(:,data_bin,:,:))),frame_con);
conv_stim_frafix=conv2(stim_bkg_ret,frame_con);
RF_bkg_ret=(conv_resp_frafix./conv_stim_frafix)'/12;
RF_bkg_ret=RF_bkg_ret(11:51,(11:71));
show_RF(RF_bkg_ret,color_min_bin_cand(data_bin),color_max_bin_cand(data_bin),P_CHIV_cell(data_bin,3),Entropy_cell(data_bin,3),'bkg_ret',fs,flip)
%%
%draw color bar
subplot('position',[0.695 0.445 0.03 0.1])
axis off
col_diff=10*(color_max_bin_cand(data_bin)-color_min_bin_cand(data_bin));
col_min=color_min_bin_cand(data_bin)*10;
col_max=color_max_bin_cand(data_bin)*10;
switch coord
    case 1,colbar_leg=[0 10 20];
    case 2,colbar_leg=[4 8 12];
end
colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],... 
    'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
%%  
%draw time course of Inforamtion value (entropy) of the exemplified neuron
%Figure 4/5 B
subplot('position',[0.3 0.62 0.4 0.1])
P_ab_thre=zeros(36,3);
P_ab_thre(Suc5(((P_CHIV_cell(:,1)<=p_thre)&(P_CHIV_dotsh_cell(:,1)<=p_thre))'),1)=1;
P_ab_thre(Suc5(((P_CHIV_cell(:,2)<=p_thre)&(P_CHIV_dotsh_cell(:,2)<=p_thre))'),2)=1;
P_ab_thre(Suc5((P_CHIV_cell(:,3)<=p_thre)'),3)=1;
Entropy_cell(:,1)=Entropy_cell(:,1)-(nanmean(I_SHUFF_cell(:,:,1))'.*P_ab_thre(:,2));
plot(1:36,Entropy_cell(:,1),'m:',1:36,Entropy_cell(:,2),'b:',1:36,Entropy_cell(:,3),'g:','LineWidth',1)
hold on
P_ab_thre(P_ab_thre==0)=NaN; 
plot(1:36,Entropy_cell(:,1).*P_ab_thre(:,1),'m-',1:36,Entropy_cell(:,2).*P_ab_thre(:,2),'b-',1:36,Entropy_cell(:,3).*P_ab_thre(:,3),'g-','LineWidth',2)
switch coord
    case 1
         legend_ent_1=5*2;
         legend_ent_2=10*2;
         legend_ent_3=15*2;
         lim_max=18*2;
set(gca,'ytick',[0 legend_ent_1 legend_ent_2 legend_ent_3],'yticklabel',{0,legend_ent_1,legend_ent_2 legend_ent_3},'TickDir','out','box','off')
    case 2
        legend_ent_1=0.5*2;
        legend_ent_2=1*2;
        lim_max=1.2*2;
set(gca,'ytick',[0 legend_ent_1 legend_ent_2],'yticklabel',{0,legend_ent_1,legend_ent_2},'TickDir','out','box','off')
end
set(gca,'xtick',[1 11 21 ],'xticklabel',{'','','','',''},'TickDir','out','box','off')
set(gca,'Linewidth',2,'Fontsize',fs,'FontName','Arial narrow')
xlim([-1 time_end])
ylim([0 lim_max])
plot([1,1],[0 lim_max],'k:','Linewidth',1)
plot([6,6],[0 lim_max],'k--','Linewidth',1)
set(gca,'xtick',[1 6 11 16 21 26],'xticklabel',{'0','50','100','150','200','250','300'},'TickDir','out','Linewidth',1.5)
ylh2=ylabel(['Information (bits /s)'],'Fontsize',fs,'FontName','Arial narrow');
%%
%draw time course of Inforamtion value (entropy) of all neurons which
%represented significant infomation
%Figure4/5 F
Entropy_thre_tmp(:,:,1:2)=(P_CHIV(:,:,1:2)<=p_thre).*(P_CHIV_dotsh(:,:,1:2)<=p_thre).*Entropy(:,:,1:2);
Entropy_thre_tmp(:,:,3)=(P_CHIV(:,:,3)<=p_thre).*Entropy(:,:,3);
for i=1:size(Entropy_thre_tmp,1)
Entropy_thre(i,:,1)=Suc5(Entropy_thre_tmp(i,:,1));
Entropy_thre(i,:,2)=Suc5(Entropy_thre_tmp(i,:,2));
Entropy_thre(i,:,3)=Suc5(Entropy_thre_tmp(i,:,3));
end
N_dotret=sum(sum(Entropy_thre(:,:,2),2)>0);
subplot('position',[0.3 0.28 0.4 0.1])
amp=10;
switch coord
    case 1
        show_SumEntropy(Entropy_thre(:,:,2)/N_dotret,1,time_end,fs,'b')
        plot([1,1],[0 amp*max(sum(Entropy_thre(:,:,2)))],'k:','Linewidth',1)
        plot([6,6],[0 amp*max(sum(Entropy_thre(:,:,2)))],'k--','Linewidth',1)
        set(gca,'ytick',[0,100 200],'yticklabel',{'0','100','200'},'TickDir','out','box','off','Linewidth',1.5)
        set(gca,'ytick',[0 2 4],'yticklabel',{'0','2','4'},'TickDir','out','box','off')
        
    case 2
        ent_shuff_mean=squeeze(nanmean(I_SHUFF,2));
        entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0).*(Entropy_thre(:,:,2)>0);
        entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
        entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
        N_dotbkg=sum(sum(entropy_thre_minus_spread,2)>0);

        show_SumEntropy(entropy_thre_minus_spread/N_dotbkg,1,time_end,fs,'m')
        plot([1,1],[0 amp*max(sum(entropy_thre_minus_spread))],'k:','Linewidth',1)
        plot([6,6],[0 amp*max(sum(entropy_thre_minus_spread))],'k--','Linewidth',1)
        set(gca,'ytick',[0 5 10 15],'yticklabel',{'0','5','10','15'},'TickDir','out','box','off','Linewidth',1.5)
        set(gca,'ytick',[0 .3 .6],'yticklabel',{'0','0.3','0.6'},'TickDir','out','box','off')
end
set(gca,'xtick',[1 6 11 16 21 26],'xticklabel',{'0','50','100','150','200','250','300'},'TickDir','out')
ylh3=ylabel(['Mean Inforamtion (bits /s)'],'Fontsize',fs,'FontName','Arial narrow');
xlh=xlabel('Time from Dot onset (ms)');

%%
%calculate half maximum and draw (dotted line in Figure4/5 C)
switch coord
case 1, information_sum=sum(squeeze(Entropy_thre(:,:,2)))/N_dotret;
[half_max]=calc_half_max(information_sum,2);
case 2, information_sum=sum(entropy_thre_minus_spread)/N_dotbkg;half=2;
max_information=max(information_sum);
half_period_lobical=information_sum<max_information/half;half_period_lobical(14)=0;
rising=find(diff(half_period_lobical)==-1);
falling=find(diff(half_period_lobical)==1);
rise_inf=information_sum([rising rising+1]);
fall_inf=information_sum([falling falling+1]);
half_max(1)=rising(1)+(max_information/half-rise_inf(1))/diff(rise_inf);
half_max(2)=falling(1)+(max_information/half-fall_inf(1))/diff(fall_inf);
end
plot([half_max(1) half_max(2)],[1*max(information_sum)/2 1*max(information_sum)/2],'k:','LineWidth',3)

%%
if 0
cd figs
switch coord
    case 1, print('fig4.ai', '-dpdf', '-painters','-bestfit')
    case 2, print('fig5.ai', '-dpdf', '-painters','-bestfit')
end
cd ../../Scripts_forFigure/
end