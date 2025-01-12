%%
%this file draw the resuls of neurons which represent all of the dot/ret, dot/bkg,and
%bkg/ret neurons (Figure7)

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
%%
coordinate='triple';
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
iCell=114;data_bin1=8;data_bin2=15;%file='m4c104r1s1n3'
load ../Data_open/data.mat

dat=data_bkg1(iCell);
file=dat.file;
mon=str2double(file(2));
cite=str2double(file(4:6));
dat.spike_correct=[];
tk=0;
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
trial_spike=get_trial_spike(dat);

%%
%draw raster plot　Figure 7A
subplot('position',[0.3 0.78+0.05+0.04 0.296774193548387 0.02])
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
        spike_peri_ras(spike_peri_ras<-20)=[];%for delete minus for matrix
        raster_mat(i_fig,floor(spike_peri_ras)+20+1)=1;%matrix
        hold on
    end
end
fil_dot=ones(20,1);
raster_mat_fil = conv2(raster_mat,fil_dot);
A=sum(raster_mat,2);
imagesc(raster_mat_fil==0,[0 1])
colormap(gca,gray)
xlim([1 (time_end-7)*10])
ylim([0 trial_number*12])

set(gca,'ytick',[1 400*12 ],'yticklabel',{'1','400'})
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow','tickdir','out','box','off')
ylh=ylabel('Trials ','Fontsize',fs,'FontName','Arial narrow');
%%
%calculate psth
edges_kthdot=-0.02:0.001:0.34;
spike_kthdot=cell(1,12);
N_use_spike=0;
for i=1:length(dat.spike_correct)
N_use_spike=N_use_spike+1;
spike_raster{i}=(dat.spike_correct{i}-dat.photo_correct{i}(1)); % 0-500
for k=1:12
spike_kthdot{k}=[spike_kthdot{k} ;dat.spike_correct{i}-dat.photo_correct{i}(k)];
end
end
for k=1:12
psth_kthdot(k,:)=histcounts( spike_kthdot{k},edges_kthdot);
end
psth_conv=conv(nanmean(psth_kthdot),con);
psth_conv=psth_conv(11:310);
psth_conv=psth_conv*1000/length(trial_spike);
%%
%draw psth Figure7 A
subplot('position',[0.3 0.78+0.05 0.296774193548387 0.04])
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
xlim([1 (time_end+1)*10])
ylim([-2*trial_number 0])
xlim([1 (time_end-7)*10])
hold on
legend_psth_1=10;
legend_psth_2=20;
legend_psth_1_norm=-2*trial_number+2*trial_number*(legend_psth_1)/max(psth_conv)/1.2;
legend_psth_2_norm=-2*trial_number+2*trial_number*(legend_psth_2)/max(psth_conv)/1.2;
legend_psth_0_norm=-2*trial_number+2*trial_number*(0)/max(psth_conv)/1.2;
set(gca,'xtick',[20 70 170 220],'xticklabel',{''},'TickDir','out','box','off')
set(gca,'ytick',[legend_psth_0_norm,legend_psth_1_norm,legend_psth_2_norm],'yticklabel',{0,legend_psth_1,legend_psth_2})
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
ylh=ylabel('Sp/s','Fontsize',fs,'FontName','Arial narrow');

%%
%Preparation to draw Recptive Fields
for bin=1:length(RF_dot_bkg)
color_min_bin_cand(bin)=nanmin(nanmin(nanmin([RF_dot_bkg{bin} RF_dot_retino{bin}])));
color_max_bin_cand(bin)=nanmax(nanmax(nanmax([RF_dot_bkg{bin} RF_dot_retino{bin}])));
end
%c
%draw receptrive field 
flip=0; %because monkey==4
%draw receptrive field for dot/ret at the earlier timing Figure7 C
subplot('position',[0.3 0.48+0.05 0.12 0.09])
show_RF(RF_dot_retino{data_bin1}',color_min_bin_cand(data_bin1),color_max_bin_cand(data_bin1),P_CHIV_cell(data_bin1,2),Entropy_cell(data_bin1,2),'retino',fs,flip)

%draw receptrive field for dot/bkg at the earlier timing Figure7 D
subplot('position',[0.57 0.48+0.05 0.12 0.09])
show_RF(RF_dot_bkg{data_bin1}',color_min_bin_cand(data_bin1),color_max_bin_cand(data_bin1),P_CHIV_cell(data_bin1,1),Entropy_cell(data_bin1,1),'bkg',fs,flip)

%draw receptrive field for bkg/ret at the earlier timing Figure7 E
subplot('position',[0.435 0.48+0.05 0.12 0.09])
%convolution with actual bkg shape
conv_resp_frafix=conv2(squeeze(nansum(spikes_bkg_ret(:,data_bin1,:,:))),frame_con);
conv_stim_frafix=conv2(stim_bkg_ret,frame_con);
RF_bkg_ret=(conv_resp_frafix./conv_stim_frafix)'/12;
RF_bkg_ret=RF_bkg_ret(11:51,(11:71));
show_RF(RF_bkg_ret,color_min_bin_cand(data_bin1),color_max_bin_cand(data_bin1),P_CHIV_cell(data_bin1,3),Entropy_cell(data_bin1,3),'bkg_ret',fs,flip)

%%color bar at the earlier timing
subplot('position',[0.695 0.48+0.05 0.03 0.1])
axis off
col_diff=10*(color_max_bin_cand(data_bin1)-color_min_bin_cand(data_bin1));
col_min=color_min_bin_cand(data_bin1)*10;
col_max=color_max_bin_cand(data_bin1)*10;
colbar_leg=[2 10 18];
colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],... 
    'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
%%
%draw receptrive field for dot/ret at the later timing Figure7 F
subplot('position',[0.3 0.34+0.05 0.12 0.09])
show_RF(RF_dot_retino{data_bin2}',color_min_bin_cand(data_bin2),color_max_bin_cand(data_bin2),P_CHIV_cell(data_bin2,2),Entropy_cell(data_bin2,2),'retino',fs,flip)

%draw receptrive field for dot/bkg at the later timing Figure7 G
subplot('position',[0.57 0.34+0.05 0.12 0.09])
show_RF(RF_dot_bkg{data_bin2}',color_min_bin_cand(data_bin2),color_max_bin_cand(data_bin2),P_CHIV_cell(data_bin2,1),Entropy_cell(data_bin2,1),'bkg',fs,flip)

%draw receptrive field for blg/ret at the later timing Figure7 H
subplot('position',[0.435 0.34+0.05 0.12 0.09])
conv_resp_frafix=conv2(squeeze(nanmean(spikes_bkg_ret(:,data_bin2,:,:))),frame_con);
conv_stim_frafix=conv2(stim_bkg_ret,frame_con);
RF_bkg_ret=(conv_resp_frafix./conv_stim_frafix)';
RF_bkg_ret=RF_bkg_ret(11:51,(11:71));
show_RF(RF_bkg_ret,color_min_bin_cand(data_bin2),color_max_bin_cand(data_bin2),P_CHIV_cell(data_bin2,3),Entropy_cell(data_bin2,3),'bkg_ret',fs,flip)

%%color bar at the later timing
subplot('position',[0.695 0.34+0.05 0.05 0.1])
axis off
col_diff=10*(color_max_bin_cand(data_bin2)-color_min_bin_cand(data_bin2));
col_min=color_min_bin_cand(data_bin2)*10;
col_max=color_max_bin_cand(data_bin2)*10;
colbar_leg=[0 2 4];
colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],...
    'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
%%
%draw time course of Inforamtion value (entropy) of the exemplified neuron
%Figure7 B
subplot('position',[0.3 0.65+0.05 0.4 0.1])
P_ab_thre=zeros(36,3);
P_ab_thre(Suc5(((P_CHIV_cell(:,1)<=p_thre)&(P_CHIV_dotsh_cell(:,1)<=p_thre))'),1)=1;
P_ab_thre(Suc5(((P_CHIV_cell(:,2)<=p_thre)&(P_CHIV_dotsh_cell(:,2)<=p_thre))'),2)=1;
P_ab_thre(Suc5((P_CHIV_cell(:,3)<=p_thre)'),3)=1;
Entropy_cell(:,1)=Entropy_cell(:,1)-(nanmean(I_SHUFF_cell(:,:,1))'.*P_ab_thre(:,2));
plot(1:36,Entropy_cell(:,1),'m:',1:36,Entropy_cell(:,2),'b:',1:36,10*Entropy_cell(:,3),'g:','LineWidth',1)
hold on
P_ab_thre(P_ab_thre==0)=NaN; 
plot(1:36,Entropy_cell(:,2).*P_ab_thre(:,2),'b-',1:36,10*Entropy_cell(:,3).*P_ab_thre(:,3),'g-',1:36,Entropy_cell(:,1).*P_ab_thre(:,1),'m-','LineWidth',2) %%bkg_ret 10 times
legend_ent_1=2;
legend_ent_2=4;
legend_ent_1=2.5*2;
legend_ent_2=5*2;
set(gca,'ytick',[0 legend_ent_1 legend_ent_2],'yticklabel',{0,legend_ent_1,legend_ent_2},'TickDir','out','box','off')
set(gca,'xtick',[1 11 21 ],'xticklabel',{'','','','',''},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
xlim([-1 time_end])
ylim([0 legend_ent_2])
plot([1,1],[0 legend_ent_2],'k:','Linewidth',1)
plot([6,6],[0 legend_ent_2],'k--','Linewidth',1)
set(gca,'xtick',[1 6 11 16 21 26],'xticklabel',{'0','50','100','150','200','250','300'},'TickDir','out')
ylh2=ylabel(['Information (bits/s)'],'Fontsize',fs,'FontName','Arial narrow');

%%
%draw time course of Inforamtion value (entropy) of all neurons which
%represented significant infomation Figure7 I
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

triple_neuron=(sum(entropy_thre_minus_spread,2)>0)&(sum(Entropy_thre(:,:,2),2)>0)&(sum(Entropy_thre(:,:,3),2)>0);
Entropy_thre(:,:,1)=Entropy_thre(:,:,1).*repmat(triple_neuron,1,36);
entropy_thre_minus_spread=entropy_thre_minus_spread.*repmat(triple_neuron,1,36);
Entropy_thre(:,:,2)=Entropy_thre(:,:,2).*repmat(triple_neuron,1,36);
Entropy_thre(:,:,3)=Entropy_thre(:,:,3).*repmat(triple_neuron,1,36);
N_triple=sum(triple_neuron);
amp=10;
%dot/bkg Figure7 I
subplot('position',[0.3 0.20+0.05-0.02 0.4 0.1])
show_SumEntropy(entropy_thre_minus_spread/N_triple,1,time_end,fs,'m')
set(gca,'ytick',[0 .4 .8],'yticklabel',{'0','0.4','0.8'},'TickDir','out','box','off')
half=2;
information_mean_dotbkg=sum(entropy_thre_minus_spread/N_triple);
plot([1,1],[0 amp*max(information_mean_dotbkg)],'k:','Linewidth',1)
plot([6,6],[0 amp*max(information_mean_dotbkg)],'k--','Linewidth',1)
[half_max_dot_bkg]=calc_half_max(information_mean_dotbkg,half);
plot([half_max_dot_bkg(1) half_max_dot_bkg(2)],[1*max(information_mean_dotbkg)/half...
    1*max(information_mean_dotbkg)/half],'k:','LineWidth',3)
ylim([0 .8])

%dot/ret Figure7 J
subplot('position',[0.3 0.11+0.05-0.02 0.4 0.07])
show_SumEntropy(Entropy_thre(:,:,2)/N_triple,1,time_end,fs,'b')
plot([1,1],[0 amp*max(sum(Entropy_thre(:,:,2)))],'k:','Linewidth',1)
plot([6,6],[0 amp*max(sum(Entropy_thre(:,:,2)))],'k--','Linewidth',1)
information_mean_dotret=sum(squeeze(Entropy_thre(:,:,2))/N_triple);
ylim([0 6.5])
set(gca,'ytick',[0 3 6 ],'yticklabel',{'0','3','6'},'TickDir','out','box','off')
[half_max_dot_ret]=calc_half_max(information_mean_dotret,half);
plot([half_max_dot_ret(1) half_max_dot_ret(2)],[1*max(information_mean_dotret)/half...
    1*max(information_mean_dotret)/half],'k:','LineWidth',3)

%bkg/ret Figure7 K
subplot('position',[0.3 0.02+0.05-0.02 0.4 0.07])
show_SumEntropy(Entropy_thre(:,:,3)/N_triple,1,time_end,fs,'g')
set(gca,'xtick',[1 6 11 16 21 26],'xticklabel',{'0','50','100','150','200','250','300'},'TickDir','out')
set(gca,'ytick',[0 1 2],'yticklabel',{'0','1','2'},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial narrow')
ylim([0 2.4])
xlim([-1 time_end])
xlh=xlabel('Time from Dot onset (ms)');

%%
%fit dot/bkg information from 2 retinal information (dot/ret info and bkg/ret)
subplot('position',[0.3 0.20+0.05-0.02 0.4 0.1])
global dot_r
global dot_b
global bkg_r
period1=1:25;

dot_b=sum(entropy_thre_minus_spread(:,period1)/N_triple);
dot_r=sum(Entropy_thre(:,period1,2))/N_triple;
bkg_r=sum(Entropy_thre(:,period1,3))/N_triple;
dot_r=[0 0 dot_r];
bkg_r=[0 0 bkg_r];
dot_r([end-1:end])=[];
bkg_r([end-1:end])=[];
p0=[.5 0];
pmin=[0 0];
pmax=[10 0];
opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',@sse_info,...
    'Aineq',[],...
    'bineq',[],...
    'x0',p0,...
    'lb',pmin,...
    'ub',pmax,'options',opts); 
ms = MultiStart('Display','iter','UseParallel','never');
[x,f] = run(ms,problem,10);
dot_b_hat=x(1)*dot_r.*(bkg_r-x(2));
plot(period1,1*dot_b_hat,'r:','LineWidth',2)
sse=sum((dot_b_hat-dot_b).^2);
sse=sse_info(x);
dc=1-sse/var(dot_b)/length(dot_b);
%%
if 0
cd figs
print('fig7.png', '-dpdf', '-painters','-bestfit')
print('fig7.ai', '-dpdf', '-painters','-bestfit')
cd ../
end