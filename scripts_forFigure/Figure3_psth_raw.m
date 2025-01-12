%%
%this file calculate the peri-stimulus histogram(PSTH) for a fixtion point
%presentation and the background (or 1st dot) presentation 
%1st figure corresponds Figure 3 A-B
%2nd figure corresponds Figure 13 E-G
%%
clear
close all
load ../Data_open/data.mat
%%
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
fs=12;
sample_trial=80;
col={[0 0 0],[1 0 1],[0 0 1],[0 1 0]};
ti=0:20;
myu=10;
sd=5;
con=exp(-(ti-myu).^2/(2*sd^2))/(sqrt(2*pi)*sd);
con=con/sum(con);
p_thre=0.004;

figure('position',[0 0 floor(100000*asp_rat_x/res_x) floor(100000*asp_rat_y/res_y)])
%%
%load neual information, Pvalue etc
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF]=get_Entropy50;
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
ent_shuff_mean=squeeze(mean(I_SHUFF,2));
entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0).*(Entropy_thre(:,:,2)>0);
entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
Entropy_thre(:,:,1)=entropy_thre_minus_spread;
Neuron_information=squeeze(sum(Entropy_thre,2)>0);
%%
%load psth_mean_baseline
load ../Data_open/psth psth psth_fp psth_dot_each ts_fp_cells ts_dot_cells psth_r
psth_raw=psth;
psth_fp_raw=psth_fp;
for iCell=1:length(data_bkg1)
if 0
psth(iCell,:)=psth(iCell,:)./nanmean(psth_fp(iCell,301:500));
psth_fixp(iCell,:)=psth_fp(iCell,:)./nanmean(psth_fp(iCell,301:500));
else
psth_fixp=psth_fp;
end
spikes(iCell).raster_dot=ts_dot_cells{iCell};
spikes(iCell).raster_fp=ts_fp_cells{iCell};
end

%%
%draw PSTH of average of all neurons  
i=1;
subplot('Position',[0.1+0.1 0.595+0.03 (0.65-0.1)*0.8 0.08])
flt=true(size(Neuron_information,1),1); 
psth_tmp=psth(flt,:);

psth_fixp_tmp=psth_fixp(flt,:);
m_PSTH=median(psth_tmp);
m_PSTH_fp=median(psth_fixp_tmp);
mean_PSTH=m_PSTH;
mean_PSTH_75=prctile(psth_tmp,75);
mean_PSTH_25=prctile(psth_tmp,25);
mean_PSTH_FP_75=prctile(psth_fixp_tmp,75);
mean_PSTH_FP_25=prctile(psth_fixp_tmp,25);

hold on
plot(401:2800,mean_PSTH(201+400:2500+500),'-','Color',col{i},'Linewidth',2)
plot(401:2800,mean_PSTH_25(201+400:2500+500),':','Color',col{i},'Linewidth',1)
plot(401:2800,mean_PSTH_75(201+400:2500+500),':','Color',col{i},'Linewidth',1)
plot([201,201],[0 40],'b:','Linewidth',1)
plot([501,501],[0 40],'b:','Linewidth',1)
plot([601,601],[0 40],'g:','Linewidth',1)

for k=0:11
    sti=801+150*k;
    plot([sti,sti],[0 40],'k:','Linewidth',1)
end
set(gca,'xtick',[201 501 601 801 951 2601],'xticklabel',{'','','','','',''},'TickDir','out','box','off','Linewidth',2)
set(gca,'ytick',[0 10 20 ],'yticklabel',{'0','10','20' },'TickDir','out','FontSize',fs,'FontName','Arial narrow')

xlim([401 2800])
ylim([.1 3.6])
ylim([0 26])
t=5;
mean_PSTH_peri_min=.4;
mean_PSTH_peri_max=1.8;

%draw average psth around dot presentation
subplot('Position',[0.1+0.7+0.02-0.085 0.595+0.03 0.08 0.08])

for k=1:10
    psth_spk_eachdot(k,:)=mean_PSTH((1001-30+(k+2-1)*150):(1000+(k+2)*150)-10);
    psth_spk_eachdot_25(k,:)=mean_PSTH_25((1001-30+(k+2-1)*150):(1000+(k+2)*150)-10);
    psth_spk_eachdot_75(k,:)=mean_PSTH_75((1001-30+(k+2-1)*150):(1000+(k+2)*150)-10);
end
plot(squeeze(nanmean(psth_spk_eachdot)),'-','Color',col{i},'Linewidth',2)
hold on
cell_PSTH_frame_min=0;cell_PSTH_frame_max=100;
plot([30 30],[min(squeeze(nanmean(psth_spk_eachdot)))*0.9,max(squeeze(nanmean(psth_spk_eachdot)))],'k:','Linewidth',1)
plot([80,80],[min(squeeze(nanmean(psth_spk_eachdot)))*0.9,max(squeeze(nanmean(psth_spk_eachdot)))],'k--','Linewidth',0.5)
plot([30 30],[0,6],'k:','Linewidth',1)
plot([80,80],[0,6],'k--','Linewidth',0.5)
set(gca,'xtick',[30 80],'xticklabel',{'','',''},'TickDir','out','box','off','Linewidth',1.5)
set(gca,'ytick',[.85 .9],'yticklabel',{'.85','.9',''})
set(gca,'ytick',[5 6],'yticklabel',{'5','6',''})

xlim([1 170])
ylim([.85 .93])
ylim([5 6])

%%
%draw PSTH of two example neurons
mean_PSTH_frame_min=0;mean_PSTH_frame_max=4;
iCell_cand=[41,75];%'m4c038r1s1n1' and 'm4c066r1s1n1' 

col_ind={[1 0 0],[0 0 1]};
col_ind={[0 0 0],[0 0 0]};
add=sample_trial;
for cell_i=1:2
iCell=iCell_cand(cell_i);   
trial_number=length( data_bkg1(iCell).photo_ts);

subplot('Position',[0.1+0.1 0.77+(cell_i-1)*0.135+0.03 (0.65-0.1)*0.8 0.04])

for i=1:floor(trial_number/sample_trial):trial_number
    spike_trial_ras_fp=ts_fp_cells{iCell}{i}*1000;
    spike_trial_ras_fp=spike_trial_ras_fp(spike_trial_ras_fp>-200&spike_trial_ras_fp<200);
    spike_trial_ras_fp=spike_trial_ras_fp+201;    
    spike_trial_ras=ts_dot_cells{iCell}{i}*1000;
    spike_trial_ras=spike_trial_ras(spike_trial_ras>-400&spike_trial_ras<2000);
    spike_trial_ras=spike_trial_ras+301+500; 
    i_fig=trial_number-i;
    if length(spike_trial_ras_fp)
        plot(spike_trial_ras_fp,i_fig,'.','Color',col_ind{cell_i},'Markersize',1)
    end
    if length(spike_trial_ras)
        plot(spike_trial_ras,i_fig,'.','Color',col_ind{cell_i},'Markersize',1)
    end
    hold on
end
xlim([1 2800])
xlim([401 2800])
ax = gca;
ax.XAxis.Visible = 'off';
if cell_i==2
    legend_ras=[0 100];
elseif cell_i==1
    legend_ras=[0 500];
end
set(gca,'ytick',legend_ras,'yticklabel',{legend_ras(1) legend_ras(2)},'box','off','TickDir','out','Linewidth',1.5)

subplot('Position',[0.1+0.1 0.69+(cell_i-1)*0.135+0.03 (0.65-0.1)*0.8 0.08])

cell_PSTH=psth_raw(iCell,:);
cell_PSTH_fp=psth_fp_raw(iCell,:);
hold on
plot(401:2800,cell_PSTH(201+400:2500+500),'-','Color',col_ind{cell_i},'Linewidth',2)
plot([201,201],[0 200],'b:','Linewidth',1)
plot([501,501],[0 200],'b:','Linewidth',1)
plot([601,601],[0 200],'g:','Linewidth',1)
for k=0:11
    sti=801+150*k;
    plot([sti,sti],[0 200],'k:','Linewidth',1)
end
xlim([1 2800])
xlim([401 2800])
ylim([0,max(cell_PSTH(201:2500))*1.1])
if cell_i==2
    legend_psth=[0 20,40];
elseif cell_i==1
    legend_psth=[0 75,150];
end
set(gca,'xtick',[201 501 601 801,951,2601],'xticklabel',{'','','','','',''},'ytick',legend_psth,...
    'yticklabel',{legend_psth(1) legend_psth(2) legend_psth(3)},'TickDir','out','box','off','Linewidth',1.5)

subplot('Position',[0.1+0.7+0.02-0.085 0.77+(cell_i-1)*0.135+0.03 0.08 0.04])
for i=1:floor(trial_number/sample_trial):trial_number
   spike_trial_ras=ts_dot_cells{iCell}{i}*1000;
   spike_trial_ras=spike_trial_ras(spike_trial_ras>(t*150-220)&spike_trial_ras<(t*150+560));
   spike_trial_ras=spike_trial_ras-(t*150-170)+1;
    i_fig=trial_number-i;
    if ~isempty(spike_trial_ras)
        plot(spike_trial_ras,i_fig,'.','Color',col_ind{cell_i},'Markersize',1)
    end
    hold on
end
ax = gca;
ax.XAxis.Visible = 'off';
set(gca,'xtick',[],'xticklabel',{},'ytick',[],'yticklabel',{},'box','off','Linewidth',1.5)
set(gca,'ytick',legend_ras,'yticklabel',{legend_ras(1) legend_ras(2)},'box','off','TickDir','out','Linewidth',1.5)
xlim([1 170]) 

%psth around dot presentaion
subplot('Position',[0.1+0.7+0.02-0.085 0.69+(cell_i-1)*0.135+0.03 0.08 0.08])
plot(squeeze(mean(psth_dot_each(iCell,3:12,:),2)),'-','Color',col_ind{cell_i},'Linewidth',2)
hold on
cell_PSTH_frame_min=0;cell_PSTH_frame_max=100;
plot([30 30],[min(squeeze(mean(psth_dot_each(iCell,3:12,:),2)))*0.9,...
    max(squeeze(mean(psth_dot_each(iCell,3:12,:),2)))*1.1],'k:','Linewidth',1)
plot([80,80],[min(squeeze(mean(psth_dot_each(iCell,3:12,:),2)))*0.9,...
    max(squeeze(mean(psth_dot_each(iCell,3:12,:),2)))*1.1],'k--','Linewidth',0.5)
set(gca,'xtick',[30 80],'xticklabel',{'','',''},'TickDir','out','box','off','Linewidth',1.5)
xlim([1 170])
if cell_i==2
    legend_psth_dot=[0 10,20];
elseif cell_i==1
    legend_psth_dot=[0 5,10];
end
set(gca,'xtick',[],'xticklabel',{},'ytick',legend_psth_dot,...
    'yticklabel',{legend_psth_dot(1) legend_psth_dot(2) legend_psth_dot(3)},'box','off','Linewidth',1.5)
end

%%
%PSTH of neurons which represent the dot/bkg(red), dot/ret(blue), or
%bkg/ret(green) information (Figure 13)

figure(2)
clf
set(gcf,'position',[0 0 floor(150000*asp_rat_x/res_x) floor(90000*asp_rat_y/res_y)])
for i=[4 3 2]
subplot('Position',[0.1+0.1 0.1+0.4 0.65-0.1 0.3])
    switch i
        case 1,flt=true(size(Neuron_information,1),1);%flt([460,914])=false;
        case 2,flt=Neuron_information(:,1);
        case 3,flt=Neuron_information(:,2);
        case 4,flt=Neuron_information(:,3);
    end
psth_tmp=psth(flt,:);
psth_fixp_tmp=psth_fixp(flt,:);
m_PSTH=median(psth_tmp);
m_PSTH_fp=median(psth_fixp_tmp);
mean_PSTH=m_PSTH;
mean_PSTH_75=prctile(psth_tmp,75);
mean_PSTH_25=prctile(psth_tmp,25);
mean_PSTH_FP_75=prctile(psth_fixp_tmp,75);
mean_PSTH_FP_25=prctile(psth_fixp_tmp,25);
plot(1:400,(m_PSTH_fp(301:700)),'-','Color',col{i},'Linewidth',2)
hold on
plot(401:2800,mean_PSTH(201+400:2500+500),'-','Color',col{i},'Linewidth',2)
plot([201,201],[0 40],'b:','Linewidth',1)
plot([501,501],[0 40],'b:','Linewidth',1)
plot([601,601],[0 40],'g:','Linewidth',1)
plot([1 400],[1 1],'k:','Linewidth',1)
plot([401 2900],[1 1],'k:','Linewidth',1)
for k=0:11
    sti=801+150*k;
    plot([sti,sti],[0 40],'k:','Linewidth',1)
end
set(gca,'xtick',[201 601 801 951 2601],'xticklabel',{'','','','',''},'TickDir','out','box','off','Linewidth',1.5)
set(gca,'ytick',[0 1 2 3],'yticklabel',{'0','1','2','3' },'TickDir','out','FontSize',fs,'FontName','Arial narrow')
set(gca,'ytick',[0 10 20],'yticklabel',{'0','10','20' },'TickDir','out','FontSize',fs,'FontName','Arial narrow')

xlim([1 2800])
xlim([401 2800])
ylim([.7 3])
ylim([0 22])
mean_PSTH_peri_min=.7;
mean_PSTH_peri_max=2;
end

%%
%draw figure 13 X(psth around dot prsentation and calculate the latency)
sub=0;
subplot('Position',[0.1+0.1 0.1 0.3 0.3])
chh=[];
for i=[3 2]
    switch i
        case 1,flt=true(size(Neuron_information,1),1);
        case 2,flt=Neuron_information(:,1);
        case 3,flt=Neuron_information(:,2);
        case 4,flt=Neuron_information(:,3);
    end
psth_tmp=psth_r(flt,:);
psth_tmp2=[];
for ii=1:3500/2
psth_tmp2(:,ii)=mean(psth_tmp(:,2*ii-1:2*ii),2);
end
psth_tmp=movmean(psth_tmp,2,2);
psth_tmp=psth_tmp(:,2:2:end);

if sub
clear Apsth_spk_eachdot
for k=1:10
    Apsth_spk_eachdot(:,k,:)=psth_tmp2(:,(501-15+(k+1)*75):(500+(k+2)*75-5));
end
B=Apsth_spk_eachdot(:,:,6:25);
BB=B(:,:);
BB=squeeze(mean(B,2));
m_psth_tmp2=mean(psth_tmp2(:,6:25),2);
sd_psth_tmp2=mean(psth_tmp2(:,6:25),2);

m_psth_tmp2=mean(BB,2);
sd_psth_tmp2=mean(BB,2);
psth_tmp2=psth_tmp2./repmat(m_psth_tmp2,1,1750);
end

m_PSTH=median(psth_tmp2);
mean_PSTH=m_PSTH;
mean_PSTH_dot=mean_PSTH(501:1500);
psth_spk_eachdot=[];
for k=1:10
    psth_spk_eachdot(k,:)=mean_PSTH((501-15+(k+1)*75):(500+(k+2)*75-5));
end

med_psth=squeeze(mean(psth_spk_eachdot));
plot(med_psth,'-','Color',col{i},'Linewidth',2)
hold on
cell_PSTH_frame_min=0;cell_PSTH_frame_max=100;
set(gca,'xtick',[15 40],'xticklabel',{'','',''},'TickDir','out','box','off','Linewidth',1.5)
set(gca,'ytick',[0 1 1.5 2],'yticklabel',{'0','1','1.5' },'TickDir','out','FontSize',fs,'FontName','Arial narrow')
set(gca,'ytick',[5 8],'yticklabel',{'5','8' },'TickDir','out','FontSize',fs,'FontName','Arial narrow')
plot([1 780],[1 1],'k:','Linewidth',1)

if i==2||i==3
if 1
med_sd=mean(med_psth(6:25))+4*std(med_psth(6:25));
else
psthA=psth_spk_eachdot(:,5:25);
med_sd=mean(psthA(:))+4*std(psthA(:));
end
plot([1 85],[med_sd med_sd],':','Color',col{i},'Linewidth',2)
ch=find((med_psth-med_sd)>0);
if ~isempty(ch)
ch(ch<15)=[];
chh=ch(1);
plot(chh,med_psth(chh),'o','Color',col{i},'Markersize',10)
else
    chh=NaN;
end
if isempty(chh)
    ch_1(iCell)=NaN;
else
   ch_1(iCell)=chh(1);
end
end
plot([15 15],[0 10],'k:','Linewidth',1)
plot([40 40],[0 10],'k--','Linewidth',0.5)
xlim([1 85])
ylim([.8 1.75])
ylim([4.5 9])
if sub 
ylim([.8 1.5])
end
end

%%
%draw figure 13 X(latencies of individua cells )
for i=2:3
switch i
    case 1,flt=true(size(Neuron_information,1),1);
    case 2,flt=Neuron_information(:,1);
    case 3,flt=Neuron_information(:,2);
    case 4,flt=Neuron_information(:,3);
end

clear ch_1 chh  psth_spk_eachdot
for iCell=find(flt)'
clear chh
disp(iCell)
ts_dot=[];ts_fp=[];
ts_fp_trial=[];ts_dot_trial=[];ts_eachdot_trial=[];
nTrials=0;
fr_spk_dot_raw_cell=[];fr_spk_dot=[];fr_spk_dot_cell=[];
%calulate time stamp all and divided with each trial sorted with dot or fixation 
for sp_i=1:length(data_bkg1(iCell).photo_ts)
    tmp_ts_dot=data_bkg1(iCell).spike_ts-data_bkg1(iCell).photo_ts{sp_i}(1);
    tmp_ts_dot(tmp_ts_dot<-1|tmp_ts_dot>2.5)=[];
    ts_dot=[ts_dot;tmp_ts_dot];
    if ~isempty(tmp_ts_dot) nTrials=nTrials+1; end
    tmp_spk_fp=data_bkg1(iCell).spike_ts-(data_bkg1(iCell).photo_ts{sp_i}(1)-data_bkg1(iCell).fp_on(sp_i));
    tmp_spk_fp(tmp_spk_fp<-0.5|tmp_spk_fp>3)=[];
    ts_fp=[ts_fp;tmp_spk_fp];
    ts_fp_trial{sp_i}=tmp_spk_fp;
    ts_dot_trial{sp_i}=tmp_ts_dot;
for iDot=1:12
     tmp_ts_eachdot=data_bkg1(iCell).spike_ts-data_bkg1(iCell).photo_ts{sp_i}(iDot);
     tmp_ts_eachdot(tmp_ts_eachdot<-.1|tmp_ts_eachdot>.3)=[];
     ts_eachdot_trial{sp_i,iDot}=tmp_ts_eachdot;
end
xe_stim=-1:0.001:2.5;
fr_spk_dot_raw_cell(sp_i,:)=histcounts(tmp_ts_dot,xe_stim);
end
psth_tmp=fr_spk_dot_raw_cell;
psth_tmp=movmean(psth_tmp,2,2);
psth_tmp=psth_tmp(:,2:2:end);

m_PSTH=mean(psth_tmp);
mean_PSTH=m_PSTH;
mean_PSTH_dot=mean_PSTH([1002:2:3000]/2);
psth_spk_eachdot=[];
for k=1:10
    psth_spk_eachdot(k,:)=mean_PSTH([(1002-30+(k+2-1)*150):2:(1002+(k+2)*150-10)]/2);
end

m_psth=squeeze(mean(psth_spk_eachdot));
if 1
m_sd=mean(m_psth([12:2:50]/2))+4*std(m_psth([12:2:50]/2));
else
A=psth_spk_eachdot(:,[12:2:50]/2);
m_sd=mean(A(:))+5*std(A(:));
end
ch=find((m_psth-m_sd)>0);
ch(ch<15)=[];
if ~isempty(ch)
chh=ch(1);
else
chh=NaN;
end
if isempty(chh)
   ch_1(iCell)=NaN;
else
   ch_1(iCell)=chh(1);
end
end
ch_1(ch_1==0)=[];
subplot('Position',[0.55 0.1-(i-3)*0.17 0.2 0.13])
histogram([ch_1-15]*2,1:2:100,'FaceColor',col{i},'EdgeColor',col{i})
set(gca,'box','off','TickDir','out','Linewidth',2)
switch i
    case 3,set(gca,'xtick',[0 10 20 30 40 50 60 70 80],'xticklabel',{})
          set(gca,'ytick',[0 2 4 6],'yticklabel',{'0','2','4','6'})
    case 2,set(gca,'xtick',[0 10 20 30 40 50 60 70 80],'xticklabel',[],'ytick',[0 2 4 6],'yticklabel',{'0','2','4','6'})
end
hold on
latencies_All{i}=ch_1;
end
%%
lat_ret=data_bkg1(Neuron_information(:,2));
latencies_ret=lat_ret(latencies_All{3}<35);
lat_bkg=data_bkg1(Neuron_information(:,1));
DD_bkg=lat_bkg(latencies_All{2}<35);
