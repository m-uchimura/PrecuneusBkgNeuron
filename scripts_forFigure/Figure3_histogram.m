%%
%this file compute mean firing rate around fp, bkg, dot presentations and
%draw histogram shown in Figure 3 D-F

%%
load ../Data_open/data.mat
%%
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
load ../Data_open/psth psth psth_fp psth_dot_each ts_fp_cells ts_dot_cells psth_r ts_eachdot_cells

%%
%calculate mean firing rate around fp, bkg, dot presentations
i=0;
for i_cell=1:length(data_bkg1)
i=i+1;
clear baseline_spikes fp_spikes bkg_spikes dot_spikes SPIKE_EXIST eachdot_spikes_up eachdot_spikes_dn
ts_dot_tmp=ts_dot_cells{i};
ts_fp_tmp=ts_fp_cells{i};
ts_eachdot_tmp=ts_eachdot_cells{i};

trial_spike=get_trial_spike(data_bkg1(i));
for iTrial=trial_spike
baseline_spikes(iTrial)= sum(ts_dot_tmp{iTrial}>-0.3&ts_dot_tmp{iTrial}<-0.2)/0.1; 
fp_spikes(iTrial)=sum(ts_fp_tmp{iTrial}>0.03&  ts_fp_tmp{iTrial}<0.1)/(0.1-0.03) ;
bkg_spikes(iTrial)=sum(ts_dot_tmp{iTrial}>-0.17&  ts_dot_tmp{iTrial}<0)/0.17 ;
dot_spikes(iTrial)=sum(ts_dot_tmp{iTrial}>0.3&  ts_dot_tmp{iTrial}<1.8)/1.5 ;
SPIKE_EXIST(iTrial)=~isempty(ts_dot_tmp{iTrial});

for iDot=1:12
eachdot_spikes_up(iTrial,iDot)=sum(ts_eachdot_tmp{iTrial,iDot}>.04&ts_eachdot_tmp{iTrial,iDot}<.115)/.075 ;
eachdot_spikes_dn(iTrial,iDot)=sum(ts_eachdot_tmp{iTrial,iDot}>-.02&ts_eachdot_tmp{iTrial,iDot}<.02)/.04 ;
end
end
%
[h_dot(i),p_dot(i)]=ttest(baseline_spikes(SPIKE_EXIST),dot_spikes(SPIKE_EXIST),'Alpha',0.005);
[h_fp(i),p_fp(i)]=ttest(baseline_spikes(SPIKE_EXIST),fp_spikes(SPIKE_EXIST),'Alpha',0.005);
[h_frame(i),p_frame(i)]=ttest(baseline_spikes(SPIKE_EXIST),bkg_spikes(SPIKE_EXIST),'Alpha',0.005);

eachdot_spikes_up_usetrial=eachdot_spikes_up(SPIKE_EXIST,:);
eachdot_spikes_dn_usetrial=eachdot_spikes_dn(SPIKE_EXIST,:);
eachdot_spikes_up_reshape=reshape(eachdot_spikes_up_usetrial,1,size(eachdot_spikes_up_usetrial,1)*size(eachdot_spikes_up_usetrial,2));
eachdot_spikes_dn_reshape=reshape(eachdot_spikes_dn_usetrial,1,size(eachdot_spikes_dn_usetrial,1)*size(eachdot_spikes_dn_usetrial,2));
eachdot_spikes_up_ave=nanmean(eachdot_spikes_up_usetrial,2);
eachdot_spikes_dn_ave=nanmean(eachdot_spikes_dn_usetrial,2);
[h_eachdot(i),p_eachdot(i)]=ttest(eachdot_spikes_up_reshape,eachdot_spikes_dn_reshape,'Alpha',0.005);
[h_eachdot_ave(i),p_eachdot_ave(i)]=ttest(eachdot_spikes_up_ave,eachdot_spikes_dn_ave,'Alpha',0.005);

mFR_bkg(i)=nanmean(bkg_spikes(SPIKE_EXIST));
mFR_fp(i)=nanmean(fp_spikes(SPIKE_EXIST));
mFR_dot(i)=nanmean(dot_spikes(SPIKE_EXIST));
mFR_base(i)=nanmean(baseline_spikes(SPIKE_EXIST));
mFR_eachdot_up(i,:)=nanmean(eachdot_spikes_up_usetrial);
mFR_eachdot_dn(i,:)=nanmean(eachdot_spikes_dn_usetrial);
end

h_dot(isnan(h_dot))=0;
h_fp(isnan(h_fp))=0;
h_frame(isnan(h_frame))=0;

mFR_fp_ratio=log10((mFR_fp)./mFR_base);
mFR_frame_ratio=log10((mFR_bkg)./mFR_base);
mFR_dot_ratio=log10((mFR_dot)./mFR_base);

h_eachdot(isnan(h_eachdot))=0;
h_eachdot_ave(isnan(h_eachdot_ave))=0;

%%
%draw histograms shown in Figure 3 D,E,F
figure(1)
ysca='linear';
yti=[10 20,30,40];
ytilabel={'10','20','30','40'};
yli=[0 30];

%Figure 3D around bkg
subplot('Position',[0.2 0.13+0.12-0.05 0.185 0.08])
hold on
histogram(mFR_frame_ratio(logical(h_frame)),-3:.02:3,'FaceColor','r','EdgeColor','r')
hist_mFR_frame_ratio=histcounts(mFR_frame_ratio,-3:.02:3);
for i=1:600
    x_i_frame(i)=-3+floor(i/2)*0.02;
    y_i_frame(i)=hist_mFR_frame_ratio(ceil(i/2));
end
plot(x_i_frame,y_i_frame,'k-','Linewidth',0.25)
xlim([-2,2])
ylim(yli)
set(gca,'xtick',[-2 -1 0 1 2],'xtickLabel',{'','','','',''},'ytick',yti,'ytickLabel',ytilabel,...
    'YMinorTick','off','yscale',ysca,'tickdir','out','box','off','Linewidth',1.5)

%Figure 3E during dot
subplot('Position',[0.415 0.13+0.12-0.05 0.185 0.08])
hist_mFR_dot_ratio=histcounts(mFR_dot_ratio,-3:.02:3);
for i=1:600
    x_i_dot(i)=-3+floor(i/2)*0.02;
    y_i_dot(i)=hist_mFR_dot_ratio(ceil(i/2));
end
hold on
histogram(mFR_dot_ratio(logical(h_dot)),-3:.02:3,'FaceColor','r','EdgeColor','r')
plot(x_i_dot,y_i_dot,'k-','Linewidth',0.25)
xlim([-2,2])
ylim([0.5,500])
ylim(yli)
set(gca,'xtick',[-2 -1 0 1 2],'xtickLabel',{'','','','',''},'ytick',yti,'ytickLabel',ytilabel,...
    'YMinorTick','off','yscale',ysca,'tickdir','out','box','off','Linewidth',1.5)

%Figure 3F around each dot
subplot('Position',[0.63 0.13+0.12-0.05 0.185 0.08])
mean_eachdot=nanmean(mFR_eachdot_up,2);%
mean_eachdot_baseline=nanmean(mFR_eachdot_dn,2);%
mFR_eachdot_ratio=log10(mean_eachdot./mean_eachdot_baseline);
hold on
histogram(mFR_eachdot_ratio(logical(h_eachdot)),-1.5:.005:1.5,'FaceColor','r','EdgeColor','r')
hist_mFR_eachdot_ratio=histcounts(mFR_eachdot_ratio,-1.5:.005:1.5);
histogram(mFR_eachdot_ratio(logical(h_eachdot)),-1.5:.005:1.5,'FaceColor','r','EdgeColor','r')

for i=1:1200
    x_i_eachdot(i)=-1.5+floor(i/2)*0.005;
    y_i_eachdot(i)=hist_mFR_eachdot_ratio(ceil(i/2));
end
plot(x_i_eachdot,y_i_eachdot,'k-','Linewidth',0.25)
xlim([-1.5,1.5])
xlim([-.6,.6])
set(gca,'xtick',[20 70],'xticklabel',{'',''},'TickDir','out')
set(gca,'xTick',[-1 0 1],'xticklabel',{'10^{-1}','1','10'},'ytick',[1 10 100],'yticklabel',{'1','10','100'},...
    'YMinorTick','off','yscale',ysca,'tickdir','out','box','off','Linewidth',1.5)
set(gca,'xtick',[-2 -1 -.5 0 .5 1 2],'xtickLabel',{'','','','','','',''},'ytick',yti,'ytickLabel',ytilabel,...
    'YMinorTick','off','yscale',ysca,'tickdir','out','box','off','Linewidth',1.5)
ylim(yli)
ylim([0 22])


%%
%check the point of examplified neuron　in Figure 3D,E,F
iCell_cand=[41,75];
D=mFR_frame_ratio(iCell_cand);
subplot('Position',[0.2 0.13+0.12-0.05 0.185 0.08])
plot(D(1),10,'r*',D(2),20,'b*')
E=mFR_dot_ratio(iCell_cand);
subplot('Position',[0.415 0.13+0.12-0.05 0.185 0.08])
plot(E(1),15,'r*',E(2),20,'b*')
F=mFR_eachdot_ratio(iCell_cand)';
subplot('Position',[0.63 0.13+0.12-0.05 0.185 0.08])
plot(F(1),15,'r*',F(2),20,'b*')

