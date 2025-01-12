%%
%this file the timecourse of draw gaze position (Figure 3C)
%%
load ../Data_open/gaze_data.mat
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
fs=12;
figure(1)
%%
monkey=0;
mon_all(mon_all==4&cite_all>215)=4.5;
mon_all=floor(mon_all);
if monkey==0
mon_flt=true(size(mon_all));
else
mon_flt=mon_all==monkey;
end
gaze_h_all=gaze_h_all(mon_flt,:);
gaze_v_all=gaze_v_all(mon_flt,:);
fp_h_all=fp_h_all(mon_flt,:);
fp_v_all=fp_v_all(mon_flt,:);
%%
%calculate relative eye position
rel_eye_pos_H = NaN(size(gaze_h_all));
rel_eye_pos_V = NaN(size(gaze_h_all));
num_trials=size(gaze_h_all,1);
for i = 1:num_trials
	rel_eye_pos_H(i,:) = gaze_h_all(i,:) - fp_h_all(i);
	rel_eye_pos_V(i,:) = gaze_v_all(i,:) - fp_v_all(i);
end
%%
%baseline average (from FP) and baseline standard deviation
baseline_h=nanmean(rel_eye_pos_H(:,140:160),2);
baseline_v=nanmean(rel_eye_pos_V(:,140:160),2);
mean_baseline=[mean(baseline_h),mean(baseline_v)];
sd_baseline=[std(baseline_h),std(baseline_v)];
EYEH=rel_eye_pos_H-repmat(baseline_h,1,700);
EYEV=rel_eye_pos_V-repmat(baseline_v,1,700);
mEYEH=mean(EYEH);
mEYEV=mean(EYEV);
sdEYEH=nanstd(EYEH);
sdEYEV=nanstd(EYEV);

%%
subplot('Position',[0.1+0.1 0.425+0.06 0.55*23/24*0.8 0.08])
plot(mEYEH,'k-','LineWidth',2)
hold on
plot(mEYEH+sdEYEH,'k-','LineWidth',1)
plot(mEYEH-sdEYEH,'k-','LineWidth',1)
plot([140 140],[-2 2],'b:','LineWidth',1) %fix in
plot([160 160],[-2 2],'g:','LineWidth',1) %bkg
plot([200 200],[-2 2],'k:','LineWidth',1) %dot
xlim([121 580])
ylim([-1.5 1.5])

for k=0:11
    sti=201+30*k;
    plot([sti,sti],[-2 2],'k:','Linewidth',1)
end
set(gca,'XTick',[140 160,200,560],'XTickLabel',{'','','',''},'TickDir','out','box','off','Linewidth',1.5)
set(gca,'ytick',[-1 -0.5 0 0.5 1],'yticklabel',{'-1','','0','','1' },'TickDir','out','FontSize',fs,'FontName','Arial narrow','box','off','Linewidth',1.5)
subplot('Position',[0.1+0.1 0.27+0.12 0.55*23/24*.8  0.08])
plot(mEYEV,'k-','LineWidth',2)
hold on
plot(mEYEV+sdEYEV,'k-','LineWidth',1)
plot(mEYEV-sdEYEV,'k-','LineWidth',1)
plot([140 140],[-2 2],'b:','LineWidth',1)
plot([160 160],[-2 2],'g:','LineWidth',1)
plot([200 200],[-2 2],'k:','LineWidth',1)
for k=0:11
    sti=201+30*k;
    plot([sti,sti],[-2 2],'k:','Linewidth',1)
end
xlim([121 580])
ylim([-1.5 1.5])
set(gca,'XTick',[140 160,200,560],'XTickLabel',{'','','','',''},'TickDir','out','box','off','Linewidth',1.5)
set(gca,'ytick',[-1 -0.5 0 0.5 1],'yticklabel',{'-1','','0','','1' },'TickDir','out','FontSize',fs,'FontName','Arial narrow','box','off','Linewidth',1.5)

