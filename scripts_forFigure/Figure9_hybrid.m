%%
%this file draw the resuls hybrid coding(figure 9)
clear
close all
addpath(genpath('../toolboxes'))
%%
load ../Data_open/hybrid_all.mat 
%%
%draw distribution Figure 9A
tc=-1:.1:.5;
col={'m','b','g','k'};
SigPeriod=permute(Entropy_thre,[3,1,2])>0;
SigPeriod(4,:,:)=sum(SigPeriod)==0;
R_mean=nanmean(R_co,4);

si=[120,50,10,1]*2;
figure(1)
for i=[1 8 15 29]
figure(i)
switch i
     case 1,set(gcf,'Position',[100 100 280 280])
     case 8,set(gcf,'Position',[100 100 520 520])
     case 15,set(gcf,'Position',[100 100 520 520])
     case 29,set(gcf,'Position',[100 100 360 360])
end 
hold on
for c=3:-1:1
switch c
     case 1,R_m1=R_mean(:,:,i).*repmat(SigPeriod(1,:,i),2,1);
     case 2,R_m1=R_mean(:,:,i).*repmat(SigPeriod(2,:,i),2,1);
     case 3,R_m1=R_mean(:,:,i).*repmat(SigPeriod(3,:,i),2,1);
     case 4,R_m1=R_mean(:,:,i).*repmat(SigPeriod(4,:,i),2,1);
end
R_m1(R_m1==0)=NaN;
if c==4,scatter(R_m1(1,:),R_m1(2,:),80,'k.'),else
scatter(R_m1(1,:),R_m1(2,:),si(c),'MarkerEdgeColor',col{c},'LineWidth',2)
end
hi(c,:)=histcounts(R_m1(1,:)-R_m1(2,:),tc);
A=R_m1(1,:)-R_m1(2,:);
N(i,c)=sum(~isnan(sum(R_m1)));
end
plot([-.6 1],[-.6 1],'k:','LineWidth',3)
plot([-.6 1],[0 0],'k:','LineWidth',2)
plot([0 0],[-.6 1],'k:','LineWidth',2)
switch i
    case 1,xlim([-.3 .4]);ylim([-.3 .4]);
    case 8,xlim([-.3 1]);ylim([-.3 1]);
    case 15,xlim([-.3 1]);ylim([-.3 1]);
    case 29,xlim([-.3 .6]);ylim([-.3 .6]);
end
axis square
set(gca,'TickDir','out','Box','off','LineWidth',2)
set(gca,'xtick',[-.3 0 .3 .6 .9],'XTickLabel',{},'ytick',[-.3 0 .3 .6 .9],'yTickLabel',{})

if 0
switch i
    case 1,print('figs/hyb_plot1.ai', '-dpdf')
    case 8,print('figs/hyb_plot8.ai', '-dpdf')
    case 15,print('figs/hyb_plot15.ai', '-dpdf')
    case 29,print('figs/hyb_plot29.ai', '-dpdf')
end
end

if i==15
figure(i+1)
set(gcf,'Position',[100 100 round(560/sqrt(2)) 150])
ed=tc;
ce=ed+.05;ce(end)=[];
hold on
bar(ce,hi(2,:),0.9,'b','FaceAlpha',.5)
bar(ce,hi(1,:),0.9,'m','FaceAlpha',.5)
set(gca,'TickDir','out','Box','off','LineWidth',2)
set(gca,'xtick',[-.8 -.4 0 .4],'XTickLabel',{},'ytick',[0 10 20 30],'yTickLabel',{})
ylim([0 20])

if 0
switch i
    case 1,print('figs/hyb_hist1.ai', '-dpdf')
    case 8,print('figs/hyb_hist8.ai', '-dpdf')
    case 15,print('figs/hyb_hist15.ai', '-dpdf')
    case 29,print('figs/hyb_hist29.ai', '-dpdf')
end
end

end
end
%%
%draw time course Figure 9B
close all
clear 
load ../Data_open/hybrid_all.mat
col={'m','b','g'};

R_mean=nanmean(R_co,4);
si=[120,50,10];
SigPeriod=permute(Entropy_thre,[3,1,2])>0;

i=15;
info_15=SigPeriod(:,:,15);
flt_bkg=info_15(1,:);flt_ret=info_15(2,:);flt_bkg_ret=info_15(3,:);

A_bkg=R_mean(1,flt_bkg,i)-R_mean(2,flt_bkg,i);
A_ret=R_mean(1,flt_ret,i)-R_mean(2,flt_ret,i);
T_bkg=find(A_bkg>0);
T_ret=find(A_ret<0);
R_bkg=R_mean(:,flt_bkg,:);
R_bkg_ret=R_mean(:,flt_bkg_ret,:);
R_ret=R_mean(:,flt_ret,:);
RR_bkg=R_bkg(:,T_bkg,:);
RR_ret=R_ret(:,T_ret,:);
%%
ho=colormap(jet);
close all
theta = linspace(0, 2*pi, 100);

for j=1:3
figure(j);
set(gcf,'Position',[100 100 340 340])
switch j
    case 1,RRR=RR_bkg;
    case 2,RRR=RR_ret;
    case 3,RRR=R_bkg_ret;
end
N_co(j)=size(RRR,2);
for i=1:35
% Plot the original data
A=RRR(:,:,i);
hold on;
plot([-.2 1],[-.2 1],'k:','LineWidth',1)
plot([-.2 1],[0 0],'k:','LineWidth',1)
plot([0 0],[-.2 1],'k:','LineWidth',1)
% Calculate mean and covariance of the data
mu = mean(A, 2);
Sigma = cov(A');
% Calculate the eigenvalues and eigenvectors of the covariance matrix
[eigenVec, eigenVal] = eig(Sigma);
ellipsePath = (eigenVec * (sqrt(eigenVal)./sqrt(size(A,2)-1))) * [cos(theta); sin(theta)] + mu;
fill(ellipsePath(1,:), ellipsePath(2,:),ho(floor(i*256/35),:),'EdgeColor','none', 'FaceAlpha', 0.3);
end
for i=1:35
plot(squeeze(mean(RRR(1,:,i:i+1),2)),squeeze(mean(RRR(2,:,i:i+1),2)),'-','Color','k','LineWidth',1)
end

ma={'square','o','d','^'};
l=0;
for i=[1 8 15 29]
l=l+1;
plot(squeeze(mean(RRR(1,:,i),2)),squeeze(mean(RRR(2,:,i),2)),['k',ma{l}],'MarkerFaceColor','w','LineWidth',1,'MarkerSize',15)
end

xlim([-.1 .61])
ylim([-.1 .61])
axis square
set(gca,'TickDir','out','Box','off','LineWidth',1)
set(gca,'xtick',[-.3 0 .3 .6 ],'XTickLabel',{},'ytick',[-.3 0 .3 .6],'yTickLabel',{})

%
if 0
switch j
    case 1,print('figs/trj_bkg.ai', '-dpdf')
    case 2,print('figs/trj_ret.ai', '-dpdf')
    case 3,print('figs/trj_bkgret.ai', '-dpdf')
end
end
end