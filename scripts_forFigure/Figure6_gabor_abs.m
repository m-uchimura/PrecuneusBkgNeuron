%%
%this file draw the results of gabor fitting in Figure 6
%%
clear
close all
addpath(genpath('../toolboxes'))
%%
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
frame_size_x=1;
frame_size_y=1;
figure('position',[0 0 floor(100000*asp_rat_x/res_x) floor(100000*asp_rat_y/res_y)])
rng(0)
load ../Data_open/gabor_bkg_abs_each_0005
load ../Data_open/gabor_retino_abs_each_0005
p_thre=0.004;
%%
load ../Data_open/data.mat
for i=1:58
    for k=1:length(data_bkg1)
    if strcmp(gabor_bkg(i).file,data_bkg1(k).file)
        A(i)=k;
        gabor_bkg(i).i=k;
    end
    end
end

for i=1:149
    for k=1:length(data_bkg1)
    if strcmp(gabor_retino(i).file,data_bkg1(k).file)
        A(i)=k;
        gabor_retino(i).i=k;
    end
    end
end
load ../Data_open/anatomy

%%
%load neual information, Pvalue etc
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
I_SHUFF(~N_bkg100,:,:,1)=20*I_SHUFF(~N_bkg100,:,:,1);
I_SHUFF(N_bkg100,:,:,1)=10*I_SHUFF(N_bkg100,:,:,1);

RF_dot_bkg =Spike_dot_bkg_filter;
RF_dot_ret =Spike_dot_ret_filter;
%%
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);

%%
%examples of dot/ret(1-4) and dot/bkg(5-8) neurons
file_cand{1}='m4c048r1s1n1';data_bin_cand(1)=14; 
file_cand{2}='m4c162r1s1n1';data_bin_cand(2)=9;
file_cand{3}='m4c166r1s1n8';data_bin_cand(3)=9;
file_cand{4}='m4c144r1s1n12';data_bin_cand(4)=13; 

file_cand{5}='m4c038r1s1n1';data_bin_cand(6-1)=10; 
file_cand{6}='m4c169r1s1n3';data_bin_cand(7-1)=13; 
file_cand{7}='m5c097r1s1n2';data_bin_cand(8-1)=9;
file_cand{8}='m5c139r1s1n6';data_bin_cand(5+3)=20;

colbar_leg_cell={[0 5 10 15]*2,[10 20 30]/2,[10 20 30],[0 2 4]...
                  [3 7 11],[2 4 6],[0 10 20],[0 2 4]};
for i=1:8
for j=1:length(data_bkg1)
    if strcmp(data_bkg1(j).file,file_cand{i})
        file_N(i)=j;
    end
end
end
%%
%draw receptie fields of data (left) and gabor fitting (right)
%for dot/ret (Figure 6A) and dot/bkg (Figure 6B)
for i=0:7
data_bin=data_bin_cand(i+1);
if i<4
    RF_example=RF_dot_ret{file_N(i+1)}{data_bin};
else
    RF_example=RF_dot_bkg{file_N(i+1)}{data_bin};
end

if file_N(i+1)>495
    RF_example=flipud(RF_example);
end
data1=RF_example';
if i<4
    N_ret=find([gabor_retino(:).i]==file_N(i+1));
    a=gabor_retino(N_ret).a;
    sigma=gabor_retino(N_ret).sigma;
    x0=gabor_retino(N_ret).x0;
    y0=gabor_retino(N_ret).y0;
    theta=gabor_retino(N_ret).theta;
    lambda= gabor_retino(N_ret).lambda;
    phase=gabor_retino(N_ret).phase;
    b=gabor_retino(N_ret).b;
    dc=gabor_retino(N_ret).dc;
    gannma=gabor_retino(N_ret).gannma;
    r=gabor_retino(N_ret).r;
    prob=gabor_retino(N_ret).prob;
    pred=gabor_retino(N_ret).pred;
else
    N_ret=find([gabor_bkg(:).i]==file_N(i+1));
    a=gabor_bkg(N_ret).a;
    sigma=gabor_bkg(N_ret).sigma;
    x0=gabor_bkg(N_ret).x0;
    y0=gabor_bkg(N_ret).y0;
    theta=gabor_bkg(N_ret).theta;
    lambda= gabor_bkg(N_ret).lambda;
    phase=gabor_bkg(N_ret).phase;
    b=gabor_bkg(N_ret).b;
    dc=gabor_bkg(N_ret).dc;
    gannma=gabor_bkg(N_ret).gannma;
    r=gabor_bkg(N_ret).r;
    prob=gabor_bkg(N_ret).prob;
    pred=gabor_bkg(N_ret).pred;
end
p=[a sigma x0 y0 theta lambda phase b gannma];
s=size(data1);
da_raw=reshape(data1,[1 s(1)*s(2)]);
data_fit=data1;
s=size(data_fit);
dat=reshape(data_fit,[1 s(1)*s(2)]);
dat=(dat-nanmean(dat))/nanstd(dat);
data_fit=reshape(dat,s);
[x, y]=meshgrid(1:length(data_fit(1,:)), 1:length(data_fit(:,1)));

center(:,i+1)=[p(3);p(4)];
for jj=1:2   
    subplot('position',[0.1+0.2*(jj-1)+(0.45+0.01)*floor(i/4) 0.8-mod(i,4)*0.135+.03 0.18 0.12]) 
    switch jj
        case 1
            color_min_bin_cand1=nanmin(nanmin(data_fit));
            color_max_bin_cand1=nanmax(nanmax(data_fit));
            imagesc(data_fit,[color_min_bin_cand1,color_max_bin_cand1])
        case 2
            color_min_bin_cand=nanmin(nanmin(pred));
            color_max_bin_cand=nanmax(nanmax(pred));
            imagesc(pred,[color_min_bin_cand,color_max_bin_cand])
            hold on
            if ismember(i,[3,7])
                text(30,5,['DC=',num2str(dc,2)],'Color','k')
            else
                text(30,5,['DC=',num2str(dc,2)],'Color','w')
            end
             draw_ellipse(p(2),p(2)/p(9),-p(5),p(3),p(4),'w')
    end
    if i<4
        rectangle('Position',[30 20 2 2] ,'EdgeColor','w','Curvature',[1 1],'LineWidth',2)
    else
        rectangle('Position',[16 11 30*frame_size_x 20*frame_size_y],'EdgeColor','w','LineWidth',2)
    end

    set(gca,'ytick',[1 10 20 30 40],'yticklabel',{'','','','',''},'TickDir','out','box','off','LineWidth',1.5)
    if mod(i,4)==3
    set(gca,'xtick',[10 30 50],'xticklabel',{'','',''},'TickDir','out','box','off','Fontsize',15,'FontAngle','normal','Linewidth',1.5)
    else
    set(gca,'xtick',[10 30 50],'xticklabel',{'','',''},'TickDir','out','box','off','Fontsize',15,'Linewidth',1.5)
    end

    if jj==2
        subplot('position',[0.1+0.2+(0.45+0.01)*floor(i/4)+0.18 0.8-mod(i,4)*0.135-0.005+.03 0.045 0.13])
        axis off
        color_max_bin_cand=color_max_bin_cand1*nanstd(da_raw)+nanmean(da_raw);
        color_min_bin_cand=color_min_bin_cand1*nanstd(da_raw)+nanmean(da_raw);
        colbar_leg=colbar_leg_cell{1+i};%i
        col_diff=10*(color_max_bin_cand-color_min_bin_cand);
        col_min=color_min_bin_cand*10;
        col_max=color_max_bin_cand*10;
        colorbar('west','Ticks',[(colbar_leg(1)-col_min)/col_diff (colbar_leg(2)-col_min)/col_diff (colbar_leg(3)-col_min)/col_diff],...
            'TickLabels',{num2str(colbar_leg(1)),num2str(colbar_leg(2)),num2str(colbar_leg(3))},'Linewidth',1.5,'Fontsize',12,'FontName','Arial narrow')
    end
end
end
%%
%extract the paremeters used in Figure 6 C-F
dc_thre=0.5;
sigma_retino=[gabor_retino(:).sigma];
gannma_retino=[gabor_retino(:).gannma];
sigma_retino_c=sigma_retino./sqrt(gannma_retino);
sigma_bkg=[gabor_bkg(:).sigma];
gannma_bkg=[gabor_bkg(:).gannma];
sigma_bkg_c=sigma_bkg./sqrt(gannma_bkg);

ecc_retino=[gabor_retino(:).ecc];
ecc_bkg=[gabor_bkg(:).ecc];
retino_reg=[ecc_retino;ones(1,length(ecc_retino))]';
bkg_reg=[ecc_bkg;ones(1,length(ecc_bkg))]';

dc_flt_ret=[gabor_retino(:).dc]>dc_thre;
dc_flt_bkg=[gabor_bkg(:).dc]>dc_thre;

iNormalize=0;
switch iNormalize 
    case 0,[~,~,R_retino,~,~]=regress(sigma_retino_c(dc_flt_ret)',retino_reg(dc_flt_ret,:));
        [~,~,R_bkg,~,~]=regress(sigma_bkg_c(dc_flt_bkg)',bkg_reg(dc_flt_bkg,:));
         R_ret=normalize(R_retino);R_bk=normalize(R_bkg);
    case 1,[~,~,R_retino,~,~]=regress(sigma_retino_c',retino_reg);
           [~,~,R_bkg,~,~]=regress(sigma_bkg_c',bkg_reg);
         R_ret=normalize(R_retino(dc_flt_ret));R_bk=normalize(R_bkg(dc_flt_bkg));
    case 2,R_ret=normalize(sigma_retino_c(dc_flt_ret));R_bk=normalize(sigma_bkg_c(dc_flt_bkg));
end
N_flt_r=false(size(dc_flt_ret));
N_flt_r(dc_flt_ret)=abs(R_ret')<3;
N_flt_b=false(size(dc_flt_bkg));
N_flt_b(dc_flt_bkg)=abs(R_bk')<3;

%%
%Draw all receptive fields of dot/ret neurons (Figure 6C)
subplot('position',[0.1 0.2-.03+.05 0.18 0.12])
dc_retino=[gabor_retino(:).dc];
x0_retino=[gabor_retino(:).x0];
y0_retino=[gabor_retino(:).y0];
sigma_retino=[gabor_retino(:).sigma];
theta_retino=[gabor_retino(:).theta];
gannma_retino=[gabor_retino(:).gannma];
ecc_retino=[gabor_retino(:).ecc];
prob_retino=[gabor_retino(:).prob];
COLO_RET={'b','c','y'};
for i=1:length(y0_retino)
    if ~isnan(sigma_retino(i))&&(dc_retino(i)>dc_thre)&&N_flt_r(i)
        draw_ellipse(sigma_retino(i),sigma_retino(i)/gannma_retino(i),theta_retino(i),x0_retino(i),(41-y0_retino(i)),'b')
        hold on
    end
end
rectangle('Position',[30 20 2 2] ,'EdgeColor','k','Curvature',[1 1],'LineWidth',2)
xlim([0 61])
ylim([0 41])
set(gca,'xtick',[10 30 50],'xticklabel',{'','',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
set(gca,'ytick',[1 10 20 30 40],'yticklabel',{'','','','',''},'TickDir','out','box','off','LineWidth',2)

%Draw all receptive fields of dot/bkg neurons(Figure 6E)
subplot('position',[0.56 0.2-.03+.05 0.18 0.12])
dc_bkg=[gabor_bkg(:).dc];
x0_bkg=[gabor_bkg(:).x0];
y0_bkg=[gabor_bkg(:).y0];
sigma_bkg=[gabor_bkg(:).sigma];
theta_bkg=[gabor_bkg(:).theta];
gannma_bkg=[gabor_bkg(:).gannma];
ecc_bkg=[gabor_bkg(:).ecc];
for i=1:length(y0_bkg)
    if ~isnan(sigma_bkg(i))&&(dc_bkg(i)>dc_thre)&&N_flt_b(i)
        draw_ellipse(sigma_bkg(i),sigma_bkg(i)/gannma_bkg(i),theta_bkg(i),x0_bkg(i),(41-y0_bkg(i)),'m')
        hold on
    end
end
rectangle('Position',[16 11 30*frame_size_x 20*frame_size_y],'EdgeColor','k','LineWidth',2)
xlim([0 61])
ylim([0 41])
set(gca,'xtick',[10 30 50],'xticklabel',{'','',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
set(gca,'ytick',[1 10 20 30 40],'yticklabel',{'','','','',''},'TickDir','out','box','off','LineWidth',2)

%%
%Draw relathinship beteen the eccentricities and size of gabor
%dot/ret neurons (Figure 6D)
subplot('position',[0.3+0.03 0.2-.03+.05 0.15 0.12])
sigma_retino_c_all=sigma_retino./sqrt(gannma_retino);
sigma_bkg_c_all=sigma_bkg./sqrt(gannma_bkg);

[R,P]=corrcoef(ecc_retino,sigma_retino_c);
plot(ecc_retino(dc_retino>dc_thre&N_flt_r),sigma_retino_c(dc_retino>dc_thre&N_flt_r),'bo')
hold on
retino_reg=[ecc_retino(dc_retino>dc_thre&N_flt_r);ones(1,length(ecc_retino(dc_retino>dc_thre&N_flt_r)))]';
[X_retino,BINT_retino,R_retino,RINT_retino,STATS_retino]=regress(sigma_retino_c(dc_retino>dc_thre&N_flt_r)',retino_reg);
y_retino=retino_reg*X_retino;
plot([ecc_retino(dc_retino>dc_thre&N_flt_r) 0],[y_retino; X_retino(2)],'b-','Linewidth',3)
plot(ecc_retino(dc_retino>dc_thre&~N_flt_r),sigma_retino_c(dc_retino>dc_thre&~N_flt_r),'b+')
find(dc_retino<dc_thre&~N_flt_r);

xlim([0 1+max(nanmax(ecc_bkg),nanmax(ecc_retino))])
ylim([0 40])
text(1,38,[' y=',num2str(X_retino(1),2),'x+',num2str(X_retino(2),2)])
text(1,32,[' R=',num2str(sqrt(STATS_retino(1)),2),' p=',num2str(STATS_retino(3),2)])
set(gca,'xtick',[0 20 40 60],'xticklabel',{'','','',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
set(gca,'ytick',[20 40],'yticklabel',{'',''})

%dot/bkg neurons(Figure 6F)
subplot('position',[0.76+0.03 0.2-.03+.05 0.15 0.12])
[R,P]=corrcoef(ecc_bkg,sigma_bkg_c);
plot(ecc_bkg(dc_bkg>dc_thre&N_flt_b),sigma_bkg_c(dc_bkg>dc_thre&N_flt_b),'mo')
hold on
frame_reg=[ecc_bkg(dc_bkg>dc_thre&N_flt_b);ones(1,length(ecc_bkg(dc_bkg>dc_thre&N_flt_b)))]';
[X_frame,BINT_frame,R_frame,RINT_frame,STATS_frame]=regress(sigma_bkg_c(dc_bkg>dc_thre&N_flt_b)',frame_reg);
y_frame=frame_reg*X_frame;
plot([ecc_bkg(dc_bkg>dc_thre&N_flt_b) 0],[y_frame;X_frame(2)] ,'m-','Linewidth',3)
plot(ecc_bkg(dc_bkg>dc_thre&~N_flt_b),sigma_bkg_c(dc_bkg>dc_thre&~N_flt_b),'m+')

set(gca,'xtick',[0 20 40 60],'xticklabel',{'','',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
set(gca,'ytick',[20 40],'yticklabel',{'',''})
xlim([0 1+max(nanmax(ecc_bkg),nanmax(ecc_retino))])
ylim([0 40])
text(1,38,[' y=',num2str(X_frame(1),2),'x+',num2str(X_frame(2),2)])
text(1,32,[' R=',num2str(sqrt(STATS_frame(1)),2),' p=',num2str(STATS_frame(3),2)])

%%
%draw the relationship between recording position and RF position(Figs 6G,I) or RF size(Figs 6H,J)
draw_gabor_Anatomy(data_bkg1,gabor_bkg,gabor_retino)

%%
if 0
cd figs
print('fig6_gabor.ai', '-dpdf', '-painters','-bestfit')
cd ../../Scripts_forFigureR1/
end