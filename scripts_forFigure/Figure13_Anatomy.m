%%
%this code draw figure 13 B
%%
clear
close all
load spine
addpath(genpath('../toolboxes'))
%%
fs=12;
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
%%
%load MRI data (resolution:1/200mm per voxel)
V2=[];
Precuneus=[];
POS_point_all=[];
load ../Data_open/data.mat
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100]=get_Entropy50;
%%
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);
%%
for i=1:length(data_bkg1)
 dep(i)= data_bkg1(i).dep;  
 ap(i)=data_bkg1(i).ap; 
 lateral(i)=abs(data_bkg1(i).xyz(1))-1; 
 monkey_vec(i)=data_bkg1(i).monkey; 
end
lateral= lateral+1;
p_thre=0.004;
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
load ../Data_open/anatomy
for i=1:length(data_bkg1)
    data_bkg1(i).dep=Anat.dep(i);
end
N_Occi=Anat.area'==1;
N_V6a=Anat.area'==2;
N_PEc=Anat.area'==3;
N_7M=Anat.area'==4;
N_MIP=Anat.area'==5;
%%
close all
for info =1:4
switch info
    case 1,Neuron_flt=Neuron_information(:,1);
    case 2,Neuron_flt=Neuron_information(:,2);
    case 3,Neuron_flt=Neuron_information(:,3);
    case 4,Neuron_flt=true(size(Neuron_information,1),1);
end

load ../Data_MRI/M1R
 for lr=2:2
     fig_1=figure('position',[100 100 floor(16000*asp_rat_x/res_x) floor(1.5*16000*asp_rat_y/res_y)]);
     image_use=squeeze(image_use_all(lr,:,:));
     imagesc(image_use,[2000,16000]);
     colormap(map)
     xli=[850 2850];yli=[500 3500];
     xlim(xli)
     ylim(yli)
     hold on
     set(gca,'xtick',[],'ytick',[],'box','off')
for monkey=[4 4.5 2 5]
mon=monkey;
save_V4_precuneus=1;
if save_V4_precuneus
        dura_point=squeeze(dura_point_all(lr,:,:));
        image_use=squeeze(image_use_all(lr,:,:));
        hold on
        switch monkey
            case 2,theta=pi*(55/180);
            case 4,theta=pi*(62/180);
            case 4.5,theta=pi*(62/180);
            case 5,theta=pi*(59/180);
        end
        for i=1:length(data_bkg1)
            if data_bkg1(i).monkey==monkey&&Neuron_flt(i)%&&(lateral(i)==lr)
                theta(i)=pi*data_bkg1(i).angle/180;
                cell_xy=floor(20*rand(1,2))-10+dura_point(data_bkg1(i).ap+8+1,:)+data_bkg1(i).dep*200*[(-1)*cos(theta(i)) sin(theta(i))];
                switch info
                    case 1, plot(cell_xy(1),cell_xy(2),'o','MarkerEdgeColor','m','Markersize',2)
                    case 2, plot(cell_xy(1),cell_xy(2),'x','MarkerEdgeColor','b','Markersize',2)
                    case 3, plot(cell_xy(1),cell_xy(2),'^','MarkerEdgeColor','g','Markersize',2)
                    case 4, plot(cell_xy(1),cell_xy(2),'.','MarkerEdgeColor','r','Markersize',2)
                end
            end
        end
    end
end
end
end
pause(2)
drawnow

if save_V4_precuneus&&0
cd figs
        figure(1),print('sag_dot_bkg.ai', '-dpdf', '-painters');    
        figure(1+1),print('sag_dot_ret.ai', '-dpdf', '-painters');
        figure(2+1),print('sag_bkg_ret.ai', '-dpdf', '-painters');
        figure(3+1),print('sag_all.ai', '-dpdf', '-painters');
cd ../
end

%%
%12mm*15MM
fig_2=figure('position',[100 100 floor(1.2*16000*asp_rat_x/res_x) floor(1.5*16000*asp_rat_y/res_y)]);
lr=2;
theta=pi*(62/180);
dura_point=squeeze(dura_point_all(lr,:,:));
image_use=squeeze(image_use_all(lr,:,:));
imagesc(image_use,[2000,16000]);
colormap(map)
xlim(xli)
ylim(yli)
hold on
set(gca,'xtick',[],'ytick',[],'box','off')

for i=1:length(data_bkg1)
   ap(i)=data_bkg1(i).ap; 
   monkey_vec(i)=data_bkg1(i).monkey; 
end

dep_edge=-0.0001:1:20; 
for i=1:8
dep_count(9-i,:)=histcounts(dep(ap==-1*i),dep_edge);
end

for j=1:8
    for i_dep=1:10
        cell_cite=squeeze(dura_point(j,:))+200*(i_dep-1)*[(-1)*cos(theta) sin(theta)];
        if dep_count(j,i_dep)>2
            plot(cell_cite(1),cell_cite(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',log(dep_count(j,i_dep)))
        elseif (dep_count(j,i_dep)<=2)&&(dep_count(j,i_dep)>0)
            plot(cell_cite(1),cell_cite(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',1)
        end
    end
end

lo=[10 20 50 100];
for i=1:4
plot(500+500*i+1,600,'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',log(lo(i)))
end  

plot([800,800],[3400,600],'k-','LineWidth',1)
plot([650,800],[1400,1400],'k-','LineWidth',1)
plot([650,800],[3400,3400],'k-','LineWidth',1)

for i=-2:4
plot([725,800],[1400+i*400,1400+i*400],'k-','LineWidth',1)
end  

if 0
cd  figs
print('fig1_sag_all.ai', '-dpdf', '-painters');
cd ../
end

%%
fig_3=figure('position',[100 100 floor(16000*asp_rat_x/res_x)*4 floor(1.5*16000*asp_rat_y/res_y)]);

N_dot_ret=Neuron_information(:,2);
N_dot_bkg=Neuron_information(:,1);
N_bkg_ret=Neuron_information(:,3);
N_Non=sum(Neuron_information,2)==0;
Num_allArea=[sum(N_Non),sum(N_bkg_ret),sum(N_dot_ret),sum(N_dot_bkg)];

subplot(1,4,1)
Num_Occi=[sum(N_Occi'&N_Non),sum(N_Occi'&N_bkg_ret),sum(N_Occi'&N_dot_ret),sum(N_Occi'&N_dot_bkg)]/sum(N_Occi);
b=barh(Num_Occi,'FaceColor','w','Linewidth',1.5);
b.FaceColor = 'flat';b.CData(1,:) = [1 1 1];b.CData(2,:) = [0 1 0];b.CData(3,:) = [0 0 1];b.CData(4,:) = [1 0 1];
set(gca,'ytick',[],'xtick',[0 0.5 1],'xticklabel',{'0','50','100'},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial')
xlim([0 1])

subplot(1,4,2)
Num_V6a=[sum(N_V6a'&N_Non),sum(N_V6a'&N_bkg_ret),sum(N_V6a'&N_dot_ret),sum(N_V6a'&N_dot_bkg)]/sum(N_V6a);
b=barh(Num_V6a,'FaceColor','w','Linewidth',1.5);
b.FaceColor = 'flat';b.CData(1,:) = [1 1 1];b.CData(2,:) = [0 1 0];b.CData(3,:) = [0 0 1];b.CData(4,:) = [1 0 1];
set(gca,'ytick',[],'xtick',[0 0.5 1],'xticklabel',{'0','50','100'},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial')
xlim([0 1])

subplot(1,4,3)
Num_PEc=[sum(N_PEc'&N_Non),sum(N_PEc'&N_bkg_ret),sum(N_PEc'&N_dot_ret),sum(N_PEc'&N_dot_bkg)]/sum(N_PEc);
b=barh(Num_PEc,'FaceColor','w','Linewidth',1.5);
b.FaceColor = 'flat';b.CData(1,:) = [1 1 1];b.CData(2,:) = [0 1 0];b.CData(3,:) = [0 0 1];b.CData(4,:) = [1 0 1];
set(gca,'ytick',[],'xtick',[0 0.5 1],'xticklabel',{'0','50','100'},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial')
xlim([0 1])

subplot(1,4,4)
Num_7M=[sum(N_7M'&N_Non),sum(N_7M'&N_bkg_ret),sum(N_7M'&N_dot_ret),sum(N_7M'&N_dot_bkg)]/sum(N_7M);
b=barh(Num_7M,'FaceColor','w','Linewidth',1.5);
b.FaceColor = 'flat';b.CData(1,:) = [1 1 1];b.CData(2,:) = [0 1 0];b.CData(3,:) = [0 0 1];b.CData(4,:) = [1 0 1];
set(gca,'ytick',[],'xtick',[0 0.5 1],'xticklabel',{'0','50','100'},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial')
xlim([0 1])

%%
if 0
cd figs
print('bar_area_infor.ai', '-dpdf', '-painters');    
cd ../../Scripts_forFigureR1/
end

