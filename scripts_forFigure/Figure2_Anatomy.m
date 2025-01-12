%%
%this code draw figure 1B(recording location) and Figure 2
%%
clear
close all
load spine
fs=12;
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
addpath(genpath('../toolboxes'))
col={'k','b','g','m','y'};
mark={'^','x','d','o','*'};
%%
load ../Data_open/data.mat
%%
%get position of recording 
V2=[];
Precuneus=[];
POS_point_all=[];
x=reshape([data_bkg1(:).xyz],3,942);
c=[data_bkg1(:).ap]>=-3&x(1,:)==3&x(3,:)<-1.9;
for i=1:length(data_bkg1)
    if c(i)
    data_bkg1(i).xyz(1)=2;
    end
end
for i=1:942
 lateral(i)=abs(data_bkg1(i).xyz(1))-1; 
 monkey_vec(i)=data_bkg1(i).monkey; 
end
%%
%load label of brain area of each cell
load ../Data_open/anatomy
ANAT=Anat.area;
lateral=Anat.X;
%%
%load MRI data (1/200mm per voxel) 
%choose monkey MRI to draw
monkey=2;
switch monkey
    case 2
        load ../Data_MRI/M2L
    case 4
        load ../Data_MRI/M1R
    case 4.5
        load ../Data_MRI/M1L
    case 5
        load ../Data_MRI/M3L
end
%%
for lr=1:size(dura_point_all,1)
    fig_1=figure('position',[100 100 floor(16000*asp_rat_x/res_x) floor(1.5*16000*asp_rat_y/res_y)]);
    dura_point=squeeze(dura_point_all(lr,:,:));
    image_use=squeeze(image_use_all(lr,:,:));
    imagesc(image_use,[2000,16000]);
    colormap(map)
    hold on
    switch monkey
        case 2
            imagesc(image_use,[150,380]);
            xlim([1550 3550])
            ylim([1300 4300])
            theta=pi*(55/180);
        case 4
            xlim([850 2850])
            ylim([500 3500])
            theta=pi*(62/180);
        case 4.5
            xlim([1250 3250])
            ylim([500 3500])
            theta=pi*(62/180);
        case 5
            xlim([1600 3600])
            ylim([0 3000])
            theta=pi*(59/180);
    end
end
file_tmp=0;
for i=1:length(data_bkg1)
    file_N=str2num(data_bkg1(i).file(4:6));
    if data_bkg1(i).monkey==monkey
        if file_tmp==file_N
        else
            switch monkey
                case 4, lr=lateral(i);
                case 4.5,lr=lateral(i)+1;%Image exist freom X=1mm
                case 5,lr=lateral(i)+1;%Image exist freom X=1mm
                case 2,lr=lateral(i);
            end
            figure(lr)
            dura_point=squeeze(dura_point_all(lr,:,:));
            theta(i)=pi*data_bkg1(i).angle/180;
            cell_xy=dura_point(data_bkg1(i).ap+8+1,:)+(data_bkg1(i).dep+rand(1)*0.1)*200*[(-1)*cos(theta(i)) sin(theta(i))]+15*rand(1,2)-7.5;
            plot(cell_xy(1),cell_xy(2),mark{ANAT(i)},'MarkerEdgeColor',col{ANAT(i)},'Markersize',3)
        end
    end
end
%%
%save figure
for lr=1:size(dura_point_all,1)
    figure(lr)
    set(gca,'xtick',[],'ytick',[],'box','of')
    if 0
        print(['figs/MRIs/figMRI',num2str(monkey*10),num2str(lr),'.ai'], '-dpdf', '-painters')
    end
end

