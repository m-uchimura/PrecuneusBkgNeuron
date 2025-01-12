%%
%this file draw the resuls of neurons which represent all of the dot/ret, dot/bkg,
%and bkg/ret neurons for Figure12
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
%%
%load neual information, Pvalue etc
[Entropy,P_CHIV,P_CHIV_dotsh,Spike_dot_bkg_filter,Spike_dot_ret_filter,spike_peridot_full,bkgpos_map,I_SHUFF,N_dotbkg_100,N_dotret_100,N_bkgret_100,N_dotbkg_50,N_dotret_50]=get_Entropy50;
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
P_CHIV_all=P_CHIV;
P_CHIV_dotsh_all=P_CHIV_dotsh;
%%
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);
%%
%set entropy zerop if p-value was abobe the threshold 

Entropy_thre_tmp(:,:,1:2)=(P_CHIV(:,:,1:2)<=p_thre).*(P_CHIV_dotsh(:,:,1:2)<=p_thre).*Entropy(:,:,1:2);
Entropy_thre_tmp(:,:,3)=(P_CHIV(:,:,3)<=p_thre).*Entropy(:,:,3);
for i=1:size(Entropy_thre_tmp,1)
    Entropy_thre(i,:,1)=Suc5(Entropy_thre_tmp(i,:,1));
    Entropy_thre(i,:,2)=Suc5(Entropy_thre_tmp(i,:,2));
    Entropy_thre(i,:,3)=Suc5(Entropy_thre_tmp(i,:,3));
end

amp=10;
ent_shuff_mean=squeeze(nanmean(I_SHUFF,2));
entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0).*(Entropy_thre(:,:,2)>0);
entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
Entropy_thre(:,:,1)=entropy_thre_minus_spread;
%%
use_co=1;  %1:dot/bkg 2;dot/ret
use_bkg=sum(Entropy_thre(:,:,use_co),2)>0;%dot/ret
sum(use_bkg)
N_dotbkg50=N_dotbkg_50;
dotbkg_50=ismember(N_dotbkg50,find(use_bkg));
N_bkg_50=N_dotbkg50(dotbkg_50);
N_bkg_100=find(N_bkg100.*use_bkg);

use_co=2;  %1:dot/bkg 2;dot/ret
use_ret=sum(Entropy_thre(:,:,use_co),2)>0;%dot/ret
sum(use_ret)
N_dotret50=N_dotret_50;

dotret_50=ismember(N_dotret50,find(use_ret));
N_ret_50=N_dotret50(dotret_50);
N_ret_100=find(N_ret100.*use_ret);


%%
%draw the timecourse of information which calulated by only 1st, 7th and
%12th dot (Figure 12 A,C)
figure('position',[0 0 floor(100000*asp_rat_x/res_x) floor(100000*asp_rat_y/res_y)])
for coord=1:2
switch coord
case 1
load ../Data_open/Entropy_data_binstart_0dotbkg_100
Entropy100=Entropy;
I_SHUFF100=I_SHUFF;
load ../Data_open/Entropy_data_binstart_0dotbkg_50
I_SHUFF50=I_SHUFF;
Entropy50=Entropy;
case 2
load ../Data_open/Entropy_data_binstart_0dotret_100
Entropy100=Entropy;
load ../Data_open/Entropy_data_binstart_0dotret_50
Entropy50=Entropy;
end
clear Entropy
Entropy=zeros(size(Entropy50));
I_SHUFF=zeros(size(I_SHUFF));
switch coord
    case 1,N_all=unique([N_bkg_100' N_bkg_50]);
        N_50=ismember(N_all,N_bkg_50);
        N_100=ismember(N_all,N_bkg_100);
        Entropy(N_50,:,:,:)=20*Entropy50(N_50,:,:,:);
        Entropy(N_100,:,:,:)=10*Entropy100(N_100,:,:,:);
        I_SHUFF(N_50,:,:,:,:)=20*I_SHUFF50(N_50,:,:,:,:);
        I_SHUFF(N_100,:,:,:,:)=10*I_SHUFF100(N_100,:,:,:,:);
        use_neurons=use_bkg;
    case 2, N_all=unique([N_ret_100' N_ret_50]);
        N_50=ismember(N_all,N_ret_50);
        N_100=ismember(N_all,N_ret_100);
        Entropy(N_50,:,:,:)=20*Entropy50(N_50,:,:,:);
        Entropy(N_100,:,:,:)=10*Entropy100(N_100,:,:,:);
        use_neurons=use_ret;
end

%for dot/bkg neurons
Neu=use_neurons;
P_CHIV=P_CHIV_all(Neu,:,:);
P_CHIV_dotsh=P_CHIV_dotsh_all(Neu,:,:);

figs=[1 7 12];
i_fig=0;
clear Entropy_thre_tmp  Entropy_thre
for i_bin=1:size(Entropy,4)
if ismember(i_bin,figs)
fi=true;
i_fig=i_fig+1;
else
fi=false;
end
%draw time course of Inforamtion value (entropy) of all neurons which
%represented significant infomation
Entropy_thre_tmp(:,:,1:2)=(P_CHIV(:,:,1:2)<=p_thre).*(P_CHIV_dotsh(:,:,1:2)<=p_thre).*Entropy(:,:,1:2,i_bin);
Entropy_thre_tmp(:,:,3)=(P_CHIV(:,:,3)<=p_thre).*Entropy(:,:,3,i_bin);
for i=1:size(Entropy_thre_tmp,1)
    Entropy_thre(i,:,1)=Suc5(Entropy_thre_tmp(i,:,1));
    Entropy_thre(i,:,2)=Suc5(Entropy_thre_tmp(i,:,2));
    Entropy_thre(i,:,3)=Suc5(Entropy_thre_tmp(i,:,3));
end

amp=10;
if coord==1
%dot/bkg
ent_shuff_mean=squeeze(nanmean(I_SHUFF,2));
entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0);%.*(Entropy_thre(:,:,2)>0);
entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
if fi
subplot('position',[0.1+(i_fig-1)*0.25 0.7 0.2 0.1])
show_SumEntropy(entropy_thre_minus_spread/size(entropy_thre_minus_spread,1),1,time_end,fs,'m')
set(gca,'xtick',[1 6 11 16,21],'xticklabel',{'0','50','100','150','200'},'TickDir','out','box','off','LineWidth',1.5)
if i_fig==1
set(gca,'ytick',[0 .5 1 ],'yticklabel',{'0','0.5','1'},'TickDir','out','box','off')
else
set(gca,'ytick',[0 .5 1 ],'yticklabel',{},'TickDir','out','box','off')
end
plot([1,1],[0 28],'k:','Linewidth',1)
plot([6,6],[0 28],'k--','Linewidth',1)
ylim([0 1.4])
xlim([-1 24])
end
%0-180
sumI_dotbkg(i_bin,1)=sum(sum(entropy_thre_minus_spread(:,1:24),2));
else
%dot/ret
if fi
subplot('position',[0.1+(i_fig-1)*0.25 0.3 0.2 0.1])
show_SumEntropy(Entropy_thre(:,:,2)/size(Entropy_thre,1),1,time_end,fs,'b')
set(gca,'xtick',[1 6 11 16,21],'xticklabel',{'0','50','100','150','200'},'TickDir','out','box','off','LineWidth',1.5)
if i_fig==1
set(gca,'ytick',[0 3 6],'yticklabel',{'0','3','6'},'TickDir','out','box','off')
else
set(gca,'ytick',[0 3 6],'yticklabel',{},'TickDir','out','box','off')
end

plot([1,1],[0 28],'k:','Linewidth',1)
plot([6,6],[0 28],'k--','Linewidth',1)
ylim([0 8])
xlim([-1 24])
end
sumI_dotret(i_bin,1)=sum(sum(Entropy_thre(:,1:24,2),2));
end
end
end

%%
%draw the average of information which calulated by each dot (dot/bkg, Figfure 12B)
subplot('position',[0.25 0.54 0.4 0.1])
plot(sumI_dotbkg(:,1)/sum(use_bkg)/24,'mo','MarkerFaceColor','m')
hold on
xlim([0 13])
ylim([0 .5])
set(gca,'xtick',1:12,'xticklabel',{},'TickDir','out','box','off','LineWidth',1.5)
set(gca,'ytick',[0 0.2 0.4],'yticklabel',{'0','0.2','0.4'},'TickDir','out','box','off')

%%
%draw the average of information which calulated by each dot (dot/ret, Figfure 12D)
subplot('position',[0.25 0.14 0.4 0.1])
plot(sumI_dotret(:,1)/size(Entropy_thre,1)/24,'bo','MarkerFaceColor','b')
hold on
xlim([0 13])
ylim([0 4])
set(gca,'xtick',1:12,'xticklabel',{},'TickDir','out','box','off','LineWidth',1.5)
set(gca,'ytick',[0 2 4],'yticklabel',{'0','2','4'},'TickDir','out','box','off')
%%
if 0
cd figs/%material
print('fig12_eachdot.ai', '-dpdf', '-painters','-bestfit')
cd ../../Scripts_forFigureR1/
end
