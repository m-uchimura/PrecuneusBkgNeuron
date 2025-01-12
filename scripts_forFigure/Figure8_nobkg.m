%%
%this file draw the resuls of the comparison between Inforamtion for the dot/ret with bkg and
%without bkg (Figure8)
%%
clear 
close all
addpath(genpath('../toolboxes'))
%%
res_x=2560;
res_y=1600;
asp_rat_x=16;
asp_rat_y=10;
time_end=36;
fs=8;
load ../Data_open/data.mat
load ../Data_open/data_nobkg.mat

%%
%load neual information, Pvalue etc WITH background
p_thre=.004;
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

%%
%remove trials in which number of spikes are less than 5
load ../Data_open/meanS.mat
Entropy(:,:,1:2)=Entropy(:,:,1:2).*double(meanS>4);

data_bkg1_01_cor=nan(size(data_bkg0_1));
for i=1:length(data_bkg0_1)
    for j=1:length(data_bkg1)
        if strcmp(data_bkg0_1(i).file,data_bkg1(j).file)
        data_bkg1_01_cor(i)=j;
        end
    end
end
flt_rem=isnan(data_bkg1_01_cor);
data_bkg0_1(flt_rem)=[];
data_bkg1_01_cor(flt_rem)=[];
%%
Ent_dot_retino_tmp=entropy_above_thres(Entropy(:,:,2),P_CHIV_dotsh(:,:,2),p_thre);
Ent_dot_retino_bkg1=Ent_dot_retino_tmp(data_bkg1_01_cor,:);
Ent_bkg_retino_tmp=entropy_above_thres(Entropy(:,:,3),P_CHIV(:,:,3),p_thre);
Ent_bkg_retino_bkg1=Ent_bkg_retino_tmp(data_bkg1_01_cor,:);

load ../Data_open/Entropy_data_nobkg

Entropy_thre_tmp_nobkg_50(:,:,1:2)=(P_CHIV_nobkg_50(:,:,1:2)<=p_thre).*(P_CHIV_dotsh_nobkg_50(:,:,1:2)<=p_thre).*Entropy_nobkg_50(:,:,1:2);
Entropy_thre_tmp_nobkg_50(:,:,3)=(P_CHIV_nobkg_50(:,:,3)<=p_thre).*Entropy_nobkg_50(:,:,3);
for i=1:size(Entropy_thre_tmp_nobkg_50,1)
Entropy_thre_nobkg_50(i,:,1)=Suc5(Entropy_thre_tmp_nobkg_50(i,:,1));
Entropy_thre_nobkg_50(i,:,2)=Suc5(Entropy_thre_tmp_nobkg_50(i,:,2));
Entropy_thre_nobkg_50(i,:,3)=Suc5(Entropy_thre_tmp_nobkg_50(i,:,3));
end
N_info_nobkg_50=squeeze(sum(Entropy_thre_nobkg_50,2)>0);

Entropy_thre_tmp_nobkg_100(:,:,1:2)=(P_CHIV_nobkg_100(:,:,1:2)<=p_thre).*(P_CHIV_dotsh_nobkg_100(:,:,1:2)<=p_thre).*Entropy_nobkg_100(:,:,1:2);
for i=1:size(Entropy_thre_tmp_nobkg_100,1)
Entropy_thre_nobkg_100(i,:,1)=Suc5(Entropy_thre_tmp_nobkg_100(i,:,1));
Entropy_thre_nobkg_100(i,:,2)=Suc5(Entropy_thre_tmp_nobkg_100(i,:,2));
end
N_info_nobkg_100=squeeze(sum(Entropy_thre_nobkg_100,2)>0);
N_dotbkg_nobkg_100=find((N_info_nobkg_100(:,1)-N_info_nobkg_50(:,1))==1);
N_dotret_nobkg_100=find((N_info_nobkg_100(:,2)-N_info_nobkg_50(:,2))==1);

Entropy_nobkg=20*Entropy_nobkg_50;
P_CHIV_nobkg=P_CHIV_nobkg_50;
P_CHIV_dotsh_nobkg=P_CHIV_dotsh_nobkg_50;
% 
Entropy_nobkg(N_dotbkg_nobkg_100,:,1)=Entropy_nobkg_100(N_dotbkg_nobkg_100,:,1);
P_CHIV_nobkg(N_dotbkg_nobkg_100,:,1)=P_CHIV_nobkg_100(N_dotbkg_nobkg_100,:,1);
P_CHIV_dotsh_nobkg(N_dotbkg_nobkg_100,:,1)=P_CHIV_dotsh_nobkg_100(N_dotbkg_nobkg_100,:,1);

Entropy_nobkg(N_dotret_nobkg_100,:,2)=10*Entropy_nobkg_100(N_dotret_nobkg_100,:,2);
P_CHIV_nobkg(N_dotret_nobkg_100,:,2)=P_CHIV_nobkg_100(N_dotret_nobkg_100,:,2);
P_CHIV_dotsh_nobkg(N_dotret_nobkg_100,:,2)=P_CHIV_dotsh_nobkg_100(N_dotret_nobkg_100,:,2);

%%
load ../Data_open/meanSNo.mat
Entropy_nobkg(:,:,1:2)=Entropy_nobkg(:,:,1:2).*double(meanSNo>4);
Entropy_nobkg(flt_rem,:,:)=[];
P_CHIV_dotsh_nobkg(flt_rem,:,:)=[];
%%
Ent_dot_retino_bkg0=entropy_above_thres(Entropy_nobkg(:,:,2),P_CHIV_dotsh_nobkg(:,:,2),p_thre);

%%
figure('position',[0 0 floor(100000*asp_rat_x/res_x) floor(100000*asp_rat_y/res_y)])
%draw the time course of  dot/ret Information with(real line) and without (dashued line) bkg
%Figure 8B
subplot('position',[0.3 0.3 0.4 0.1])
Retino_bkg1=sum(Ent_dot_retino_bkg1,2)>0;
Retino_bkg0=sum(Ent_dot_retino_bkg0,2)>0;
flt=true(1,length(Retino_bkg0));
show_SumEntropy_nobkg(Ent_dot_retino_bkg1(flt,:)/sum(Retino_bkg1|Retino_bkg0),Ent_dot_retino_bkg0(flt,:)/sum(Retino_bkg1|Retino_bkg0),1,time_end,fs)
set(gca,'ytick',[0 1 2 ],'yticklabel',{'0','1','2'},'TickDir','out','box','off')
ylim([0 2.3])

%draw the time course of difference of the dot/ret Information
%Figure 8C
subplot('position',[0.3 0.21 0.4 0.07])
show_SumEntropy_nobkg([],Ent_dot_retino_bkg0(flt,:)/sum(Retino_bkg1|Retino_bkg0)-Ent_dot_retino_bkg1(flt,:)/sum(Retino_bkg1|Retino_bkg0),1,time_end,fs)
set(gca,'xtick',[1 6 11 16 21 26],'xticklabel',{'0','50','100','150','200','250','300'},'TickDir','out')
set(gca,'ytick',[0 .5 1],'yticklabel',{'0','0.5','1'},'TickDir','out','box','off')
ymax=35.98;
ylh3=ylabel(['Mean Inforamtion (bits/s)'],'Fontsize',fs,'FontName','Arial');
ylim([0 1.2])

if 0
%draw the time course of bkg/ret Information if needed (not show in Figure 8)
subplot('position',[0.3 0.12 0.4 0.07])
Ent_dot_retino_bkg0(flt,:)/sum(Retino_bkg1|Retino_bkg0)-Ent_dot_retino_bkg1(flt,:)/sum(Retino_bkg1|Retino_bkg0)
show_SumEntropy(Ent_bkg_retino_bkg1(Retino_bkg1|Retino_bkg0,:)/sum(sum(Ent_bkg_retino_bkg1(Retino_bkg1|Retino_bkg0,:),2)>0),1,time_end,fs,'g')
set(gca,'ytick',[0 1 2],'yticklabel',{'0','1','2'},'TickDir','out','box','off')
set(gca,'xtick',[1 6 11 16 21 26],'xticklabel',{'0','50','100','150','200','250','300'},'TickDir','out')
xlh=xlabel('Time from dot onset (ms)');
ylim([0 2])
end
sum(Retino_bkg1|Retino_bkg0)
sum(Retino_bkg1)
sum(Retino_bkg0)

%%
if 0
%save time course of dot/ret inforamtion for Figure S1
save output/result_nobkg Ent_dot_retino_bkg0 Ent_dot_retino_bkg1 Ent_bkg_retino_bkg1 flt
cd figs
print('fig8_nobkg.ai', '-dpdf', '-painters','-bestfit')
cd ..
end


