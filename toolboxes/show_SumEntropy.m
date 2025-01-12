%%
%this fucntion draw the timecourse of summation of information value (Entropy) 
%input: ENT_R1:time course matrix of information value (Entropy)
%l:magnification rate for entropy (usually 10 to adjust the duration of time bin (0.1sed))
%bin_n: nubmer of bin to display
%fs:fontsize
%col:line color

%%
function show_SumEntropy(ENT_R1,l,bin_n,fs,col)
switch col
    case 'r', light_color=[1,0.8,0.8];
    case 'b', light_color=[0.8,0.8,1];
    case 'g', light_color=[0.8,1,0.8];
    case 'm', light_color=[1,0.8,1];
end
ENT_R1_thre_cumsum_retino=cumsum(ENT_R1);
hold on
for i=1:size(ENT_R1_thre_cumsum_retino,1)
    if sum(ENT_R1(i,1:bin_n))
        plot(ENT_R1_thre_cumsum_retino(i,1:bin_n)*l,[col,'-'],'Linewidth',0.5,'Color',light_color)
    end
end
plot(sum(ENT_R1(:,1:bin_n))*l,[col,'-'],'Linewidth',2)

set(gca,'xtick',[1 6 11 16 21 26 31],'xticklabel',{'','','','',''},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial')
xlim([-1 bin_n])
ymax=max(sum(ENT_R1(:,1:bin_n))*l);
ylim([0 ymax*1.1])