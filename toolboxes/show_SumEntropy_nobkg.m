%%
%this fucntion draw the timecourse of summation of information value(Entropy) 
% for nobkg experimetn (Figure 7)
%input: ENT_bkg:time course matrix of information value (Entropy) with bkg
%input: ENT_nobkg:time course matrix of information value (Entropy) without bkg
%l:magnification rate for entropy (usually 10 to adjust the duration of time bin (0.1sed))
%bin_n: nubmer of bin to display
%fs:fontsize

%%
function show_SumEntropy_nobkg(ENT_withbkg,ENT_nobkg,l,bin_n,fs)
ENT_R1_thre_cumsum_nobkg=cumsum(ENT_nobkg);
%draw the timecourse of entropy without bgk
plot(sum(ENT_nobkg(:,1:bin_n))*l,'b--','Linewidth',2)
hold on
for i=1:size(ENT_R1_thre_cumsum_nobkg,1)
    if sum(ENT_nobkg(i,1:bin_n))
        plot(ENT_R1_thre_cumsum_nobkg(i,1:bin_n)*l,'b--','Linewidth',0.3,'Color',[0.8,0.8,0.9])
    end
end

%draw the timecourse of entropy with bgk
if ~isempty(ENT_withbkg)
ENT_R1_thre_cumsum_bkg=cumsum(ENT_withbkg);
plot(sum(ENT_withbkg(:,1:bin_n))*l,'b-','Linewidth',2)
hold on
for i=1:size(ENT_R1_thre_cumsum_bkg,1)
    if sum(ENT_withbkg(i,1:bin_n))
        plot(ENT_R1_thre_cumsum_bkg(i,1:bin_n,1)*l,'b-','Linewidth',0.3,'Color',[0.8,0.8,0.9])
    end
end
plot(sum(ENT_withbkg(:,1:bin_n))*l,'b-','Linewidth',2)
end
plot(sum(ENT_nobkg(:,1:bin_n)*l),'b--','Linewidth',2)

set(gca,'ytick',[0  20 40 60 100 150],'yticklabel',{'0','20','40','60','15'},'TickDir','out','box','off')
set(gca,'xtick',[1 11 21],'xticklabel',{'','','','',''},'TickDir','out','box','off')
set(gca,'Linewidth',1.5,'Fontsize',fs,'FontName','Arial')
if isempty(ENT_withbkg)
    ymax=max(sum(ENT_nobkg(:,1:bin_n))*l);
else
    ymax=max([sum(ENT_withbkg(:,1:bin_n))*l sum(ENT_nobkg(:,1:bin_n))*l]);
end
xlim([-1 bin_n])
ylim([0 ymax*1.1])