%%
%this fucntion draw the receptive fields 
%input: RF:2d matrix which represent spike number per stimulus
%min,max:limitation value of receptive image 
%P: P-value which is used to display
%ent: entropy which is used to display
%coordinate: coordiantion of RF
%fs:fontsize
%flip:which determine whter flip or not (to adjust the recorded hemisphre)
%%
function show_RF(RF,min,max,P,ent,coordinate,fs,flip)
if flip
    imagesc(fliplr(RF),[min,max])
else
    imagesc(RF,[min,max])
end
hold on
switch coordinate
    case 'retino', rectangle('Position',[30 20 2 2] ,'EdgeColor','w','Curvature',[1 1],'LineWidth',2)
    case 'bkg',rectangle('Position',[16 11 30 20],'EdgeColor','w','LineWidth',2)
    case 'bkg_ret',rectangle('Position',[30 20 2 2] ,'EdgeColor','w','Curvature',[1 1],'LineWidth',2)
end
t2_1= text(5,-3,['P=',num2str(P,2),'  I=',num2str(ent,2),'bits/s'],'Fontsize',fs);
if strcmp(coordinate,'retino')
set(gca,'ytick',[1 11 21 31 41],'yticklabel',{'20°','','0°','','-20°'},'TickDir','out','box','off','LineWidth',1.5)
else
set(gca,'ytick',[11 21 31],'yticklabel',{'','',''},'TickDir','out','box','off','LineWidth',1.5)    
end
set(gca,'xtick',[11 31 51],'xticklabel',{'-20°','0°','20°'},'TickDir','out','box','off','Fontsize',fs)