
function conv_resp=draw_BkgRet_pos(dat,quad)

bkg_x=dat.frame(:,1);
bkg_y=dat.frame(:,2);
fix_x=dat.fixation(:,1);
fix_y=dat.fixation(:,2);
trial_number=length(bkg_x);
t=0;
for i=1:trial_number
   if sum(( dat.spike_ts>dat.photo_ts{i}(1)-1)&( dat.spike_ts<dat.photo_ts{i}(1)+2.5))
         if ~isnan(fix_x(i))
           t=t+1;
           trial_spike(t)=i;
           %ts_cell{i}=dat.spike_ts((dat.spike_ts>dat.photo_ts{i}(1)-1)&( dat.spike_ts<dat.photo_ts{i}(1)+2));
         end
   end
end
%Spikes= ts_cell(trial_spike);
%Photo=dat.photo_ts(trial_spike);
%%
switch quad
    case 1,quad_g=1;
    case 2,quad_g=0;
    case 3,quad_g=2;
    case 4,quad_g=3;
end
gain_x=bkg_x-fix_x>=0;
gain_y=bkg_y-fix_y>=0;
gain=gain_x+gain_y*2;
g=gain(trial_spike);
bkg_pos_x=bkg_x(trial_spike)-fix_x(trial_spike);
bkg_pos_y=bkg_y(trial_spike)-fix_y(trial_spike);
bkg_pos_x=bkg_pos_x(g==quad_g);
bkg_pos_y=bkg_pos_y(g==quad_g);
bkgpos_map=zeros(51,41);
 for i=1:length(bkg_pos_x)
 bkgpos_map((26+bkg_pos_x(i)),21+bkg_pos_y(i))=bkgpos_map((26+bkg_pos_x(i)),21+bkg_pos_y(i))+1;
 end
frame_con=zeros(30,20);
frame_con(1,:)=1;
frame_con(30,:)=1;
frame_con(:,1)=1;
frame_con(:,20)=1;
conv_resp=conv2(bkgpos_map,frame_con);
conv_resp=conv_resp(11:71,11:51);
imagesc(conv_resp',[0 28])
hold on
rectangle('Position',[30 20 2 2] ,'EdgeColor','r','Curvature',[1 1],'LineWidth',2)
set(gca,'xtick',[11 31 51],'xticklabel',{'','',''},'ytick',[11 21 31],'yticklabel',{'','',''},'TickDir','out','box','off','LineWidth',1.5)    
