%%
%draw the relationship between recording position and RF position(6G,I) or RF size(H,J)  6 G-J
function draw_gabor_Anatomy(data_bkg1,gabor_bkg,gabor_retino)
load ../Data_open/anatomy
%%
dc_thre=.5;
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

data_retino=data_bkg1([gabor_retino(:).i]);
data_bkg=data_bkg1([gabor_bkg(:).i]);
%%
load ../Data_MRI/M2L dura_point_all
dura_point_allN{1}=dura_point_all;
load ../Data_MRI/M1R dura_point_all
dura_point_allN{2}=dura_point_all;
load ../Data_MRI/M1L dura_point_all
dura_point_allN{3}=dura_point_all;
load ../Data_MRI/M3L dura_point_all
dura_point_allN{4}=dura_point_all;

%%
for i=1:length(data_bkg)
for j=1:942
    if strcmp(data_bkg1(j).file,data_bkg(i).file)
        i_new_bkg(i)=j;
    end
end
end
for i=1:length(data_retino)
for j=1:942
    if strcmp(data_bkg1(j).file, data_retino(i).file)
        i_new_ret(i)=j;
    end
end
end
lateral=Anat.X;
lr_ret=lateral(i_new_ret);
lr_bkg=lateral(i_new_bkg);
%%
col={'r','b','k','m','c'};
for use_co=0:1
if use_co
gab=gabor_retino;
dat1=data_retino;
lr1=lr_ret;
else
gab=gabor_bkg;
dat1=data_bkg;
lr1=lr_bkg;
end
ti={'XCenter','YCenter','Sigma'};
    save_V4_precuneus=1;
    if save_V4_precuneus
        hold on
        for i=1:length(dat1)
            lr=lr1(i);
            switch dat1(i).monkey
                case 2, theta=pi*(55/180);
                    dura_point=squeeze(dura_point_allN{1}(lr,:,:));
                    base=squeeze(dura_point_allN{1}(2,-5+8+1,:));
                case 4,theta=pi*(62/180);
                    dura_point=squeeze(dura_point_allN{2}(lr,:,:));
                    base=squeeze(dura_point_allN{2}(2,-5+8+1,:));
                case 4.5,theta=pi*(62/180);
                    dura_point=squeeze(dura_point_allN{3}(lr+1,:,:));
                    base=squeeze(dura_point_allN{3}(2+1,-5+8+1,:));
                case 5,theta=pi*(59/180);
                    dura_point=squeeze(dura_point_allN{4}(lr+1,:,:));
                    base=squeeze(dura_point_allN{4}(2+1,-5+8+1,:));
            end
            theta(i)=pi*dat1(i).angle/180;
            cell_xy=(floor(5*rand(1,2))-2)+dura_point(dat1(i).ap+8+1,:)+dat1(i).dep*200*[(-1)*cos(theta(i)) sin(theta(i))];
            if use_co
            yz_retino(i,:)=(dura_point(dat1(i).ap+8+1,:)+dat1(i).dep*200*[(-1)*cos(theta(i)) sin(theta(i))]-base')/200;
            else
            yz_bkg(i,:)=(dura_point(dat1(i).ap+8+1,:)+dat1(i).dep*200*[(-1)*cos(theta(i)) sin(theta(i))]-base')/200;
            end
end
end
end
%%
xyz=reshape([data_retino(:).xyz],3,149);
xyz(2:3,:)=yz_retino';
xyz(1,1:104)=-abs(xyz(1,1:104));
XYZ_retino=xyz;
gannma_retino=[gabor_retino(:).gannma];
sigma_retino=[gabor_retino(:).sigma];
y0_retino=40-[gabor_retino(:).y0];
x0_retino=[gabor_retino(:).x0];
ec_retinoc=[gabor_retino(:).ecc];
dc_retino=[gabor_retino(:).dc];
xyz=reshape([data_bkg(:).xyz],3,58);
xyz(2:3,:)=yz_bkg';
xyz(1,1:49)=-abs(xyz(1,1:49));
XYZ_bkg=xyz;
gannma_bkg=[gabor_bkg(:).gannma];
sigma_bkg=[gabor_bkg(:).sigma];
ecc_bkg=[gabor_bkg(:).ecc];
y0_bkg=40-[gabor_bkg(:).y0];
x0_bkg=[gabor_bkg(:).x0];
dc_bkg=[gabor_bkg(:).dc];

%%
Z_retino=-XYZ_retino(3,:)+30;
Z_bkg=-XYZ_bkg(3,:)+30;
subplot('position',[0.3+0.03-0.015 0.2-.04-0.15 0.165 0.12])
sigma_retino_c_all=sigma_retino./sqrt(gannma_retino);
sigma_bkg_c_all=sigma_bkg./sqrt(gannma_bkg);
[R,P]=corrcoef(Z_retino,sigma_retino_c);
plot(Z_retino(dc_retino>dc_thre&N_flt_r),sigma_retino_c(dc_retino>dc_thre&N_flt_r),'bo')
hold on
retino_reg=[Z_retino(dc_retino>dc_thre&N_flt_r);ones(1,length(Z_retino(dc_retino>dc_thre&N_flt_r)))]';
[X_retino,BINT_retino,R_retino,RINT_retino,STATS_retino]=regress(sigma_retino_c(dc_retino>dc_thre&N_flt_r)',retino_reg);
y_retino=retino_reg*X_retino;
plot([Z_retino(dc_retino>dc_thre&N_flt_r) 0],[y_retino; X_retino(2)],'b-','Linewidth',3)
plot(Z_retino(dc_retino>dc_thre&~N_flt_r),sigma_retino_c(dc_retino>dc_thre&~N_flt_r),'b+')
find(dc_retino<dc_thre&~N_flt_r); %Empty

xlim([-.5+min(nanmin(Z_bkg),nanmin(Z_retino)) .5+max(nanmax(Z_bkg),nanmax(Z_retino))])
ylim([0 1+max(nanmax(sigma_bkg_c),nanmax(sigma_retino_c))])
ylim([0 40])
text(-5+30,38,[' y=',num2str(X_retino(1),2),'x+',num2str(X_retino(2),2)])
text(-5+30,32,[' R=',num2str(sqrt(STATS_retino(1)),2),' p=',num2str(STATS_retino(3),2)])
set(gca,'xtick',[-5 0]+30,'xticklabel',{'',''},'ytick',[20 40],'yticklabel',{'',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
set(gca)
%%
subplot('position',[0.76+0.03-0.015 0.2-.04-0.15 0.16 0.12])
[R,P]=corrcoef(Z_bkg,sigma_bkg_c);
plot(Z_bkg(dc_bkg>dc_thre&N_flt_b),sigma_bkg_c(dc_bkg>dc_thre&N_flt_b),'mo')
hold on
bkg_reg=[Z_bkg(dc_bkg>dc_thre&N_flt_b);ones(1,length(Z_bkg(dc_bkg>dc_thre&N_flt_b)))]';
[X_bkg,BINT_bkg,R_bkg,RINT_bkg,STATS_bkg]=regress(sigma_bkg_c(dc_bkg>dc_thre&N_flt_b)',bkg_reg);
y_bkg=bkg_reg*X_bkg;
plot([Z_bkg(dc_bkg>dc_thre&N_flt_b) 0],[y_bkg; X_bkg(2)],'m-','Linewidth',3)
plot(Z_bkg(dc_bkg>dc_thre&~N_flt_b),sigma_bkg_c(dc_bkg>dc_thre&~N_flt_b),'m+')
find(dc_bkg<dc_thre&~N_flt_b);
xlim([-.5+min(nanmin(Z_bkg),nanmin(Z_retino)) .5+max(nanmax(Z_bkg),nanmax(Z_retino))])
ylim([0 1+max(nanmax(sigma_bkg_c),nanmax(sigma_retino_c))])
ylim([0 40])
text(-5+30,38,[' y=',num2str(X_bkg(1),2),'x+',num2str(X_bkg(2),2)])
text(-5+30,32,[' R=',num2str(sqrt(STATS_bkg(1)),2),' p=',num2str(STATS_bkg(3),2)])
set(gca,'xtick',[-5 0]+30,'xticklabel',{'',''},'ytick',[20 40],'yticklabel',{'',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
%%
subplot('position',[0.1 0.2-.04-0.15 0.16 0.12])
[R,P]=corrcoef(Z_retino,y0_retino);
plot(Z_retino(dc_retino>dc_thre),y0_retino(dc_retino>dc_thre),'bo')
hold on
retino_reg=[Z_retino(dc_retino>dc_thre);ones(1,length(Z_retino(dc_retino>dc_thre)))]';
[X_retino,BINT_retino,R_retino,RINT_retino,STATS_retino]=regress(y0_retino(dc_retino>dc_thre)',retino_reg);
y_retino=retino_reg*X_retino;
plot([Z_retino(dc_retino>dc_thre) 0],[y_retino; X_retino(2)],'b-','Linewidth',3)
xlim([-.5+min(nanmin(Z_bkg),nanmin(Z_retino)) .5+max(nanmax(Z_bkg),nanmax(Z_retino))])

ylim([0 1+max(nanmax(sigma_bkg_c),nanmax(y0_retino))])
ylim([0 40])
text(-5+30,38,[' y=',num2str(X_retino(1),2),'x+',num2str(X_retino(2),2)])
text(-5+30,32,[' R=',num2str(sqrt(STATS_retino(1)),2),' p=',num2str(STATS_retino(3),2)])
set(gca,'xtick',[-5 0]+30,'xticklabel',{'',''},'ytick',[20 40],'yticklabel',{'',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)

%%
subplot('position',[0.56 0.2-.04-0.15 0.16 0.12])
[R,P]=corrcoef(Z_bkg,y0_bkg);
plot(Z_bkg(dc_bkg>dc_thre),y0_bkg(dc_bkg>dc_thre),'mo')
hold on
bkg_reg=[Z_bkg(dc_bkg>dc_thre);ones(1,length(Z_bkg(dc_bkg>dc_thre)))]';
[X_bkg,BINT_bkg,R_bkg,RINT_bkg,STATS_bkg]=regress(y0_bkg(dc_bkg>dc_thre)',bkg_reg);

y_bkg=bkg_reg*X_bkg;
plot([Z_bkg(dc_bkg>dc_thre) 0],[y_bkg; X_bkg(2)],'m-','Linewidth',3)
xlim([-.5+min(nanmin(Z_bkg),nanmin(Z_retino)) .5+max(nanmax(Z_bkg),nanmax(Z_retino))])
ylim([0 40])

text(-5+30,38,[' y=',num2str(X_bkg(1),2),'x+',num2str(X_bkg(2),2)])
text(-5+30,32,[' R=',num2str(sqrt(STATS_bkg(1)),2),' p=',num2str(STATS_bkg(3),2)])
set(gca,'xtick',[-5 0]+30,'xticklabel',{'',''},'ytick',[20 40],'yticklabel',{'',''},'TickDir','out','box','off','Linewidth',1.5,'FontName','Arial narrow','Fontsize',12)
