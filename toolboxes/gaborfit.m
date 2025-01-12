%%
%Fit the receptive fields (data1) with gabor function 
%output parameters and statisitics for the fitting
%%
function [parafit,a,sigma,x0,y0,theta,lambda,phase,b,dc,r,prob,gannma]=gaborfit(data1)
global data_fit 
data_fit=data1;
s=size(data_fit);
dat=reshape(data_fit,[1 s(1)*s(2)]);
dat=(dat-nanmean(dat))/nanstd(dat);
amax=nanmax(nanmax(abs(dat)));
data_fit=reshape(dat,s);

%initial parameters and limitations
p0=[amax s(1)/2 s(1)/2 s(2)/2 0 s(1)/2 0 0 0.5];
pmin=[0 3 1 1 0 3 0 -1 0.01];
pmax=[2*amax s(1)*5 s(2) s(1) pi s(1)*5 pi 1 100];

%fitting
opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',@ssegab,...
    'Aineq',[],...
    'bineq',[],...
    'x0',p0,...
    'lb',pmin,...
    'ub',pmax,'options',opts); 
ms = MultiStart('Display','iter','UseParallel','never');
tic
[x,f] = run(ms,problem,100)
toc
%draw the fitted image
p=x;
[x, y]=meshgrid(1:length(data_fit(1,:)), 1:length(data_fit(:,1)));
subplot(2,2,3)
imagesc(data_fit)
axis equal
axis off
subplot(2,2,4)
imagesc(gaborf(x,y,p))
axis equal
axis off

dc=1-ssegab(p)/var(dat)/(s(1)*s(2));
[r, prob]=corrcoef(reshape(gaborf(x,y,p),[1 s(1)*s(2)]),dat);
parafit=p;
a=p(1);
sigma=p(2);
x0=p(3);
y0=p(4);
theta=p(5);
lambda=p(6);
phase=p(7);
b=p(8);
gannma=p(9);


