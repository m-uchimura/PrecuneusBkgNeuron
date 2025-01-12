%%
%calculate of sum of square errors fot gabor fitting
%%
function sser = ssegab(param)
global data_fit
[x, y]=meshgrid(1:length(data_fit(1,:)), 1:length(data_fit(:,1)));
z=gaborf(x, y, param);
e=z-data_fit;
se=e.*e;
sser=sum(sum(se));
end
