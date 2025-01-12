%%
%this file retrun ectropy vector for which the P value has been below the threshold 
% for 5 consecutive bins and set to 0 for less than 4
%%
function ent_thre=entropy_above_thres(Entropy,P,thresh)
Entropy_thre=Entropy;
P_thre_tmp=P<=thresh;
for i=1:size(P,1)
P_thre(i,:)=Suc5(P_thre_tmp(i,:));
end
ent_thre=Entropy_thre.*P_thre;