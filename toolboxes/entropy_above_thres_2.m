%%
%this file retrun ectropy vector for which both of 2 P value (P1 and P2) has been below the threshold 
% for 5 consecutive bins and set to 0 for less than 4
%%
function ent_thre=entropy_above_thres_2(Entropy,P1,P2,thresh)
Entropy_thre=Entropy;
P_thre_tmp=(P1<=thresh)&(P2<=thresh);
for i=1:size(P1,1)
P_thre(i,:)=Suc5(P_thre_tmp(i,:));
end
ent_thre=Entropy_thre.*P_thre;