%calculate mutualI, chi2 with 2*2 compartment
function [mutualI,kai_value,chi2ans]=calc_chi2_2_rotate(forward,sti,center_x,center_y)
x=2;y=2;
kai_map_raw=zeros(x,y);%%
kai_sti=zeros(x,y);%%
%[size_x,size_y]=size(forward);

kai_map_raw(1,1)=nansum(forward(find(forward(:,2)<=center_x&forward(:,3)<=center_y)));
kai_map_raw(2,1)=nansum(forward(find(forward(:,2)>center_x&forward(:,3)<=center_y)));
kai_map_raw(1,2)=nansum(forward(find(forward(:,2)<=center_x&forward(:,3)>center_y)));
kai_map_raw(2,2)=nansum(forward(find(forward(:,2)>center_x&forward(:,3)>center_y)));
kai_sti(1,1)=nansum(sti(find(sti(:,2)<=center_x&sti(:,3)<=center_y)));
kai_sti(2,1)=nansum(sti(find(sti(:,2)>center_x&sti(:,3)<=center_y)));
kai_sti(1,2)=nansum(sti(find(sti(:,2)<=center_x&sti(:,3)>center_y)));
kai_sti(2,2)=nansum(sti(find(sti(:,2)>center_x&sti(:,3)>center_y)));

% kai_map_raw
% kai_sti
Ave_raw=nanmean(reshape(kai_map_raw,x*y,1));
Sum_raw=nansum(reshape(kai_map_raw,x*y,1));
Ave_sti=nanmean(reshape(kai_sti,x*y,1));
Sum_sti=nansum(reshape(kai_sti,x*y,1));
E_sti=kai_sti*Ave_raw/Ave_sti;
information_matrix=-(kai_sti/Sum_sti).*log2(kai_sti/Sum_sti)-(-(kai_map_raw/Sum_raw).*log2(kai_map_raw/Sum_raw));

mutualI=nansum(reshape(information_matrix,x*y,1));

kai_value_matrix=((kai_map_raw-E_sti).^2)./E_sti;
kai_value=nansum(nansum(kai_value_matrix));

chi2ans=chi2cdf(kai_value,x*y-1,'upper');

