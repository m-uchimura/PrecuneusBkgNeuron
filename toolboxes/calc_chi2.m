%calculate mutualI, chi2 with 4*3 compartment

function [mutualI,kai_value,chi2ans]=calc_chi2(forward,sti)
kai_map_raw=zeros(4,3);%%
kai_sti=zeros(4,3);%%
for x=1:4
    for y=1:3
        raw_dis=forward(15*x-14:15*x,13*y-12:13*y);
        kai_map_raw(x,y)=nansum(nansum(raw_dis));
        sti_dis=sti(15*x-14:15*x,13*y-12:13*y);
        kai_sti(x,y)=nansum(nansum(sti_dis));
    end
end

Ave_raw=nanmean(reshape(kai_map_raw,x*y,1));
Sum_raw=nansum(reshape(kai_map_raw,x*y,1));
Ave_sti=nanmean(reshape(kai_sti,x*y,1));
Sum_sti=nansum(reshape(kai_sti,x*y,1));
E_sti=kai_sti*Ave_raw/Ave_sti;
information_matrix_sti=-(kai_sti/Sum_sti).*log2(kai_sti/Sum_sti);
information_matrix_spike=(-(kai_map_raw/Sum_raw).*log2(kai_map_raw/Sum_raw));
information_matrix_sti(isnan(information_matrix_sti))=0;
information_matrix_spike(isnan(information_matrix_spike))=0;
information_matrix=information_matrix_sti-information_matrix_spike;
mutualI=nansum(reshape(information_matrix,x*y,1));

kai_value_matrix=((kai_map_raw-E_sti).^2)./E_sti;
kai_value=nansum(nansum(kai_value_matrix));

chi2ans=chi2cdf(kai_value,x*y-1,'upper');

