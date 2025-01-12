function [kai_value,chi2ans]=calc_chi2_gain(Ao,Bo,A,B)
kai_map_raw0=zeros(4,3);%%
kai_sti0=zeros(4,3);%%
kai_map_raw=zeros(4,3);%%
kai_sti=zeros(4,3);%%
for x=1:4
    for y=1:3
        raw_dis0=Ao(15*x-14:15*x,13*y-12:13*y);
        kai_map_raw0(x,y)=nansum(nansum(raw_dis0));
        sti_dis0=Bo(15*x-14:15*x,13*y-12:13*y);
        kai_sti0(x,y)=nansum(nansum(sti_dis0));
    end
end
E_sti=kai_map_raw0./kai_sti0;
for x=1:4
    for y=1:3
        raw_dis=A(15*x-14:15*x,13*y-12:13*y);
        kai_map_raw(x,y)=nansum(nansum(raw_dis));
        sti_dis=B(15*x-14:15*x,13*y-12:13*y);
        kai_sti(x,y)=nansum(nansum(sti_dis));
    end
end

E_sti=kai_map_raw0./kai_sti0.*kai_sti;
E_sti=E_sti*(sum(kai_map_raw(:))/sum(E_sti(:)));
kai_value_matrix=((kai_map_raw-E_sti).^2)./E_sti;
kai_value=nansum(nansum(kai_value_matrix));
chi2ans=chi2cdf(kai_value,x*y-1,'upper');
end