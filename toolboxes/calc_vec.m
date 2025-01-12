function vec=calc_vec(sti,forward)
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
vec_temp=kai_map_raw./kai_sti;
vec=vec_temp(:);