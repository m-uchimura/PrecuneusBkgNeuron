%%
%this file Leave the values other than 0 in 5 consecutive bins, 
% and set to 0 for bins in 4 or less.
%%
function [data_suc5,vlength_frame]=Suc5(data_suc)
data_suc5=data_suc;
a = find([0 data_suc5] ~= 0);
if isempty(a)
    vlength_frame=0;
else
    vs = a(find(diff([0 a]) > 1))-1;
    ve = [a(find(diff(a) > 1)) a(end)]-1;
    for n=1:length(vs)
        vres{n} = data_suc5(vs(n):ve(n));
        if ve(n)-vs(n)<4
            data_suc5(vs(n):ve(n))=0;
        end
    end
    vlength_frame = cellfun('length', vres);
end