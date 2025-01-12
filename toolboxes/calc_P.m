function [P_I]= calc_P(maxI,entropy,rep_num)
Ord_I=max(find(sort(maxI)<entropy));
if isempty( Ord_I)
    Ord_I=0;
end
P_I=(rep_num-Ord_I)/rep_num;