%%
%this file calculate half maximum time of information_sum vecter
%%
function [half_max]=calc_half_max(information_sum,half)
max_information=max(information_sum);
half_period_lobical=information_sum<max_information/half;
rising=find(diff(half_period_lobical)==-1);
falling=find(diff(half_period_lobical)==1);
rise_inf=information_sum([rising rising+1]);
fall_inf=information_sum([falling falling+1]);
half_max(1)=rising(1)+(max_information/half-rise_inf(1))/diff(rise_inf);
half_max(2)=falling(1)+(max_information/half-fall_inf(1))/diff(fall_inf);