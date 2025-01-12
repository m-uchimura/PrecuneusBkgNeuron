%%
%calculate sum of square errors (sse) between the timecourse of dot/bkg and
%product of dot/ret and bkg/ret with fitted paramters 
%%
function sse=sse_info(params)
global dot_r
global dot_b
global bkg_r

sse=sum((dot_b-params(1)*dot_r.*(bkg_r-params(2))).^2);