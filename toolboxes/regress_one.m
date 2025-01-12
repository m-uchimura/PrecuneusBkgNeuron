function [b,stats,x_fit,y_fit]=regress_one(x,y)
%x,y both tate vector
if size(x,1)==1,x=x';end
if size(y,1)==1,y=y';end

X=[ones(length(x),1) x];
if (max(x)-min(x))<.02
x_fit=min(x):.001:max(x);
else
x_fit=min(x):.01:max(x);
end
[b ,bint,r,rint,stats]= regress(y,X);
y_fit=b(1)+b(2)*x_fit;