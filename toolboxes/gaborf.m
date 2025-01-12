%%
%gabor funtion with parameters as param
%%
function z = gaborf(x, y, param)
%parameters for Gaussian
a=param(1);
sigma=param(2);
x0=param(3);
y0=param(4);

%paramiters for sinusoidal wave
theta=param(5);
lambda=param(6);
phase=param(7);
b=param(8);
g=param(9);

xp = (x-x0)*cos(theta) + (y-y0)*sin(theta);
yp =-(x-x0)*sin(theta) + (y-y0)*cos(theta);

%The Gabor is defined as follows:
z = a*exp(-(xp.^2+g.^2*yp.^2)/2/sigma^2).*cos(2*pi*xp/lambda + phase) + b;

