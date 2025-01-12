%%
%draw ellipse
%%
function draw_ellipse(long_axis,short_axis,angle,origin_x,origin_y,col)
div = 100;
theta = (0:div)/div*2*pi;
P = [long_axis*cos(theta); short_axis*sin(theta)];
%rotate the angle by linear transformation
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
R = [cos(angle) sin(angle); -sin(angle) cos(angle)];

Q = R*P;
%shift the cetner
Q(1,:) = Q(1,:) + origin_x;
Q(2,:) = Q(2,:) + origin_y;
plot(Q(1,:),Q(2,:),'-','Color',col);