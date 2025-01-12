%%
%this file make Figure 1G,H,I 
%%
clear all
close all

addpath(genpath('../toolboxes'))
%%
% size andresolution of the monitor 
monitor_width_cm = 166; % cm
monitor_height_cm = 166 * 9 / 16; % cm (16:9a aspect ratio)
monitor_resolution_x = 1920; % pixels
monitor_resolution_y = 1080; % pixels

% distance from eye to monitor
distance_to_eye_cm = 42; % cm

% view angle of a rectangle
rect_width_deg = 60; % horizontal
rect_height_deg = 40; % vertical

num_div_x = 4;
num_div_y = 3;

rect_div_width_deg = rect_width_deg / num_div_x;
rect_div_height_deg = rect_height_deg / num_div_y;

x_deg=zeros(1,20);
y_deg=zeros(1,20);

 rect_x_start_deg =  - rect_width_deg / 2;
 rect_y_start_deg =  - rect_height_deg / 2;

for i=1:(num_div_x+1)
    for j=1:(num_div_y+1)
        x_deg(j+ (i-1)*4)=rect_x_start_deg + (i-1) * rect_div_width_deg;
        y_deg(j+ (i-1)*4)=rect_y_start_deg + (j-1) * rect_div_height_deg;
    end
end

%%
%draw figur
figure('position', [488.0000  196.2  1021  465.8])
subplot(1, 3, 1)
hold on;
draw_rects(x_deg, y_deg, 4, 5, 'k-', 2, 1)
eye_x=0; eye_y=0; col='r-'; lw=2; 
for i=1:length(x_deg)
    [x_real(i), y_real(i)]=pseudo_deg2real_deg(x_deg(i), y_deg(i), eye_x, eye_y);
end

draw_rects(x_real, y_real, 4, 5, col, lw, 1)
axis equal;
set(gca, 'xtick', -40:10:40, 'ytick', -30:10:30)
set(gca,'box','off','tickdir','out','fontname','arial narrow','fontsize',12,'linewidth',1.5)
xlabel('X (deg)');
ylabel('Y (deg)');
xlim([-40 40]);
ylim([-30 30]);
text(-60, 40, 'G','FontName', 'Arial', 'FontSize', 17)
put_number(1:12);

subplot(1, 3, 3)
hold on
for s=5:(-1):1
    switch s
        case 1, eye_x=0; eye_y=0; col='r-'; lw=2; %eye_x, eye_y is gaze position in linear angle coordinate
        case 2, eye_x=5; eye_y=5; col='m-'; lw=1;
        case 3, eye_x=-5; eye_y=5; col='c-'; lw=1;
        case 4, eye_x=-5; eye_y=-5; col='g-'; lw=1;
        case 5, eye_x=5; eye_y=-5; col='b-'; lw=1;
    end
for i=1:length(x_deg)
    [x_real(i), y_real(i)]=pseudo_deg2real_deg(x_deg(i), y_deg(i), eye_x, eye_y);
end

draw_rects(x_real, y_real, 4, 5, col, lw, 0)
end

axis equal;
set(gca, 'xtick', -40:10:40, 'ytick', -30:10:30)
set(gca,'box','off','tickdir','out','fontname','arial narrow','fontsize',12,'linewidth',1.5)
xlabel('X (deg)');
ylabel('Y (deg)');

xlim([-40 40]);
ylim([-30 30]);

text(-60, 40, 'I','FontName', 'Arial', 'FontSize', 16)
put_number(1:12)

subplot(1, 3, 2)
hold on;
for s=5:(-1):1
    switch s
        case 1, eye_x=0; eye_y=0; col='r+'; lw=2;%eye_x, eye_y is gaze position in linear angle coordinate
        case 2, eye_x=5; eye_y=5; col='m+'; lw=2;
        case 3, eye_x=-5; eye_y=5; col='c+'; lw=2;
        case 4, eye_x=-5; eye_y=-5; col='g+'; lw=2;
        case 5, eye_x=5; eye_y=-5; col='b+'; lw=2;
    end
plot(eye_x, eye_y, col, 'LineWidth', lw, 'MarkerSize', 12)
end
axis equal;
set(gca, 'xtick', -40:10:40, 'ytick', -30:10:30)
set(gca,'box','off','tickdir','out','fontname','arial narrow','fontsize',12,'linewidth',1.5)
xlabel('X (deg)');
ylabel('Y (deg)');
xlim([-40 40]);
ylim([-30 30]);
text(-60, 40, 'H','FontName', 'Arial', 'FontSize', 16)

cd figs/
print('fig1_ghi.ai', '-dpdf', '-painters')
cd ../

%%

function put_number(n)
l=length(n);
for j=1:l
    i=n(j);
    x=-24+15*mod(i-1,4);
    y=13-13*floor((i-1)/4);
    text(x, y, num2str(i), 'FontName', 'Arial Narrow','FontSize', 12)
end
end

function draw_rects(x, y, n, m, col, lw, out_flag)

%draws edges of rectangles when n x m vertices are given 
%out_flag: 1  draws all edges
%out_flag: 0 no circumferences
if out_flag==1
% vertical edges
for j=1:m
for i=1:(n-1)
    x1=x(n*(j-1)+i); x2=x(n*(j-1)+i+1);
    y1=y(n*(j-1)+i); y2=y(n*(j-1)+i+1);
    plot([x1, x2],[y1, y2], col, 'LineWidth', lw)
    hold on
end
end

% horizontal edges
for j=1:n
for i=1:(m-1)
    x1=x(n*(j-1)+i); x2=x(n*j+i);
    y1=y(n*(j-1)+i); y2=y(n*j+i);
    plot([x1, x2],[y1, y2], col, 'LineWidth', lw)
    hold on
end
end

else

for j=2:(m-1)
for i=1:(n-1)
    x1=x(n*(j-1)+i); x2=x(n*(j-1)+i+1);
    y1=y(n*(j-1)+i); y2=y(n*(j-1)+i+1);
    plot([x1, x2],[y1, y2], col, 'LineWidth', lw)
    hold on
end
end

% horizontal edges
for j=1:n
for i=2:(m-2)
    x1=x(n*(j-1)+i); x2=x(n*j+i);
    y1=y(n*(j-1)+i); y2=y(n*j+i);
    plot([x1, x2],[y1, y2], col, 'LineWidth', lw)
    hold on
end
end
end
end


