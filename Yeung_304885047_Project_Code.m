%% Clear Cache 
clear all; 
close all;
clc;

%% Initial Conditions
run = input('Enter the case number (1-8) you want:');
if (run ~= 1 && run ~= 2 &&run ~= 3 && run ~= 4 && run ~= 5 && run ~= 6 && run ~= 7 && run ~= 8)
   fprintf('Error! ');
end

switch (run)
    
    case 1
        a0=2;
        m = 3;
        k = 200;
        c = 2;

    case 2
        a0=2;
        m =4;
        k = 50;
        c = 45;

    case 3
        a0=2;
        m = 5;
        k = 125;
        c= 50;
    case 4
        a0=2;
        m = 8;
        k = 25;
        c = 35;
    case 5
        a0=2;
        m = 10;
        k=100;
        c = 10;
    case 6 %over
        a0=10;
        m = 1;
        k=20;
        c = 20;
    case 7 %crit
        a0=10;
        m=1;
        k = 100;
        c = 20;
    case 8 % under
        a0=1;
        m = 20;
        k = 10;
        c = 20;
end

% type = input('Press (1) for forward euler or (2) for runge kutta: ' );
%% Pre-Defined Variables
   t0 = 0;
    t_f = 10; 
    dt = 0.005;
    nt = (t_f-t0)/dt;
    t=linspace(t0,t_f,nt); 
    w = sqrt(k/m);
    epilson = 0.5*c*sqrt(1/(m*k));
    x0=[1,0];
    x1 = zeros(nt,2);
    x2 = zeros(nt,2);
    x3 = zeros(nt,2);

%% Coefficient Change
omega = sqrt(k/m);
s = c/2*sqrt(1/(m*k));
dt = 0.005;
v=0;
%% Forward Euler numeric scheme
position = zeros(1,nt);
position(1) = 1;
velocity = zeros(1,nt);


    [x] = VibrationPosition(x0,a0,t,nt,epilson,m,k,w,dt,'forwardeuler',position,velocity);

%% 2nd Order RK-Scheme
for i = 1:1:2
    if i == 1
        a0 = 0; % w/ homogenous
        j = 1;
    else
        j = 2;
    end
figure(1)
    subplot(3,2,j);
    a = plot(t,x(1,:));
    set(a,'LineWidth',5)
    hold on
    xlabel('Time(s)')
    ylabel('Position(m)')
    xlim([0,10])
    ylim([-1 1])
    if i == 1
    title('Forward Homogenous','FontSize',24)
    else
        title('Forward InHomogenous','FontSize',24)
    end
    set(gcf,'Position',[50 50 1000 800])
    set(gca,'LineWidth',3,'FontSize',20)

    [x] = VibrationPosition(x0,a0,t,nt,epilson,m,k,w,dt,'2nd order RK',position,velocity);
    subplot(3,2,j+2);
    a = plot(t,x(1,:));
    set(a,'LineWidth',5)
    hold on
    
    xlabel('Time(s)')
    ylabel('Position(m)')
    xlim([0,10])
    ylim([-1 1])
    if i == 1
        title('2nd RK Homogenous','FontSize',24)
    else
        title('2nd RK InHomogenous','FontSize',24)
    end    
    set(gcf,'Position',[50 50 1000 800])
    set(gca,'LineWidth',3,'FontSize',20)
       
   %% 4th Order RK-Scheme     
    [x] = VibrationPosition(x0,a0,t,nt,epilson,m,k,w,dt,'4th order RK',position,velocity);
    subplot(3,2,j+4);    
    a = plot(t,x(1,:));
    set(a,'LineWidth',5)
    hold on
    hold on
    xlabel('Time(s)')
    ylabel('Position(m)')
    xlim([0,10])
    ylim([-1 1])
    if i == 1
        title('4thRK Homogenous','FontSize',24)
    else
        title('4th RK InHomogenous','FontSize',24)
    end
    set(gcf,'Position',[50 50 1000 800])
    set(gca,'LineWidth',3,'FontSize',20)

end    
    
%% Bode Plots for Frequency and Gain
%initializations
m_new = m;
k_new = k;
c_new = c;
w_n = sqrt(k_new/m_new);
w_out = [(1/10):dt:100];
epilson = 0.01*(c_new/2)*sqrt(1/(m_new*k_new));
g = zeros(1,length(w_out));
l = zeros(1,length(w_out));
io = zeros(1,length(w_out));

for t = 1:1:length(w_out)
   l(t) = w_out(t)/w_n;
   g(t) = 1/(sqrt((1-l(t)^2)^2)+(2*epilson*l(t))^2);
   io(t) = -atan((2*epilson*l(t))/(1-l(t)^2));
end

%plots for absolute gain and phase shift
io_new = rad2deg(io)-90;
figure (2)
subplot(1,2,1)
plot(l,g)
xlim([0.01 10]);
ylim([0.01 10])
xlabel('Normal frequency')
ylabel('Absolute gain')
title('Gain plot')

subplot(1,2,2)
plot(l,io_new)
xlim([0.01 10]);
xlabel('Normal frequency')
ylabel('Phase Shift')
title('Phase plot')

%% Animation of the Time Evolution

vidfile = VideoWriter('testmovie8.mp4','MPEG-4');
% vidfile = VideoWriter('ExpMovie.avi','Motion JPEG AVI');
vidfile.FrameRate = 30;
open(vidfile);

for i = 2:1:length(x(1,:))
    
    
figure (3)
subplot(1,2,1)
a0 = 0;
title('Homogenous 4th Order RK')
plot(0,x(1,i-1),'ws','MarkerSize',20,'MarkerFaceColor','w')
plot(0,x(1,i),'rs','MarkerSize',20,'MarkerFaceColor','b')
axis([-1,2,-1,1.5])
hold on
plot(0,position(1),'rs','MarkerSize',20)

subplot(1,2,2)
title('Inhomogenous 4th Order RK')
plot(0,x(1,i-1),'ws','MarkerSize',20,'MarkerFaceColor','w')
plot(0,x(1,i),'rs','MarkerSize',20,'MarkerFaceColor','b')
axis([-1,2,-1,1.5])
hold on
plot(0,position(1),'rs','MarkerSize',20)
Fr(k+1) = getframe(gcf);
writeVideo(vidfile,Fr(k+1));
end






%% Functions

% Name:function [x].This function would approximate the position of a mass-spring-damper
% system from the second partial differential equation for an individual
% time-step. 
function [x] = VibrationPosition(x0,a0,t,nt,epilson,m,k,w,dt,type,position,velocity)

for k = 1:1:(nt - 1)
    f1 = a0*sin(t(k)/(2*pi))/m; % forcing function
    f2 = a0*sin((t(k)+0.5*dt)/(2*pi))/m;
    f3 = a0*sin(t(k+1)/(2*pi))/m;
    switch type 
        case'forwardeuler'
            position(k+1) = velocity(k)*dt + position(k);
            velocity(k+1) = velocity(k)+dt*(-2*epilson*w*velocity(k)-w^2*position(k)+f1);           
        case '2nd order RK'
            cx1=dt*velocity(k);
            cv1 = dt*(-2*epilson*w*velocity(k)-w^2*position(k)+(f1));
            cx2= dt*(velocity(k)+0.5*cv1);
            cv2 = dt*(-2*epilson*w*(velocity(k)+0.5*cv1)-w^2*(position(k)+0.5*cx1)+f2);
            position(k+1) = position(k) + cx2;
            velocity(k+1) = velocity(k) + cv2;
        case '4th order RK'
            cx1=dt*velocity(k);
            cv1 = dt*(-2*epilson*w*velocity(k)-w^2*position(k)+f1);
            cx2= dt*(velocity(k)+0.5*cv1);
            cv2 = dt*(-2*epilson*w*(velocity(k)+0.5*cv1)-w^2*(position(k)+0.5*cx1)+f2);
            cx3 = dt*(velocity(k) + 0.5*cv2);
            cv3 = dt*(-2*epilson*w*(velocity(k)+0.5*cv2)-w^2*(position(k)+0.5*cx2)+f2);
            cx4 = dt*(velocity(k) + 0.5*cv3);
            cv4 = dt*(-2*epilson*w*(velocity(k)+cv3)-w^2*(position(k)+cx3)+f3);
            position(k+1) = position(k)+ 1/6*(cx1+2*cx2+2*cx3+cx4);
            velocity(k+1) = velocity(k) + 1/6*(cv1+2*cv2+2*cv3+cv4);
    end
end
x = [position;velocity];
end