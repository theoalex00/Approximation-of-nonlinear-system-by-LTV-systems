close all;
clear all;

%% Parameters
J1 = 45;
J2 = 60;
J3 = 75;
d1 = 5*1e-2;
d2 = 2*1e-2;
d3 = -5*1e-2;

kp1 = 10;
kp2 = 10;
kp3 = 10;
kI1 = 5;
kI2 = 5;
kI3 = 5;
k1 = 0.00;
k2 = 0.00;
k3 = 0.00;
omega1_ref = 0.01;
omega2_ref = -0.02;
omega3_ref = 0.05;

M = diag([J1,J2,J3,1,1,1]);

b = [d1 ; d2 ; d3 ; sqrt(kI1)*omega1_ref ; sqrt(kI2)*omega2_ref ; sqrt(kI3)*omega3_ref];
tmax = 100;
x0 = [0; 0; 0; 0; 0; 0];

%% ---------------------- Nonlinear system ----------------------

A_fn = @(x) [
      -kp1       ,   -J3*x(3)   ,   J2*x(2)   ,   sqrt(kI1) ,     0        ,    0       ;
      J3*x(3)    ,   -kp2       ,  -J1*x(1)   ,     0       ,   sqrt(kI2)  ,    0       ;           
     -J2*x(2)    ,  J1*x(1)     ,   -kp3      ,     0       ,     0        ,  sqrt(kI3) ;
     -sqrt(kI1)  ,     0        ,     0       ,     0       ,     0        ,    0       ;
         0       ,  -sqrt(kI2)  ,     0       ,     0       ,     0        ,    0       ;
         0       ,     0        ,  -sqrt(kI3) ,     0       ,     0        ,    0       ; 
];

odefun = @(t, x) M \ ( A_fn(x) * x + b);
tspan = [0; tmax;];
opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',0.1);
[t, x] = ode23tb(odefun, tspan, x0, opts);

u1 = (-kp1*x(:,1)+sqrt(kI1)*x(:,4));
u2 = (-kp2*x(:,2)+sqrt(kI2)*x(:,5));
u3 = (-kp2*x(:,3)+sqrt(kI2)*x(:,6));
tplot = tmax;

linewidth = 1.7;

hfig = figure(1);
hold on;
plot(t,x(:,1),'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$\omega_{1}$ $(rad/s)$','interpreter','latex');
xlabel('$Time$ $(s)$','interpreter','latex');
plot_template;

hfig = figure(2);
hold on;
plot(t,x(:,2),'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$\omega_{2}$ $(rad/s)$','interpreter','latex');
xlabel('$Time$ $(s)$','interpreter','latex');
plot_template;

hfig = figure(3);
hold on;
plot(t,x(:,3),'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$\omega_{3}$ $(rad/s)$','interpreter','latex');
yticks([0,0.05]);
xlabel('$Time$ $(s)$','interpreter','latex');
plot_template;

hfig = figure(4);
hold on;
plot(t,x(:,4),'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$z_{1}$','interpreter','latex');
plot_template;

hfig = figure(5);
hold on;
plot(t,x(:,5),'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$z_{2}$','interpreter','latex');
plot_template;

hfig = figure(6);
hold on;
plot(t,x(:,6),'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$z_{3}$','interpreter','latex');
plot_template;

%% ---------------------- Linear Time-varying systems ----------------------
t = zeros(tmax+1,1);
for i = 1:tmax+1
    t(i) = i-1;
end

% Initialization
x = zeros(tmax+1,6);
x(:,1) = ones(tmax+1,1)*x0(1);
x(:,2) = ones(tmax+1,1)*x0(2);
x(:,3) = ones(tmax+1,1)*x0(3);
x(:,4) = ones(tmax+1,1)*x0(4);
x(:,5) = ones(tmax+1,1)*x0(5);
x(:,6) = ones(tmax+1,1)*x0(6);

for j = 1:2   % number of iterations
    data    = [t, x];
    t_data  = data(:,1);        
    d_data  = data(:,2:end);    

    x1 = @(tt) interp1(t_data, d_data(:,1), tt, 'pchip', 'extrap');
    x2 = @(tt) interp1(t_data, d_data(:,2), tt, 'pchip', 'extrap');
    x3 = @(tt) interp1(t_data, d_data(:,3), tt, 'pchip', 'extrap');
    x4 = @(tt) interp1(t_data, d_data(:,4), tt, 'pchip', 'extrap');
    x5 = @(tt) interp1(t_data, d_data(:,5), tt, 'pchip', 'extrap');
    x6 = @(tt) interp1(t_data, d_data(:,6), tt, 'pchip', 'extrap');

A_fn = @(t,x) [
      -kp1       ,  -J3*x3(t) ,  J2*x2(t)   , sqrt(kI1) ,    0        ,   0       ;
      J3*x3(t)   ,  -kp2      , -J1*x1(t)   ,   0       ,   sqrt(kI2) ,   0       ;           
     -J2*x2(t)   , J1*x1(t)   ,  -kp3       ,   0       ,    0        , sqrt(kI3) ;
      -sqrt(kI1) ,    0       ,    0        ,   0       ,    0        ,    0      ;
         0       , -sqrt(kI2) ,    0        ,   0       ,    0        ,    0      ;
         0       ,    0       ,  -sqrt(kI3) ,   0       ,    0        ,    0      ;
];

    odefun = @(t, x) M \ ( A_fn(t, x) * x + b);
    tspan = [t_data(1), t_data(end)];
    opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',0.1);
    [t, x] = ode23tb(odefun, tspan, x0, opts);

    disp(j)

    %% Plots
    for i=1:6
        figure(i);
        hold on;
        if j==1
            plot(t,x(:,i),'Linewidth',linewidth,'LineStyle','--');
        elseif j==2
            plot(t,x(:,i),'Linewidth',linewidth,'LineStyle','-.');
        end
    end
end

for i=1:6
    figure(i);
    hold on;
    legend('nonlinear','1st','2nd','interpreter','latex','Box','off');
    hold off;
end
