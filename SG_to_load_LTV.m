close all;
clear all;

%% Parameters
P = 2;
Lq = 4.8*1e-3;
Ld = 5*1e-3;
Lf = 576.92*1e-3;
kf = 40*1e-3*sqrt(3/2);
kD = 8.4*1e-4;
kQ = 8*1e-4;
LQ = 5.1618*1e-4;
LD = 35.3685*1e-4;
LfD = 20*1e-3*sqrt(3/2);
Rs = 0.0031;
Rf = 0.715;
RD = 1.1421;%0.9523;%;1.1421;
RQ = 0.9523;
b = 0.0001;
J = 27547;
Rload = 4;
omega_ref = 100*pi;
If_ref = 500;

kpm = 20;
kIm = 5;
kpf = 3;
kIf = 3;
k1 = 0.00;
k2 = 0.00;

tmax = 20000;

M  =  [Ld kf  kD  0  0     0     0   0;
       kf Lf  LfD 0  0     0     0   0;
       kD LfD LD  0  0     0     0   0;
       0  0   0   Lq kQ    0     0   0;
       0  0   0   kQ LQ    0     0   0;
       0  0   0   0  0  4*J/P^2  0   0;
       0  0   0   0  0     0     1   0;
       0  0   0   0  0     0     0   1
      ]; 

%% ---------------------- Nonlinear system ------------------

A_fn = @(x) [
    -(Rs+Rload) ,     0       ,  0   ,          0                 ,  0   ,  -(Lq*x(4)+kQ*x(5))       ,   0        ,   0;
        0       , -(Rf+kpf)   ,  0   ,          0                 ,  0   ,           0               , sqrt(kIf)  ,   0;
        0       ,     0       , -RD  ,          0                 ,  0   ,           0               ,   0        ,   0;
        0       ,     0       ,  0   ,      -(Rs+Rload)           ,  0   , Ld*x(1)+kf*x(2)+kD*x(3)   ,   0        ,   0;
        0       ,     0       ,  0   ,          0                 , -RQ  ,           0               ,   0        ,   0;
 Lq*x(4)+kQ*x(5),     0       ,  0   , -(Ld*x(1)+kf*x(2)+kD*x(3)) ,  0   ,      -(4*b/P^2+kpm)       ,   0        , sqrt(kIm);
        0       , -sqrt(kIf)  ,  0   ,          0                 ,  0   ,           0               ,   0        ,   0;
        0       ,     0       ,  0   ,          0                 ,  0   ,        -sqrt(kIm)         ,   0        ,   0;
];
c = [0; 0; 0; 0; 0; 0; sqrt(kIf)*If_ref; sqrt(kIm)*omega_ref;];

odefun = @(t, x) M \ ( A_fn(x) * x + c);
tspan = [0; tmax;];
x0    = [0; 0; 0; 0; 0; 0; 0; 0;];
opts = odeset('RelTol',1e-4,'AbsTol',1e-4 );
[t, x] = ode23tb(odefun,tspan,x0,opts);

linewidth = 1.7;

hfig = figure(1);
hold on;
plot(t,x(:,1),'Linewidth',linewidth);
xlim([0,1500]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$I_{d}$ $(A)$','interpreter','latex');
plot_template;
hold off;

hfig = figure(2);
hold on;
plot(t,x(:,2),'Linewidth',linewidth);
xlim([0,15]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$I_{f}$ $(A)$','interpreter','latex');
yticks([0 200 400 500 600]);
plot_template;
hold off;

hfig = figure(3);
hold on;
plot(t,x(:,3),'Linewidth',linewidth);
xlim([0,15]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$I_{D}$ $(A)$','interpreter','latex');
plot_template;
hold off;

hfig = figure(4);
hold on;
plot(t,x(:,4),'Linewidth',linewidth);
xlim([0,1500]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$I_{q}$ $(A)$','interpreter','latex');
plot_template;
hold off;

hfig = figure(5);
hold on;
plot(t,x(:,5),'Linewidth',linewidth);
xlim([0,1500]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$I_{Q}$ $(A)$','interpreter','latex');
plot_template;
hold off;

hfig = figure(6);
hold on;
plot(t,x(:,6),'Linewidth',linewidth);
xlim([0,2200]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$\omega_{e}$ $(rad/s)$','interpreter','latex');
yticks([0 200 314 400 600]);
plot_template;
hold off;

hfig = figure(7);
hold on;
plot(t,x(:,7),'Linewidth',linewidth);
xlim([0,15]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$z_{1}$','interpreter','latex');
plot_template;
hold off;

hfig = figure(8);
hold on;
plot(t,x(:,8),'Linewidth',linewidth);
xlim([0,1500]);
xlabel('$Time$ $(s)$','interpreter','latex');
ylabel('$z_{2}$','interpreter','latex');
plot_template;
hold off;

%% ---------------------- Linear Time-varying systems ------------------

t = zeros(tmax+1,1);
for i = 1:tmax+1
    t(i) = i-1;
end

% Initialization
x = zeros(tmax+1,8);
x(:,1) = ones(tmax+1,1)*x0(1);
x(:,2) = ones(tmax+1,1)*x0(2);
x(:,3) = ones(tmax+1,1)*x0(3);
x(:,4) = ones(tmax+1,1)*x0(4);
x(:,5) = ones(tmax+1,1)*x0(5);
x(:,6) = ones(tmax+1,1)*x0(6);
x(:,7) = ones(tmax+1,1)*x0(7);
x(:,8) = ones(tmax+1,1)*x0(8);

for j = 1:4  % number of iterations
    data    = [t, x];
    t_data  = data(:,1);        
    d_data  = data(:,2:end);    

    d1_fun = @(tt) interp1(t_data, d_data(:,1), tt, 'pchip', 'extrap');
    d2_fun = @(tt) interp1(t_data, d_data(:,2), tt, 'pchip', 'extrap');
    d3_fun = @(tt) interp1(t_data, d_data(:,3), tt, 'pchip', 'extrap');
    d4_fun = @(tt) interp1(t_data, d_data(:,4), tt, 'pchip', 'extrap');
    d5_fun = @(tt) interp1(t_data, d_data(:,5), tt, 'pchip', 'extrap');
    d6_fun = @(tt) interp1(t_data, d_data(:,6), tt, 'pchip', 'extrap');
    d7_fun = @(tt) interp1(t_data, d_data(:,7), tt, 'pchip', 'extrap');
    d8_fun = @(tt) interp1(t_data, d_data(:,8), tt, 'pchip', 'extrap');

    A_fn = @(t, x) [
    -(Rs+Rload)           ,     0       ,  0   ,          0                                ,  0   ,  -(Lq*d4_fun(t)+kQ*d5_fun(t))          ,  0        ,    0;
        0                 , -(Rf+kpf)   ,  0   ,          0                                ,  0   ,           0                            , sqrt(kIf) ,    0;
        0                 ,     0       , -RD  ,          0                                ,  0   ,           0                            ,  0        ,    0;
        0                 ,     0       ,  0   ,         -(Rs+Rload)                       ,  0   , Ld*d1_fun(t)+kf*d2_fun(t)+kD*d3_fun(t) ,  0        ,    0;
        0                 ,     0       ,  0   ,          0                                , -RQ  ,           0                            ,  0        ,    0;
Lq*d4_fun(t)+kQ*d5_fun(t) ,     0       ,  0   , -(Ld*d1_fun(t)+kf*d2_fun(t)+kD*d3_fun(t)) ,  0   ,     -(4*b/P^2+kpm)                     ,  0        , sqrt(kIm);
        0                 ,  -sqrt(kIf) ,  0   ,          0                                ,  0   ,           0                            ,  0        ,    0;
        0                 ,     0       ,  0   ,          0                                ,  0   ,        -sqrt(kIm)                      ,  0        ,    0;];

    odefun = @(t, x) M \ ( A_fn(t, x) * x + c);
    opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
    [t, x] = ode23tb(odefun, tspan, x0, opts);
    
    disp(j);
    if (j==1||j==2||j==4)
    for i=1:8
        figure(i);
        hold on;
        if i==6
            if j==1
            plot(t,x(:,i),'Linewidth',linewidth,'LineStyle','--');
            elseif j==2
            plot(t,x(:,i),'Linewidth',linewidth,'LineStyle','-.');
            elseif j==4
            plot(t,x(:,i),'Linewidth',linewidth,'LineStyle',':');
            end
        else
        if j==1
        plot(t,x(:,i),'Linewidth',linewidth,'LineStyle','--');
        elseif j==2
        plot(t,x(:,i),'Linewidth',linewidth,'LineStyle','-.');
        elseif j==4
        plot(t,x(:,i),'Linewidth',linewidth,'LineStyle',':');   
        hold off;
        end
        end
    end
    end
end 

for i=1:8
    figure(i);
    hold on;
    legend('nonlinear','1st','2nd','4th','interpreter','latex','Box','off');
    hold off;
end



       



