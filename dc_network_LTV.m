close all;
clear all;

%% Parameters
G1 = 1/7;
G2 = 1/6;
G3 = 1/8;
G12 = 1/0.4;
G23 = 1/0.5;
G13 = 1/0.8;
G22 = G1+G12+G13;
G44 = G2+G12+G23;
G77 = G3+G13+G23;

Rc1 = 0.4;
Lc1 = 0.004;
C1 = 1*1e-3;
Vin1 = 80;
kp1 = 1;
kI1 = 5;
c1 = 100;
Vref1 = 100;

Rc2 = 0.3;
Lc2 = 0.001;
C2 = 2*1e-3;
Vin2 = 80;
kp2 = 1;
kI2 = 5;
c2 = 100;
Vref2 = 100;

R = 0.5;
L = 0.01;
C3 = 3*1e-3;
omega = 100*pi;
Vd = 10;
Vq = 10;
kpd = 3;
kpq = 3;
kId = 5;
kIq = 5;
c3 = 100;
Id_ref = 6;
Iq_ref = 0;

M = diag([Lc1,C1,Lc2,C2,L,L,C3,1,1,1,1]);
b = [Vin1-Rc1*c1 ; G13*c3 ; Vin2-Rc2*c2 ; G23*c3 ; Vd ; Vq ; -c3*(G13+G23+G3) ; sqrt(kI1)*Vref1 ; sqrt(kI2)*Vref2 ; sqrt(kId)*Id_ref ; sqrt(kIq)*Iq_ref];
tmax = 10;
rel_tol = 1e-4;
abs_tol = 1e-4;
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

%% ---------------------- Nonlinear system ----------------------

no_states = 11;

u1 = @(x)(-kp1*x(2)+sqrt(kI1)*x(8))/c1;
u2 = @(x)(-kp2*x(4)+sqrt(kI2)*x(9))/c2;
md = @(x)-(-kpd*x(5)+sqrt(kId)*x(10))/c3;
mq = @(x)-(-kpq*x(6)+sqrt(kIq)*x(11))/c3;

A_fn = @(x) [
      -Rc1     ,     -u1(x)        ,    0         ,    0       ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
      u1(x)    ,   -(G22+kp1)      ,    0         ,   G12      ,      0            ,   0        ,    G13       ,  sqrt(kI1)   ,       0       ,    0        ,    0       ;
         0     ,        0          ,  -Rc2        ,  -u2(x)    ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
         0     ,       G12         ,   u2(x)      , -(G44+kp2) ,      0            ,   0        ,    G23       ,     0        ,    sqrt(kI2)  ,    0        ,    0       ;
         0     ,        0          ,    0         ,    0       ,    -(R+kpd)       , -L*omega   ,   -md(x)     ,     0        ,       0       ,  sqrt(kId)  ,    0       ;
         0     ,        0          ,    0         ,    0       ,    L*omega        , -(R+kpq)   ,   -mq(x)     ,     0        ,       0       ,    0        ,  sqrt(kIq) ;
         0     ,       G13         ,    0         ,   G23      ,      md(x)        ,   mq(x)    ,   -G77       ,     0        ,       0       ,    0        ,    0       ;
         0     ,     -sqrt(kI1)    ,    0         ,    0       ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
         0     ,        0          ,    0         , -sqrt(kI2) ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;        
         0     ,        0          ,    0         ,    0       ,   -sqrt(kId)      ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
         0     ,        0          ,    0         ,    0       ,      0            , -sqrt(kIq) ,     0        ,     0        ,       0       ,    0        ,    0       ;
];

odefun = @(t,x) M \ ( A_fn(x) * x + b);
tspan = [0; tmax;];
opts = odeset('RelTol',rel_tol,'AbsTol',abs_tol);
[t, x] = ode23tb(odefun, tspan, x0, opts);

%% Plots
tplot = tmax;
linewidth = 1.7;
for i=1:11
    hfig = figure(i);
    hold on;
    if i==1
        plot(t,x(:,i)+c1,'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$I_{c1}$ $(A)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==2
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$V_{1}$ $(V)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==3
        plot(t,x(:,i)+c2,'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$I_{c2}$ $(A)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==4
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$V_{2}$ $(V)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==5
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$I_{d}$ $(A)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==6
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$I_{q}$ $(A)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==7
        plot(t,x(:,i)+c3,'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$V_{3}$ $(V)$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==8
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$z_{1}$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==9
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$z_{2}$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==10
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$z_{d}$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    elseif i==11
        plot(t,x(:,i),'LineWidth',linewidth);
        xlim([0,tplot]);
        ylabel('$z_{q}$','interpreter','latex');
        xlabel('$Time$ $(s)$','Interpreter','latex');
        plot_template;
    end
end

u1 = (-kp1*x(:,2)+sqrt(kI1)*x(:,8))/c1;
u2 = (-kp2*x(:,4)+sqrt(kI2)*x(:,9))/c2;
md = -(-kpd*x(:,5)+sqrt(kId)*x(:,10))/c3;
mq = -(-kpq*x(:,6)+sqrt(kIq)*x(:,11))/c3;

hfig = figure(12);
plot(t,u1,'LineWidth',linewidth);
hold on;
plot(t,u2,'LineWidth',linewidth);
plot(t,md,'LineWidth',linewidth);
plot(t,mq,'LineWidth',linewidth);
xlim([0,tplot]);
ylabel('$Duty$ $ratios$','interpreter','latex');
xlabel('$Time$ $(s)$','Interpreter','latex');
legend('$u_{1}$','$u_{2}$','$m_{d}$','$m_{q}$','Interpreter','latex')
plot_template;
hold off;

%% ---------------------- Linear Time-varying systems ----------------------

no_iter = 15;
t = zeros(tmax+1,1);
for i = 1:tmax+1
    t(i) = i-1;
end

% Initialization
x = zeros(tmax+1,11);
x(:,1) = ones(tmax+1,1)*x0(1);
x(:,2) = ones(tmax+1,1)*x0(2);
x(:,3) = ones(tmax+1,1)*x0(3);
x(:,4) = ones(tmax+1,1)*x0(4);
x(:,5) = ones(tmax+1,1)*x0(5);
x(:,6) = ones(tmax+1,1)*x0(6);
x(:,7) = ones(tmax+1,1)*x0(7);
x(:,8) = ones(tmax+1,1)*x0(8);
x(:,9) = ones(tmax+1,1)*x0(9);
x(:,10) = ones(tmax+1,1)*x0(10);
x(:,11) = ones(tmax+1,1)*x0(11);

for j = 1:no_iter
    data = [t, x];
    t_data = data(:,1);        
    d_data = data(:,2:end);    

    x1 = @(t) interp1(t_data, d_data(:,1), t, 'spline', 'exap');
    x2 = @(t) interp1(t_data, d_data(:,2), t, 'spline', 'extrap');
    x3 = @(t) interp1(t_data, d_data(:,3), t, 'spline', 'extrap');
    x4 = @(t) interp1(t_data, d_data(:,4), t, 'spline', 'extrap');
    x5 = @(t) interp1(t_data, d_data(:,5), t, 'spline', 'extrap');
    x6 = @(t) interp1(t_data, d_data(:,6), t, 'spline', 'extrap');
    x7 = @(t) interp1(t_data, d_data(:,7), t, 'spline', 'extrap');
    x8 = @(t) interp1(t_data, d_data(:,8), t, 'spline', 'extrap');
    x9 = @(t) interp1(t_data, d_data(:,9), t, 'spline', 'extrap');
    x10 = @(t) interp1(t_data, d_data(:,10), t, 'spline', 'extrap');
    x11 = @(t) interp1(t_data, d_data(:,11), t, 'spline', 'extrap');

    u1 = @(t)(-kp1*x2(t)+sqrt(kI1)*x8(t))/c1;
    u2 = @(t)(-kp2*x4(t)+sqrt(kI2)*x9(t))/c2;
    md = @(t)-(-kpd*x5(t)+sqrt(kId)*x10(t))/c3;
    mq = @(t)-(-kpd*x6(t)+sqrt(kId)*x11(t))/c3;

A_fn = @(t) [
      -Rc1     ,     -u1(t)        ,    0         ,    0       ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
      u1(t)    ,   -(G22+kp1)      ,    0         ,   G12      ,      0            ,   0        ,    G13       ,  sqrt(kI1)   ,       0       ,    0        ,    0       ;
         0     ,        0          ,  -Rc2        ,  -u2(t)    ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
         0     ,       G12         ,   u2(t)      , -(G44+kp2) ,      0            ,   0        ,    G23       ,     0        ,    sqrt(kI2)  ,    0        ,    0       ;
         0     ,        0          ,    0         ,    0       ,    -(R+kpd)       , -L*omega   ,   -md(t)     ,     0        ,       0       ,  sqrt(kId)  ,    0       ;
         0     ,        0          ,    0         ,    0       ,    L*omega        , -(R+kpq)   ,   -mq(t)     ,     0        ,       0       ,    0        ,  sqrt(kIq) ;
         0     ,       G13         ,    0         ,   G23      ,      md(t)        ,   mq(t)    ,   -G77       ,     0        ,       0       ,    0        ,    0       ;
         0     ,     -sqrt(kI1)    ,    0         ,    0       ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
         0     ,        0          ,    0         , -sqrt(kI2) ,      0            ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;        
         0     ,        0          ,    0         ,    0       ,   -sqrt(kId)      ,   0        ,     0        ,     0        ,       0       ,    0        ,    0       ;
         0     ,        0          ,    0         ,    0       ,      0            , -sqrt(kIq) ,     0        ,     0        ,       0       ,    0        ,    0       ;
];

    %% RHS for ODE solver: ẋ = M⁻¹·(A(t,x)·x)
    odefun = @(t, x) M\(A_fn(t)*x + b);

    %% Integrate
    tspan = [t_data(1), t_data(end)];
    opts = odeset('RelTol',rel_tol,'AbsTol',abs_tol);
    [t, x] = ode23tb(odefun, tspan, x0, opts);

    if (j==1||j==5||j==13)
    for i=1:11
        figure(i);
        if i==1
            if j==1
                plot(t,x(:,i)+c1,'LineWidth',linewidth,'LineStyle','--');
            elseif j==5
                plot(t,x(:,i)+c1,'LineWidth',linewidth,'LineStyle','-.');
            elseif j==13
                plot(t,x(:,i)+c1,'LineWidth',linewidth,'LineStyle',':');
            end
        elseif i==3
            if j==1
                plot(t,x(:,i)+c2,'LineWidth',linewidth,'LineStyle','--');
            elseif j==5
                plot(t,x(:,i)+c2,'LineWidth',linewidth,'LineStyle','-.');
            elseif j==13
                plot(t,x(:,i)+c2,'LineWidth',linewidth,'LineStyle',':');
            end
        elseif i==7
            if j==1
                plot(t,x(:,i)+c3,'LineWidth',linewidth,'LineStyle','--');
            elseif j==5
                plot(t,x(:,i)+c3,'LineWidth',linewidth,'LineStyle','-.');
            elseif j==13
                plot(t,x(:,i)+c3,'LineWidth',linewidth,'LineStyle',':');
            end
        else
            if j==1
                plot(t,x(:,i),'LineWidth',linewidth,'LineStyle','--');
            elseif j==5
                plot(t,x(:,i),'LineWidth',linewidth,'LineStyle','-.');
            elseif j==13
                plot(t,x(:,i),'LineWidth',linewidth,'LineStyle',':');
            end
        if i==5
            yticks([-10,0,6,10]);
        end
        end
    end
    end
end

for i=1:11     
    figure(i);
    hold on;
    legend('nonlinear','1st','5th','13th','interpreter','latex','Box','off');
    hold off;
end