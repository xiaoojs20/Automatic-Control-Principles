%% Initialization
clear, clc;
format long G;
global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p m Vt0 xdSigma xdSigma_p;
parameter();
%% Task A1: Stability analysis of open-loop systems
% Assignment2 transfer function(Vf------>Vt)
% sys4 =
%
%    0.07552 s^2 + 0.00151 s + 0.009872
%   ------------------------------------
%   s^3 + 0.572 s^2 + 0.1879 s + 0.04846

s = tf('s');
fprintf("~~~~~~~> transfer function <~~~~~~~\n");
sys_generator = (0.07552*s^2 + 0.00151*s + 0.009872)/(s^3 + 0.572*s^2 + 0.1879*s + 0.04846);
Routhtable = Routh_table(sys_generator.Denominator{1});

%% Task A2: Stability performance analysis of open-loop systems
% figure;
[dVt,t] = step(sys_generator,100);
for i=1:3
    subplot(1,3,i);
    plot(t,i*0.05*dVt+Vt0(2),'LineWidth',2);
    ylim([0.975 1.015]);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    title(['$\Delta V_{f} = $' num2str(i*0.05) 'V'], 'Interpreter','latex');
    sgtitle(['The effect of different' ' $\Delta V_{f}$ ' 'on' ' $V_{t}$ '...
        '(transfer function)'],'Interpreter','Latex');
end

%% Task A2: Linear Model
dVf=0.05;
WhetherToPlotPowerAngleCurve=0;
tspan=[0,100];

workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
system_matrix=Linearization(workingpoint);
A=system_matrix(:,1:3);B=system_matrix(:,4);C=system_matrix(:,5);
[t_linear,y]=ode45(@(t,y) gen_linear(A,B,t,y,dVf),tspan,[0;0;0]);
delta_t=y(:,1)+workingpoint(1);w_t=y(:,2)+workingpoint(2);
Eq_p_t=y(:,3)+workingpoint(3);Vt_t=y*C+workingpoint(4);
x_linear=[delta_t,w_t,Eq_p_t,Vt_t];

% figure;
for Vf_factor = 1:3
    workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
    system_matrix=Linearization(workingpoint);
    A=system_matrix(:,1:3);B=system_matrix(:,4);C=system_matrix(:,5);
    [t_linear,y]=ode45(@(t,y) gen_linear(A,B,t,y,dVf*Vf_factor),tspan,[0;0;0]);
    delta_t=y(:,1)+workingpoint(1);w_t=y(:,2)+workingpoint(2);
    Eq_p_t=y(:,3)+workingpoint(3);Vt_t=y*C+workingpoint(4);
    x_linear=[delta_t,w_t,Eq_p_t,Vt_t];

    subplot(1,3,Vf_factor);
    plot(t_linear,x_linear(:,4),'LineWidth',2);
    grid on;
    ylim([0.975 1.015]);
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    title(['$\Delta V_{f} = $' num2str(Vf_factor*0.05) 'V'], 'Interpreter','latex');
    sgtitle(['The effect of different' ' $\Delta V_{f}$ ' 'on' ' $V_{t}$ '...
        '(linear model)'],'Interpreter','Latex');
end

%% Task A2: Nonlinear Model
tspan=[0,100];
dVf=0.05;
WhetherToPlotPowerAngleCurve=0;
workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
y0=workingpoint(1:3);
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear(t,y,Vf+dVf),tspan,y0);
Vt_t=Vt_observer(y(:,1),y(:,3));
x_nonlinear=[y,Vt_t];

% figure;
for Vf_factor = 1:3
    workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
    y0=workingpoint(1:3);
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear(t,y,Vf+dVf*Vf_factor),tspan,y0);
    Vt_t=Vt_observer(y(:,1),y(:,3));
    x_nonlinear=[y,Vt_t];
    subplot(1,3,Vf_factor);
    plot(t_nonlinear,x_nonlinear(:,4),'LineWidth',2);
    grid on;
    ylim([0.975 1.015]);
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    title(['$\Delta V_{f} = $' num2str(Vf_factor*0.05) 'V'], 'Interpreter','latex');
    sgtitle(['The effect of different' ' $\Delta V_{f}$ ' 'on' ' $V_{t}$ '...
        '(nonlinear model)'],'Interpreter','Latex');
end

%% Task A3: Dynamic performance analysis of open-loop systems
% Analysis overshoot, damping_factor, stable_time
%% Task A3: Linear Model
dVf=1;
tspan=[0,100];

[t_linear,y]=ode45(@(t,y) gen_linear(A,B,t,y,dVf),tspan,[0;0;0]);
delta_t=y(:,1)+workingpoint(1);w_t=y(:,2)+workingpoint(2);
Eq_p_t=y(:,3)+workingpoint(3);Vt_t=y*C+workingpoint(4);
x_linear=[delta_t,w_t,Eq_p_t,Vt_t];
x_linear_no_initial=[y(:,1) y(:,2) y(:,3) y*C];

% [overshoot_linear, damping_factor_linear, stable_time_linear] = ...
%     dynamic_performance(t_linear,x_linear_no_initial(:,4),0);

%% Task A3: Nonlinear Model
tspan=[0,300];
y0=workingpoint(1:3);
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear(t,y,Vf+dVf),tspan,y0);
Vt_t=Vt_observer(y(:,1),y(:,3));
x_nonlinear=[y,Vt_t];
x_nonlinear_no_initial = [y,Vt_t-workingpoint(4)];

% [overshoot_nonlinear, damping_factor_nonlinear, stable_time_nonlinear] = ...
%     dynamic_performance(t_nonlinear,x_nonlinear_no_initial(:,4),1);

%% Task B1: Closed-loop system modeling
dVf=0;
global K dU_ref U_ref V_disturb;
V_disturb = 0;
K_list = [0.1 0.5 1 5 20 100];

%% Task B12: Linear Model closed-loop
tspan=[0,100];
system_matrix=Linearization(workingpoint);
A=system_matrix(:,1:3);B=system_matrix(:,4);C=system_matrix(:,5);
[~,y]=ode45(@(t,y) gen_linear(A,B,t,y,dVf),tspan,[0;0;0]);
delta_t=y(:,1)+workingpoint(1);w_t=y(:,2)+workingpoint(2);
Eq_p_t=y(:,3)+workingpoint(3);Vt_t=y*C+workingpoint(4);
x_linear=[delta_t,w_t,Eq_p_t,Vt_t];

dU_ref = 0.05;
% figure;
for K = K_list
    U_ref = workingpoint(4)+Vf/K;
    fprintf("~~~~~~~> K=%f, U_ref=%f, dU_ref=%f <~~~~~~~\n",K,U_ref,dU_ref);
    sys_generator;
    sys_gain = K;
    sys_gain_generator =  series(sys_gain,sys_generator);
    sys_gain_generator_feedback = feedback(sys_gain_generator,1);
    [Vt_linear_feedback_step,t_linear] = step(sys_gain_generator_feedback,200);
    
    subplot(2,3,find(K_list==K));
    plot(t_linear,dU_ref*Vt_linear_feedback_step+workingpoint(4),'LineWidth',2);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Linear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of different K on' ' $V_{t}$ '...
        '(linear model)'],'Interpreter','Latex');
%     RouthTable_linear_closedloop = Routh_table(sys_gain_generator_feedback.Denominator{1});
end

%% Task B2: Linear Model open-loop
% figure;
[dVt,t] = step(sys_generator,100);
plot(t,0.05*dVt+workingpoint(4),'LineWidth',2);
% ylim([0.975 1.015]);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
title(['$\Delta V_{f} = 0.05V$' '(linear model open-loop)'],...
    'Interpreter','latex','FontSize',16);

%% Task B12: Nonlinear Model
tspan=[0,100];
% figure;
for K = K_list
    U_ref = workingpoint(4)+Vf/K;
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
    Vt_t=Vt_observer(y(:,1),y(:,3));
    
    subplot(2,3,find(K_list==K));
    plot(t_nonlinear,Vt_t,'LineWidth',2);
%     ylim([0.977 1.03]);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of different K on' ' $V_{t}$ '...
        '(nonlinear model)'],'Interpreter','Latex');
end

%% Task B2: Nonlinear Model open-loop
dVf=0.05;
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear(t,y,Vf+dVf),tspan,y0);
Vt_t=Vt_observer(y(:,1),y(:,3));
x_nonlinear=[y,Vt_t];

% figure;
plot(t_nonlinear,x_nonlinear(:,4),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
title(['$\Delta V_{f} = 0.05V$' '(nonlinear model open-loop)'],...
    'Interpreter','latex','FontSize',16);

%% Task B3: Linear Model closed-loop
tspan=[0,100];
K = 5;
dU_ref_list = [0.05 0.10 0.15 0.5];
% figure;
for dU_ref = dU_ref_list
    U_ref = workingpoint(4)+Vf/K;
    fprintf("~~~~~~~> K=%f, U_ref=%f, dU_ref=%f <~~~~~~~\n",K,U_ref,dU_ref);
    sys_generator;
    sys_gain = K;
    sys_gain_generator =  series(sys_gain,sys_generator);
    sys_gain_generator_feedback = feedback(sys_gain_generator,1);
    [Vt_linear_feedback_step,t_linear] = step(sys_gain_generator_feedback,100);

    G0=dcgain(sys_gain_generator);
    error_steady_state = dU_ref*1/(1+G0)
    
    subplot(2,2,find(dU_ref_list==dU_ref));
    plot(t_linear,dU_ref*Vt_linear_feedback_step+workingpoint(4),'LineWidth',2);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Linear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of different ' '$\Delta U_{ref}$' ' on' ' $V_{t}$ '...
        '(linear model)'],'Interpreter','Latex');
    RouthTable_linear_closedloop = Routh_table(sys_gain_generator_feedback.Denominator{1});
end


%% Task B3: Linear Model open-loop
dVfs = [0.05 0.10 0.15 0.5];
% figure;
for dVf = dVfs
    [dVt,t] = step(sys_generator,100);
    subplot(2,2,find(dVf==dVfs));
    plot(t,dVf*dVt+workingpoint(4),'LineWidth',2);
    % ylim([0.975 1.015]);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    title(['$\Delta V_{f} = $' num2str(dVf) 'V (linear model open-loop)'],...
        'Interpreter','latex');
    sgtitle(['The effect of different ' '$\Delta V_{f}$' ' on' ' $V_{t}$ '...
        '(linear model open-loop)'],'Interpreter','Latex');
end

%% Task B3: Nonlinear Model closed-loop
tspan=[0,100];
% figure;
for dU_ref = dU_ref_list
    U_ref = workingpoint(4)+Vf/K;
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
    Vt_t=Vt_observer(y(:,1),y(:,3));    
    subplot(2,2,find(dU_ref_list==dU_ref));
    plot(t_nonlinear,Vt_t,'LineWidth',2);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of different ' '$\Delta U_{ref}$' ' on' ' $V_{t}$ '...
        '(nonlinear model)'],'Interpreter','Latex');
end

%% Task B4: Anti-interference analysis
%% Task B4: Linear Model closed-loop
K = 5;
U_ref = workingpoint(4)+Vf/K;
dU_ref = 0;
% figure;
for dVf = dVfs
    sys_generator;
    sys_gain = K;
    sys_gain_generator =  series(sys_gain,sys_generator);
    sys_gain_generator_feedback = feedback(sys_gain_generator,1);
    sys_gain_generator_feedback_disturb = sys_gain_generator_feedback/sys_gain;

    [Vt_linear_feedback_disturb_step,t_linear] = step(sys_gain_generator_feedback_disturb,100);
    
    subplot(2,2,find(dVf==dVfs));
    plot(t_linear,dVf*Vt_linear_feedback_disturb_step+workingpoint(4),'LineWidth',2);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Linear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta V_{f}=$' num2str(dVf) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of disturbance ' '$\Delta V_{f}$' ' on' ' $V_{t}$ '...
        '(linear model)'],'Interpreter','Latex');
end

%% Task B4: Nonlinear Model closed-loop
% figure;
for V_disturb = dVfs
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
    Vt_t=Vt_observer(y(:,1),y(:,3));    
    subplot(2,2,find(V_disturb==dVfs));
    plot(t_nonlinear,Vt_t,'LineWidth',2);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta V_{f}=$' num2str(V_disturb) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of disturbance ' '$\Delta V_{f}$' ' on' ' $V_{t}$ '...
        '(nonlinear model)'],'Interpreter','Latex');
end

V_disturb = 0;

%% Task B5: Dynamic performance analysis of AVR closed-loop systems

%% Task B5: Linear Model
K = 5;
dU_ref = 1;
U_ref = workingpoint(4)+Vf/K;

sys_gain = K;
sys_gain_generator =  series(sys_gain,sys_generator);
sys_gain_generator_feedback = feedback(sys_gain_generator,1);
[Vt_linear_feedback_step,t_linear] = step(sys_gain_generator_feedback,200);

% [overshoot_linear_closedloop, damping_factor_linear_closedloop, stable_time_linear_closedloop]=...
%     dynamic_performance(t_linear,Vt_linear_feedback_step,0);

%% Task B5: Nonlinear Model
tspan=[0,300];
K = 5;
U_ref = workingpoint(4)+Vf/K;
dU_ref = 1;
V_disturb = 0;
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
Vt_nonlinear_feedback_step=Vt_observer(y(:,1),y(:,3));

% [overshoot_nonlinear_closedloop, damping_factor_nonlinear_closedloop, stable_time_nonlinear_closedloop]=...
%     dynamic_performance(t_nonlinear,Vt_nonlinear_feedback_step-workingpoint(4),1);

%% Project 4
%% Task A: AVR closed-loop Output Vt 
%% Task A12: Root tracks(Output Vt)
zero_sys_generator = zero(sys_generator);
pole_sys_generator = pole(sys_generator);
% figure;
% rlocus(sys_generator);
% rlocfind(sys_generator);
%% Task A3: Time domain simulation, nonlinear closed-loop
U_ref = workingpoint(4)+Vf/K;
dU_ref = 0.05;
V_disturb = 0;

%% Task A3: Vt for different K, nonlinear closed-loop
tspan=[0,200];
figure;
for K = K_list
    U_ref = workingpoint(4)+Vf/K;
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
    Vt_t=Vt_observer(y(:,1),y(:,3));
    
    subplot(2,3,find(K_list==K));
    plot(t_nonlinear,Vt_t,'LineWidth',2);
%     ylim([0.977 1.03]);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$','Interpreter','latex','Location','southeast');
    title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of different K on' ' $V_{t}$ '...
        '(nonlinear model closed loop)'],'Interpreter','Latex');
end

%% Task A3: omega for different K, nonlinear closed-loop
tspan=[0,300];
figure;
for K = K_list
    U_ref = workingpoint(4)+Vf/K;
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
    Vt_t=Vt_observer(y(:,1),y(:,3));
    subplot(2,3,find(K_list==K));
    plot(t_nonlinear,y(:,2),'LineWidth',2);
    ylim([0.95 1.05]);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$\omega /(rad/s)$','FontSize',12,'Interpreter','Latex');
    legend('$\omega$','Interpreter','latex','Location','southeast');
    title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
        ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
        'FontSize',10,'Interpreter','latex');
    sgtitle(['The effect of different K on' ' $\omega$ '...
        '(nonlinear model closed loop)'],'Interpreter','Latex');
end

%% Task B1: delta, omega, Eq (closed-loop)
% C1 = [1 0 0];
% C2 = [0 1 0];
% C3 = [0 0 1];
% sys_generator_delta = tf(ss(A,B,C1,0));
% sys_generator_omega = tf(ss(A,B,C2,0));
% sys_generator_Eq_ = tf(ss(A,B,C3,0));
% % sys_generator_Vt = tf(ss(A,B,C',0));
% syms delta omega Eq_p
% syms ddelta domega dEq_p
% Pem=m*(Eq_p*Vs)*sin(delta)/xdSigma_p;
% dPem=diff(Pem,delta)*ddelta+diff(Pem,omega)*domega+diff(Pem,Eq_p)*dEq_p;
% C_Pem=[diff(dPem,ddelta),diff(dPem,domega),diff(dPem,dEq_p)];
% C_Pem=double(subs(C_Pem,[delta Eq_p],[workingpoint(1) workingpoint(3)]));
% sys_generator_Pem = tf(ss(A,B,C_Pem,0));
% 
% sys_AVR=feedback(sys_gain,sys_generator);
% sys_AVR_delta=series(sys_AVR,sys_generator_delta);
% sys_AVR_delta_feedback=feedback(sys_AVR_delta,1);
% sys_AVR_omega=series(sys_AVR,sys_generator_omega);
% sys_AVR_omega_feedback=feedback(sys_AVR_omega,1);
% sys_AVR_Pem=series(sys_AVR,sys_generator_Pem);
% sys_AVR_Pem_feedback=feedback(sys_AVR_Pem,1);
% 
% % figure;
% for K=K_list
%     sys_gain=K;
%     sys_AVR=feedback(sys_gain,sys_generator);
%     sys_AVR_delta=series(sys_AVR,sys_generator_delta);
%     sys_AVR_delta_feedback=feedback(sys_AVR_delta,1);
%     subplot(2,3,find(K_list==K));
%     rlocus(sys_AVR_delta);
%     title(['The root trajectory of $\delta$, $K=$' num2str(K)],'FontSize',12,'Interpreter','latex');
%     sgtitle('The root trajectory of $\delta$','FontSize',16,'Interpreter','latex');
% end
% % figure;
% for K=K_list
%     sys_gain=K;
%     sys_AVR=feedback(sys_gain,sys_generator);
%     sys_AVR_omega=series(sys_AVR,sys_generator_omega);
%     sys_AVR_omega_feedback=feedback(sys_AVR_omega,1);
%     subplot(2,3,find(K_list==K));
%     rlocus(sys_AVR_omega);
%     title(['The root trajectory of $\omega$, $K=$' num2str(K)],'FontSize',12,'Interpreter','latex');
%     sgtitle('The root trajectory of $\omega$','FontSize',16,'Interpreter','latex');
% end
% % figure;
% for K=K_list
%     sys_gain=K;
%     sys_AVR=feedback(sys_gain,sys_generator);
%     sys_AVR_Pem=series(sys_AVR,sys_generator_Pem);
%     sys_AVR_Pem_feedback=feedback(sys_AVR_Pem,1);
%     subplot(2,3,find(K_list==K));
%     rlocus(sys_AVR_Pem);
%     title(['The root trajectory of $P_{em}$, $K=$' num2str(K)],'FontSize',12,'Interpreter','latex');
%     sgtitle('The root trajectory of $P_{em}$','FontSize',16,'Interpreter','latex');
% end

%% Task B2: delta, omega, Eq' for different K (Nonlinear model closed loop)
%% Task B2: delta for different K (Nonlinear model closed loop)
% tspan=[0,300];
% % figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
%     subplot(2,3,find(K_list==K));
%     plot(t_nonlinear,y(:,1),'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\delta /rad$','FontSize',12,'Interpreter','Latex');
%     legend('$\delta$','Interpreter','latex','Location','southeast');
%     title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
%         ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $\delta$ '...
%         '(nonlinear model closed loop)'],'Interpreter','Latex');
% end
  
%% Task B2: omega for different K (Nonlinear model closed loop)
% tspan=[0,300];
% % figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);    
%     subplot(2,3,find(K_list==K));
%     plot(t_nonlinear,y(:,2),'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\omega /(rad/s)$','FontSize',12,'Interpreter','Latex');
%     legend('$\omega$','Interpreter','latex','Location','southeast');
%     title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
%         ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $\omega$ '...
%         '(nonlinear model closed loop)'],'Interpreter','Latex');
% end

%% Task B2: Pem for different K (Nonlinear model closed loop)
% tspan=[0,300];
% % figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
%     Pem=Pem_observer(y(:,1),y(:,3));
%     
%     subplot(2,3,find(K_list==K));
%     plot(t_nonlinear,Pem,'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$P_{em} /W$','FontSize',12,'Interpreter','Latex');
%     legend('$P_{em}$','Interpreter','latex','Location','southeast');
%     title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
%         ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $P_{em}$ '...
%         '(nonlinear model closed loop)'],'Interpreter','Latex');
% end

%% Task B3: Root trajectory of the system with different outputs
% figure;
% rlocus(sys_AVR_delta);
% title("The root trajectory of $\delta$",'FontSize',16,'Interpreter','latex');
% figure;
% rlocus(sys_AVR_omega);
% title("The root trajectory of $\omega$",   'FontSize',16,'Interpreter','latex');
% figure;
% rlocus(sys_AVR_Pem);
% title("The root trajectory of $P_{em}$",'FontSize',16,'Interpreter','latex');

%% Task B3: Time domain simulation (Nonlinear model closed-loop double feedback)
%% Task B3: delta (Nonlinear model closed-loop double feedback)
% global K2
% tspan=[0,200];
% K=5;
% dU_ref=0.05;
% K2=5;
% figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(1);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,1),tspan,y0);
%     Vt_t=Vt_observer(y(:,1),y(:,3));
%     subplot(2,3,find(K_list==K));
%     plot(t_nonlinear,y(:,1),'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\delta/rad$','FontSize',12,'Interpreter','Latex');
%     legend('$\delta$','Interpreter','latex','Location','southeast');
%     title(['$U_{ref}=$ ' num2str(U_ref) '$V$,' ' $\Delta U_{ref}=$' num2str(dU_ref) '$V$' ...
%         ' $K=$ ' num2str(K) ' $K_2=$ ' num2str(K2)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $\delta$ '...
%         '(nonlinear model closed loop double feedback)'],'Interpreter','Latex');
% end

%% Task B3: omega (Nonlinear model closed-loop double feedback)
% K_list=[0.1 0.5 1 5];
% % figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(2);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,2),tspan,y0);
%     Vt_t=Vt_observer(y(:,1),y(:,3));
%     subplot(2,3,find(K_list==K));
%     plot(t_nonlinear,y(:,2),'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\omega/(rad/s)$','FontSize',12,'Interpreter','Latex');
%     legend('$\omega$','Interpreter','latex','Location','southeast');
%     title(['$U_{ref}=$ ' num2str(U_ref) '$V$,' ' $\Delta U_{ref}=$' num2str(dU_ref) '$V$' ...
%         ' $K=$ ' num2str(K) ' $K_2=$ ' num2str(K2)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $\omega$ '...
%         '(nonlinear model closed loop double feedback)'],'Interpreter','Latex');
% end

%% Task B3: Pem (Nonlinear model closed-loop double feedback)
% tspan=[0 300];
% K=0.1;
% K_list=[0.1 0.5 1 5 20 100];
% K2_list=[0.1 1 5 10 50 100];
% figure;
% for K2 = K2_list
%     Pem_workingpoint=Pem_observer(workingpoint(1),workingpoint(3));
%     U_ref = workingpoint(4)+Vf/K+K2.*Pem_workingpoint;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,3),tspan,y0);
%     Vt_t=Vt_observer(y(:,1),y(:,3));
%     Pem=Pem_observer(y(:,1),y(:,3));
%     subplot(2,3,find(K2_list==K2));
%     plot(t_nonlinear,Pem,'LineWidth',2);
%     ylim([0.799 0.8015]);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$P_{em}/W$','FontSize',12,'Interpreter','Latex');
%     legend('$P_{em}$','Interpreter','latex','Location','southeast');
%     title(['$U_{ref}=$ ' num2str(U_ref) '$V$,' ' $\Delta U_{ref}=$' num2str(dU_ref) '$V$' ...
%         ' $K=$ ' num2str(K) ' $K_2=$ ' num2str(K2)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $P_{em}$ '...
%         '(nonlinear model closed loop double feedback)'],'Interpreter','Latex');
% end


