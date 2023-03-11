%% Initialization
clear, clc;
format long G;
global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p m Vt0 xdSigma xdSigma_p;
parameter();
%% Task A1: Linear delta, omega, Eq (closed-loop)
global K K2 dU_ref U_ref V_disturb;

dVf=0.05;
WhetherToPlotPowerAngleCurve=0;
tspan=[0,100];
V_disturb=0;

workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
system_matrix=Linearization(workingpoint);
A=system_matrix(:,1:3);B=system_matrix(:,4);C=system_matrix(:,5);
[t_linear,y]=ode45(@(t,y) gen_linear(A,B,t,y,dVf),tspan,[0;0;0]);
delta_t=y(:,1)+workingpoint(1);w_t=y(:,2)+workingpoint(2);
Eq_p_t=y(:,3)+workingpoint(3);Vt_t=y*C+workingpoint(4);
x_linear=[delta_t,w_t,Eq_p_t,Vt_t];

C1 = [1 0 0];
C2 = [0 1 0];
C3 = [0 0 1];
syms delta omega Eq_p
syms ddelta domega dEq_p
Pem=m*(Eq_p*Vs)*sin(delta)/xdSigma_p;
dPem=diff(Pem,delta)*ddelta+diff(Pem,omega)*domega+diff(Pem,Eq_p)*dEq_p;
C_Pem=[diff(dPem,ddelta),diff(dPem,domega),diff(dPem,dEq_p)];
C_Pem=double(subs(C_Pem,[delta Eq_p],[workingpoint(1) workingpoint(3)]));

%% Initial calculation of transfer function
s=tf('s');
sys_generator_Vt=(0.07552*s^2 + 0.00151*s + 0.009872)/(s^3 + 0.572*s^2 + 0.1879*s + 0.04846);
sys_generator_delta = tf(ss(A,B,C1,0));
sys_generator_omega = tf(ss(A,B,C2,0));
sys_generator_Eq_ = tf(ss(A,B,C3,0));
sys_generator_Pem = tf(ss(A,B,C_Pem,0));

sys_generator_num=cell2mat(sys_generator_Vt.num);
sys_generator_den=cell2mat(sys_generator_Vt.den);
sys_generator_omega_num=cell2mat(sys_generator_omega.num);
sys_generator_delta_num=cell2mat(sys_generator_delta.num);
sys_generator_Pem_num=cell2mat(sys_generator_Pem.num);

K=100;
K2=1;
sys_gain=K;
sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);
sys_AVR_omega_zero=zero(sys_AVR_omega);
sys_AVR_omega_pole=pole(sys_AVR_omega);
sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
sys_AVR_delta_zero=zero(sys_AVR_delta);
sys_AVR_delta_pole=pole(sys_AVR_delta);
sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
sys_AVR_Pem_zero=zero(sys_AVR_Pem);
sys_AVR_Pem_pole=pole(sys_AVR_Pem);

%% Task 1 Plotting Bode diagrams of "lead-lag link controller and PID controller"
% s=tf('s');
% G_PID=((1/0.7)*s+1)*((1/2)*s+1)/(s);
% G_PID_zero=zero(G_PID);
% G_PID_pole=pole(G_PID);
% G_lead_lag=(((1/2)*s+1)/((1/100)*s+1))*(((1/0.7)*s+1)/((1/0.014)*s+1));
% G_lead_lag_zero=zero(G_lead_lag);
% G_lead_lag_pole=pole(G_lead_lag);
% G_amplify=100*s/s;
% % Bode diagram
% figure;
% bode(G_lead_lag);
% grid on;
% hold on;
% bode(G_PID);
% bode(G_amplify);
% hold off;
% legend("超前-滞后环节控制器","PID控制器","放大环节控制器");
% title("The Bode diagram of $G_{lead-lag}(s)$, $G_{PID}(s)$,  $G_{amplify}(s)$",'FontSize',16,'Interpreter','latex');

%% Task 2 lead-lag correction + AVR --> f domain character, magnitude margin, phase angle margin
% figure;
margin(sys_AVR_omega);
grid on;
zpk(sys_AVR_omega);
% 配置串联超前滞后校正
omega_turn=[0.3576 8.098];
omega_c=0.6;
beta_T2=1/(omega_c/5);
beta=50;
T2=beta_T2*beta;
alpha_T1=1/omega_turn(1);
T1=alpha_T1/beta;
s=tf('s');
corr_lag=(beta_T2*s+1)/(T2*s+1);
corr_lead=(alpha_T1*s+1)/(T1*s+1);
% k_corr=1.03;
k_corr=1.0278;
G_corr=k_corr*series(corr_lag,corr_lead);
zpk(G_corr);
% figure;
margin(G_corr);
grid on;
sys_AVR_corr_omega = series(sys_AVR_omega,G_corr);
zpk(sys_AVR_corr_omega);
% figure;
margin(sys_AVR_corr_omega);
grid on;

%% 校正前的时域仿真
% before corr nonlinear model --> w and Vt
K=5;K2=1;
U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(2);
y0=workingpoint(1:3);
tspan=[0 300];
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_PSS(t,y,"omega"),tspan,y0);
% figure;
plot(t_nonlinear,y(:,2),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
legend('$\omega$','Interpreter','latex','Location','southeast');
title(['$\omega$' '(Noninear, before corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_PSS(t,y,"omega"),tspan,y0);
Vt_t=Vt_observer(y(:,1),y(:,3));
% figure;
plot(t_nonlinear,Vt_t,'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$V_{t}$','FontSize',12,'Interpreter','Latex');
title(['$V_{t}$' '(Noninear, before corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

% before corr linear model --> w and Vt
K=100;K2=1;dU_ref=0.05;
sys_AVR_PSS_omega=feedback(sys_AVR_omega,1);
[domega,t] = step(sys_AVR_PSS_omega,300);
% figure;
plot(t,dU_ref*domega+workingpoint(2),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
title(['$\omega$' '(Linear, before corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

sys_AVR_PSS_Vt=(K*sys_generator_Vt)/(1+K*sys_generator_Vt+K*sys_generator_omega);
[dVt,t] = step(sys_AVR_PSS_Vt,300);
% figure;
plot(t,dU_ref*dVt,'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
title(['$\Delta V_{t}$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

%% 校正后的时域仿真
% % after corr linear model --> w and Vt
% 
% % 施加0.05阶跃, 看系统频率的变化量 +-0.2Hz
% % figure;
% % nyquist(sys_AVR_corr_omega);
% % title("The Nyquist diagram of $G_{0\omega}(s)$",'FontSize',16,'Interpreter','latex');
% sys_AVR_PSS_corr_omega=feedback(sys_AVR_corr_omega,1);
% dU_ref=0.05;
% K=100;
% [domega,t] = step(sys_AVR_PSS_corr_omega,1000);
% % figure;
% plot(t,dU_ref*domega+workingpoint(2),'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
% title(['$\omega$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);
% 
% % figure;
% plot(t,50*(dU_ref*domega+workingpoint(2)),'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$f/Hz$','FontSize',12,'Interpreter','Latex');
% title(['$f$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);
% 
% % 施加0.05阶跃, 看系统机端电压变化量
% sys_AVR_PSS_corr_Vt=(K*G_corr*sys_generator_Vt)/(1+K*sys_generator_Vt+G_corr*K*sys_generator_omega);
% zpk(sys_AVR_PSS_corr_Vt);
% [dVt,t] = step(sys_AVR_PSS_corr_Vt,3000);
% % figure;
% plot(t,dU_ref*dVt+workingpoint(4),'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
% title(['$V_{t}$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);
% 
% % figure;
% plot(t,dU_ref*dVt,'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
% title(['$\Delta V_{t}$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

%% Task 3 Optimize parameters 稳准快
% lead-lag correction + AVR --> f domain character, magnitude margin, phase angle margin
% 配置串联超前滞后校正
sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);
omega_turn=[0.3576 8.098];
% omega_c=0.357;
% beta_T2=1/(omega_c/5);
beta_T2=1/0.5;
beta=1/50;
T2=beta_T2/beta;
% alpha=1/beta;
alpha=60;
alpha_T1=1/0.6;
T1=alpha_T1/alpha;
s=tf('s');
corr_lag=(beta_T2*s+1)/(T2*s+1);
corr_lead=(alpha_T1*s+1)/(T1*s+1);
% corr_lead=(alpha_T1*s+1)/1;
% k_corr=1.0278;
k_corr=1.045;
G_corr=k_corr*series(corr_lag,corr_lead);
% figure;
margin(G_corr);
grid on;
sys_AVR_corr_omega = series(sys_AVR_omega,G_corr);
% figure;
margin(sys_AVR_corr_omega);
grid on;

%% 校正后的时域仿真
% after corr linear model --> w and Vt

% 施加0.05阶跃, 看系统频率的变化量 +-0.2Hz
sys_AVR_PSS_corr_omega=feedback(sys_AVR_corr_omega,1);
dU_ref=0.05;
K=100;
[domega,t] = step(sys_AVR_PSS_corr_omega,1000);
% figure;
% plot(t,50*(dU_ref*domega+workingpoint(2)),'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$f/Hz$','FontSize',12,'Interpreter','Latex');
% title(['$f$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

% figure;
plot(t,dU_ref*domega+workingpoint(2),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
title(['$\omega$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

% 施加0.05阶跃, 看系统机端电压变化量
sys_AVR_PSS_corr_Vt=(K*G_corr*sys_generator_Vt)/(1+K*sys_generator_Vt+G_corr*K*sys_generator_omega);
[dVt,t] = step(sys_AVR_PSS_corr_Vt,1000);
% figure;
% plot(t,dU_ref*dVt,'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
% title(['$\Delta V_{t}$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

figure;
plot(t,dU_ref*dVt+workingpoint(4),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
title(['$V_{t}$' '(Linear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

% after corr nonlinear model --> w and Vt
K=100;K2=1;
U_ref = (workingpoint(4)+Vf/K)/dcgain(G_corr)+K2.*workingpoint(2); % 注意修改
y0=workingpoint(1:3);
tspan=[0 1000];
[A_corr,B_corr,C_corr,D_corr] = tf2ss(G_corr.num{1}, G_corr.den{1});
x0 =-inv(A_corr)*B_corr*(U_ref-workingpoint(2));
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_PSS_corr(G_corr,t,y,"omega"),tspan,[y0';x0]);
figure;
plot(t_nonlinear,y(:,2),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
legend('$\omega$','Interpreter','latex','Location','southeast');
title(['$\omega$' '(Noninear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);

Vt_t=Vt_observer(y(:,1),y(:,3));
figure;
plot(t_nonlinear,Vt_t,'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$V_{t}$','FontSize',12,'Interpreter','Latex');
title(['$V_{t}$' '(Noninear, after corr, when' '$\Delta U_{ref} = $' num2str(dU_ref) ')'], 'Interpreter','latex','FontSize',16);


