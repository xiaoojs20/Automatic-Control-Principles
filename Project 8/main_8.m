%% Initialization
clear, clc;
format long G;
global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p m Vt0 xdSigma xdSigma_p workingpoint;
parameter();

%% Task A1: Linear delta, omega, Eq (closed-loop)
global K K2 dU_ref U_ref V_disturb;

K=1;
dU_ref=0.05;
dVf=0.05;
K2=1;
V_disturb=0;

WhetherToPlotPowerAngleCurve=0;
tspan=[0,100];
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
syms ddelta dVT dEq_p
Pem=m*(Eq_p*Vs)*sin(delta)/xdSigma_p;
dPem=diff(Pem,delta)*ddelta+diff(Pem,omega)*dVT+diff(Pem,Eq_p)*dEq_p;
C_Pem=[diff(dPem,ddelta),diff(dPem,dVT),diff(dPem,dEq_p)];
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

%% Task 1：延时环节
%% Task 1.1: 一阶惯性环节
taus=[0.1 1 10 100];
% Nyquist图判断稳定性
% figure;
% for i=1:4
%     tau=taus(i);
%     sys_delay=1/(tau*s+1);
%     sys_delay_Vt=100*sys_delay*sys_generator_Vt;
%     subplot(2,2,i);
%     nyquist(sys_delay_Vt);
%     title(['$\tau_0=$' num2str(tau)], 'Interpreter','latex');
%     sgtitle("The Nyquist diagram of $G_{V_t}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end
% % 时域仿真验证稳定性
% tlim=300;
% for i=1:2
%     tau=taus(i);
%     sys_delay=1/(tau*s+1);
%     sys_delay_Vt=100*sys_delay*sys_generator_Vt;
%     sys_AVR_delay_Vt=feedback(sys_delay_Vt,1);
%     [dVt,t] = step(sys_AVR_delay_Vt,tlim);
%     figure;
%     plot(t,dU_ref*dVt,'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
%     title(['$\Delta V_{t}$' '(Linear, when' ' Timedelay = ' num2str(tau) ')'], 'Interpreter','latex','FontSize',16);
% end

% Timedelay=100, 用根轨迹找MAX K
K=329;
tau=100;
sys_delay=1/(tau*s+1);
sys_delay_Vt=K*sys_delay*sys_generator_Vt;
zpk(sys_delay_Vt);
% figure;
% rlocus(sys_delay_Vt);

% % 时域仿真MAX K~=300
% tlim=5000;
% sys_AVR_delay_Vt=feedback(sys_delay_Vt,1);
% [dVt,t] = step(sys_AVR_delay_Vt,tlim);
% figure;
% plot(t,dU_ref*dVt,'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
% title(['$\Delta V_{t}$' '(Linear, when' ' K = ' num2str(K) ')'], 'Interpreter','latex','FontSize',16);

%% Task 1.2: 纯延时环节
K=100;
sys_puredelay=exp(-tau*s);
sys_puredelay_Vt=K*sys_puredelay*sys_generator_Vt;
sys_AVR_puredelay_Vt=feedback(sys_puredelay_Vt,1);

% Nyquist图判断稳定性
% taus=[0.1 1 10 100];
% figure;
% for i=1:2
%     tau=taus(i);
%     sys_puredelay=exp(-tau*s);
%     sys_puredelay_Vt=K*sys_puredelay*sys_generator_Vt;
%     subplot(1,2,i);
%     nyquist(sys_puredelay_Vt);
%     title(['$\tau_0=$' num2str(tau)], 'Interpreter','latex');
%     sgtitle("The Nyquist diagram of $G_{V_t}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end
% % 时域仿真验证稳定性
% tlim=300;
% for i=1:2
%     tau=taus(i);
%     sys_puredelay=exp(-tau*s);
%     sys_puredelay_Vt=K*sys_puredelay*sys_generator_Vt;
%     sys_AVR_puredelay_Vt=feedback(sys_puredelay_Vt,1);
%     [dVt,t] = step(sys_AVR_puredelay_Vt,tlim);
%     figure;
%     plot(t,dU_ref*dVt,'LineWidth',2);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
%     title(['$\Delta V_{t}$' '(Linear, when' ' Timedelay = ' num2str(tau) ')'], 'Interpreter','latex','FontSize',16);
% end

% Timedelay=1, 用根轨迹找MAX K
tau=1;
sys_puredelay=exp(-tau*s);
sys_puredelay_Vt=pade(sys_puredelay)*sys_generator_Vt;
% figure;
% rlocus(sys_puredelay_Vt);
% figure;
% margin(25.46*sys_puredelay*sys_generator_Vt)

% % 时域仿真MAX K~=300
% K=25;
% tlim=300;
% sys_AVR_puredelay_Vt=feedback(K*sys_puredelay*sys_generator_Vt,1);
% [dVt,t] = step(sys_AVR_puredelay_Vt,tlim);
% 
% figure;
% plot(t,dU_ref*dVt,'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\Delta V_{t}$','FontSize',12,'Interpreter','Latex');
% title(['$\Delta V_{t}$' '(Linear, when' ' K = ' num2str(K) ')'], 'Interpreter','latex','FontSize',16);

%% Task 2: 未知扰动
%% Task 2.1: 开环系统的传递函数±10%的波动 -->Nyquist
% % AVR
% K=3.3;
% sys_AVR=K*sys_generator_Vt;
% Ks=[0.9,1,1.1];
% figure;
% for i=1:3
%     subplot(1,3,i);
%     nyquist(Ks(i)*sys_AVR);
%     title(['$amp=$' num2str(Ks(i))], 'Interpreter','latex');
%     sgtitle("The Nyquist diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end

% AVR+PSS
K=100;
% 配置串联超前滞后校正
beta_T2=1/0.5;
beta=1/50;
T2=beta_T2/beta;
alpha=60;
alpha_T1=1/0.6;
T1=alpha_T1/alpha;
s=tf('s');
corr_lag=(beta_T2*s+1)/(T2*s+1);
corr_lead=(alpha_T1*s+1)/(T1*s+1);
k_corr=1.045;
G_corr=k_corr*series(corr_lag,corr_lead);

sys_AVR_Vt=feedback(K,sys_generator_Vt);
sys_AVR_PSS_omega=G_corr*sys_AVR_Vt*sys_generator_omega;

Ks=[0.9,1,1.1];
figure;
for i=1:3
    subplot(1,3,i);
    nyquist(Ks(i)*sys_AVR_PSS_omega);
    title(['$amp=$' num2str(Ks(i))], 'Interpreter','latex');
    sgtitle("The Nyquist diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
end
figure;
margin(sys_AVR_PSS_omega);

%% Task 2.2: 偏差的幅频特性的大小不超过 0.4 -->Bode
% % AVR
% K=3.3;
% amp=[10^(-0.4/20), 1, 10^(0.4/20)];
% sys_AVR=K*sys_generator_Vt;
% figure;
% for i=1:3
%     subplot(1,3,i);
%     margin(amp(i)*sys_AVR);
%     title(['$amp=$' num2str(amp(i))], 'Interpreter','latex');
%     sgtitle("The Bode diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
%     % nyquist(Ac*sys_AVR_omega);
% end
% figure;
% margin(sys_AVR);

% AVR+PSS
K=100;
amp=[10^(-0.4/20), 1, 10^(0.4/20)];
sys_AVR_PSS_omega=G_corr*sys_AVR_omega;
figure;
for i=1:3
    subplot(1,3,i);
    margin(amp(i)*sys_AVR_PSS_omega);
    title(['$amp=$' num2str(amp(i))], 'Interpreter','latex');
    sgtitle("The Bode diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
    % nyquist(Ac*sys_AVR_PSS_omega);
end
figure;
margin(sys_AVR_PSS_omega);

%% Task 2.3: 偏差的幅频特性的大小不超过原系统的 20% -->Bode
% Ac=1.2;
% K=3.3;
% sys_AVR_omega=K*sys_generator_omega;
% figure;
% margin(Ac*sys_AVR_omega);
% % nyquist(Ac*sys_AVR_omega);
% % AVR+PSS
% K=100;
% sys_AVR_PSS_omega=K*G_corr*sys_AVR_Vt*sys_generator_omega;
% figure;
% margin(Ac*sys_AVR_PSS_omega);
% % nyquist(Ac*sys_AVR_PSS_omega);









