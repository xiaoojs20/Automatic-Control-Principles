%% Initialization
clear, clc;
format long G;
global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p m Vt0 xdSigma xdSigma_p workingpoint;
parameter();

%% Task A1: Linear delta, omega, Eq (closed-loop)
global K K2 dU_ref U_ref V_disturb;

K=100;
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

s=tf('s');
sys_gain=K;
sys_generator_Vt=(0.07552*s^2 + 0.00151*s + 0.009872)/(s^3 + 0.572*s^2 + 0.1879*s + 0.04846);
sys_generator_delta = tf(ss(A,B,C1,0));
sys_generator_omega = tf(ss(A,B,C2,0));
sys_generator_Eq_ = tf(ss(A,B,C3,0));
sys_generator_Pem = tf(ss(A,B,C_Pem,0));
sys_AVR=feedback(sys_gain*sys_generator_Vt,1);

%% Initial calculation of transfer function
% s=tf('s');
% sys_generator_Vt=(0.07552*s^2 + 0.00151*s + 0.009872)/(s^3 + 0.572*s^2 + 0.1879*s + 0.04846);
% sys_generator_delta = tf(ss(A,B,C1,0));
% sys_generator_omega = tf(ss(A,B,C2,0));
% sys_generator_Eq_ = tf(ss(A,B,C3,0));
% sys_generator_Pem = tf(ss(A,B,C_Pem,0));
% 
% sys_generator_num=cell2mat(sys_generator_Vt.num);
% sys_generator_den=cell2mat(sys_generator_Vt.den);
% sys_generator_omega_num=cell2mat(sys_generator_omega.num);
% sys_generator_delta_num=cell2mat(sys_generator_delta.num);
% sys_generator_Pem_num=cell2mat(sys_generator_Pem.num);
% 
% K=100;
% K2=1;
% sys_gain=K;
% sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);
% sys_AVR_omega_zero=zero(sys_AVR_omega);
% sys_AVR_omega_pole=pole(sys_AVR_omega);
% sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
% sys_AVR_delta_zero=zero(sys_AVR_delta);
% sys_AVR_delta_pole=pole(sys_AVR_delta);
% sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
% sys_AVR_Pem_zero=zero(sys_AVR_Pem);
% sys_AVR_Pem_pole=pole(sys_AVR_Pem);

%% Task 1：基于线性模型，分析系统的能控性与能观性
A_AVR=A-K*B*C';
B_AVR=K*B;
C_AVR=C';
D_AVR = 0;
% 能控性
S=ctrb(A_AVR,B_AVR);
rank_S=rank(S);
if rank_S==3
    fprintf("~~~~~~~> rank(S)=%d, 系统能控 <~~~~~~~\n",rank_S);
else
    fprintf("~~~~~~~> rank(S)=%d, 系统不能控 <~~~~~~~\n",rank_S);
end
% 能观性
V=obsv(A_AVR,C_AVR);
rank_V=rank(V);
if rank_V==3
    fprintf("~~~~~~~> rank(V)=%d, 系统能观 <~~~~~~~\n",rank_V);
else
    fprintf("~~~~~~~> rank(V)=%d, 系统不能观 <~~~~~~~\n",rank_V);
end

%% Task 2：设计一个状态观测器
sys_Vt = ss(A_AVR,B_AVR,C_AVR,D_AVR);
[dVt,t] = step(sys_Vt,100);
figure;
dU_ref=0.05;
plot(t,dU_ref*dVt+workingpoint(4),'LineWidth',2);
grid on;
% hold on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$V_{t}$','FontSize',12,'Interpreter','Latex');
title(['$V_{t}$' '(Linear)'], 'Interpreter','latex','FontSize',16);

pole_obsv=[-5 -2+1j -2-1j];
G=place(A_AVR',C_AVR',pole_obsv)';
[num,den]=ss2tf(A_AVR,B_AVR,C_AVR,D_AVR);
G_AVR_Vt=ss2tf_fraction_from_NUM_DEN(num,den,s);
% vpa(G_AVR_Vt)

A_obsv=A_AVR-G*C_AVR;

%% Linear obsv
U_ref = workingpoint(4)+Vf/K;
tspan=[0 500];
y0_linear=zeros(6,1);
[t_linear,y_linear]=ode45(@(t,y) gen_linear_AVR_obsv(t,y,A_AVR,B_AVR,C_AVR,G),tspan,y0_linear);

%% Nonlinear obsv
K=100;
U_ref = workingpoint(4)+Vf/K;
tspan=[0 500];
A_obsv=A_AVR-G*C_AVR;
y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
[t_,y_]=ode45(@(t,y) gen_linear(A_AVR,B_AVR,t,y,dU_ref),tspan,[0 0 0]);
[t_nonlinear,y_nonlinear]=ode45(@(t,y) gen_nonlinear_AVR_obsv(t,y,A_obsv,B_AVR,G),tspan,y0_nonlinear);

%% xd changes 但是只针对被控系统，不针对观测器
xd=1.998-0.2;
xdSigma = xd + xtl;
xdSigma_p = xd_p + xtl;
Td_p=Td0_p*(xd_p+xtl)/(xd+xtl);
K=100;
dVf=0.05;
K2=1;
V_disturb=0;

workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
system_matrix=Linearization(workingpoint);
A=system_matrix(:,1:3);B=system_matrix(:,4);C=system_matrix(:,5);
A_AVR=A-K*B*C';B_AVR=K*B;C_AVR=C';D_AVR=0;
G=place(A_AVR',C_AVR',pole_obsv)';
A_obsv=A_AVR-G*C_AVR;
K=100;
U_ref = workingpoint(4)+Vf/K;
y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
[t_linear_xd,y_xd]=ode45(@(t,y) gen_linear_AVR_obsv(t,y,A_AVR,B_AVR,C_AVR,G),tspan,y0_linear);
[t_nonlinear_xd,y_nonlinear_xd]=ode45(@(t,y) gen_nonlinear_AVR_obsv(t,y,A_obsv,B_AVR,G),tspan,y0_nonlinear);

%% plot
% delta
figure;
plot(t_linear_xd,y_xd(:,1),'LineWidth',2);
hold on;
plot(t_linear,y_linear(:,4),'LineWidth',2);
hold off;
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
% legend('$\delta$ origin','$\delta$ obsv error','Interpreter','latex');
legend(['$\delta$ origin, xd=' num2str(xd)],'$\delta$ obsv error','Interpreter','latex');
title(['$\delta$ ' '(Linear, obsv)'], 'Interpreter','latex','FontSize',16);
% omega
figure;
plot(t_linear_xd,y_xd(:,2),'LineWidth',2);
hold on;
plot(t_linear,y_linear(:,5),'LineWidth',2);
hold off;
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
% legend('$\omega$ origin','$\omega$ obsv error','Interpreter','latex');
legend(['$\omega$ origin, xd=' num2str(xd)],'$\omega$ obsv error','Interpreter','latex');
title(['$\omega$ ' '(Linear, obsv)'], 'Interpreter','latex','FontSize',16);
% Eq'
figure;
plot(t_linear_xd,y_xd(:,3),'LineWidth',2);
hold on;
plot(t_linear,y_linear(:,6),'LineWidth',2);
hold off;
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$E_{q}$','FontSize',12,'Interpreter','Latex');
% legend('$E_{q}$ origin','$E_{q}$ obsv error','Interpreter','latex');
legend(['$E_{q}$ origin, xd=' num2str(xd)],'$E_{q}$ obsv error','Interpreter','latex');
title(['$E_{q}$ ' '(Linear, obsv)'], 'Interpreter','latex','FontSize',16);
% K_list=[1,10,100,500];
% plot_linear_AVR_obsv_K(K_list,A,B,C,G);

% pole_obsv_1=[-5 -2+1j -2-1j];
% pole_obsv_2=[-5 -2+8j -2-8j];
% pole_obsv_3=[-20 -2+1j -2-1j];
% pole_obsv_4=[-5 0+1j 0-1j];
% pole_list={pole_obsv_1,pole_obsv_2,pole_obsv_3,pole_obsv_4};
% plot_linear_AVR_obsv_pole(pole_list,A,B,C,G);

%% Nonlinear obsv
% K=100;
% U_ref = workingpoint(4)+Vf/K;
% tspan=[0 500];
% A_obsv=A_AVR-G*C_AVR;
% y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
% [t_,y_]=ode45(@(t,y) gen_linear(A_AVR,B_AVR,t,y,dU_ref),tspan,[0 0 0]);
% [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_obsv(t,y,A_obsv,B_AVR,G),tspan,y0_nonlinear);

% delta
figure;
plot(t_nonlinear_xd,y_nonlinear_xd(:,1)-workingpoint(1),'LineWidth',2);
hold on;
plot(t_,y_(:,1),'LineWidth',2);
plot(t_nonlinear,y_nonlinear(:,4),'LineWidth',2);
hold off;
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
legend(['$\delta$ origin, xd=' num2str(xd)],'$\delta$ obsv','$\delta$ obsv error','Interpreter','latex');
title(['$\delta$ ' '(Nonlinear, obsv)'], 'Interpreter','latex','FontSize',16);
% omega
figure;
plot(t_nonlinear_xd,y_nonlinear_xd(:,2)-workingpoint(2),'LineWidth',2);
hold on;
plot(t_,y_(:,2),'LineWidth',2);
plot(t_nonlinear,y_nonlinear(:,5),'LineWidth',2);
hold off;
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
legend(['$\omega$ origin, xd=' num2str(xd)],'$\omega$ obsv','$\omega$ obsv error','Interpreter','latex');
title(['$\omega$ ' '(Nonlinear, obsv)'], 'Interpreter','latex','FontSize',16);
% Eq'
figure;
plot(t_nonlinear_xd,y_nonlinear_xd(:,3)-workingpoint(3),'LineWidth',2);
hold on;
plot(t_,y_(:,3),'LineWidth',2);
plot(t_nonlinear,y_nonlinear(:,6),'LineWidth',2);
hold off;
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$E_{q}$','FontSize',12,'Interpreter','Latex');
legend(['$E_{q}$ origin, xd=' num2str(xd)],'$E_{q}$ obsv', '$E_{q}$ obsv error','Interpreter','latex');
title(['$E_{q}$ ' '(Nonlinear, obsv)'], 'Interpreter','latex','FontSize',16);

% K_list=[1,10,100,500];
% plot_nonlinear_AVR_obsv_K(K_list,A,B,C,G);

% pole_obsv_1=[-5 -2+1j -2-1j];
% pole_obsv_2=[-5 -2+8j -2-8j];
% pole_obsv_3=[-20 -2+1j -2-1j];
% pole_obsv_4=[-5 0+1j 0-1j];
% pole_list={pole_obsv_1,pole_obsv_2,pole_obsv_3,pole_obsv_4};
% plot_nonlinear_AVR_obsv_pole(pole_list,A,B,C,G);

%% Task 3: 含状态观测器的状态反馈控制器
% pole_ctrb=[-1.63 -0.71+0.25j -0.71-0.25j];
pole_ctrb=[-5 -0.44+0.4j -0.44-0.4j];
% pole_ctrb=[-0.5 -0.1+0.34j -0.1-0.34j];

K_ctrb=place(A_AVR,B_AVR,pole_ctrb);

%% 增加K状态反馈
%% Linear obsv+ctrb
% U_ref = workingpoint(4)+Vf/K;
% tspan=[0 50];
% y0_linear=zeros(6,1);
% [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv_ctrb(t,y,A_AVR,B_AVR,C_AVR,G,K_ctrb),tspan,y0_linear);
% % Vt
% figure;
% Vt_t=Vt_observer(y(:,1),y(:,3));
% plot(t_linear,Vt_t-Vt_t(1),'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$V_{t}$','FontSize',12,'Interpreter','Latex');
% % legend('$V_{t}$ origin','Interpreter','latex');
% title(['$V_{t}$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
% 
% % delta
% figure;
% plot(t_linear,y(:,1),'LineWidth',2);
% % hold on;
% % plot(t_linear,y(:,4),'LineWidth',2);
% % hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
% % legend('$\delta$ origin','$\delta$ obsv error','Interpreter','latex');
% title(['$\delta$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
% % omega
% figure;
% plot(t_linear,y(:,2),'LineWidth',2);
% % hold on;
% % plot(t_linear,y(:,5),'LineWidth',2);
% % hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
% % legend('$\omega$ origin','$\omega$ obsv error','Interpreter','latex');
% title(['$\omega$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
% % Eq'
% figure;
% plot(t_linear,y(:,3),'LineWidth',2);
% % hold on;
% % plot(t_linear,y(:,6),'LineWidth',2);
% % hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$E_{q}$','FontSize',12,'Interpreter','Latex');
% % legend('$E_{q}$ origin','$E_{q}$ obsv error','Interpreter','latex');
% title(['$E_{q}$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);



% pole_obsv=[-5 -2+1j -2-1j];
% G=place(A_AVR',C_AVR',pole_obsv)';
% 
% pole_obsv_ctrb_1=[-1000 -500+300j -500-300j];
% pole_obsv_ctrb_2=[-50 -50+5j -50-5j];
% pole_obsv_ctrb_3=[-10 -10+0.3j -10-0.3j];
% pole_obsv_ctrb_4=[-100 -100+300j -100-300j];
% 
% pole_list={pole_obsv_ctrb_1,pole_obsv_ctrb_2,pole_obsv_ctrb_3,pole_obsv_ctrb_4};
% plot_linear_AVR_obsv_ctrb_pole(pole_obsv,pole_list,A,B,C,G);

%% Noninear obsv+ctrb
% K=100;
% U_ref = workingpoint(4)+Vf/K;
% tspan=[0 100];
% A_obsv=A_AVR-G*C_AVR;
% y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
% [t_,y_]=ode45(@(t,y) gen_linear(A_AVR,B_AVR,t,y,dU_ref),tspan,[0 0 0]);
% [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_obsv_ctrb(t,y,A_obsv,B_AVR,G,K_ctrb),tspan,y0_nonlinear);
% % Vt
% figure;
% Vt_t=Vt_observer(y(:,1),y(:,3));
% plot(t_nonlinear,Vt_t-workingpoint(4),'LineWidth',2);
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$V_{t}$','FontSize',12,'Interpreter','Latex');
% % legend('$V_{t}$ origin','Interpreter','latex');
% title(['$V_{t}$ ' '(Nonlinear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
% 
% % delta
% figure;
% plot(t_nonlinear,y(:,1)-workingpoint(1),'LineWidth',2);
% % hold on;
% % plot(t_,y_(:,1),'LineWidth',2);
% % plot(t_nonlinear,y(:,4),'LineWidth',2);
% % hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
% % legend('$\delta$ origin','$\delta$ obsv','$\delta$ obsv error','Interpreter','latex');
% title(['$\delta$ ' '(Nonlinear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
% % omega
% figure;
% plot(t_nonlinear,y(:,2)-workingpoint(2),'LineWidth',2);
% % hold on;
% % plot(t_,y_(:,2),'LineWidth',2);
% % plot(t_nonlinear,y(:,5),'LineWidth',2);
% % hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
% % legend('$\omega$ origin','$\omega$ obsv','$\omega$ obsv error','Interpreter','latex');
% title(['$\omega$ ' '(Nonlinear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
% % Eq'
% figure;
% plot(t_nonlinear,y(:,3)-workingpoint(3),'LineWidth',2);
% % hold on;
% % plot(t_,y_(:,3),'LineWidth',2);
% % plot(t_nonlinear,y(:,6),'LineWidth',2);
% % hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$E_{q}$','FontSize',12,'Interpreter','Latex');
% % legend('$E_{q}$ origin','$E_{q}$ obsv', '$E_{q}$ obsv error','Interpreter','latex');
% title(['$E_{q}$ ' '(Nonlinear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);












