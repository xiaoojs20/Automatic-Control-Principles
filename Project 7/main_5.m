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
sys_generator=(0.07552*s^2 + 0.00151*s + 0.009872)/(s^3 + 0.572*s^2 + 0.1879*s + 0.04846);
sys_generator_delta = tf(ss(A,B,C1,0));
sys_generator_omega = tf(ss(A,B,C2,0));
sys_generator_Eq_ = tf(ss(A,B,C3,0));
sys_generator_Pem = tf(ss(A,B,C_Pem,0));

sys_generator_num=cell2mat(sys_generator.num);
sys_generator_den=cell2mat(sys_generator.den);
sys_generator_omega_num=cell2mat(sys_generator_omega.num);
sys_generator_delta_num=cell2mat(sys_generator_delta.num);
sys_generator_Pem_num=cell2mat(sys_generator_Pem.num);

K=5;
K2=1;
sys_gain=K;
sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num)
sys_AVR_omega_zero=zero(sys_AVR_omega);
sys_AVR_omega_pole=pole(sys_AVR_omega);
sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
sys_AVR_delta_zero=zero(sys_AVR_delta);
sys_AVR_delta_pole=pole(sys_AVR_delta);
sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
sys_AVR_Pem_zero=zero(sys_AVR_Pem);
sys_AVR_Pem_pole=pole(sys_AVR_Pem);

%% Task
w_pos=logspace(-3,3,1000);
w_neg=-logspace(-3,3,1000);
w_neg(1,:)=w_neg(1,end:-1:1);
w=[w_neg,w_pos];

%% Task A2: frecscan nyquist
% abs_angle_freqscan(sys_AVR_omega,"omega",w);
% nyquist_freqscan(sys_AVR_omega,"omega",w);
% abs_angle_freqscan(sys_AVR_delta,"delta",w);
% nyquist_freqscan(sys_AVR_delta,"delta",w);
% abs_angle_freqscan(sys_AVR_Pem,"Pem",w);
% nyquist_freqscan(sys_AVR_Pem,"Pem",w);

%% Task A3: nyquist
% figure;
% nyquist(sys_AVR_omega,'b.-');
% title("The Nyquist diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
% figure;
% nyquist(sys_AVR_delta,'b.-');
% title("The Nyquist diagram of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
% figure;
% nyquist(sys_AVR_Pem,'b.-');
% title("The Nyquist diagram of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');

%% Task A4: different K --> nyquist
% K_list = [0.1 0.5 1 5 10 100 500 1000];
% figure;
% for K = K_list
%     sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);
%     subplot(2,4,find(K_list==K));
%     nyquist(sys_AVR_omega,'b.-');
%     title(['K = ', num2str(K)],'FontSize',12,'Interpreter','latex')
%     sgtitle("The Nyquist diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end
% figure;
% for K = K_list
%     sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
%     subplot(2,4,find(K_list==K));
%     nyquist(sys_AVR_delta,'b.-');
%     title(['K = ', num2str(K)],'FontSize',12,'Interpreter','latex')
%     sgtitle("The Nyquist diagram of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end
% figure;
% for K = K_list
%     sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
%     subplot(2,4,find(K_list==K));
%     nyquist(sys_AVR_Pem,'b.-');
%     title(['K = ', num2str(K)],'FontSize',12,'Interpreter','latex')
%     sgtitle("The Nyquist diagram of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end

%% Assignment: 4 task: a3 b1
tspan=[0,300];
dVf=0;
WhetherToPlotPowerAngleCurve=0;
workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
y0=workingpoint(1:3);
[t_nonlinear,y]=ode45(@(t,y) gen_nonlinear(t,y,Vf+dVf),tspan,y0);
Vt_t=Vt_observer(y(:,1),y(:,3));
x_nonlinear=[y,Vt_t];
K_list = [0.1 0.5 1 5 10 100 500 1000];
%% 没有加状态量反馈
% figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
%     subplot(2,4,find(K_list==K));
%     plot(t_nonlinear,y(:,2),'LineWidth',2);
%     ylim([0.95 1.05]);
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
% figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
%     subplot(2,4,find(K_list==K));
%     plot(t_nonlinear,y(:,1),'LineWidth',2);
%     ylim([0.5 0.75]);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
%     legend('$\delta$','Interpreter','latex','Location','southeast');
%     title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
%         ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $\delta$ '...
%         '(nonlinear model closed loop)'],'Interpreter','Latex');
% end
% figure;
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop(t,y),tspan,y0);
%     Pem=Pem_observer(y(:,1),y(:,3));
%     subplot(2,4,find(K_list==K));
%     plot(t_nonlinear,Pem,'LineWidth',2);
%     ylim([0.69 0.91]);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$P_{e}$','FontSize',12,'Interpreter','Latex');
%     legend('$P_{e}$','Interpreter','latex','Location','southeast');
%     title(['Nonlinear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
%         ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $P_{e}$ '...
%         '(nonlinear model closed loop)'],'Interpreter','Latex');
% end

%% 加了状态量反馈
% tspan=[0,200];
% K=5;
% dU_ref=0.05;
% K2=1;
% figure;
% K_list = [0.1 0.5 1 5];
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(1);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,"omega"),tspan,y0);
%     subplot(2,2,find(K_list==K));
%     plot(t_nonlinear,y(:,2),'LineWidth',2);
%     ylim([1 8.5]);
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
% figure;
% K_list = [0.1 0.5 1 5 10 100 500 1000];
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(1);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,"delta"),tspan,y0);
%     subplot(2,4,find(K_list==K));
%     plot(t_nonlinear,y(:,1),'LineWidth',2);
%     ylim([0, 0.8]);
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
% figure;
% K_list = [0.1 0.5 1 5 10 100 500 1000];
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(1);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,"Pem"),tspan,y0);
%     Pem=Pem_observer(y(:,1),y(:,3));
%     subplot(2,4,find(K_list==K));
%     plot(t_nonlinear,Pem,'LineWidth',2);
%     ylim([0.788, 0.81]);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$P_{em}$','FontSize',12,'Interpreter','Latex');
%     legend('$P_{em}$','Interpreter','latex','Location','southeast');
%     title(['$U_{ref}=$ ' num2str(U_ref) '$V$,' ' $\Delta U_{ref}=$' num2str(dU_ref) '$V$' ...
%         ' $K=$ ' num2str(K) ' $K_2=$ ' num2str(K2)],...
%         'FontSize',10,'Interpreter','latex');
%     sgtitle(['The effect of different K on' ' $P_{em}$ '...
%         '(nonlinear model closed loop double feedback)'],'Interpreter','Latex');
% end

%% Task B2: frecscan bode
% bode_freqscan(sys_AVR_omega,"omega",w);
% bode_freqscan(sys_AVR_delta,"delta",w);
% bode_freqscan(sys_AVR_Pem,"Pem",w);


%% Task B3: bode
% zpk(sys_AVR_omega);
% zpk(sys_AVR_delta);
% zpk(sys_AVR_Pem);
% figure;
% bode(sys_AVR_omega,'b.-');
% grid on;
% title("The Bode diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
% figure;
% bode(sys_AVR_delta,'b.-');
% grid on;
% title("The Bode diagram of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
% figure;
% bode(sys_AVR_Pem,'b.-');
% grid on;
% title("The Bode diagram of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');

%% Task B4: different K --> bode
K_list = [0.1 0.5 1 5 10 100 500 1000];
%% subplot 分开画
% figure;
% for K = K_list
%     sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);    
%     subplot(2,4,find(K_list==K));
%     bode(sys_AVR_omega,'b.-');
%     grid on;
%     title(['K = ', num2str(K)],'FontSize',12,'Interpreter','latex')
%     sgtitle("The Bode diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end
% figure;
% for K = K_list
%     sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
%     subplot(2,4,find(K_list==K));
%     bode(sys_AVR_delta,'b.-');
%     grid on;
%     title(['K = ', num2str(K)],'FontSize',12,'Interpreter','latex')
%     sgtitle("The Bode diagram of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end
% figure;
% for K = K_list
%     sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
%     subplot(2,4,find(K_list==K));
%     bode(sys_AVR_Pem,'b.-');
%     grid on;
%     title(['K = ', num2str(K)],'FontSize',12,'Interpreter','latex')
%     sgtitle("The Bode diagram of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');
% end

%% 合在一张图比较 Bode图 加开环
% s=tf('s');
% sys=20/(s*(s+1)*(0.1*s+1))
% sys=(-0.3089*s)/(s^3 + 1.327 *s^2 + 0.203 *s + 0.1472)
% [Gm,Pma,Wg,Wc] = margin(sys)

% K_list = [0.1 0.5 1 5 10];
% figure;
% for K = K_list
%     sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);
%     bode(sys_AVR_omega);
%     fprintf(['omega', 'K = ', num2str(K)]);
%     [G_margin,P_margin,Wg,Wc]=margin(sys_AVR_omega);
%     G_margin=20*log10(G_margin);
%     P_margin;
%     hold on;
% end
% bode(sys_generator_omega);
% bode(5*sys_generator_omega,1);
% hold off;
% grid on;
% % legend('K=0.1','K=0.5','K=1','K=5','K=10','K=500','K=1000','Interpreter','latex');
% legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');
% sgtitle("The Bode diagram of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
% 
% figure;
% for K = K_list
%     sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
%     bode(sys_AVR_delta);
%     fprintf(['delta', 'K = ', num2str(K)]);
%     [G_margin,P_margin,Wg,Wc]=margin(sys_AVR_delta);
%     G_margin=20*log10(G_margin);
%     P_margin;
%     hold on;
% end
% bode(sys_generator_delta);
% bode(5*sys_generator_delta);
% hold off;
% grid on;
% % legend('K=0.1','K=0.5','K=1','K=5','K=10','K=500','K=1000','Interpreter','latex');
% legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');
% sgtitle("The Bode diagram of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
% 
% figure;
% for K = K_list
%     sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
%     bode(sys_AVR_Pem);
%     fprintf(['Pem', 'K = ', num2str(K)]);
%     [G_margin,P_margin,Wg,Wc]=margin(sys_AVR_Pem);
%     G_margin=20*log10(G_margin);
%     P_margin;
%     hold on;
% end
% bode(sys_generator_Pem);
% bode(5*sys_generator_Pem);
% hold off;
% grid on;
% % legend('K=0.1','K=0.5','K=1','K=5','K=10','K=500','K=1000','Interpreter','latex');
% legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');
% sgtitle("The Bode diagram of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');

%% 合在一张图比较 非线性时域仿真 加开环
% tspan=[0,200];
% dU_ref=0.05;
% figure;
% K_list = [0.1 0.5 1 5];
% % K_list = [0.1 0.5 1];
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(2);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,"omega"),tspan,y0);
%     plot(t_nonlinear,y(:,2),'LineWidth',2);
%     hold on;
% end
% K_list = [1 5];
% for K = K_list
%     U_ref = Vf/K+K2.*workingpoint(2);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_state_feedback(t,y,"omega"),tspan,y0);
%     plot(t_nonlinear,y(:,2),'LineWidth',2);
%     hold on;
% end
% hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\omega/(rad/s)$','FontSize',12,'Interpreter','Latex');
% legend('K=0.1','K=0.5','K=1','K=5','open-loop K=1','open-loop K=5','Interpreter','latex');
% % legend('K=0.1','K=0.5','K=1','open-loop K=1','Interpreter','latex');
% sgtitle(['The effect of different K on' ' $\omega$ ', 'and open-loop'],'Interpreter','Latex');
% 
% figure;
% K_list = [0.1 0.5 1 5 10];
% for K = K_list
%     U_ref = workingpoint(4)+Vf/K+K2.*workingpoint(1);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,"delta"),tspan,y0);
%     plot(t_nonlinear,y(:,1),'LineWidth',2);
%     hold on;
% end
% K_list = [1 5];
% for K = K_list
%     U_ref = Vf/K+K2.*workingpoint(1);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_state_feedback(t,y,"delta"),tspan,y0);
%     plot(t_nonlinear,y(:,1),'LineWidth',2);
%     hold on;
% end
% hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$\delta/rad$','FontSize',12,'Interpreter','Latex');
% legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');
% sgtitle(['The effect of different K on' ' $\delta$ ', 'and open-loop'],'Interpreter','Latex');
% 
% figure;
% K_list = [0.1 0.5 1 5 10];
% for K = K_list
%     Pem_workingpoint=Pem_observer(workingpoint(1),workingpoint(3));
%     U_ref = workingpoint(4)+Vf/K+K2.*Pem_workingpoint;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_closedloop_double_feedback(t,y,"Pem"),tspan,y0);
%     Pem=Pem_observer(y(:,1),y(:,3));
%     plot(t_nonlinear,Pem,'LineWidth',2);
%     hold on;
% end
% K_list = [1 5];
% for K = K_list
%     U_ref = Vf/K+K2.*Pem_workingpoint;
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_state_feedback(t,y,"Pem"),tspan,y0);
%     Pem=Pem_observer(y(:,1),y(:,3));
%     plot(t_nonlinear,Pem,'LineWidth',2);
%     hold on;
% end
% hold off;
% grid on;
% xlabel('time/s','FontSize',12,'Interpreter','Latex');
% ylabel('$P_{em}$','FontSize',12,'Interpreter','Latex');
% legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');
% sgtitle(['The effect of different K on' ' $P_{em}$ ', 'and open-loop'],'Interpreter','Latex');


%% 分析闭环系统的频带宽度和谐振峰值 用AVR和反馈闭环的传递函数画幅频特性
K_list = [0.1 0.5 1 5];
figure;
for K = K_list
    sys_AVR_omega=tf(K*K2*sys_generator_omega_num,sys_generator_den+K*sys_generator_num);
    sys_AVR_PSS_omega=feedback(sys_AVR_omega,1);
    abs_angle_freqscan(sys_AVR_PSS_omega,"omega",w_pos);
    fprintf(['omega', 'K = ', num2str(K)]);
    hold on;
end
K_list = [1 5];
for K = K_list
    sys_PSS_omega=feedback(K*sys_generator_omega,1);
    abs_angle_freqscan(sys_PSS_omega,"omega",w_pos);
end
hold off;
grid on;
legend('K=0.1','K=0.5','K=1','K=5','open-loop K=1','open-loop K=5','Interpreter','latex');

K_list = [0.1 0.5 1 5 10];
figure;
for K = K_list
    sys_AVR_delta=tf(K*K2*sys_generator_delta_num,sys_generator_den+K*sys_generator_num);
    sys_AVR_PSS_delta=feedback(sys_AVR_delta,1);
    abs_angle_freqscan(sys_AVR_PSS_delta,"delta",w_pos);
    fprintf(['delta', 'K = ', num2str(K)]);
    hold on;
end
K_list = [1 5];
for K = K_list
    sys_PSS_delta=feedback(K*sys_generator_delta,1);
    abs_angle_freqscan(sys_PSS_delta,"delta",w_pos);
end
hold off;
grid on;
legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');

K_list = [0.1 0.5 1 5 10];
figure;
for K = K_list
    sys_AVR_Pem=tf(K*K2*sys_generator_Pem_num,sys_generator_den+K*sys_generator_num);
    sys_AVR_PSS_Pem=feedback(sys_AVR_Pem,1);
    abs_angle_freqscan(sys_AVR_PSS_Pem,"Pem",w_pos);
    fprintf(['Pem', 'K = ', num2str(K)]);
    hold on;
end
K_list = [1 5];
for K = K_list
    sys_PSS_Pem=feedback(K*sys_generator_Pem,1);
    abs_angle_freqscan(sys_PSS_Pem,"Pem",w_pos);
end
hold off;
grid on;
legend('K=0.1','K=0.5','K=1','K=5','K=10','open-loop K=1','open-loop K=5','Interpreter','latex');







