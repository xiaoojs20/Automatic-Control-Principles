function dy=plot_linear_AVR_obsv_ctrb_pole(pole_obsv,pole_list,A,B,C,G)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0 workingpoint;
global K dU_ref U_ref V_disturb K2;

tspan=[0 1];
y0_linear=zeros(6,1);
figure;
for i = 1:4
    pole_ctrb=pole_list{1,i};
    A_AVR=A-K*B*C';
    B_AVR=K*B;
    C_AVR=C';
    D_AVR=0;
    K_ctrb=place(A_AVR,B_AVR,pole_ctrb);
    U_ref = workingpoint(4)+Vf/K;
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv_ctrb(t,y,A_AVR,B_AVR,C_AVR,G,K_ctrb),tspan,y0_linear);
    % Vt
    subplot(2,2,i);
    Vt_t=Vt_observer(y(:,1),y(:,3));
    plot(t_linear,Vt_t,'LineWidth',2);
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$V_{t}$','FontSize',12,'Interpreter','Latex');
    legend('$V_{t}$ origin','Interpreter','latex');
    title(['pole ctrb= ' num2str(pole_ctrb)]);
    sgtitle(['$V_{t}$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
end

figure;
for i = 1:4
    pole_ctrb=pole_list{1,i};
    A_AVR=A-K*B*C';
    B_AVR=K*B;
    C_AVR=C';
    D_AVR=0;
    K_ctrb=place(A_AVR,B_AVR,pole_ctrb);
    U_ref = workingpoint(4)+Vf/K;
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv_ctrb(t,y,A_AVR,B_AVR,C_AVR,G,K_ctrb),tspan,y0_linear);
    % delta
    subplot(2,2,i);
    plot(t_linear,y(:,1),'LineWidth',2);
    hold on;
    plot(t_linear,y(:,4),'LineWidth',2);
    hold off;
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
    legend('$\delta$ origin','$\delta$ obsv error','Interpreter','latex');
    title(['pole ctrb= ' num2str(pole_ctrb)]);
    sgtitle(['$\delta$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
end
figure;
for i = 1:4
    pole_ctrb=pole_list{1,i};
    A_AVR=A-K*B*C';
    B_AVR=K*B;
    C_AVR=C';
    D_AVR=0;
    K_ctrb=place(A_AVR,B_AVR,pole_ctrb);
    U_ref = workingpoint(4)+Vf/K;
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv_ctrb(t,y,A_AVR,B_AVR,C_AVR,G,K_ctrb),tspan,y0_linear);
    % omega
    subplot(2,2,i);
    plot(t_linear,y(:,2),'LineWidth',2);
    hold on;
    plot(t_linear,y(:,5),'LineWidth',2);
    hold off;
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
    legend('$\omega$ origin','$\omega$ obsv error','Interpreter','latex');
    title(['pole ctrb= ' num2str(pole_ctrb)]);
    sgtitle(['$\omega$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
end
figure;
for i = 1:4
    pole_ctrb=pole_list{1,i};
    A_AVR=A-K*B*C';B_AVR=K*B;C_AVR=C';D_AVR=0;
    K_ctrb=place(A_AVR,B_AVR,pole_ctrb);
    U_ref = workingpoint(4)+Vf/K;
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv_ctrb(t,y,A_AVR,B_AVR,C_AVR,G,K_ctrb),tspan,y0_linear);
    % Eq'
    subplot(2,2,i);
    plot(t_linear,y(:,3),'LineWidth',2);
    hold on;
    plot(t_linear,y(:,6),'LineWidth',2);
    hold off;
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$E_{q}$','FontSize',12,'Interpreter','Latex');
    legend('$E_{q}$ origin','$E_{q}$ obsv error','Interpreter','latex');
    title(['pole ctrb= ' num2str(pole_ctrb)]);
    sgtitle(['$E_{q}$ ' '(Linear, obsv ctrb)'], 'Interpreter','latex','FontSize',16);
end
   
end

