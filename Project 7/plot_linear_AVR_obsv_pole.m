function dy=plot_linear_AVR_obsv_pole(pole_list,A,B,C,G)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0 workingpoint;
global K dU_ref U_ref V_disturb K2;

figure;
for i = 1:4
    pole_obsv=pole_list{1,i};
    A_AVR=A-K*B*C';
    B_AVR=K*B;
    C_AVR=C';
    D_AVR=0;
    G=place(A_AVR',C_AVR',pole_obsv)';
    U_ref = workingpoint(4)+Vf/K;
    tspan=[0 500];
    y0_linear=zeros(6,1);
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv(t,y,A_AVR,B_AVR,C_AVR,G),tspan,y0_linear);
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
    title(['pole obsv= ' num2str(pole_obsv)]);
    sgtitle(['$\delta$ ' '(Linear)'], 'Interpreter','latex','FontSize',16);
end
figure;
for i = 1:4
    pole_obsv=pole_list{1,i};
    A_AVR=A-K*B*C';
    B_AVR=K*B;
    C_AVR=C';
    D_AVR=0;
    G=place(A_AVR',C_AVR',pole_obsv)';
    U_ref = workingpoint(4)+Vf/K;
    tspan=[0 500];
    y0_linear=zeros(6,1);
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv(t,y,A_AVR,B_AVR,C_AVR,G),tspan,y0_linear);
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
    title(['pole obsv= ' num2str(pole_obsv)]);
    sgtitle(['$\omega$ ' '(Linear)'], 'Interpreter','latex','FontSize',16);
end
figure;
for i = 1:4
    pole_obsv=pole_list{1,i};
    A_AVR=A-K*B*C';
    B_AVR=K*B;
    C_AVR=C';
    D_AVR=0;
    G=place(A_AVR',C_AVR',pole_obsv)';
    U_ref = workingpoint(4)+Vf/K;
    tspan=[0 500];
    y0_linear=zeros(6,1);
    [t_linear,y]=ode45(@(t,y) gen_linear_AVR_obsv(t,y,A_AVR,B_AVR,C_AVR,G),tspan,y0_linear);
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
    title(['pole obsv= ' num2str(pole_obsv)]);
    sgtitle(['$E_{q}$ ' '(Linear)'], 'Interpreter','latex','FontSize',16);
end
   
end

