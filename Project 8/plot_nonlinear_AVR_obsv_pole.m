function dy=plot_nonlinear_AVR_obsv_pole(pole_list,A,B,C,G)
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
    A_obsv=A_AVR-G*C_AVR;
    y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
    [t_,y_]=ode45(@(t,y) gen_linear(A_AVR,B_AVR,t,y,dU_ref),tspan,[0 0 0]);
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_obsv(t,y,A_obsv,B_AVR,G),tspan,y0_nonlinear);
    % delta
    subplot(2,2,i);
    plot(t_nonlinear,y(:,1)-workingpoint(1),'LineWidth',2);
    hold on;
    plot(t_,y_(:,1),'LineWidth',2);
    plot(t_nonlinear,y(:,4),'LineWidth',2);
    hold off;
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$\delta$','FontSize',12,'Interpreter','Latex');
    legend('$\delta$ origin','$\delta$ obsv','$\delta$ obsv error','Interpreter','latex');
    title(['pole obsv= ' num2str(pole_obsv)]);
    sgtitle(['$\delta$ ' '(Nonlinear)'], 'Interpreter','latex','FontSize',16);
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
    A_obsv=A_AVR-G*C_AVR;
    y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
    [t_,y_]=ode45(@(t,y) gen_linear(A_AVR,B_AVR,t,y,dU_ref),tspan,[0 0 0]);
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_obsv(t,y,A_obsv,B_AVR,G),tspan,y0_nonlinear);
     % omega
    subplot(2,2,i);
    plot(t_nonlinear,y(:,2)-workingpoint(2),'LineWidth',2);
    hold on;
    plot(t_,y_(:,2),'LineWidth',2);
    plot(t_nonlinear,y(:,5),'LineWidth',2);
    hold off;
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$\omega$','FontSize',12,'Interpreter','Latex');
    legend('$\omega$ origin','$\omega$ obsv','$\omega$ obsv error','Interpreter','latex');
    title(['pole obsv= ' num2str(pole_obsv)]);
    sgtitle(['$\omega$ ' '(Nonlinear)'], 'Interpreter','latex','FontSize',16);
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
    A_obsv=A_AVR-G*C_AVR;
    y0_nonlinear=[workingpoint(1:3) zeros(1,3)];
    [t_,y_]=ode45(@(t,y) gen_linear(A_AVR,B_AVR,t,y,dU_ref),tspan,[0 0 0]);
    [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear_AVR_obsv(t,y,A_obsv,B_AVR,G),tspan,y0_nonlinear);
    % Eq'
    subplot(2,2,i);
    plot(t_nonlinear,y(:,3)-workingpoint(3),'LineWidth',2);
    hold on;
    plot(t_,y_(:,3),'LineWidth',2);
    plot(t_nonlinear,y(:,6),'LineWidth',2);
    hold off;
    grid on;
    xlabel('time/s','FontSize',12,'Interpreter','Latex');
    ylabel('$E_{q}$','FontSize',12,'Interpreter','Latex');
    legend('$E_{q}$ origin','$E_{q}$ obsv', '$E_{q}$ obsv error','Interpreter','latex');
    title(['pole obsv= ' num2str(pole_obsv)]);
    sgtitle(['$E_{q}$ ' '(Nonlinear)'], 'Interpreter','latex','FontSize',16);
end