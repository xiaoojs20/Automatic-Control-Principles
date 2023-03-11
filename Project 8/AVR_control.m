% %% Initialization
% clear, clc;
% global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p Vt0
% parameter();
% %% Task A1: Stability analysis of open-loop systems
% % Assignment2 transfer function(Vf------>Vt)
% % sys4 =
% %  
% %    0.07552 s^2 + 0.00151 s + 0.009872
% %   ------------------------------------
% %   s^3 + 0.572 s^2 + 0.1879 s + 0.04846
% 
% s = tf('s');
% fprintf("~~~~~~~> transfer function <~~~~~~~\n");
% sys = (0.07552*s^2 + 0.00151*s + 0.009872)/(s^3 + 0.572*s^2 + 0.1879*s + 0.04846);
% sys
% % create Routh table
% Routhtable = Routh_table([1 0.572 0.1879 0.04846]);
% % Routhtable_test = Routh_table([1 2 2 4 1 1]);
% 
% %% Task A2: Stability performance analysis of open-loop systems
% [dVt,t] = step(sys,100);
% for i=1:3
%     subplot(1,3,i);
%     plot(t,i*0.05*dVt+Vt0(2),'LineWidth',2);
%     ylim([0.975 1.015]);
%     grid on;
%     xlabel('time/s','FontSize',12,'Interpreter','Latex');
%     ylabel('$V_t/V$','FontSize',12,'Interpreter','Latex');
%     title(['$\Delta V_{f} = $' num2str(i*0.05) 'V'], 'Interpreter','latex') 
%     sgtitle(['The effect of different' ' $\Delta V_{f}$ ' 'on' ' $V_{t}$ '...
%         '(transfer function)'],'Interpreter','Latex');
% end
% 
% %% Assignment 2 Linear Model Validation
% dVf=0.05;
% ylabel0={'$\delta$(p.u.)','$\omega$(p.u.)','$E^{''}_q$(p.u.)','$V_t$(p.u.)'};
% WhetherToPlotPowerAngleCurve=0;
% tspan=[0,100];
% figure;
% for Vf_factor = 1:3
%     workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
%     tspan=[0,100];y0=workingpoint(1:3);
%     [t_nonlinear,y]=ode45(@(t,y) gen_nonlinear(t,y,Vf+dVf*Vf_factor),tspan,y0);
%     Vt_t=Vt_observer(y(:,1),y(:,3));
%     x_nonlinear=[y,Vt_t];
%     subplot(1,3,Vf_factor);
%     plot(t_nonlinear,x_nonlinear(:,4),'LineWidth',2);
%     grid on;
%     ylim([0.975 1.015]);
%     xlabel('time/s','FontSize',16,'Interpreter','Latex');
%     ylabel('$V_{t}$','FontSize',16,'Interpreter','Latex');
%     sgtitle(['The effect of different' ' $\Delta V_{f}$ ' 'on' ' $V_{t}$ '...
%         '(nonlinear model)'],'Interpreter','Latex');
% end
% 
% %% Assignment 2 Nonlinear Model Validation
% figure;
% for Vf_factor = 1:3
%     workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
%     system_matrix=Linearization(workingpoint);
%     A=system_matrix(:,1:3);B=system_matrix(:,4);C=system_matrix(:,5);
%     [t_linear,y]=ode45(@(t,y) gen_linear(A,B,t,y,dVf*Vf_factor),tspan,[0;0;0]);
%     delta_t=y(:,1)+workingpoint(1);w_t=y(:,2)+workingpoint(2);
%     Eq_p_t=y(:,3)+workingpoint(3);Vt_t=y*C+workingpoint(4);
%     x_linear=[delta_t,w_t,Eq_p_t,Vt_t];
% 
%     subplot(1,3,Vf_factor);
%     plot(t_linear,x_linear(:,4),'LineWidth',2);
%     grid on;
%     ylim([0.975 1.015]);
%     xlabel('time/s','FontSize',16,'Interpreter','Latex');
%     ylabel('$V_{t}$','FontSize',16,'Interpreter','Latex');
%     sgtitle(['The effect of different' ' $\Delta V_{f}$ ' 'on' ' $V_{t}$ '...
%         '(linear model)'],'Interpreter','Latex');
% end
% 
% %% Task A3: Dynamic performance analysis of open-loop systems
% 
% 
% 
% 
% 
% 
% %% Closed-loop system modeling
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
