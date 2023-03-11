function y=stability_performance_linear(sys,maxtime)
global K dU_ref U_ref;
[Vt_t,t_linear] = step(sys,maxtime);
figure;
WhetherToPlotPowerAngleCurve=0;
workingpoint=stablepoint(WhetherToPlotPowerAngleCurve);
plot(t_linear,dU_ref*Vt_t+workingpoint(4),'LineWidth',2);
grid on;
xlabel('time/s','FontSize',12,'Interpreter','Latex');
ylabel('$V_{t}/V$','FontSize',12,'Interpreter','Latex');
legend('$V_{t}$','Interpreter','latex','Location','southeast');
title(['Linear closed-loop system ' '$U_{ref}=$ ' num2str(U_ref) '$V$,' ...
    ' $\Delta U_{ref}=$' num2str(dU_ref) '$V, K=$ ' num2str(K)],...
    'FontSize',10,'Interpreter','latex');
RouthTable_linear_closedloop = Routh_table(sys.Denominator{1});
y=0;
end