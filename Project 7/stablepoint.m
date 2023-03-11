function y=stablepoint(WhetherToPlotPowerAngleCurve)
    global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p
    % 1.calculate steady-state delta
    Td_p=Td0_p*(xd_p+xtl)/(xd+xtl);
    w=w0;
    fun=@(delta)Td_p/Td0_p*(Vf+(xd-xd_p)/(xd_p+xtl)*Vs*cos(delta))*Vs/(xd_p+xtl).*sin(delta)-Pm;
    point=fsolve(fun,[0,pi/2]);

    % 2.plot power-angle curve
    if WhetherToPlotPowerAngleCurve==1
        delta_vector=0:0.001:pi/2;
        Pm_vector=Pm-Td_p/Td0_p*(Vf+(xd-xd_p)/(xd_p+xtl)*Vs*cos(delta_vector))*Vs/(xd_p+xtl).*sin(delta_vector);
        plot(delta_vector,Pm_vector,'LineWidth',2);
        hold on;grid on;
        xlabel('$\delta$','FontSize',16,'Interpreter','Latex');
        ylabel('$\Delta P$','FontSize',16,'Interpreter','Latex');
        plot([0,pi/2],[0,0],'--r','LineWidth',2);
    end

    % 3.judge the stability of the steady-state delta
    stable_delta=[];
    for i=1:length(point)
        delta=point(i);
        d_delta=1*10^(-6);
        dPm=Td_p/Td0_p*(Vf+(xd-xd_p)/(xd_p+xtl)*Vs*cos(delta+d_delta))*Vs/(xd_p+xtl).*sin(delta+d_delta)-Td_p/Td0_p*(Vf+(xd-xd_p)/(xd_p+xtl)*Vs*cos(delta))*Vs/(xd_p+xtl).*sin(delta);
        if dPm>0
            stable_delta=[stable_delta,delta];
        end
    end

    % 4.calculate other state variable
    delta=stable_delta;
    w=w0*ones(1,length(delta));
    Eq_p=Td_p/Td0_p*(Vf+(xd-xd_p)/(xd_p+xtl)*Vs*cos(delta));
    Vt=Vt_observer(delta,Eq_p);
    y=[delta,w,Eq_p,Vt];
    
end