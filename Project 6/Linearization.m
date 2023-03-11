function y=Linearization(workingpoint)
global xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0
delta_0=workingpoint(1);Eq_p_0=workingpoint(3);Vt_0=workingpoint(4);

A=[0,1,0;
    -w0/2/H*Eq_p_0*Vs/(xd_p+xtl)*cos(delta_0),-D/2/H,-w0/2/H*Vs/(xd_p+xtl)*sin(delta_0);
    -1/Td0_p*(xd-xd_p)/(xd_p+xtl)*Vs*sin(delta_0),0,-1/Td_p];
B=[0;0;1/Td0_p];

lamda1=1/2/Vt_0*(2*(Eq_p_0-(Eq_p_0-Vs*cos(delta_0))...
            /(xd_p+xtl)*xd_p)*Vs*(-sin(delta_0))/(xd_p+xtl)*xd_p+2*xq*Vs/(xd+xtl)*...
            sin(delta_0)*xq*Vs/(xd+xtl)*cos(delta_0));
lamda2=1/2/Vt_0*2*(Eq_p_0-(Eq_p_0-Vs*cos(delta_0))/(xd_p+xtl)*xd_p)...
            *(1-xd_p/(xd_p+xtl));

y=[A,B,[lamda1;0;lamda2]];
end