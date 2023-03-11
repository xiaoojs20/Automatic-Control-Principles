function dy=gen_nonlinear_AVR_obsv(t,y,A_obsv,B,G)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0 workingpoint;
global K dU_ref U_ref V_disturb K2;

Vt_t=Vt_observer(y(1),y(3));
y_obsv=[y(4);y(5);y(6)];
Vf_in = K.*(U_ref + dU_ref - Vt_t) + V_disturb;
dy1=y(2)-w0;
dy2=w0/2/H*Pm-D/2/H*(y(2)-w0)-w0/2/H*y(3)*Vs/(xd_p+xtl).*sin(y(1));
dy3=-1/Td_p*y(3)+1/Td0_p*(xd-xd_p)/(xd_p+xtl)*Vs*cos(y(1))+1/Td0_p*Vf_in;
dy_obsv=A_obsv*y_obsv+B*(dU_ref)+G*(Vt_t-workingpoint(4));

dy=[dy1;dy2;dy3;dy_obsv];
end

