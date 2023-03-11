function dy=gen_nonlinear_AVR(t,y)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
global K dU_ref U_ref V_disturb;

Vt_t=Vt_observer(y(1),y(3));
Pem=Pem_observer(y(1),y(3));
Vf_in = K.*(U_ref + dU_ref - Vt_t) + V_disturb;
dy1=y(2)-w0;
dy2=w0/2/H*Pm-D/2/H*(y(2)-w0)-w0/2/H*y(3)*Vs/(xd_p+xtl).*sin(y(1));
dy3=-1/Td_p*y(3)+1/Td0_p*(xd-xd_p)/(xd_p+xtl)*Vs*cos(y(1))+1/Td0_p*Vf_in;
dy=[dy1;dy2;dy3;];
end

