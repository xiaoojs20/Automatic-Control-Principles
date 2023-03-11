function dy=gen_nonlinear_state_feedback(t,y,mode)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
global K dU_ref U_ref V_disturb K2;
% mode=1,delta,Pem

Vt_t=Vt_observer(y(1),y(3));
Pem=Pem_observer(y(1),y(3));
if mode=="delta"
    Vf_in = K.*(U_ref + dU_ref - K2.*y(1)) + V_disturb;
elseif mode=="omega"
    Vf_in = K.*(U_ref + dU_ref - K2.*y(2)) + V_disturb;
elseif mode=="Pem"
    Vf_in = K.*(U_ref + dU_ref - K2.*Pem) + V_disturb;
end
dy1=y(2)-w0;
dy2=w0/2/H*Pm-D/2/H*(y(2)-w0)-w0/2/H*y(3)*Vs/(xd_p+xtl).*sin(y(1));
dy3=-1/Td_p*y(3)+1/Td0_p*(xd-xd_p)/(xd_p+xtl)*Vs*cos(y(1))+1/Td0_p*Vf_in;
dy=[dy1;dy2;dy3;];
end

