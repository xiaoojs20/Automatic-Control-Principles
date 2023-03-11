function dy=gen_nonlinear_AVR_PSS(t,y,mode)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
global K dU_ref U_ref V_disturb K2;
% mode=1,delta,Pem

Vt_t=Vt_observer(y(1),y(3));
Pem=Pem_observer(y(1),y(3));
if mode=="delta"
    Vf_in = K.*(U_ref + dU_ref - K2.*y(1) - Vt_t) + V_disturb;
elseif mode=="omega"
    Vf_in = K.*(U_ref + dU_ref - K2.*y(2) - Vt_t) + V_disturb;
elseif mode=="Pem"
    Vf_in = K.*(U_ref + dU_ref - K2.*Pem - Vt_t) + V_disturb;
end
dy=gen_nonlinear(t,y,Vf_in);
end

