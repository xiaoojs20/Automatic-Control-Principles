function dy=gen_nonlinear_AVR(t,y)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
global K dU_ref U_ref V_disturb;

Vt_t=Vt_observer(y(1),y(3));
Pem=Pem_observer(y(1),y(3));
Vf_in = K.*(U_ref + dU_ref - Vt_t) + V_disturb;
dy=gen_nonlinear(t,y,Vf_in);
end

