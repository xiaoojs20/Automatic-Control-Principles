function dy=gen_nonlinear_AVR_PSS_corr(Gc,t,y,mode)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
global K dU_ref U_ref V_disturb K2;
% mode=1,delta,Pem
Vt_t=Vt_observer(y(1),y(3));
[A_corr,B_corr,C_corr,D_corr] = tf2ss(Gc.num{1}, Gc.den{1});
Gc_x_p=A_corr*[y(4);y(5)]+B_corr*(U_ref+dU_ref-y(2));
Gc_y=C_corr*[y(4);y(5)]+D_corr*(U_ref+dU_ref-y(2));
% y(4),y(5)只是用于解方程代表x'=Ax+Bu中的x,(U_ref+dU_ref-y(2))代表着输入量u
if mode=="delta"
    Vf_in = K.*(Gc_y-Vt_t) + V_disturb;
elseif mode=="omega"
    Vf_in = K.*(Gc_y-Vt_t) + V_disturb;
elseif mode=="Pem"
    Vf_in = K.*(Gc_y-Vt_t) + V_disturb;
end

dstate=gen_nonlinear(t,[y(1);y(2);y(3)],Vf_in);
dy=[dstate;Gc_x_p;];

end
