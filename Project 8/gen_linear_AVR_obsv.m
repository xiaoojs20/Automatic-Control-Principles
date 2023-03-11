function dy=gen_linear_AVR_obsv(t,y,A,B,C,G)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
global K dU_ref U_ref V_disturb K2;

Vt_t=Vt_observer(y(1),y(3));
u=U_ref;

y_avr=[y(1);y(2);y(3)];
y_obsv=[y(4);y(5);y(6)];
dy_avr=A*y_avr+B*u;
dVt=C*y_avr;
y=dVt;
A_obsv=A-G*C;
dy_obsv=A_obsv*y_obsv+B*u+G*y;
dy=[dy_avr;dy_obsv];
end

