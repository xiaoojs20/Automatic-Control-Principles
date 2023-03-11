function dy=gen_nonlinear(t,y,Vf_in)
global  xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0;
dy1=y(2)-w0;
dy2=w0/2/H*Pm-D/2/H*(y(2)-w0)-w0/2/H*y(3)*Vs/(xd_p+xtl).*sin(y(1));
dy3=-1/Td_p*y(3)+1/Td0_p*(xd-xd_p)/(xd_p+xtl)*Vs*cos(y(1))+1/Td0_p*Vf_in;
dy=[dy1;dy2;dy3;];
end