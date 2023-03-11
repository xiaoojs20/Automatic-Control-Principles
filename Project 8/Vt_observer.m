function y=Vt_observer(delta,Eq_p,lamda1,lamda2)
    global xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0
    
    if nargin==2
        id=(Eq_p-Vs*cos(delta))/(xd_p+xtl);
        Vtq=Eq_p-id*xd_p;
        iq=Vs/(xd+xtl)*sin(delta);
        Vtd=iq*xq;
        y=sqrt(Vtd.^2+Vtq.^2);
        
    elseif nargin==4
        y=lamda1*delta+lamda2*Eq_p;
    end
    
end