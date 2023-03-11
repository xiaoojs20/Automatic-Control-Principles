function y=Pem_observer(delta,Eq_p,lamda1,lamda2)
    global xd xq xd_p H D Td_p Td0_p xtl Vf Vs Pm w0 m xdSigma_p K2
    
    if nargin==2
        y=m.*(Eq_p*Vs).*sin(delta)./xdSigma_p;      
    elseif nargin==4
        y=lamda1*delta+lamda2*Eq_p;
    end
    
end