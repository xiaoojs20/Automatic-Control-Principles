function parameter()
    global xd xq xd_p H D Td0_p xtl Vf Vs Pm w0 Td_p m Vt0 xdSigma xdSigma_p
    xd=1.998;
    xq=1.998;
    xd_p=0.311;
    H=2.5;
    D=0.1;
    Td0_p=6.11;
    xtl=0.4;
    Vf=1.1;
    Vs=1.0;
    Pm=0.8;
    w0=1;
    xdSigma = xd + xtl;
    xdSigma_p = xd_p + xtl;
    Td_p=Td0_p*(xd_p+xtl)/(xd+xtl);
    m=1;
    Vt0 = [0.966440716284484 0.977046711563641 0.986973031776744];
end