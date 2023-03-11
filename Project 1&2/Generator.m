%% Modeling parameters
clear;
clc;
xd = 1.998;
xq = 1.998;
xd2 = 0.311;
H = 2.5;
D = 0.1;
Td02 = 6.11;
xtl = 0.4;
Vf = 1.1;
Vs = 1.0;
Pm = 0.8;
omega0 = 1.0;
xdSigma = xd + xtl;
xdSigma2 = xd2 + xtl;
Td2 = Td02*(xd2 + xtl)/(xd + xtl);
format long
%% Variable is delta
% % 这里第一次使用了单变量方法，精度差
% % syms delta;
% % omega = omega0;
% % Eq2 = (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(delta));
% % Eq2 = vpa(Eq2);
% % omegad= @(delta) (omega0/(2*H))*Pm - ((omega0/(2*H))*((((Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(delta)))*Vs)/xdSigma2))*sin(delta);
% % x1 = fsolve(omegad, 0.7)
% % x2 = fsolve(omegad, 1.2)
% 
% 
syms delta;
omega = omega0;
F = @(x) [x(2) - (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(x(1)));...
          (omega0/(2*H))*Pm - (omega0/(2*H))*((x(2)*Vs)/xdSigma2)*sin(x(1))];

format long
x1 = vpa(fsolve(F, [0; 0]));
x2 = vpa(fsolve(F, [2; 0]));

deltasolve = [x1(1) x2(1)];
Eq2solve = [x1(2) x2(2)];

%% 给功角加一个小扰动
syms Eq2;
disturb =0.001; %0.0573degree
deltad = omega - omega0;
Eq2dist = (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(delta));
omegad = (omega0/(2*H))*Pm - (D/(2*H))*(omega-omega0) - (omega0/(2*H))*((Eq2dist*Vs)/xdSigma2)*sin(delta);
Eq2d = -Eq2/Td2 + (1/Td02)*((xd-xd2)/xdSigma2)*Vs*cos(delta) + Vf/Td02;
% 为什么这里disturb需要进行微调
% 因为上面的解其实不能保证完全精确，所以调整至disturb增加，Eq2d减小的程度，才能保证这个扰动分别让其左偏和右偏

% Eq2dbalancetest = subs(Eq2d,delta,deltasolve(1)); % shoule be zero
Eq2dPosDisturb1 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(1)+disturb,Eq2solve(1)}));
omegadPosDisturb1 = vpa(subs(omegad,{delta,Eq2},{deltasolve(1)+disturb,Eq2solve(1)}));
Eq2dNegDisturb1 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(1)-disturb,Eq2solve(1)}));
omegadNegDisturb1 = vpa(subs(omegad,{delta,Eq2},{deltasolve(1)-disturb,Eq2solve(1)}));

Eq2dPosDisturb2 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(2)+disturb,Eq2solve(2)}));
omegadPosDisturb2 = vpa(subs(omegad,{delta,Eq2},{deltasolve(2)+disturb,Eq2solve(2)}));
Eq2dNegDisturb2 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(2)-disturb,Eq2solve(2)}));
omegadNegDisturb2 = vpa(subs(omegad,{delta,Eq2},{deltasolve(2)-disturb,Eq2solve(2)}));

Eq2dDisturb1 = [Eq2dPosDisturb1 Eq2dNegDisturb1];
omegadDisturb1 = [omegadPosDisturb1 omegadNegDisturb1];
Eq2dDisturb2 = [Eq2dPosDisturb2 Eq2dNegDisturb2];
omegadDisturb2 = [omegadPosDisturb2 omegadNegDisturb2];
Eq2dDisturb = [Eq2dDisturb1; Eq2dDisturb2];
omegadDisturb = [omegadDisturb1; omegadDisturb2];


if(omegadDisturb1(1)<0 && omegadDisturb1(2)>0)
    stable(1) = 1;
else 
    stable(1) = 0;
end
if(omegadDisturb2(1)<0 && omegadDisturb2(2)>0)
    stable(2) = 1;
else 
    stable(2) = 0;
end
%% 线性化矩阵，李雅普诺夫间接法（something wrong, ignore）
% delta1 = 0.735373558459363;
% delta2 = 1.122194080447383;
% % A对吗?
% A = [0 1 0;...
%     (-omega0/(2*H))*(Vs/xdSigma2).*(Eq2.*cos(delta)+sin(delta)*-(1/Td02)*((xd-xd2)/xdSigma2)*Vs.*sin(delta)) -D/(2*H) 0 ; ...
%     -(1/Td02)*((xd-xd2)/xdSigma2)*Vs.*sin(delta) 0 -1/Td2];
% 
% Anew1 = subs(A,delta,0.735373558459363);
% e1 = eig(Anew1);
% Anew2 = subs(A,delta,1.122194080447383);
% e2 = eig(Anew2);

%% 构造在稳态运行点处线性化状态函数方程
% 平衡点参数
deltabalance = deltasolve(1);
omegabalance = 1.0;
Eq2balance = Eq2solve(1);

% 构造线性化方程
syms delta omega Eq2 dVf ;
syms ddelta domega dEq2; % 代表线性化的变量
id = (Eq2-Vs*cos(delta))/xdSigma2;
Vtq = Eq2-id*xd2;
iq = Vs*sin(delta)/xdSigma;
Vtd = iq*xq;
Vt = sqrt(Vtd^2 + Vtq^2);
dVt = diff(Vt,delta)*ddelta + diff(Vt,omega)*domega + diff(Vt,Eq2)*dEq2;
dVt = vpa(subs(dVt,[delta omega Eq2],[deltabalance omegabalance Eq2balance]));

deltad = omega - omega0;
omegad = (omega0/(2*H))*Pm - D*(omega-omega0)/(2*H) - (omega0/(2*H))*((Eq2*Vs)/xdSigma2)*sin(delta);
Eq2d = -Eq2/Td2 + (1/Td02)*((xd-xd2)/xdSigma2)*Vs*cos(delta) + Vf/Td02;

ddeltad = diff(deltad,delta)*ddelta + diff(deltad,omega)*domega + diff(deltad,Eq2)*dEq2;
domegad = diff(omegad,delta)*ddelta + diff(omegad,omega)*domega + diff(omegad,Eq2)*dEq2;
dEq2d = diff(Eq2d,delta)*ddelta + diff(Eq2d,omega)*domega + diff(Eq2d,Eq2)*dEq2;

ddeltad = vpa(subs(ddeltad,[delta omega Eq2],[deltabalance omegabalance Eq2balance]));
domegad = vpa(subs(domegad,[delta omega Eq2],[deltabalance omegabalance Eq2balance]));
dEq2d = vpa(subs(dEq2d + (1/Td02)*dVf,[delta omega Eq2],[deltabalance omegabalance Eq2balance]));

dd = [ddeltad domegad dEq2d]';
d = [ddelta domega dEq2 dVf];
A = double([diff(dd,d(1)) diff(dd,d(2)) diff(dd,d(3))]);
B = double([diff(dd,d(4))]);
C = double([diff(dVt,d(1)) diff(dVt,d(2)) diff(dVt,d(3))]);



writematrix(double(A),"linear_A.csv");
writematrix(double(B),"linear_B.csv");
writematrix(double(C),"linear_C.csv");
