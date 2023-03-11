%% Modeling parameters
clear;
clc;
xd = 1.998;
% xd = 2.5;
xq = xd;
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
DeltaVf = 0.05;
% DeltaVf = 1;
format long

%% 计算不同稳态初值
syms delta;
%% 1.1+0.05（任务一）
omega = omega0;
F = @(x) [x(2) - (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(x(1)));...
          (omega0/(2*H))*Pm - (omega0/(2*H))*((x(2)*Vs)/xdSigma2)*sin(x(1))];

format long
x1 = vpa(fsolve(F, [0; 0]));
x2 = vpa(fsolve(F, [2; 0]));

deltaSolve = [x1(1) x2(1)];
Eq2Solve = [x1(2) x2(2)];

%% 1.15+0.05 Upper Vf
Vf = 1.15;
omega = omega0;
F = @(x) [x(2) - (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(x(1)));...
          (omega0/(2*H))*Pm - (omega0/(2*H))*((x(2)*Vs)/xdSigma2)*sin(x(1))];

format long
x1 = vpa(fsolve(F, [0; 0]));
x2 = vpa(fsolve(F, [2; 0]));

deltaUpperSolve = [x1(1) x2(1)];
Eq2UpperSolve = [x1(2) x2(2)];

%% 1.05+0.05 Lower Vf
Vf = 1.05;
omega = omega0;
F = @(x) [x(2) - (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(x(1)));...
          (omega0/(2*H))*Pm - (omega0/(2*H))*((x(2)*Vs)/xdSigma2)*sin(x(1))];

format long
x1 = vpa(fsolve(F, [0; 0]));
x2 = vpa(fsolve(F, [2; 0]));

deltaLowerSolve = [x1(1) x2(1)];
Eq2LowerSolve = [x1(2) x2(2)];

%% 检验另外两个稳态初值的稳定性
% 运行这段代码的时候记得把上面的deltaLowerSolve、Eq2LowerSolve改成deltasolve、Eq2solve
% syms Eq2;
% disturb =0.001; %0.0573degree
% deltad = omega - omega0;
% Eq2dist = (Td2/Td02)*(Vf + ((xd-xd2)/xdSigma2)*Vs*cos(delta));
% omegad = (omega0/(2*H))*Pm - (D/(2*H))*(omega-omega0) - (omega0/(2*H))*((Eq2dist*Vs)/xdSigma2)*sin(delta);
% Eq2d = -Eq2/Td2 + (1/Td02)*((xd-xd2)/xdSigma2)*Vs*cos(delta) + Vf/Td02;
% vpa(omegad);
% 
% Eq2dPosDisturb1 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(1)+disturb,Eq2solve(1)}));
% omegadPosDisturb1 = vpa(subs(omegad,{delta,Eq2},{deltasolve(1)+disturb,Eq2solve(1)}));
% Eq2dNegDisturb1 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(1)-disturb,Eq2solve(1)}));
% omegadNegDisturb1 = vpa(subs(omegad,{delta,Eq2},{deltasolve(1)-disturb,Eq2solve(1)}));
% 
% Eq2dPosDisturb2 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(2)+disturb,Eq2solve(2)}));
% omegadPosDisturb2 = vpa(subs(omegad,{delta,Eq2},{deltasolve(2)+disturb,Eq2solve(2)}));
% Eq2dNegDisturb2 = vpa(subs(Eq2d,{delta,Eq2},{deltasolve(2)-disturb,Eq2solve(2)}));
% omegadNegDisturb2 = vpa(subs(omegad,{delta,Eq2},{deltasolve(2)-disturb,Eq2solve(2)}));
% 
% Eq2dDisturb1 = [Eq2dPosDisturb1 Eq2dNegDisturb1];
% omegadDisturb1 = [omegadPosDisturb1 omegadNegDisturb1];
% Eq2dDisturb2 = [Eq2dPosDisturb2 Eq2dNegDisturb2];
% omegadDisturb2 = [omegadPosDisturb2 omegadNegDisturb2];
% Eq2dDisturb = [Eq2dDisturb1; Eq2dDisturb2]
% omegadDisturb = [omegadDisturb1; omegadDisturb2]

%% 三阶状态方程,非线性动态过程
% syms delta omega Eq2 dVf ;
% deltad = omega - omega0;
% omegad = (omega0/(2*H))*Pm - D*(omega-omega0)/(2*H) - (omega0/(2*H))*((Eq2*Vs)/xdSigma2)*sin(delta);
% Eq2d = -Eq2/Td2 + (1/Td02)*((xd-xd2)/xdSigma2)*Vs*cos(delta) + Vf/Td02;

% 123列分别为Vf,upperVf,lowerVf的初值
% 用double函数否则"错误使用 odearguments输入必须为单精度或双精度浮点值"

x0 = double([[deltaLowerSolve(1);omega0;Eq2LowerSolve(1)]...
    [deltaSolve(1);omega0;Eq2Solve(1)]...
    [deltaUpperSolve(1);omega0;Eq2UpperSolve(1)]]);
Vfs = [1.05;1.1; 1.15];
tspan = [0 100];
for i=1:3 
stateFunction = @(t,x)[x(2)-omega0;...
    (omega0/(2*H))*Pm - D*(x(2)-omega0)/(2*H) - (omega0/(2*H))*((x(3)*Vs)/xdSigma2)*sin(x(1));...
    -x(3)/Td2 + (1/Td02)*((xd-xd2)/xdSigma2)*Vs*cos(x(1)) + (Vfs(i)+DeltaVf)/Td02];
    [t,y] = ode45(stateFunction,tspan,x0(:,i));
    
    delta = y(:,1);
    Eq2 = y(:,3);
    id = (Eq2-Vs.*cos(delta))./xdSigma2;
    Vtq = Eq2-id.*xd2;
    iq = Vs.*sin(delta)/xdSigma;
    Vtd = iq.*xq;
    Vt = sqrt(Vtd.^2 + Vtq.^2);
    y(:,4) = Vt;

    subplot(3,3,i);
    plot(t,y);
    titles = ["task1 Vf=1.05V" "task1 Vf=1.1V" "task1 Vf=1.15V"];
    title(titles(i)); 
    legend('\delta','\omega','E_{q}’','V_{t}');
end

%% 线性模型动态过程
syms delta omega Eq2 dVf ;
syms ddelta domega dEq2; % 代表线性化的变量
x = [delta; omega; Eq2];
x0 = double([[deltaLowerSolve(1);omega0;Eq2LowerSolve(1)]...
    [deltaSolve(1);omega0;Eq2Solve(1)]...
    [deltaUpperSolve(1);omega0;Eq2UpperSolve(1)]]);

id = (Eq2-Vs.*cos(delta))./xdSigma2;
Vtq = Eq2-id.*xd2;
iq = Vs.*sin(delta)/xdSigma;
Vtd = iq.*xq;
Vt = sqrt(Vtd.^2 + Vtq.^2);

deltad = omega - omega0;
omegad = (omega0./(2.*H)).*Pm - D.*(omega-omega0)/(2.*H) - (omega0/(2.*H)).*((Eq2.*Vs)/xdSigma2).*sin(delta);
Eq2d = -Eq2./Td2 + (1./Td02).*((xd-xd2)./xdSigma2).*Vs.*cos(delta) + Vf./Td02;

dVt = diff(Vt,delta).*ddelta + diff(Vt,omega).*domega + diff(Vt,Eq2).*dEq2;
ddeltad = diff(deltad,delta).*ddelta + diff(deltad,omega).*domega + diff(deltad,Eq2).*dEq2;
domegad = diff(omegad,delta).*ddelta + diff(omegad,omega).*domega + diff(omegad,Eq2).*dEq2;
dEq2d = diff(Eq2d,delta).*ddelta + diff(Eq2d,omega).*domega + diff(Eq2d,Eq2).*dEq2;

C1 = [1 0 0];
C2 = [0 1 0];
C3 = [0 0 1];

for i=1:3
    dVts = vpa(subs(dVt,x,x0(:,i)));
    ddeltads = vpa(subs(ddeltad,x,x0(:,i)));
    domegads = vpa(subs(domegad,x,x0(:,i)));
    dEq2ds = vpa(subs(dEq2d + (1/Td02).*dVf,x,x0(:,i)));
    Vt0(i) = double(subs(Vt,[delta,Eq2], [x0(1,i) x0(3,i)]));
        
    dd = [ddeltads domegads dEq2ds]';
    d = [ddelta domega dEq2 dVf];
    A = double([diff(dd,d(1)) diff(dd,d(2)) diff(dd,d(3))]);
    B = double([diff(dd,d(4))]);
    C = double([diff(dVts,d(1)) diff(dVts,d(2)) diff(dVts,d(3))]);
   
    linearFuntion = @(t,linear_x)(A*linear_x+B*DeltaVf);
    [t,y] = ode45(linearFuntion,tspan,[0;0;0]);
    y(:,4) = C(1)*y(:,1)+C(2)*y(:,2)+C(3)*y(:,3);
    for j = 1:3
    y(:,j) = y(:,j) + x0(j,i);
    end
    y(:,4) = y(:,4) + Vt0(i);

    subplot(3,3,i+3);
    plot(t,y);
    legend('\delta','\omega','E_{q}’','V_{t}');
    titles = ["task2 Vf=1.05V" "task2 Vf=1.1V" "task2 Vf=1.15V"];
    title(titles(i)); 

%% 系统传递函数
    ltiSys1 = ss(A,B,C1,0);
    ltiSys2 = ss(A,B,C2,0);
    ltiSys3 = ss(A,B,C3,0);
    ltiSys4 = ss(A,B,C,0);

    sys1 = tf (ltiSys1)
    sys2 = tf (ltiSys2)
    sys3 = tf (ltiSys3)
    sys4 = tf (ltiSys4)

    [ddeltat,t1] = step(sys1,100);
    [domegat,t2] = step(sys2,100);
    [dEq2t,t3] = step(sys3,100);
    [dVtt,t4] = step(sys4,100);
    
    subplot(3,3,i+6);
    % 0.05 step + 稳态初值
    if DeltaVf == 0.05
        plot(t1,0.05*ddeltat+x0(1,i)); 
        hold on
        plot(t2,0.05*domegat+x0(2,i));
        plot(t3,0.05*dEq2t+x0(3,i));
        plot(t4,0.05*dVtt+Vt0(i));
        hold off
    end
    % 1 step + 稳态初值
    if DeltaVf == 1.0
        plot(t1,ddeltat+x0(1,i)); 
        hold on
        plot(t2,domegat+x0(2,i));
        plot(t3,dEq2t+x0(3,i));
        plot(t4,dVtt+Vt0(i));
        hold off
    end
    legend('\delta','\omega','E_{q}’','V_{t}');
    titles = ["task3 Vf=1.05V" "task3 Vf=1.1V" "task3 Vf=1.15V"];
    title(titles(i)); 

% 输出方程

end

%% 修改电抗参数

A
B
C
C1
C2
C3


