    clc
  clear
  %%Initializing 
  A =[-1 0; 1 0];
   B = [0.9; 0];
   C =[0 0.5];
    D = 0;
  L=place(A',C',[-10 -15])';
  eig(A-L*C)
  eig(A)
  x=[-1;1]; % initial state
  xhat=[0;0]; % initial estimate
  XX=x;
  XXhat=xhat;
  T=40;
  UU=ones(1,T); % input signal

for k=1:T,
u=UU(k);
y(k)=C*x+D*u;
yhat(k)=C*xhat+D*u;
x=A*x+B*u;
error(k+1)=(y(k)-yhat(k));
xhat=A*xhat+B*u+L*(y(k)-yhat(k));
XX=[XX,x];
XXhat=[XXhat,xhat];
end

figure 
plot(1:T,yhat);
hold on 
plot(1:T,y);
hold off

figure 
plot(error);