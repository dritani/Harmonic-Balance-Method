clc;
clear;



f=inline('t*cos(y)','t','y');
[t,y]=ode45(f,[0,1000],-1);

plot(t,y,'-b');

xlabel('t');
ylabel('y');