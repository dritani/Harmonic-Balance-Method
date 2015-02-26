clc;
clear;%f(t,y)=cosy;
a=0;%start time
b=10;%end time 
j=-1;
h=10^j;
N=(b-a)/h;
t=a;
y0=-1;
w=y0;
v=y0;

A=zeros(N,2); 
A(1,1)=a;
A(1,2)=y0; 

B=zeros(N,2); 
B(1,1)=a;
B(1,2)=y0; 
for i=1:1:N;
    w=w+h*t*cos(w);%euler
    
    K1=h*t*cos(v);%Runge-Kutta
    K2=h*(t+h/2)*cos(v+K1/2);
    K3=h*(t+h/2)*cos(v+K2/2);
    K4=h*(t+h)*cos(v+K3);
    v=v+(K1+2*K2+2*K3+K4)/6;
    
    t=a+i*h; 
    
    A(i+1,1)=t;
    A(i+1,2)=w; 
    
    B(i+1,1)=t;
    B(i+1,2)=v;
    
end

%plot
[T,Y]=meshgrid(0:.4:10,-3:.24:3);
dY=T.*cos(Y);
dT=ones(size(dY));
dYu=dY./sqrt(dT.^2+dY.^2);
dTu=dT./sqrt(dT.^2+dY.^2);
quiver(T,Y,dTu,dYu,'k');
hold on
plot(A(:,1),A(:,2),'-ob'); %plot euler method
plot(B(:,1),B(:,2),'-*r'); %plot runge-kutta
hold off
xlabel('t');
ylabel('y');
    