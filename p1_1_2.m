clear;
%f(t,y)=cosy;
a=0;%start time
b=10;%end time 
j=0;
h=10^j;%time steps
N=(b-a)/h;% number of steps

%initial conditions
t=a;
y0=-1;
%initial condition for Euler's method
w=y0; 
%initial condition for Runge-Kutta method
v=y0;

A=zeros(N,2); 
A(1,1)=a;
A(1,2)=y0; 

B=zeros(N,2); 
B(1,1)=a;
B(1,2)=y0; 
for i=1:1:N;
    w=w+h*cos(w);%euler
    
    K1=h*cos(v);%Runge-Kutta
    K2=h*cos(v+K1/2);
    K3=h*cos(v+K2/2);
    K4=h*cos(v+K3);
    v=v+(K1+2*K2+2*K3+K4)/6;
    
    t=a+i*h; 
    
    %output of euler
    A(i+1,1)=t;
    A(i+1,2)=w; 
   
    %output of runge-kutta
    B(i+1,1)=t;
    B(i+1,2)=v;
    
end

%plot
%plot vector field
[T,Y]=meshgrid(0:.4:10,-3:.24:3);
dY=cos(Y);
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