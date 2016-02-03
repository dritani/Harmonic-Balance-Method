clc;
clear;%f(t,y)=cosy;
a=0;%start time
b=10;%end time changed from 10 to 100
b_n=100;
j=-3;
h=10^j;%time step
N=(b-a)/h;% numbers of the steps
N_n=(b_n-a)/h;
%initial condition
t=a;
y0=-1;
v=y0;

B=zeros(N,2); 
B(1,1)=a;
B(1,2)=y0; 


for i=1:1:N;
        
    %Runge-Kutta
    K1=h*t*cos(v);
    K2=h*(t+h/2)*cos(v+K1/2);
    K3=h*(t+h/2)*cos(v+K2/2);
    K4=h*(t+h)*cos(v+K3);
    v=v+(K1+2*K2+2*K3+K4)/6;
    
    t=a+i*h; 
   
    B(i+1,1)=t;
    B(i+1,2)=v;
    
end

t=a;
y0=-1;
v=y0;

B_n=zeros(N_n,2);
B_n(1,1)=a;
B_n(1,2)=y0; 

for i=1:1:N_n;
        
    %Runge-Kutta
    K1=h*t*cos(v);
    K2=h*(t+h/2)*cos(v+K1/2);
    K3=h*(t+h/2)*cos(v+K2/2);
    K4=h*(t+h)*cos(v+K3);
    v=v+(K1+2*K2+2*K3+K4)/6;
    
    t=a+i*h; 
   
    B_n(i+1,1)=t;
    B_n(i+1,2)=v;
    
end

%plot
%plot vector field

%subplot(2,1,1);
plot(B_n(:,1),B_n(:,2),'-*b'); %plot runge-kutta
hold on
plot(B(:,1),B(:,2),'-*r');
hold off
xlabel('t');
ylabel('y');

    