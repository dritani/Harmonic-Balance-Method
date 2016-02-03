clc;
clear;
a=0;%start time
b=10;%end time
format long;
max_error=zeros(6,5);
y_exact=@(t)2*atan(tanh(0.25*(t.^2-4*atanh(tan(0.5)))));
f=@(t,y)t*cos(y);
%
for j=0:1:6
    h=10^-j;
    N=(b-a)/h;
    y0=-1;
    t=zeros(N,1);
    %initial values
    t(1)=a;
    w=y0;
    v=y0;
    
    A=zeros(N,8);
    A(1,1)=a;
    A(1,2)=w;
    A(1,3)=y0;
    
    A(1,4)=abs(y0-w);
    A(1,5)=((y0-w)/h);%local truncation error of euler method
    
    A(1,6)=v;
    A(1,7)=abs(y0-v);%global truncation error of euler method
    A(1,8)=((y0-v)/h);%local truncation error of euler method    
   
    
    %p_y=y0;
    p_w=y0;
    p_v=y0;

    for i=2:1:N
        %numerical solution after each steps
        w=p_w+h*f(t(i-1),p_w);
          
        %global truncation error of euler method
        A(i,4)=abs(y_exact(t(i-1))-p_w);
        
        A(i,5)=abs(y_exact(t(i-1))-w)/h;
        %A(i+1,5)=((y-p_w)/h-cos(p_w));%local truncation error of euler method    
       
        
        K1=h*f(t(i-1),p_v);%Runge-Kutta
        K2=h*f((t(i-1)+h/2),(p_v+K1/2));
        K3=h*f((t(i-1)+h/2),(p_v+K2/2));
        K4=h*f((t(i-1)+h),(p_v+K3));
        phi=(K1+2*K2+2*K3+K4)/(6*h);
        v=p_v+h*phi;
           
         
        t(i)=a+i*h;  
        
        A(i,1)=t(i);
        A(i,2)=p_w;
        A(i,3)=y_exact(t(i-1)); %exact value
        
        %column 4,5,7,8 of matrix A are local and trucation errors of two
        %different method
        %Global is defined as the difference between approximation solution
        %and exact solution
             
        
        A(i,6)=p_v;
        A(i,7)=abs(y_exact(t(i-1))-p_v);%global truncation error of runge-kutta method
        
        
        %A(i+1,8)=((y-p_v)/h-phi);%local truncation error of runge-kutta method
        A(i,8)=abs(y_exact(t(i-1))-p_v)/h;
        %p_y=y;
        p_w=w;
        p_v=v;

    end
    %A
    
    
    max_error(j+1,1)= h;
    max_error(j+1,2)= max(A(:,4)); %global euler error
    max_error(j+1,3)= max(A(:,5)); %local truncation error--euler
    max_error(j+1,4)= max(A(:,7)); %global runge-kutta error
    max_error(j+1,5)= max(A(:,8)); %local truncation error--runge-kutta
       
end
%A
%max_error
loglog(max_error(:,1),max_error(:,2),max_error(:,1),max_error(:,3),max_error(:,1),max_error(:,4),max_error(:,1),max_error(:,5));
hleg1=legend('Euler-Global Error','Euler-Local Error','Runge-Kutta-Global Error','Runge-Kutta-Local Error');
set(hleg1,'Location','NorthEast')
set(hleg1,'Interpreter','none')
xlabel('time-steps h');
ylabel('Maximum Truncation Errors from t=0 to t=10');
set(gca,'XDir','reverse');
hold off
