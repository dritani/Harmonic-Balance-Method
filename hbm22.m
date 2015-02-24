function [] = hbm22 () % HBM2.2 Solves the 8 equations
clc
syms x1 x2 a1 a2 a3 a4 b1 b2 b3 b4 t
w = 1.1;
% As instructed, x1(t) and x2(t) are expressed as a truncated
% Fourier series
x1=a1*cos(w*t)+a2*sin(w*t)+a3*cos(3*w*t)+a4*sin(3*w*t);
x2=b1*cos(w*t)+b2*sin(w*t)+b3*cos(3*w*t)+b4*sin(3*w*t);
% Define the first and second derivative for x1 and x2
x1_first=diff(x1,t);
x1_second=diff(x1,t,2);
x2_first=diff(x2,t);
x2_second=diff(x2,t,2);
% Equation 1 and 2 
eq1=x1_second+0.5*(x1_first-x2_first)+2*x1-x2-0.1*cos(w*t)+0.5*(x1^3);
eq2=x2_second+0.5*(x2_first-x1_first)+2*x2-x1;
% Define the limits of integration
a=0;
b=(2*pi()/w);
% Integrate equation 1 and 2 as instructed to obtain 8 nonlinear
% algebraic equations in the 8 unknowns
fcn(1)=int(eq1*cos(w*t), t, a, b);
fcn(2)=int(eq1*sin(w*t), t, a, b);
fcn(3)=int(eq1*cos(3*w*t), t, a, b);
fcn(4)=int(eq1*sin(3*w*t), t, a, b);
%
fcn(5)=int(eq2*cos(w*t), t, a, b);
fcn(6)=int(eq2*sin(w*t), t, a, b);
fcn(7)=int(eq2*cos(3*w*t), t, a, b);
fcn(8)=int(eq2*sin(3*w*t), t, a, b);
% The Jacobian is calculated the easy way. 'fcn' contains the 8
% above equations, while 'v' contains the variables we will differentiate
% with respect to
v = [a1 a2 a3 a4 b1 b2 b3 b4];
B = jacobian(fcn,v);
% The tolerance is set, as well as the initial guess for the 8 unknowns.
tol = 10^(-8);
LeNorme=1;
m = 1.0;
a1 = m; a2 = m; a3 = m; a4 = m;
b1 = m; b2 = m; b3 = m; b4 = m;
% All the 'fcn' functions, are symbolic. To evaluate them, the 'eval'
% command is used. The values are then stored in an array.
    for i=1:8
        f(i,1)=eval(fcn(i));
    end
% The same thing is done for the Jacobian matrix B
    B=eval(B);
% The initial guess for the 8 unknowns is stored in the array x. This 
% is also the array where the converged results will be stored
    x = [a1 a2 a3 a4 b1 b2 b3 b4]';
% Broyden's method with 50 max iterations and 10^(-8) tolerance
for k = 1 : 50
    s = B\(-f);
    x = x + s;
    a1 = x(1); a2 = x(2); a3 = x(3); a4 = x(4);
    b1 = x(5); b2 = x(6); b3 = x(7); b4 = x(8);
    for i = 1 : 8
        fnew(i,1) = eval(fcn(i));
    end
    y = fnew-f;
    LeNorme=norm(x);
    if LeNorme < tol
        break
    end
    f = fnew;
    B = B + ((y-B*s)*s')/(s'*s);
end
% Plots the results in the time domain.
t = linspace(0,100,1000);
x1 = a1*cos(w*t) + a2*sin(w * t) + a3*cos(3*w*t) + a4*sin(3*w*t);
x2 = b1*cos(w*t) + b2*sin(w*t)+ b3*cos(3*w*t) + b4*sin(3*w*t);
%
x
subplot(2,1,1)
plot(t,x1)
xlabel('Time(s)');
ylabel('Amplitude');
title('x1(t)');
%
subplot(2,1,2)
plot(t,x2)
xlabel('Time(s)');
ylabel('Amplitude');
title('x2(t)');
end