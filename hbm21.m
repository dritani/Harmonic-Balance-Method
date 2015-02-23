function [] = hbm21 () 
clc
% This code is an adaptation of the algorithm posted on the 
% MECH 309 webCT page in the folder "week 11".
%
a = 0;                  % start time
b = 200;                % end time
N = 2000;              % # of steps
h= (b-a)/N;             % step size
%
t = a; 
z(1) = 1;                 % initial value of x1 at t = a
z(2) = 0;                 % initial value of x1' at t = a
z(3) = 1;                 % initial value of x2 at t = a
z(4) = 0;                 % initial value of x2' at t = a
% Defines a matrix that will store each iteration for x1 and x2.
% It does not store the values of x1' or x2'
P=zeros(N,2);
% Runge Kutta Algorithm as discussed in class notes.
for i = 1:1:N
    K1 = h*ODE(t,z);
    K2 = h*ODE(t + h/2, z + K1/2);
    K3 = h*ODE(t + h/2, z + K2/2);
    K4 = h*ODE(t + h, z + K3);
%
    z = z + (K1 + 2*K2 + 2*K3 + K4)/6;
% Each iteration of x1 and x2 (i.e. z(1) and z(3)) are stored in P.
    P(i,1)=z(1);
    P(i,2)=z(3);
    t = a + i*h;
end
% Plots x1 and x2
t=linspace(a,b,N);
subplot(2,1,1)
plot(t,P(:,1))
xlabel('Time(s)');
ylabel('Amplitude');
title('x1(t)');
%
subplot(2,1,2)
plot(t,P(:,2))
xlabel('Time(s)');
ylabel('Amplitude');
title('x2(t)');
%
% Defines the system of two 2nd order ODEs as a system of
% four 1st order ODEs using the substitution z1=x1, z2=x1',
% z3=x2, z4=x2'
%
function z_prime = ODE(t, z)
%
m=1; k=1; d=0.5; e=0.5; w=1.1; F=0.1;
z_prime(1) = z(2);
z_prime(2) = (F*cos(w*t)-e*(z(1)^3)-d*z(2)+d*z(4)-2*k*z(1)+k*z(3))/m;
z_prime(3) = z(4);
z_prime(4) = (-d*z(4)+d*z(2)-2*k*z(3)+k*z(1))/m;
%
end
end