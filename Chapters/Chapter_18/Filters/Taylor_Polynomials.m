%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the DG method with LGL
%with 2nd or 3rd Order RK.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
N=4;

x=[-1:0.1:1]';
[xmax,lmax]=size(x);

T=zeros(N,xmax);

for i=1:N+1
    for j=1:xmax
        T(i,j)=1.0/(factorial(i-1))*x(j)^(i-1);
    end
end

%PLOT Polynomials
figure;
hold on;
%C = {'k','b','r','g','y','m'};
C=jet(N+1);
for i=1:N+1
    %plot(x(:),T(i,:),'color',C(i,:),'LineWidth',2);
    plot(x(:),T(i,:),'LineWidth',2);
end
legend('T_0','T_1','T_2','T_3','T_4');

xlabel('x','FontSize',18);
ylabel('T_k(x)','FontSize',18);
set(gca, 'FontSize', 18);