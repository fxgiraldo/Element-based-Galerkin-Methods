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
N=8;
s=round(2*N/3);

%Quad Filter
alpha_quad=1.0;
F_quad=2.0;

%Exponential Filter
alpha_exp=17.0;
F_exp=18.0;

%ERFC Filter
alpha_erf=1.0;
F_erf=12.0;

%Compute Filter Matrix
k=[0:0.03333:N]';
[kmax,lmax]=size(k);

%Quadratic Filter
alpha=alpha_quad;
F=F_quad;
sigma_quad=Quadratic_Filter(N,alpha,F,k,s,kmax);

%Exponential Filter
alpha=alpha_exp;
F=F_exp;
sigma_exp=Exponential_Filter(N,alpha,F,k,s,kmax);

%ERFC-Log Filter
alpha=alpha_erf;
F=F_erf;
erfc_log=1.0;
sigma_erfc_log=ERFC_Filter(N,alpha,F,k,s,kmax,erfc_log);

%ERFC Filter
alpha=alpha_erf;
F=F_erf;
erfc_log=0;
sigma_erfc=ERFC_Filter(N,alpha,F,k,s,kmax,erfc_log);

%PLOT Filters
figure;
hold on;
% plot(k,sigma_quad,'r-','LineWidth',2);
% plot(k,sigma_exp,'b-','LineWidth',2);
plot(k,sigma_erfc_log,'k','LineWidth',2);
plot(k,sigma_erfc,'m-','LineWidth',2);
% legend('quadratic','exponential','erfc-log','erfc');
legend('erfc-log','erfc');

xlabel('k','FontSize',18);
ylabel('\sigma_k(x)','FontSize',18);
set(gca, 'FontSize', 18);
hold on;

%Plot 2/3N line
sigma_s=[0:0.1:1]';
k_s=s*ones(11,1);
plot(k_s,sigma_s,'c:','LineWidth',2);