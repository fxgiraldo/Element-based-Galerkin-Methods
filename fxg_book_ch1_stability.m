%---------------------------------------------------------------------%
%This code plots the Stability of for the 1D Wave Equation.
%Written by F.X. Giraldo on 10/97
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%Variables:
%---------------------------------------------------------------------%
clear; clear all; clf; close all;

method=input(' 1=FE & Upwind \n 2=FE & Centered \n 3=FE & Downwind \n 4=BE & Upwind \n 5=BE & Centered \n 6=BE & Downwind \n 7=TR & Upwind \n 8=TR & Centered \n 9=TR & Downwind \n Enter Method: ');

% m=100;
n=100;
dtheta=2*pi/n;
theta=0:dtheta:2*pi;
if (method <= 3) 
    cmax=1.0
elseif (method >= 4)
    cmax=10.0
end
m=10*cmax;
dc=cmax/m;
c=0:dc:cmax;
em1=exp(-i*theta);
ep1=exp(+i*theta);

%Plot Stability of Specific Equations
if (method == 1) %FE + Upwind
    for j=2:m+1
        z=1 - c(j)*(1 - em1); %1st order upwind (explicit 1st order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 2) %FE + Centered
    for j=2:m+1
        z=1 - 0.5*c(j)*(ep1 - em1); %2nd order centered (explicit 1st order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 3) %FE + Downwind
    for j=2:m+1
        z=1 - c(j)*(ep1 - 1); %1st order downwind (explicit 1st order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 4) %BE + Upwind
    for j=2:m+1
        z=1./(1 + c(j)*(1 - em1)); %1st order uwpind, (implicit 1st order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 5) %BE + Centered
    for j=2:m+1
        z=1./(1 + 0.5*c(j)*(ep1 - em1)); %2nd order centered, (implicit 1st order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 6) %BE + Downwind
    for j=2:m+1
        z=1./(1 + c(j)*(ep1 - 1)); %1st order downwind, (implicit 1st order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 7) %TR + Upwind
    for j=2:m+1
        z=(1 - 0.5*c(j)*(1 - em1))./(1 + 0.5*c(j)*(1 - em1)); %1st order upwind (Trapezoidal 2nd order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 8) %TR + Centered
    for j=2:m+1
        z=(1 - 0.25*c(j)*(ep1 - em1))./(1 + 0.25*c(j)*(ep1 - em1)); %2nd order centered (TR 2nd order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
elseif (method == 9) %TR + Downwind
    for j=2:m+1
        z=(1 - 0.5*c(j)*(ep1 - 1))./(1 + 0.5*c(j)*(ep1 - 1)); %1st order downwind (TR 2nd order in time)
        plot_handle=plot(z,'--b');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
end

%Plot Unit Circle
z=ep1;
plot_handle=plot(z,'-r');
set(plot_handle,'LineWidth',2);
hold on;
    
xlabel('Re(q)','FontSize',18); 
ylabel('Im(q)','FontSize',18); 
set(gca,'FontSize',18);
axis equal
% text(0,0,'\sigma=0.25','FontSize',18);
% text(0,0.25,'\sigma=0.5','FontSize',18);
% text(0,0.5,'\sigma=0.75','FontSize',18);
% text(0,1.0,'\sigma=1.0','FontSize',18);