%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the 
%HDG method with 2 Order RK.
%This version constructs the Global Matrices which are good for 
%comparing CG and DG.
%Written by F.X. Giraldo on July 3, 2012
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nelem=16; %Number of Elements
nop=8;    %Interpolation Order

kstages=3; %RK2, RK3, RK34
dt=1e-3; %time-step, fraction of one revolution
Courant_max=0.3;
time_final=1.0; %final time in revolutions
nplots=50; %plotting variable - Ignore
iplot_solution=0; %Switch to Plot or Not.
iplot_matrices=0;
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
space_method_type=3; %1=CG, 2=DG, 3=HDG
cgdg_method=2; %1=separate or 2=unified
ti_method=5; %1=SDIRK1, 2=SDIRK2, 3=SDIRK3; 4=SDIRK4; 5=SDIRK5

icase=1; %case number: 1 is a Gaussian, 2 is a square wave, 3 is
         %Gaussian with source and 4 is square wave with source.
xmu=0.0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=100000000; %time-step frequency that the filter is applied.
diss=1;

%Constants
ngl=nop + 1;
    
%Time-Step Array
ti_method_array=[1 2 3 4 5];
ti_method_begin=1;
ti_method_end=5;
error0=1e1;
dt_array=[1e0 2.5e-1 1e-1 5e-2 2.5e-2 1e-2 5e-3 2.5e-3 1e-3];
idt_begin=1;
idt_end=7;
for ti_method_loop=ti_method_begin:ti_method_end
    ti_method=ti_method_array(ti_method_loop);
    icount=0;
    
    for idt=idt_begin:idt_end
        dt=dt_array(idt);
        icount=icount + 1;
        l2_norm(icount,ti_method_loop)=10*10^(ti_method-1)*dt^(ti_method);
        dt_total(icount,ti_method_loop)=dt;
    end %dt
end %ti_method_loop

h=figure;
figure(h);
for ti_method_loop=ti_method_begin:ti_method_end
    switch ti_method_loop
    case {1}
        plot_handle=loglog(dt_total(:,ti_method_loop),l2_norm(:,ti_method_loop),'r-');
    case (2)
        plot_handle=loglog(dt_total(:,ti_method_loop),l2_norm(:,ti_method_loop),'b-');
    case (3)
        plot_handle=loglog(dt_total(:,ti_method_loop),l2_norm(:,ti_method_loop),'g-');
    case (4)
        plot_handle=loglog(dt_total(:,ti_method_loop),l2_norm(:,ti_method_loop),'k-');
    case (5)
        plot_handle=loglog(dt_total(:,ti_method_loop),l2_norm(:,ti_method_loop),'m-');
    end %switch
    set(plot_handle,'LineWidth',2);
    hold on
end %for
title_text=['Theoretical Convergence Rates'];
title([title_text],'FontSize',18);
xlabel('\Delta t','FontSize',18);
ylabel('Normalized L^2 Error','FontSize',18);
legend('K=1','K=2','K=3','K=4','K=5');
set(gca, 'FontSize', 18);
axis([1e-3 1e0 1e-8 1e0]);
