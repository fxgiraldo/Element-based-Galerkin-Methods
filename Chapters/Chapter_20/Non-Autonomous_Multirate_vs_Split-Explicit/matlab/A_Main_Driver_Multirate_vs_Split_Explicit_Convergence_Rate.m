%---------------------------------------------------------------------%
%This code compares various Time-Integrators for a simple two-rate ODE
%using various RK methods and Gragg Extrapolation.
%Written by F.X. Giraldo on 9/6/2019
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all;
close all;

plot_solution = 0;
plot_convergence = 1;

c=-100; %speed of fast wave
u0=1; %initial condition
time_final=0.5*pi; %final time in revolutions
ti_method=15;   
                %ti_method==1 is RK2,
                %ti_method==2 is SSP-RK,
                %ti_method==3 is Wicker-Skamarock RK3
                %ti_method==4 is WS RK General
                %ti_method==5 is WS RK General in Butcher form
                %ti_method==6 is Multistep WS RK General in Butcher form  
                %ti_method==7 is Multistep WS RK General in Butcher form
                %and substepping
                %ti_method==8 is Split-Explicit Wicker-Skamarock RK3
                %ti_method==9 is Split-Explicit Wicker-Skamarock General RK        
                %ti_method==10 is Split-Explicit Wicker-Skamarock RK3-RK2
                %ti_method==11 is Split-Explicit Wicker-Skamarock RK3-RK3
                %ti_method==12 is Split-Explicit Wicker-Skamarock RK3-ARK1
                %ti_method==13 is ERK in Butcher tableau
                %ti_method==14 is Multirate with ERK in Butcher tableau
                %ti_method==15 is Multirate with ERK and substepping
                %ti_method==16 is LSRK 5-stage or 14-stage 4th-order methods                                
                %ti_method==17 is Multirate LSRK(5,4) and (14,4)
                %ti_method==18 is Multirate LSRK(5,4) and (14,4) and
                %substepping
                %ti_method==19 is Multirate LSRK(5,4) and (14,4) and
                %adaptive substepping
                %ti_method==20 is LSRK in Butcher form
                %ti_method==21 is Multirate LSRK in Butcher form
                %ti_method==22 is Multirate LSRK in Butcher form and
                %substepping
                
Istages=3;
Jstages=3;
Msteps=10;
Msteps_max=0;
DT_Ratio=10; %ideal slow to fast ratio
if (ti_method == 1) %RK2
    Istages=min(3,Istages);
elseif (ti_method == 2) %SSP
    if (Istages ~= 2 && Istages ~=3)
        disp(['Error with Istages for this method']);
        disp(['Options are: 2 and 3. Exiting...']);
        return
    end
elseif (ti_method >= 5 && ti_method <= 7) %WSRK in Butcher form
    [alphaI,betaI] = construct_WSRK_coefficients(Istages);
    [alphaJ,betaJ] = construct_WSRK_coefficients(Jstages);
elseif (ti_method >= 8 && ti_method <= 12) %Split-Explicit WS-RK
    if (mod(Msteps,6) ~= 0)
        disp(['Error with Msteps for this method. ']);
        disp(['Msteps must be divisible by 6. Making it 6 ']);
        Msteps=6;
    end
elseif (ti_method >= 13 && ti_method <= 15) %ERK in Butcher form
    Istages=Istages;
    Jstages=Jstages;
    if (Istages == 2 || Istages == 3 || Istages == 4 || Istages == 6 || Istages == 8 || Istages == 10)
        [alphaI,betaI] = construct_ERK_coefficients(Istages);
    else
        disp(['Error with IStage value for this Method. Exiting... ']);
        disp(['Options are: 2, 3, 4, 6, 8, 10. Exiting... ']);
        return
    end
    if (Jstages == 2 || Jstages == 3 || Jstages == 4 || Jstages == 6 || Jstages == 8 || Jstages == 10)
        [alphaJ,betaJ] = construct_ERK_coefficients(Jstages);
    else
        disp(['Error with JStage value for this Method. Exiting... ']);
        disp(['Options are: 2, 3, 4, 6, 8, 10. Exiting... ']);
        return
    end
elseif (ti_method >= 16 && ti_method <= 19) %LSRK
    [alphaI,betaI,cI] = construct_LSRK_coefficients(Istages);
    [alphaJ,betaJ,cJ] = construct_LSRK_coefficients(Jstages);
elseif (ti_method >= 20 && ti_method <= 22) %LSRK in Butcher Form
    [A_vecI,B_vecI,C_vecI] = construct_LSRK_coefficients(Istages);
    [alphaI,betaI,cI] = convert_LSRK_to_Butcher(A_vecI,B_vecI);
    [A_vecJ,B_vecJ,C_vecJ] = construct_LSRK_coefficients(Jstages);
    [alphaJ,betaJ,cJ] = convert_LSRK_to_Butcher(A_vecJ,B_vecJ);
end

%Store time values
ntime_vector=[25 50 100 200 400 800 1600 3200];
ntime_begin=4;
ntime_end=8;

ntime_counter=0;
for ntime_loop = ntime_begin:ntime_end
    ntime=ntime_vector(ntime_loop);
    dt=time_final/(ntime-1);
    dt_fast_ideal=dt/DT_Ratio;
    ntime_counter=ntime_counter + 1;
    dt_vector(ntime_counter)=dt;
    
    %Compute Exact Solution
    time=0;
    time_vector=linspace(time,time_final,ntime);
    qe=zeros(ntime,1);
    qn=zeros(ntime,1);
    for i=1:ntime
        t=time_vector(i);
        qe(i) = u0*exp(c*t) + ( exp(c*t) - c*sin(t) - cos(t) )/(1 + c^2);
    end
    
    %Initialize State Vector
    qn(1)=u0;
    
    %Time Integration
    tic
    for itime=1:ntime-1
        ti=time;        
        tf=time+ dt;
        
        qq0=qn(itime);       
        if ti_method == 1 %RK2 Midpoint
            qqp = rk2(qq0,c,dt,time);
            method_text=['RK2'];
        elseif ti_method == 2 %SSP-RK
            qqp = ssp_rk(qq0,c,Istages,dt,time);
            method_text=['SSP-RK',num2str(Istages)];
        elseif ti_method == 3 %Wicker-Skamarock RK3
            qqp = WS_rk3(qq0,c,dt,time);
            method_text=['WS-RK3'];
        elseif ti_method == 4 %Wicker-Skamarock General RK
            qqp = WS_general_rk(qq0,c,Istages,dt,time);
            method_text=['General WS-RK',num2str(Istages)];
        elseif ti_method == 5 %Wicker-Skamarock General RK in Butcher Form
            qqp = erk_butcher(qq0,c,alphaI,betaI,Istages,dt,time);
            method_text=['WS-RK Butcher',num2str(Istages)];
        elseif ti_method == 6 %Multirate WSRK in Butcher form
            qqp = multirate_erk_butcher(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,dt,time);
            method_text=['Multirate WS Butcher I=',num2str(Istages),' J=',num2str(Jstages)];
        elseif ti_method == 7 %Multirate WSRK in Butcher form and Substepping
            qqp = multirate_erk_butcher_substepping(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,Msteps,dt,time);
            method_text=['Multirate WS Butcher I=',num2str(Istages),' J=',num2str(Jstages),' Msteps = ', num2str(Msteps)];
       elseif ti_method == 8 %Split-Explicit Wicker-Skamarock RK3
            %qqp = split_explicit_WS_rk3(qq0,c,Msteps,dt,time); %works in NUMA but poor convergence
            qqp = split_explicit_WS_rk3_original(qq0,c,Msteps,dt,time); %doesnt work in NUMA but 1st order convergence
            method_text=['Split-Explicit WS-RK3',' Msteps = ', num2str(Msteps)];
        elseif ti_method == 9 %Split-Explicit Wicker-Skamarock General RK
            %qqp = split_explicit_WS_rk3(qq0,c,Msteps,dt,time); 
            qqp = split_explicit_WS_general_rk(qq0,c,Istages,Msteps,dt,time); 
            method_text=['Split-Explicit WS-RK',num2str(Istages),' Msteps = ',num2str(Msteps)];
        elseif ti_method == 10 %Split-Explicit Wicker-Skamarock RK3-RK2
            qqp = split_explicit_WS_rk3_rk2(qq0,c,Msteps,dt,time); 
            method_text=['Split-Explicit WS-RK3-RK2',' Msteps = ', num2str(Msteps)]; 
        elseif ti_method == 11 %Split-Explicit Wicker-Skamarock RK3-RK3
            qqp = split_explicit_WS_rk3_rk3(qq0,c,Msteps,dt,time); 
            method_text=['Split-Explicit WS-RK3-RK3',' Msteps = ', num2str(Msteps)]; 
        elseif ti_method == 12 %Split-Explicit Wicker-Skamarock RK3-ARK1
            qqp = split_explicit_WS_rk3_ark1(qq0,c,Msteps,dt,time); 
            method_text=['Split-Explicit WS-RK3-ARK1',' Msteps = ', num2str(Msteps)]; 
        elseif ti_method == 13 %ERK in Butcher tableau form
            qqp = erk_butcher(qq0,c,alphaI,betaI,Istages,dt,time);
            method_text=['ERK',num2str(Istages)];
        elseif ti_method == 14 %Multirate with ERK in Butcher tableau form
            qqp = multirate_erk_butcher(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,dt,time);
            %qqp = multirate_erk_butcher_v2(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,dt,time);
            %Not working: can handle non-montonic time-interval
            method_text=['Multirate ERK I=',num2str(Istages),' J=',num2str(Jstages)];
        elseif ti_method == 15 %Multirate with ERK and Substepping
            qqp = multirate_erk_butcher_substepping(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,Msteps,dt,time);
            method_text=['Multirate ERK I=',num2str(Istages),' J=',num2str(Jstages),' Msteps = ', num2str(Msteps)];
        elseif ti_method == 16 %LSRK
            qqp = lsrk(qq0,c,alphaI,betaI,cI,dt,time);
            method_text=['LSRK I=',num2str(Istages)];
        elseif ti_method == 17 %Multirate LSRK
            qqp = multirate_lsrk(qq0,c,alphaI,betaI,cI,alphaJ,betaJ,cJ,dt,time);
            method_text=['Multirate LSRK I=',num2str(Istages),' J=',num2str(Jstages)];
        elseif ti_method == 18 %Multirate LSRK and substepping
            qqp = multirate_lsrk_substepping(qq0,c,alphaI,betaI,cI,alphaJ,betaJ,cJ,Msteps,dt,time);
            method_text=['Multirate LSRK I=',num2str(Istages),' J=',num2str(Jstages),' Msteps = ', num2str(Msteps)];
        elseif ti_method == 19 %Multirate LSRK and adaptive substepping
            [qqp,Msteps_max] = multirate_lsrk_substepping_adaptive(qq0,c,alphaI,betaI,cI,alphaJ,betaJ,cJ,dt,time,DT_Ratio,Msteps_max);
            method_text=['Multirate LSRK I=',num2str(Istages),' J=',num2str(Jstages),' with adaptive Msteps = ', num2str(Msteps_max)];
        elseif ti_method == 20 %LSRK in Butcher tableau form
            qqp = erk_butcher(qq0,c,alphaI,betaI,Istages,dt,time);
            method_text=['LSRK-Butcher',num2str(Istages)];
        elseif ti_method == 21 %Multirate LSRK in Butcher form
            qqp = multirate_erk_butcher(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,dt,time);
            method_text=['Multirate ERK I=',num2str(Istages),' J=',num2str(Jstages)];
        elseif ti_method == 22 %Multirate LSRK in Butcher form and Substepping
            qqp = multirate_erk_butcher_substepping(qq0,c,alphaI,betaI,Istages,alphaJ,betaJ,Jstages,Msteps,dt,time);
            method_text=['Multirate ERK I=',num2str(Istages),' J=',num2str(Jstages),' Msteps = ', num2str(Msteps)];
        end %if
        qn(itime+1)=qqp(1);
        time=time + dt;
    end %itime
    toc
    
    %Compute Norm
    top=0;
    bot=0;
    
    for i=1:ntime
        top=top + (qn(i)-qe(i))^2;
        bot=bot + qe(i)^2;
    end
    l2_norm(ntime_counter)=sqrt( top/bot );
    
    dt
    l2_norm
    method_text
    
end %ntime_loop

%Plot Solution
if (plot_solution == 1)
    h=figure;
    figure(h);
    plot_handle=plot(time_vector,qn,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(time_vector,qe,'b--');
    set(plot_handle,'LineWidth',2);
    xlabel('t','FontSize',18);
    ylabel('y(t)','FontSize',18);
    title_text=[method_text,', L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
end

%Compute Convergence Rate
rate=0;
for i=1:ntime_counter-1
    rate = rate + log(l2_norm(i+1)/l2_norm(i)) / log(dt_vector(i+1)/dt_vector(i));
end
rate=rate/(ntime_counter-1);

%Plot Convergence Rate
if (plot_convergence == 1)
    h=figure;
    figure(h);
    plot_handle=loglog(dt_vector,l2_norm,'r-x');
    set(plot_handle,'LineWidth',2);
    xlabel('dt','FontSize',18);
    ylabel('|| q||_{L^2}','FontSize',18);
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    title_text=[method_text,', Convergence Rate = ',num2str(rate)];
    title([title_text],'FontSize',18);
end
