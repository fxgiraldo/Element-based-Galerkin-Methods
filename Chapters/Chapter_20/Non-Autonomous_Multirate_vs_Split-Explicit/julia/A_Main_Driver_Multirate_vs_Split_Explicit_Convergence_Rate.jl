#---------------------------------------------------------------------#
# This code compares various Time-Integrators for a simple two-rate ODE
# using various RK methods and Gragg Extrapolation.
# Written by F.X. Giraldo / P.R. Mugg on 9/20/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#

# Modules/Imports/Packages/Files
include("construct_ERK_coefficients.jl")
include("rhs_function.jl")
include("rk2.jl")
include("ssp_rk.jl")
include("WS_rk3.jl")
include("WS_general_rk.jl")
include("split_explicit_WS_rk3.jl")
#include("split_explicit_WS_rk3_original.jl")    #use ti_method==5
include("split_explicit_WS_general_rk.jl")
include("erk_butcher.jl")
include("multirate_erk_butcher.jl")
include("multirate_erk_butcher_substepping.jl")
using Plots; pyplot()
using LinearAlgebra

# Plot Setup (0-false; 1-true)
plot_solution = 1
plot_convergence = 1

c = -100          # speed of fast wave
u0 = 1            # initial condition
time_final = π/2  # final time in revolutions
ti_method = 7
                  # ti_method==1 is RK2,
                  # ti_method==2 is SSP-RK2 or SSP-RK3 (choose using Istages)
                  # ti_method==3 is Wicker-Skamarock RK3
                  # ti_method==4 is Wicker-Skamarock General RK
                  # ti_method==5 is Split-Explicit Wicker-Skamarock RK3
                  # ti_method==6 is Split-Explicit Wicker-Skamarock General RK
                  # ti_method==7 is ERK in Butcher tableau
                  # ti_method==8 is Multirate with ERK in Butcher tableau
                  # ti_method==9 is Multirate with ERK and substepping

Istages = 2
Jstages = 2
Msteps = 1
if (ti_method == 1)
    Istages = min(3,Istages)
elseif (ti_method == 5 || ti_method == 6)
    if ((Msteps % 6) != 0)
        println("Error with Msteps for this method. \n
                 Msteps must be divisible by 6. Making Msteps = 6 \n")
        Msteps = 6;
    end
elseif (ti_method >= 7)
    Istages += 1
    Jstages += 1
    if (Istages == 2 || Istages == 3 || Istages == 4 || Istages == 6 || Istages == 8 || Istages == 10)
        αI, βI = construct_ERK_coefficients(Istages)
    else
        println("Error with IStage value for this Method. Exiting... \n
                 Options are: 2, 3, 4, 6, 8, 10. Exiting... \n")
        return
    end
    if (Jstages == 2 || Jstages == 3 || Jstages == 4 || Jstages == 6 || Jstages == 8 || Jstages == 10)
        αJ, βJ = construct_ERK_coefficients(Jstages)
    else
        println("Error with JStage value for this Method. Exiting... \n
                 Options are: 2, 3, 4, 6, 8, 10. Exiting... \n")
        return
    end
end

#Store time values
ntime_vector = [25 50 100 200 400 800 1600 3200]
ntime_begin = 4
ntime_end = 8
Δt_vector = zeros(ntime_end - (ntime_begin-1))

l2_norm = []
ntime_counter = 0
for ntime_loop = ntime_begin:ntime_end
    ntime = ntime_vector[ntime_loop]
    Δt = time_final/(ntime-1)
    global ntime_counter += 1
    Δt_vector[ntime_counter] = Δt

    #Compute Exact Solution
    time = 0
    global time_vector = collect(range(time,length=ntime,stop=time_final))
    global qe = zeros(ntime,1)
    global qn = zeros(ntime,1)
    for i = 1:ntime
        t = time_vector[i]
        qe[i] = u0*exp(c*t) + (exp(c*t) - c*sin(t) - cos(t))/(1 + c^2)
    end

    #Initialize State Vector
    qn[1] = u0

    #Time Integration
    t1 = time_ns()
    for itime = 1:ntime-1
        qq0 = qn[itime]
        if ti_method == 1 #RK2 Midpoint
            qqp = rk2(qq0,c,Δt,time)
            global method_text = "RK2"
        elseif ti_method == 2 #SSP-RK2 or SSP-RK3 (choose Istages = 2 or 3)
            qqp = ssp_rk(qq0,c,Istages,Δt,time)
            global method_text = "SSP-RK$Istages"
        elseif ti_method == 3 #Wicker-Skamarock RK3
            qqp = WS_rk3(qq0,c,Δt,time)
            global method_text = "WS-RK3"
        elseif ti_method == 4 #Wicker-Skamarock General RK
            qqp = WS_general_rk(qq0,c,Istages,Δt,time)
            global method_text = "General WS-RK$Istages"
        elseif ti_method == 5 #Split-Explicit Wicker-Skamarock RK3
            qqp = split_explicit_WS_rk3(qq0,c,Msteps,Δt,time) #works in NUMA but poor convergence
            #qqp = split_explicit_WS_rk3_original(qq0,c,Msteps,Δt,time) #doesnt work in NUMA but 1st order convergence
            global method_text = "Split-Explicit WS-RK3, Msteps = $Msteps"
        elseif ti_method == 6 #Split-Explicit Wicker-Skamarock General RK
            qqp = split_explicit_WS_general_rk(qq0,c,Istages,Msteps,Δt,time)
            global method_text = "Split-Explicit General WS-RK$Istages, Msteps = $Msteps"
        elseif ti_method == 7 #ERK in Butcher tableau form
            qqp = erk_butcher(qq0,c,αI,βI,Istages,Δt,time)
            global method_text = "ERK$(Istages-1)"
        elseif ti_method == 8 #Multirate with ERK in Butcher tableau form
            qqp = multirate_erk_butcher(qq0,c,αI,βI,Istages,αJ,βJ,Jstages,Δt,time)
            global method_text = "Multirate ERK$(Istages-1)$(Jstages-1)"
        elseif ti_method == 9 #Multirate with ERK and Substepping
            qqp = multirate_erk_butcher_substepping(qq0,c,αI,βI,Istages,αJ,βJ,Jstages,Msteps,Δt,time)
            global method_text = "Multirate ERK$(Istages-1)$(Jstages-1), Msteps = $Msteps"
        end #if
        qn[itime+1] = qqp[1]
        time += Δt
    end #itime
    t2 = time_ns()
    total_time = (t2-t1)*1*10^(-9)

    #Compute Norm and store in l2_norm array
    push!(l2_norm, norm((qn-qe),2))

    @show method_text
    @show Δt
    @show l2_norm
    @show total_time; println()

end #ntime_loop

#Plot Solution
if (plot_solution == 1)
    plot_s = plot(time_vector,qn,color=:red,linestyle=:dash,lw=2,label="Num Sol")
    plot!(plot_s,time_vector,qe,color=:blue,linestyle=:dot,lw=2,label="Exact Sol")
    xlabel!(plot_s, "t")
    ylabel!(plot_s, "y(t)")
    title!(plot_s, "$method_text, \n L2 Norm (min) = $(minimum(l2_norm))")
    display(plot_s)
end

#Compute Convergence Rate
rate = 0
for i=1:ntime_counter-1
    global rate += log(l2_norm[i+1]/l2_norm[i]) / log(Δt_vector[i+1]/Δt_vector[i])
end
rate = rate/(ntime_counter-1)

#Plot Convergence Rate
if (plot_convergence == 1)
    plot_c = plot(Δt_vector,l2_norm,color=:red,xaxis=:log,yaxis=:log,lw=2,legend=false)
    scatter!(plot_c,Δt_vector,l2_norm,color=:red,markersize=4,markershape=:xcross)
    xlabel!(plot_c, "Δt")
    ylabel!(plot_c, "L2 Norm")
    title!(plot_c, "$method_text,\n Convergence Rate = $rate")
    display(plot_c)
end
