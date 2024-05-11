#=
-------------------------------------------------------------------------------------------------------------
This file runs the 1D Integration using 
The interpolation points used are the following:
ipoints=1: Lobatto
ipoints=2: Legendre
ipoints=3: Chebyshev
ipoints=4: Equi-spaced

This is part of Project 1 described in Sec. 4.4 (Eq. 4.13) in F.X. Giraldo's Introduction to Element-based Galerkin Methods using 
Tensor-Product Bases: Analysis, Algorithms, and Applications.

Written by F.X. Giraldo on July 6, 2021.
Department of Applied Mathematics
Naval Postgraduate School
Monterey, CA 93943
-------------------------------------------------------------------------------------------------------------
=#

using Plots, LinearAlgebra, FastGaussQuadrature

include("QuadraturePoints.jl")

#Some Constants
DFloat = Float64
Quadrature_type = "fxg"
#Quadrature_type = "julia"
Nmin=1
Nmax=64
Ns=101
c=π/2
iplot_interp=1
Npoints=4
machine_zero=eps(DFloat)

#Allocate Arrays
Narray=zeros(Int64,Nmax)
l1_norm_integration=zeros(DFloat,Nmax,Npoints)
l2_norm_integration=zeros(DFloat,Nmax,Npoints)

function main()

    @show(DFloat,Quadrature_type)

    #Loop through polynomial orders
    for ipoints=1:Npoints
        inop=0
        for N=Nmin:Nmax
            Q=N
            Np=N+1
            inop+=1
            Narray[inop]=N

            #Select Points
            if Quadrature_type == "fxg"
                if ipoints == 1
                    (ξ,ω) = QuadraturePoints.lobatto_gauss(Np) #Lobatto
                    ξq=copy(ξ); ωq=copy(ω)
                elseif ipoints == 2
                    (ξ,ω) = QuadraturePoints.legendre_gauss(Np) #Legendre
                    ξq=copy(ξ); ωq=copy(ω)
                elseif ipoints == 3
                    (ξ,ω) = QuadraturePoints.chebyshev_gauss(Np) #Chebyshev
                    (ξq,ωq) = QuadraturePoints.legendre_gauss(Np) #Legendre
                elseif ipoints == 4
                    (ξ,ω) = QuadraturePoints.equispaced_points(Np) #Equi-spaced
                    (ξq,ωq) = QuadraturePoints.legendre_gauss(Np) #Legendre
                end
            elseif Quadrature_type == "julia"
                if ipoints == 1
                    (ξ,ω) = gausslobatto(Np) #Lobatto
                    ξq=copy(ξ); ωq=copy(ω)
                elseif ipoints == 2
                    (ξ,ω) = gausslegendre(Np) #Legendre
                    ξq=copy(ξ); ωq=copy(ω)
                elseif ipoints == 3
                    (ξ,ω) = gausschebyshev(Np) #Chebyshev
                    (ξq,ωq) = gausslegendre(Np) #Legendre
                elseif ipoints == 4
                    (ξ,ω) = QuadraturePoints.equispaced_points(Np) #Equi-spaced
                    (ξq,ωq) = gausslegendre(Np) #Legendre
                end
            end

            #--------------------------------------------------#
            #Integration
            #--------------------------------------------------#
            #Compute Sample Space
            
            #=-------------------------------------------------
            Students Include the Lagrange Basis Function
            (ψ,dψ) = QuadraturePoints.lagrange_basis(Np,Ns,ξ,ξq)
            ----------------------------------------------------=#

            #Compute Expansion Coefficients
            q_coeff=zeros(DFloat,Np)
            for i=1:Np
                x=ξ[i]
                q_coeff[i]=cos(c*x)
            end #i

            #Compute Nth Order Integral
            qn=0
            for i=1:Np
                qsum=0
                for j=1:Np
                    qsum=qsum + ψ[j,i]*q_coeff[j]
                end #j
                qn=qn + ωq[i]*qsum
            end #i

            #Compute Exact Integral
            qe=4/π

            #Compute L1; L2; & L8 Norm
            l1_norm_integration[inop,ipoints]=abs( (qn - qe)/qe )
            l2_norm_integration[inop,ipoints]=max(sqrt( (qn - qe)^2/qe^2 ),machine_zero)
        end #N

    end #ipoints

    plot_handle=plot(Narray,l2_norm_integration,xlabel="N",ylabel="Error Norm",legend=true,lw=3,yaxis=:log,label=["Lobatto" "Legendre" "Chebyshev" "Equispaced"],title="L2 Integration Error")
    display(plot_handle)

    #Plot Interpolation
    println("Done") #output

end

#----------------------------------#
# Run the main function
#----------------------------------#
main()

