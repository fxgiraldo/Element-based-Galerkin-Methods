using Plots, FastGaussQuadrature

include("QuadraturePoints.jl")

#Some Constants
DFloat = Float64
Quadrature_type = "julia" #fxg or julia
Nmin=1
Nmax=64
Ns=51
c=π/2
iplot_interp=1
Npoints=4
machine_zero=eps(DFloat)

#ipoints=1: Lobatto
#ipoints=2: Legendre
#ipoints=3: Chebyshev
#ipoints=4: Equi-spaced

#Allocate Arrays
Narray=zeros(Int64,Nmax)
l1_norm_interpolation=zeros(DFloat,Nmax,Npoints)
l2_norm_interpolation=zeros(DFloat,Nmax,Npoints)
l8_norm_interpolation=zeros(DFloat,Nmax,Npoints)

#{{{ Main
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
                elseif ipoints == 2
                    (ξ,ω) = QuadraturePoints.legendre_gauss(Np) #Legendre
                elseif ipoints == 3
                    (ξ,ω) = QuadraturePoints.chebyshev_gauss(Np) #Chebyshev
                elseif ipoints == 4
                    (ξ,ω) = QuadraturePoints.equispaced_points(Np) #Equi-spaced
                end
            elseif Quadrature_type == "julia"
                if ipoints == 1
                    (ξ,ω) = gausslobatto(Np) #Lobatto
                elseif ipoints == 2
                    (ξ,ω) = gausslegendre(Np) #Legendre
                elseif ipoints == 3
                    (ξ,ω) = gausschebyshev(Np) #Chebyshev
                elseif ipoints == 4
                    (ξ,ω) = QuadraturePoints.equispaced_points(Np) #Equi-spaced
                end
            end

            #--------------------------------------------------#
            #Interpolation
            #--------------------------------------------------#
            #Compute Sample Space
            xs=zeros(Ns)
            xs=range(-1,length=Ns,stop=1)
            (ψ,dψ) = QuadraturePoints.lagrange_basis(Np,Ns,ξ,xs)

            #Compute Expansion Coefficients
            q_coeff=zeros(DFloat,Np)
            for i=1:Np
                x=ξ[i]
                q_coeff[i]=cos(c*x)
            end #i

            #Compute Nth Order Interpolant
            qn=zeros(DFloat,Ns)
            for i=1:Ns
                qsum=0
                for j=1:Np
                    qsum=qsum + dψ[j,i]*q_coeff[j]
                end #j
                qn[i]=qsum
            end #i

            #Compute Exact Solution
            qe=zeros(DFloat,Ns)
            for i=1:Ns
                x=xs[i]
                qe[i]=-c*sin(c*x)
            end #i

            #Compute L1; L2; & L8 Norm
            l1_top=0; l1_bot=0
            l2_top=0; l2_bot=0
            l8_top=-1000; l8_bot=-1000
            for i=1:Ns
                l1_top=l1_top + abs(qn[i]-qe[i])
                l1_bot=l1_bot + abs(qe[i])
                l2_top=l2_top + (qn[i]-qe[i])^2
                l2_bot=l2_bot + qe[i]^2
                l8_top=max(l8_top, abs(qn[i]-qe[i]))
                l8_bot=max(l8_bot, abs(qe[i]))
            end
            l1_norm_interpolation[inop,ipoints]=( l1_top/l1_bot )
            l2_norm_interpolation[inop,ipoints]=max(sqrt( l2_top/l2_bot ), machine_zero)
            l8_norm_interpolation[inop,ipoints]=( l8_top/l8_bot )
        end #N

    end #ipoints

    closeall
    plot_handle=plot(Narray,l2_norm_interpolation,xlabel="N",ylabel="Error Norm",legend=true,lw=3,yaxis=:log,label=["Lobatto" "Legendre" "Chebyshev" "Equispaced"],title="L2 Derivative Error")
    display(plot_handle)

    #Plot Interpolation
    println("Done") #output

end
#}}} Main

#----------------------------------#
# Run the main function
#----------------------------------#
main()

