#=
-----------------------------------------------------------------------------------------------
This program  constructs points and weights using:
 1) Chebyshev, 
 2) Legendre, 
 3) Lobatto, and
 4) Equi-spaced points

 This constitutes Project 1 and the results are described in Secs. 3.4 and 4.4 in the textbook.
 Written by F.X. Giraldo
            Department of Applied Mathematics
            Naval Postgraduate School
            Monterey, California 93943-5216
-----------------------------------------------------------------------------------------------
=#

using Plots
include("QuadraturePoints.jl")

#Some Constants
DFloat = Float64
#----------------------------------------Only Change this---------------------------------------
N=4
#-----------------------------------------------------------------------------------------------

function main()

    #Allocate arrays
    ξarray=zeros(DFloat,N+1,4)
    ωarray=zeros(DFloat,N+1,4)

    #Some Output
    @show (N)

    #Chebyshev Points
    (ξ,ω) = QuadraturePoints.chebyshev_gauss(N+1)
    println("Chebyshev Points ") #output
    @show (ξ,ω)
    ξarray[:,1]=ξ
    ωarray[:,1]=ω

    #Legendre Points
    (ξ,ω) = QuadraturePoints.legendre_gauss(N+1)
    println("Legendre Points ") #output
    @show (ξ,ω)
    ξarray[:,2]=ξ
    ωarray[:,2]=ω

    #Lobatto Points
    (ξ,ω) = QuadraturePoints.lobatto_gauss(N+1)
    println("Lobatto Points ") #output
    @show (ξ,ω)
    ξarray[:,3]=ξ
    ωarray[:,3]=ω

    #Equispaced Points
    (ξ,ω) = QuadraturePoints.equispaced_points(N+1)
    println("Equispaced Points ") #output
    @show (ξ,ω)
    ξarray[:,4]=ξ
    ωarray[:,4]=ω

    #Plot Example
    plot_handle=plot(ξarray,ωarray,xlabel="Root",ylabel="Weight",legend=true,lw=3,label=["Chebyshev" "Legendre" "Lobatto" "Equispaced"],title="Roots and Weights")
    display(plot_handle)
    savefig(plot_handle,"QuadraturePoints.png")

    #Plot Interpolation
    println("Done") #output

end

#Call Main Program
main()
