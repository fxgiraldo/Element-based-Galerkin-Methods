using Plots

include("QuadraturePoints_students.jl")

N=1

function main()

    #Some Constants
    DFloat = Float64
    N=4

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
    #    plot_handle=plot(ξarray,ωarray,xlabel="Root",ylabel="Weight",legend=true,lw=3,label=["Chebyshev" "Legendre" "Lobatto" "Equispaced"],title="Roots and Weights",seriestype=:scatter)
    plot_handle=plot(ξarray,ωarray,xlabel="Root",ylabel="Weight",legend=true,lw=3,label=["Chebyshev" "Legendre" "Lobatto" "Equispaced"],title="Roots and Weights")
    display(plot_handle)
    savefig(plot_handle,"QuadraturePoints.png")

    #Plot Interpolation
    println("Done") #output

end

#Call Main Program
main()
