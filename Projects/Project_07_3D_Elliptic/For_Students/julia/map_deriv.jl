#---------------------------------------------------------------------#
#This function computes the Derivative Mapping for General 2D Quad Grids.
#Written by F.X. Giraldo on 5/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function map_deriv(ψ,dψ,f,Np,Nq,DFloat)

    #Allocate arrays
    f_ξ=zeros(DFloat,Nq,Nq,Nq)
    f_η=zeros(DFloat,Nq,Nq,Nq)
    f_ζ=zeros(DFloat,Nq,Nq,Nq)

    #Loop through Integration Points
    for n=1:Nq, m=1:Nq, l=1:Nq
        sum_ξ=0
        sum_η=0
        sum_ζ=0
        for k=1:Np, j=1:Np, i=1:Np
            sum_ξ=sum_ξ + dψ[i,l]*ψ[j,m]*ψ[k,n]*f[i,j,k]
            sum_η=sum_η + ψ[i,l]*dψ[j,m]*ψ[k,n]*f[i,j,k]
            sum_ζ=sum_ζ + ψ[i,l]*ψ[j,m]*dψ[k,n]*f[i,j,k]
        end #i,j,k
        f_ξ[l,m,n]=sum_ξ
        f_η[l,m,n]=sum_η
        f_ζ[l,m,n]=sum_ζ
    end #l,m,n
    return (f_ξ,f_η,f_ζ)
end
