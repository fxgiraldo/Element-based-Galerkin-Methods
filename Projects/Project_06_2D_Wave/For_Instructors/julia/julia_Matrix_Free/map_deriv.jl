#---------------------------------------------------------------------#
#This function computes the Derivative Mapping for General 2D Quad Grids.
#Written by F.X. Giraldo on 5/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function map_deriv(ψ,dψ,f,Np,Nq,DFloat)

    #Allocate arrays
    f_ξ=zeros(DFloat,Nq,Nq)
    f_η=zeros(DFloat,Nq,Nq)

    #Loop through Integration Points
    for l=1:Nq, k=1:Nq
        sum_ξ=0
        sum_η=0
        for j=1:Np, i=1:Np
            sum_ξ=sum_ξ + dψ[i,k]*ψ[j,l]*f[i,j]
            sum_η=sum_η + ψ[i,k]*dψ[j,l]*f[i,j]
        end #i,j
        f_ξ[k,l]=sum_ξ
        f_η[k,l]=sum_η
    end #k,l
    return (f_ξ,f_η)
end
