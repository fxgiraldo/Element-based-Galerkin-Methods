#=
---------------------------------------------------------------------
This function computes the Flux integral contribution to the RHS vector
for a Weak Form unified CGDG method for the 2D Wave Equation.

Written by F.X. Giraldo on July 14, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function create_rhs_flux!(rhs,q,u,face,normals,jac_face,ωq,mapL,mapR,intma,periodicity,DFloat)

    #Get lengths of arrays
    Np=size(intma,1)
    Nface=size(face,2)
       
    #local arrays
    ql=zeros(DFloat,Np,1)
    qr=zeros(DFloat,Np,1)
    ul=zeros(DFloat,Np,1)
    vl=zeros(DFloat,Np,1)

    #Construct Flux Integral contribution
    for f=1:Nface

        #Store Left Side Variables
        el=face[3,f]
        if (el != -6) #periodic bc
            ilocl=face[1,f]
            for l=1:Np
                #Get Pointers
                il=mapL[1,l,ilocl]
                jl=mapL[2,l,ilocl]
                IL=intma[il,jl,el]

                #Left Element
                ql[l]=q[IL]
                ul[l]=u[1,IL]
                vl[l]=u[2,IL]
            end #l  

            #Store Right Side Variables
            er=face[4,f]
            if (er > 0) 
                ilocr=face[2,f];            
                for l=1:Np
                    #Get Pointers
                    ir=mapR[1,l,ilocr]
                    jr=mapR[2,l,ilocr]
                    IR=intma[ir,jr,er]
                    
                    #Right Element
                    qr[l]=q[IR]
                    #ur[l]=u[IR]; #Not needed since U is continuous
                    #vr[l]=v[IR]; #Not needed since V is continuous
                end #l
            end #if er

            #Do Gauss-Lobatto Integration
            for l=1:Np
                wq=ωq[l]*jac_face[l,f]
                
                #Store Normal Vectors
                nxl=normals[1,l,f]
                nyl=normals[2,l,f]
                nxr=-nxl
                nyr=-nyl

                #Interpolate onto Quadrature Points
                qlq_k=ql[l]
                qrq_k=qr[l]
                u_k=ul[l]
                v_k=vl[l]    
                ul_k=u_k
                vl_k=v_k
                ur_k=u_k
                vr_k=v_k
                        
                #Compute Rusanov flux Constant
                unl=nxl*ul_k + nyl*vl_k
                unr=nxl*ur_k + nyl*vr_k
                claml=abs(unl)
                clamr=abs(unr)
                clam=max(claml,clamr)

                #Flux Variables
                fxl=qlq_k*ul_k
                fyl=qlq_k*vl_k

                fxr=qrq_k*ur_k
                fyr=qrq_k*vr_k

                #Normal Flux Component
                flux_ql=( nxl*fxl + nyl*fyl )
                flux_qr=( nxr*fxr + nyr*fyr )
                flux_q=flux_ql - flux_qr

                #Dissipation Term
                diss_q=clam*(qrq_k - qlq_k)

                #Construct Rusanov Flux
                flux_q=0.5*(flux_q - diss_q)

                #Loop through Side Interpolation Points
                #--------------Left Side------------------#
                il=mapL[1,l,ilocl]
                jl=mapL[2,l,ilocl]
                IL=intma[il,jl,el]
                IL=periodicity[IL]
                rhs[IL] -= wq*flux_q
                
                #--------------Right Side------------------#
                if (er > 0)
                    ir=mapR[1,l,ilocr]
                    jr=mapR[2,l,ilocr]
                    IR=intma[ir,jr,er]
                    IR=periodicity[IR]
                    rhs[IR] += wq*flux_q
                end #if         
            end #l   
        end #if el      
    end #f

    return(rhs)
end
