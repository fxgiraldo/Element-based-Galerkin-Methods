#=
---------------------------------------------------------------------
This function computes the integer array ISIDE which is needed to 
create the FACE arrays.

Written by F.X. Giraldo on July 12, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function create_side(intma,bsido,Npoin,Ne,Nboun,Nface,Np,DFloat)

    #global arrays
    iside = zeros(Int64,4,Nface)
    jeside= zeros(Int64,4,Ne)

    #local arrays
    lwher = zeros(Int64,Npoin,1)
    lhowm = zeros(Int64,Npoin,1)
    icone = zeros(Int64,5*Npoin,1)
    inode = zeros(Int64,4,1)
    jnode = zeros(Int64,4,1)

    #Fix lnode
    inode[1]=1
    inode[2]=Np
    inode[3]=Np
    inode[4]=1
    jnode[1]=1
    jnode[2]=1
    jnode[3]=Np
    jnode[4]=Np

    #count how many elements own each node
    for in0=1:4
        for e=1:Ne
            I=intma[inode[in0],jnode[in0],e]
            lhowm[I]=lhowm[I] + 1
        end #ie
    end #in0

    #track elements owning each node
    lwher[1]=0
    for I=2:Npoin
        lwher[I]=lwher[I-1] + lhowm[I-1]
    end #ip

    #another tracker array
    lhowm = zeros(Int64,Npoin,1)
    for in0=1:4
        for e=1:Ne
            I=intma[inode[in0],jnode[in0],e]
            lhowm[I]=lhowm[I] + 1
            jloca=lwher[I] + lhowm[I]
            icone[jloca]=e
        end #ie
    end #in

    #LOOP OVER THE NODES
    iloca=0
    for I=1:Npoin
        iloc1=iloca
        iele=lhowm[I]
        
        if (iele != 0 ) 
            iwher=lwher[I]

            #LOOP OVER THOSE ELEMENTS SURROUNDING NODE I
            ip1=I
            for el=1:iele
                e=icone[iwher+el]

                #find out position of ip in intma
                i0=0
                for in0=1:4
                    i0=in0
                    ipt=intma[inode[in0],jnode[in0],e]
                    if (ipt == I) 
                        break
                    end
                end #in0
                in0=i0

                #Check Edge of Element E which claims I
                j=0
                for jnod = 1:2:3
                    iold=0
                    j=j+1
                    in2=in0 + jnod
                    if (in2 > 4) 
                        in2=in2-4
                    end
                    ip2=intma[inode[in2],jnode[in2],e]
                    if (ip2 >= ip1) 

                        #check whether side is old | new
                        jloca=0
                        if (iloca != iloc1) 
                            for s=iloc1+1:iloca
                                jloca=s
                                if (iside[2,s] == ip2) 
                                    iold=1
                                    break
                                end
                            end #s 
                        end #iloca   
                        
                        if (iold == 0)
                            #NEW SIDE
                            iloca=iloca + 1
                            iside[1,iloca]=ip1
                            iside[2,iloca]=ip2
                            iside[2+j,iloca]=e
                        elseif (iold == 1) 
                            #OLD SIDE
                            iside[2+j,jloca]=e
                        end #iold
                    end #ip2   
                end #jnod
            end #iel

            #Perform some Shifting to order the nodes of a side in CCW direction
            for s=iloc1+1:iloca
                if (iside[3,s] == 0) 
                    iside[3,s]=iside[4,s]
                    iside[4,s]=0
                    iside[1,s]=iside[2,s]
                    iside[2,s]=ip1
                end #iside   
            end #s
        end #if iele    
    end #ip

    if (iloca != Nface) 
        println( "Error in CREATE_SIDE. iloca = ",iloca," Nface = ",Nface)
        exit
    end 

    #RESET THE BOUNDARY MARKERS
    for s=1:Nface
        if (iside[4,s] == 0) 
        il=iside[1,s]
        ir=iside[2,s]
        e=iside[3,s]
        for ib=1:Nboun
            ibe=bsido[3,ib]
            ibc=bsido[4,ib]
            if (ibe == e) 
                ilb=bsido[1,ib]
                irb=bsido[2,ib]
                if  (ilb == il && irb == ir)
                    iside[4,s]=-ibc
                    break
                end #ilb
            end #ibe
        end #ib
        end #iside
    end #s

    #FORM ELEMENT/SIDE CONNECTIVITY ARRAY
    for s=1:Nface
        el=iside[3,s]
        er=iside[4,s]
        is1=iside[1,s]
        is2=iside[2,s]

        #LEFT SIDE
        for in0=1:4
            i1=intma[inode[in0],jnode[in0],el]
            in1=in0 + 1
            if (in1 > 4) 
                in1=1
            end #in1
            i2=intma[inode[in1],jnode[in1],el]
            if ((is1 == i1) && (is2 == i2)) 
                jeside[in0,el]=s
            end #is1
        end #in0   

        #RIGHT SIDE
        if (er > 0) 
            for in0=1:4
                i1=intma[inode[in0],jnode[in0],er]
                in1=in0 + 1
                if (in1 > 4) 
                    in1=1
                end #in1   
                i2=intma[inode[in1],jnode[in1],er]
                if ((is1 == i2) && (is2 == i1)) 
                    jeside[in0,er]=s
                end #is1   
            end #in   
        end #ier
    end #s

    return(iside,jeside)
end