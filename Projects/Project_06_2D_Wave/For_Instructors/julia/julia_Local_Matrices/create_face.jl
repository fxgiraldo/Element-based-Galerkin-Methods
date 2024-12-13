#=
---------------------------------------------------------------------
This function computes the integer arrays FACE, mapL, and mapR, which 
allows us to loop through the faces of the grid in order to perform flux integrals.

Written by F.X. Giraldo on July 14, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function create_face(iside,intma,Nface,Np)

    #global arrays
    face=zeros(Int64,4,Nface)
    mapL=zeros(Int64,2,Np,4)
    mapR=zeros(Int64,2,Np,4)
    
    #local arrays
    inode=zeros(Int64,4,1)
    jnode=zeros(Int64,4,1)
    
    #Construct Boundary Pointer
    inode[1]=1
    inode[2]=Np
    inode[3]=Np
    inode[4]=1
    jnode[1]=1
    jnode[2]=1
    jnode[3]=Np
    jnode[4]=Np
    
    #Construct IMAP arrays
    for l=1:Np
    
       #eta=-1
       mapL[1,l,1]=l
       mapL[2,l,1]=1
       mapR[1,l,1]=Np+1-l
       mapR[2,l,1]=1
    
       #ksi=+1
       mapL[1,l,2]=Np
       mapL[2,l,2]=l
       mapR[1,l,2]=Np
       mapR[2,l,2]=Np+1-l
    
       #eta=+1
       mapL[1,l,3]=Np+1-l
       mapL[2,l,3]=Np
       mapR[1,l,3]=l
       mapR[2,l,3]=Np
    
       #ksi=-1
       mapL[1,l,4]=1
       mapL[2,l,4]=Np+1-l
       mapR[1,l,4]=1
       mapR[2,l,4]=l
    end #l
    
    #loop thru the sides
    for i=1:Nface
    
       ip1=iside[1,i]
       ip2=iside[2,i]
       iel=iside[3,i]
       ier=iside[4,i]
    
       #check for position on Left Element
       for j=1:4
          j1=j
          j2=j+1
          if (j2 > 4) 
             j2=1
          end #j2   
          jp1=intma[inode[j1],jnode[j1],iel]
          jp2=intma[inode[j2],jnode[j2],iel]
    
          if (ip1 == jp1 && ip2 == jp2) 
             face[1,i]=j
             break #leave J loop
          end #ip1
       end #j
    
       #check for position on Right Element
       if (ier > 0)
          for j=1:4
             j1=j
             j2=j+1
             if (j2 > 4) 
                j2=1
             end #j2   
             jp1=intma[inode[j1],jnode[j1],ier]
             jp2=intma[inode[j2],jnode[j2],ier]
    
             if (ip1 == jp2 && ip2 == jp1) 
                face[2,i]=j
                break #leave J loop
             end #ip1
          end #j
       end #ier
       
       #Store Elements into FACE
       face[3,i]=iel
       face[4,i]=ier
    end #i

    return(face,mapL,mapR)

end