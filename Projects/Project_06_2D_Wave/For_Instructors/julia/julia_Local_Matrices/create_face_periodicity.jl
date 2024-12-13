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

function create_face_periodicity!(face,iside,coord,Nface,Nboun)

    #Constant
    tol=DFloat(1e-6)
    
    #Local arrays
    Nhull=Int(Nboun/4)
    ileft=zeros(Int64,Nhull,1)
    iright=zeros(Int64,Nhull,1)
    itop=zeros(Int64,Nhull,1)
    ibot=zeros(Int64,Nhull,1)
    
    #initialize
    nleft=0
    nright=0
    ntop=0
    nbot=0
    
    #Find Extrema of Domain
    xmax=maximum(coord[1,:])
    xmin=minimum(coord[1,:])
    ymax=maximum(coord[2,:])
    ymin=minimum(coord[2,:])
    
    #loop thru sides and extract Left; Right; Bot; & Top
    for s=1:Nface
    
        #Check for Periodicity Edges
        ier=face[4,s]
        if (ier == -6) 
    
            i1=iside[1,s]
            i2=iside[2,s]
            xm=0.5*( coord[1,i1] + coord[1,i2] )
            ym=0.5*( coord[2,i1] + coord[2,i2] )
    
            #check Grid Point
            if ( abs(xm - xmin) < tol ) #left boundary
               nleft=nleft + 1
               ileft[nleft]=s
            elseif ( abs(xm - xmax) < tol ) #right boundary
               nright=nright + 1
               iright[nright]=s
            elseif ( abs(ym - ymin) < tol ) #bottom boundary
               nbot=nbot + 1
               ibot[nbot]=s
            elseif ( abs(ym - ymax) < tol ) #top boundary
               ntop=ntop + 1
               itop[ntop]=s
            else
               println("Error in CREATE_FACE_PERIODICITY. No match in PERIODIC_BCS for: s = ",s," ier = ",ier)
                exit
            end #if
        end #ier
    end #s
    
    #-----Loop through Periodic BCs
    
    #First: Do Left & Right
    for i=1:nleft
        sl=ileft[i]
        i1=iside[1,sl]
        yl1=coord[2,i1]
    
        #Search for Corresponding Right Edge
        for j=1:nright
            sr=iright[j]
            i2=iside[2,sr]
            yr2=coord[2,i2]
            if ( abs(yl1-yr2) < tol ) #they match
               face[2,sl]=face[1,sr]
               face[4,sl]=face[3,sr]
               face[1,sr]=sl #store the face number that this face replicates on the Left side
               face[3,sr]=-6; #means skip it due to Periodicity
               iside[4,sl]=iside[3,sr]
               iside[3,sr]=-6
               break
            end #if
        end #j   
    end #i
    
    #Second: Do Top & Bottom
    for i=1:nbot
        sl=ibot[i]
        i1=iside[1,sl]
        xl1=coord[1,i1]
    
        #Search for Corresponding Top Edge
        for j=1:ntop
            sr=itop[j]
            i2=iside[2,sr]
            xr2=coord[1,i2]
            if ( abs(xl1-xr2) < tol ) #they match
               face[2,sl]=face[1,sr]
               face[4,sl]=face[3,sr]
               face[1,sr]=sl #store the face number that this face replicates on the Left side
               face[3,sr]=-6; #means skip it due to Periodicity
               iside[4,sl]=iside[3,sr]
               iside[3,sr]=-6
               break
            end #if
         end #j   
    end #i

    return(face,iside)
end
