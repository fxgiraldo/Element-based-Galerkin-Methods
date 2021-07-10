#---------------------------------------------------------------------#
#This code computes the Equispaced points & weights
#Written by F.X. Giraldo on 4/2008
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function equispaced_points(P::Integer)

    #Initialize arrays
    xgl=zeros(P)
    wgl=zeros(P)

    #Construct array of points
    xgl=range(-1, length=P, stop=+1)

    return xgl,wgl
end #function
