#---------------------------------------------------------------------#
#This code computes the Chebyshev points & weights
#Written by F.X. Giraldo on 7/2014
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function chebyshev_gauss(P::Integer,DFloat)

    #Initialize arrays
    xgl=zeros(DFloat,P)
    wgl=zeros(DFloat,P)

    #Construct arrays
    for i=1:P
        xgl[P+1-i]=cos( (2*i-1)*pi/(2*P) );
        wgl[P+1-i]=pi/P
    end
    return (xgl,wgl)
end #function
