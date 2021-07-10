#---------------------------------------------------------------------#
#This function computes the Initial & Exact Solutions.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function exact_solution(coord,Npoin,time,case,DFloat)

#Set some constants
xmin=0
xmax=+2
xc=0.5*(xmax+xmin)
xl=xmax-xmin
c=2.0

#Initialize
qe=zeros(DFloat,Npoin,2)
qb=zeros(DFloat,Npoin,1)
gravity=DFloat(1.0)
hmean=DFloat(1.0)

#Construct Solution
for I=1:Npoin
    x=coord[I]
    h=0.5*cos(c*pi*x)*cos(c*π*time)
    u=0.5*sin(c*pi*x)*sin(c*π*time)
    hb=hmean
    qe[I,1]=h
    qe[I,2]=u   
    qb[I]=hb
end #I
return (qe,qb,gravity)
end
