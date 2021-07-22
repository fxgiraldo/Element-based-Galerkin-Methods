#---------------------------------------------------------------------#
#This function warps a 2D grid
#Written by F.X. Giraldo on 5/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function warp_grid!(coord,Npoin,DFloat)

    for i=1:Npoin
        x, y, z = coord[1,i], coord[2,i], coord[3,i]
        coord[1,i]=x + sin(π * x) * sin(2 * π * y) * sin(2 * π * z) / 50
        coord[2,i]=y + sin(2 * π * x) * sin(π * y) * sin(2 * π * z) / 50
        coord[3,i]=z + sin(2 * π * x) * sin(2 * π * y) * sin(π * z) / 50
    end

end
