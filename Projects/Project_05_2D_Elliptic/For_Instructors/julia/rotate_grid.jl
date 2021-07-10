#---------------------------------------------------------------------#
#This function warps a 2D grid
#Written by F.X. Giraldo on 5/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function rotate_grid!(coord,Npoin,DFloat)

    α=π/4
    for i=1:Npoin
        x, y = coord[1,i], coord[2,i]
        coord[1,i]=cos(α)*x - sin(α)*y
        coord[2,i]=sin(α)*x + cos(α)*y
    end

end
