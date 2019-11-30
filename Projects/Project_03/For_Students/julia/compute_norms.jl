#---------------------------------------------------------------------#
#This function computes the Grid and Elements in 1D.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function compute_norms(qn,qe,npoin,DFloat)

    norms=zeros(DFloat,3)

    l1_top, l1_bot = 0, 0
    l2_top, l2_bot = 0, 0
    l8_top, l8_bot =-1000, -1000
    for i=1:npoin
        l1_top=l1_top + abs(qn[i]-qe[i])
        l1_bot=l1_bot + abs(qe[i])
        l2_top=l2_top + (qn[i]-qe[i])^2
        l2_bot=l2_bot + qe[i]^2
        l8_top=max(l8_top, abs(qn[i]-qe[i]))
        l8_bot=max(l8_bot, abs(qe[i]))
    end
    norms[1]=l1_top/l1_bot
    norms[2]=sqrt(l2_top/l2_bot)
    norms[3]=l8_top/l8_bot

    return (norms)
end
