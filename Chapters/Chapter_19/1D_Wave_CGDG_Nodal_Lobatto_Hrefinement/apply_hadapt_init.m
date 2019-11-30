%---------------------------------------------------------------------%
%This function creates the H-refinement arrays for all levels of
%refinement.
%Written by F.X. Giraldo on 3/10/2016
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [active] = apply_hadapt_init(children,active,elem_lev_pointer,coord,ngl,hadapt_lev,h_eps)

%Generate Refined grids
for iloop=1:hadapt_lev
    for e=elem_lev_pointer(1,iloop):elem_lev_pointer(2,iloop)
        %Form New Coords
        x0=coord(1,e);
        xN=coord(ngl,e);
        if (abs(x0) < h_eps && abs(xN) < h_eps)
            iparent=e;
            ichild1=children(1,iparent);
            ichild2=children(2,iparent);
            active(ichild1)=1;
            active(ichild2)=1;
            active(iparent)=0;
        end
    end %e
end %iloop

