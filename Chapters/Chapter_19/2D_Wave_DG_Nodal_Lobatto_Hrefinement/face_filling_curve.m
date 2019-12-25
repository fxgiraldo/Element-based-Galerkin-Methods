%---------------------------------------------------------------------%
%This subroutine creates the space filling curve out of the element tree
% note that nelem should be the original number of elements at the top
% level
% elmptr is the element pointer array which defines the element tree
% recursively
% nsfc is the length of sfc
%Written by M.A. Kopera on 10/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

function [ffc,nffc] = face_filling_curve(nface,facepa)

clear ffc;
m=0;
for f=1:nface
    if facepa(f)==1
        m=m+1;
        ffc(m)=f;
    end

end
%     if(facep(e)==0)
%        m=m+1;
%        ffc(m)=e;
%     else if facep(e)>0
%        el=facep(e);
%        for k=1:2
%            if (facep(el)>=0)
%             m=m+1;
%             ffc(m)=el;
%             el=el+1;
%            end
%        end
%     end
% end
nffc=m;
end