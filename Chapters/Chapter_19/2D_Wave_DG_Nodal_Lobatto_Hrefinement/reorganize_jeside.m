%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function jeside=reorganize_jeside(jeside,nelem)

    jcopy = [0 0 0 0];
    
    for ie=1:nelem
       
        jcopy = jeside(ie,:);
        jeside(ie,1) = jcopy(1);
        jeside(ie,2) = jcopy(3);
        jeside(ie,3) = jcopy(4);
        jeside(ie,4) = jcopy(2);
        
    end

end
