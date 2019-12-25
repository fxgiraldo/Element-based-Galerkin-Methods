%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function gpflag = flag_global_points(intma,sfc,nsfc,npoinT,ngl,jeside,face);

gpflag = zeros(npoinT,1);

%local pointers
ipoint = zeros(ngl,4);
jpoint = zeros(ngl,4);
for f=1:4
    switch f
        case 1
            ipoint(1:ngl,f)=1:ngl;
            jpoint(1:ngl,f)=1;
        case 2
            ipoint(1:ngl,f)=ngl;
            jpoint(1:ngl,f)=1:ngl;
        case 3
            ipoint(1:ngl,f)=1:ngl;
            jpoint(1:ngl,f)=ngl;
        case 4
            ipoint(1:ngl,f)=1;
            jpoint(1:ngl,f)=1:ngl;
    end
end

for is=1:nsfc
    ie = sfc(is);
    for i=1:ngl
        for j=1:ngl
            ip = intma(ie,i,j);
            gpflag(ip) = 1;
        end
    end 
    for f=1:4
        ifc = jeside(ie,f);
        iflag = face(ifc,9);
        if (iflag==1)
            if(face(ifc,5)~=0)&&(face(ifc,7)~=0)&&(face(ifc,3)==f)
                for i=1:ngl
                   ip = intma(ie,ipoint(i,f),jpoint(i,f));
                   gpflag(ip) = 2;
                end
                gpflag(face(ifc,1)) = 1;
                gpflag(face(ifc,2)) = 1;
            else if (face(ifc,6)~=0)&&(face(ifc,8)~=0)&&(face(ifc,4)==f)
                for i=1:ngl
                   ip = intma(ie,ipoint(i,f),jpoint(i,f));
                   gpflag(ip) = 2;
                end
                gpflag(face(ifc,1)) = 1;
                gpflag(face(ifc,2)) = 1;
                end
            end
        end
    end
end

end
