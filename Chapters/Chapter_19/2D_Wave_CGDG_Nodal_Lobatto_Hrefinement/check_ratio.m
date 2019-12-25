%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [cr,pref,npref] = check_ratio(parent,jeside,face)
    cr = 0;
    k=0;
    for f=1:4 %go around element faces
        
        ifex = jeside(parent,f); %store the number of each parent face
        flag = face(ifex,9);
        
        if flag==0
            continue
        elseif flag==1
            if face(ifex,5)==parent %determine where parent element is stored in the face structure
                ipar = 5;
                iopp = 6;
            elseif face(ifex,6)==parent
                ipar = 6;
                iopp = 5;
            end
            
            if face(ifex,ipar+2)==0
                continue
            end
            
            ngb = face(ifex,iopp);
            
            if tl(ngb)<tl(parent)
                cr=1;
                k=k+1;
                pref(k) = ngb;
            end 
        end    
    end

end
