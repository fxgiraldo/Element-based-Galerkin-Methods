%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [tm,facepa,face,jeside] = refine_elements(tc,tm,iref,nref,facepa,jeside,face)

    fe(1,1) = 1; %joins face side and children elements
    fe(1,2) = 2;
    fe(2,1) = 2;
    fe(2,2) = 3;
    fe(3,1) = 3;
    fe(3,2) = 4;
    fe(4,1) = 1;
    fe(4,2) = 4;
    
    ip(1) = 3;
    ip(2) = 4;
    ip(3) = 1;
    ip(4) = 2;

for i=1:nref
    parent = iref(i);
    
    %adjust the tree marker
    tm(parent) = 0;
    for j=1:4
        tm(tc(j,parent)) = 1;
    end
    
    
    
    %adjust the face structure
    
    %enable internal faces
    facepa(jeside(tc(1,parent),2))=1;
    facepa(jeside(tc(1,parent),3))=1;
    facepa(jeside(tc(3,parent),1))=1;
    facepa(jeside(tc(3,parent),4))=1;
    
    
    %adjust external faces
    for j=1:4
        pface = jeside(parent,j);
        if(face(pface,9)==0)
            if(face(pface,5) == parent)
               face(pface,5) = tc(fe(j,1),parent);
               face(pface,7) = tc(fe(j,2),parent);
            elseif(face(pface,6) == parent)
               face(pface,6) = tc(fe(j,1),parent);
               face(pface,8) = tc(fe(j,2),parent);               
            end
            face(pface,9) = 1;
            jeside(tc(fe(j,1),parent),j) = pface;
            jeside(tc(fe(j,2),parent),j) = pface;
            
        elseif(face(pface,9)==1)
            facepa(pface) = 0; %disable parent face
            facepa(jeside(tc(fe(j,1),parent),j)) = 1; %enable children faces
            facepa(jeside(tc(fe(j,2),parent),j)) = 1;
            
            %adjust neighbour connectivity information
            if(face(pface,5) == parent)
                jeside(face(pface,6),ip(j)) = tc(fe(j,1),parent);
                jeside(face(pface,8),ip(j)) = tc(fe(j,2),parent);
            elseif(face(pface,6) == parent)
                jeside(face(pface,5),ip(j)) = tc(fe(j,1),parent);
                jeside(face(pface,7),ip(j)) = tc(fe(j,2),parent);         
            end
            face(pface,9) = 2;
        end
            
        
    end
    
    
   
end

end
