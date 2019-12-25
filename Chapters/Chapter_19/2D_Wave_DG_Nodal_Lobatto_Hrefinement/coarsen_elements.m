%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face)

fe(1,1) = 1; %joins face side and children elements
fe(1,2) = 2;
fe(2,1) = 2;
fe(2,2) = 3;
fe(3,1) = 4;
fe(3,2) = 3;
fe(4,1) = 1;
fe(4,2) = 4;
    
ip(1) = 3;
ip(2) = 4;
ip(3) = 1;
ip(4) = 2;

for i=1:ncor
    
    if(tm(icor(i)) == 0) %skip if already coarsened
        continue
    end
    
    %adjust element structure
    parent = tp(icor(i));
    tm(parent) = 1;

    for j=1:4
        tm(tc(j,parent)) = 0;
    end
    
    %adjust faces structure
    
    %internal faces of the element (children)
    facepa(jeside(tc(1,parent),2))=0;
    facepa(jeside(tc(1,parent),3))=0;
    facepa(jeside(tc(3,parent),1))=0;
    facepa(jeside(tc(3,parent),4))=0;
    
    %external faces of the element
    

    
    for j=1:4 %go over four external parent faces
        pface = jeside(parent,j);
        child1 = tc(fe(j,1),parent);
        child2 = tc(fe(j,2),parent);
        ch1face = jeside(child1,j);
        ch2face = jeside(child2,j);
        if(face(ch1face,9)==1) %% adjust if hanging node 
            if(face(ch1face,6)==child1 || face(ch1face,6) == child2)
                face(ch1face,6) = parent;
                face(ch1face,8) = 0;
            elseif(face(ch1face,5)==child1 || face(ch1face,5) == child2)
                face(ch1face,5) = parent;
                face(ch1face,7) = 0;
            end
            face(ch1face,9)=0;
            jeside(parent,j) = ch1face;
            facepa(ch1face) = 1;
            
            %adjust neighbour if it is a domain boundary
            
            fb = face(ch1face,11);
            if (fb>0)
               face(fb,9) = 0;
                if(face(fb,6)==child1 || face(fb,6) == child2)
                    face(fb,6) = parent;
                    face(fb,8) = 0;
                elseif(face(fb,5)==child1 || face(fb,5) == child2)
                    face(fb,5) = parent;
                    face(fb,7) = 0;
                end
            end
            
        elseif(ch1face~=ch2face) %% adjust if conforming grid
            facepa(ch1face) = 0; %deactivate both conforming faces
            facepa(ch2face) = 0;
            

                if(face(ch1face,5)==child1)
                    ngb1 = face(ch1face,6);
                elseif(face(ch1face,6)==child1)
                    ngb1 = face(ch1face,5);
                end
                if(face(ch2face,5)==child2)
                    ngb2 = face(ch2face,6);
                elseif(face(ch2face,6)==child2)
                    ngb2 = face(ch2face,5);
                end
                
            
            face(pface,9) = 1; %% adjust parent face to type 1 (hanging node)
            face(pface,3) = j;
            face(pface,4) = ip(j);
            face(pface,5) = parent;
            face(pface,7) = 0;
            face(pface,6) = ngb1;
            face(pface,8) = ngb2;
            
            
            %adjust the neighbours connectivity information
            fb1 = face(ch1face,11);
            if(fb1==0)
                jeside(ngb1,ip(j)) = pface;
                jeside(ngb2,ip(j)) = pface;
            end
            
            fb = face(pface,11);
            %adjust faces if on the boundary
            if(fb1>0)
                %disable conforming faces at the other side of the boundary
                facepa(face(ch1face,11)) = 0;
                facepa(face(ch2face,11)) = 0;
                
                face(fb,9) = 1; %% adjust parent face to type 1 (hanging node)
                face(fb,4) = j;
                face(fb,3) = ip(j);
                face(fb,6) = parent;
                face(fb,8) = 0;
                face(fb,5) = ngb1;
                face(fb,7) = ngb2;

                face(fb,11) = pface;
                
                pnface = jeside(tp(ngb1),ip(j));
                
                facepa(jeside(ngb1,ip(j))) = 0;
                facepa(jeside(ngb2,ip(j))) = 0;
                jeside(ngb1,ip(j)) = pnface;
                jeside(ngb2,ip(j)) = pnface;
                facepa(pnface) = 1;
            end
        end
    end
    %%enable external faces of parent element
    for j=1:4
        facepa(jeside(parent,j)) = 1;
    end
    
end

end
