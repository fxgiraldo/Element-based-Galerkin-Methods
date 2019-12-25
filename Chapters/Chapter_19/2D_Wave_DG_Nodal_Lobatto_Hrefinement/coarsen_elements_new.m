%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [jeside,face,facepa,tm,ft,nft] = coarsen_elements_new(icor,ncor,jeside,jesideh,face,facepa,ngl,tc,tp,tm)

% disp('in coarsen');
ft = 0; %faces touched
nft = 0; %number of faces touched

%assume icor has no siblings !!!

[fe,ip,~] = get_pointers(ngl);

for i=1:ncor
    
    %store parent
    parent = tp(icor(i)); 
    
    %store children
    for j=1:4
        child(j) = tc(j,parent); 
        tm(child(j)) = 0; %disable children elements
    end
    tm(parent) = 1; %enable parent element
    
    %disable internal faces
    facepa(jeside(child(1),4)) = 0;
    facepa(jeside(child(1),2)) = 0;
    facepa(jeside(child(3),1)) = 0;
    facepa(jeside(child(3),3)) = 0;
    
    %go over external faces of children
    for f=1:4
        ch(1) = child(fe(f,1))
        ch(2) = child(fe(f,2))
        ifex  = jeside(ch(1),f)
        ifex2 = jeside(ch(2),f)
        
            
        if face(ifex,9)==1 %if face is non-conforming
           
            %find the location of parent
            if face(ifex,5)==ch(1) || face(ifex,5)==ch(2) %determine where parent element is stored in the face structure
               ipar = 5;
               iopp = 6;
            elseif face(ifex,6)==ch(1) || face(ifex,6)==ch(2)
               ipar = 6;
               iopp = 5;
            end
            
            face(ifex,ipar) = parent; %replace children by parent
            face(ifex,ipar+2) = 0;
            
            face(ifex,9) = 0; %mark face as conforming
            
            bnf = face(ifex,11);
            if bnf>0 %if boundary face update the opposite side
                %find the location of parent
                if face(bnf,5)==ch(1) || face(bnf,5)==ch(2) %determine where parent element is stored in the face structure
                   ipar = 5;
                   iopp = 6;
                elseif face(bnf,6)==ch(1) || face(bnf,6)==ch(2)
                   ipar = 6;
                   iopp = 5;
                end

                face(bnf,ipar) = parent; %replace children by parent
                face(bnf,ipar+2) = 0;

                face(bnf,9) = 0; %mark face as conforming             
            end %boundary update
            
        elseif face(ifex,9)==0 %if the face is conforming
            
            facepa(ifex)=0; %disable children faces
            facepa(ifex2)=0;
            
            ifexh = jesideh(parent,f) %historical parent face
            facepa(ifexh) = 1; %enable parent face
            
            jeside(parent,f) = ifexh; % put it to jeside
            
            if face(ifexh,5)==ch(1) || face(ifexh,5)==ch(2) %determine where parent element is stored in the face structure
               iparh = 5;
               iopph = 6;
            elseif face(ifexh,6)==ch(1) || face(ifexh,6)==ch(2)
               iparh = 6;
               iopph = 5;
            end
            ngb(1) = face(ifexh,iopph); %store neigbour elements
            ngb(2) = face(ifexh,iopph+2);
            
            face(ifexh,iparh) = parent; %update historic face
            face(ifexh,iparh+2) = 0;
            
            face(ifexh,9) = 1; %mark is as non-conforming
            
            jeside(ngb(1),ip(f)) = ifexh; % update neighbour jeside
            jeside(ngb(2),ip(f)) = ifexh;
            
            bnfa(1) = face(jeside(ch(1),f),11); %locate boundary neighbour faces
            bnfa(2) = face(jeside(ch(2),f),11);
            
            if bnfa(1)>0 %if boundary face
                for j=1:2 %go over two children faces
                
                    %locate opposite element and face
                    if face(bnfa(j),5)==ch(j) 
                       ipar = 5;
                       iopp = 6;
                    elseif face(bnfa(j),6)==ch(j)
                       ipar = 6;
                       iopp = 5;
                    end
                    
%                     ngb(j) = face(bnfa(j),iopp); %store neighbour elements
                    facepa(bnfa(j)) = 0; %disable neighbour faces HERE TOO
        
                end
                
                ngbparent = tp(ngb(1));
                ifexn = jesideh(ngbparent,ip(f)); % parent neighbour face
%                 
%                 ifexn
%                 ngb
%                 
                face(ifexh,11) = ifexn; %link parent faces
                face(ifexn,11) = ifexh;
                jeside(ngb(1),ip(f)) = ifexn; % update neighbour jeside
                jeside(ngb(2),ip(f)) = ifexn;
                
                facepa(ifexn) = 1; %enable neighbour face
                facepa(ifexh) = 1; %enable parent face
                
                %adjust neighbour parent face
                if face(ifexn,5)==ngb(1) || face(ifexn,5)==ngb(2)
                       ipar = 5;
                       iopp = 6;
                       ipars = 3;
                       iopps = 4;
                elseif face(ifexn,6)==ngb(1) || face(ifexn,6)==ngb(2)
                       ipars = 4;
                       iopps = 3;
                       ipar = 6;
                       iopp = 5;
                end
%                 face(ifexn,ipar) = ngb(1);
%                 face(ifexn,ipar+2) = ngb(2);
%                 face(ifexn,iopp) = parent;
%                 face(ifexn,iopp+2) = 0;
                face(ifexn,ipars) = f;
                face(ifexn,iopps) = ip(f);
                face(ifexn,iopp) = ngb(1);
                face(ifexn,iopp+2) = ngb(2);
                face(ifexn,ipar) = parent;
                face(ifexn,ipar+2) = 0;
                face(ifexn,9) = 1;
                
                %arrange face properly
                face = arrange_face_local(face,ifexn);
                %mark for recalculation
                nft=nft+1;
                ft(nft) = ifexn;
            end %if boundary face
            
            %arrange face properly
            face = arrange_face_local(face,ifexh);
            %mark for recalculation
            nft=nft+1;
            ft(nft) = ifexh;
        end %conforming - non-conforming face
    end
end

end
