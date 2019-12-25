%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [jeside,face,facepa,tm,ft,nft] = refine_elements_tree(iref,nref,jeside,jesideh,face,facepa,ngl,tc,tm)

% disp('in refine tree');
ft = 0; %faces touched
nft = 0; %number of faces touched

[fe,ip,~] = get_pointers(ngl);

for i=1:nref
    
    parent = iref(i); %store parent
    tm(parent) = 0; %disable parent element
    for j=1:4
        child(j) = tc(j,parent);
        tm(child(j)) = 1; %enable children elements
    end
    
    %enable internal faces
    facepa(jeside(child(1),4)) = 1;
    facepa(jeside(child(1),2)) = 1;
    facepa(jeside(child(3),1)) = 1;
    facepa(jeside(child(3),3)) = 1;
    
    for f=1:4 %go over external faces
        
        ch(1) = child(fe(f,1)); %store children elements on a face
        ch(2) = child(fe(f,2));
        
        ifexp = jeside(parent,f); %parent face
        if face(ifexp,5)==parent %determine where parent element is stored in the face structure
           ipar = 5;
           iopp = 6;
        elseif face(ifexp,6)==parent
           ipar = 6;
           iopp = 5;
        end
 
        
        %update parent face
        face(ifexp,ipar) = ch(1);
        face(ifexp,ipar+2) = ch(2);
        
        if face(ifexp,9)==0
           
            jeside(ch(1),f) = ifexp; %update jeside
            jeside(ch(2),f) = ifexp;
            face(ifexp,9)=1; %mark face as non-conforming
            
            %arrange face properly
            face = arrange_face_local(face,ifexp);
            %mark for recalculation
            nft=nft+1;
            ft(nft) = ifexp;
            
            bnf = face(ifexp,11);
            if bnf>0 %if boundary
                if face(bnf,5)==parent 
                   ipar = 5;
                   iopp = 6;
                elseif face(bnf,6)==parent
                   ipar = 6;
                   iopp = 5;
                end
                face(bnf,ipar) = ch(1);
                face(bnf,ipar+2) = ch(2);
                face(bnf,9) = 1;
                
                %arrange face properly
                face = arrange_face_local(face,bnf);
                %mark for recalculation
                nft=nft+1;
                ft(nft) = bnf;
            end
            
        elseif face(ifexp,9)==1
            ngb(1) = face(ifexp,iopp);
            ngb(2) = face(ifexp,iopp+2);
            face(ifexp,9) = 2; %mark as double-refined
            facepa(ifexp) = 0; %disable parent face
            for j=1:2
                ifexc(j) = jesideh(ch(j),f); %store children faces
                facepa(ifexc(j)) = 1; %enable children faces
                jeside(ch(j),f) = ifexc(j); %update jeside
%                 %make sure the face is set up correctly
%                 face(ifexc(j),3) = f;
%                 face(ifexc(j),4) = ip(f);
%                 face(ifexc(j),5) = ch(j);
%                 face(ifexc(j),7) = ngb(j);
                jeside(ngb(j),ip(f)) = ifexc(j); %update neighbour jeside
                
            end
            
            bnf = face(ifexp,11);
            if bnf>0 %if boundary
                if face(bnf,5)==parent 
                   ipar = 5;
                   iopp = 6;
                elseif face(bnf,6)==parent
                   ipar = 6;
                   iopp = 5;
                end
                face(bnf,ipar) = ch(1); %update old neighbour face
                face(bnf,ipar+2) = ch(2);
                face(bnf,9) = 2;
                
                ngb(1) = face(bnf,iopp); %store neighbours
                ngb(2) = face(bnf,iopp+2);
                
                facepa(bnf) = 2; %mark bnf as double-refined
                facepa(bnf) = 0; %disable old neighbour face
                for j=1:2
                    ifexn(j) = face(ifexc(j),11); %store neighbour face
                    facepa(ifexn(j)) = 1; %enable neighbour face
                    jeside(ngb(j),ip(f)) = ifexn(j); %update jeside
                    
                    face(ifexn(j),11) = ifexc(j);
                    
                end
            end
        end
    end
end
end
