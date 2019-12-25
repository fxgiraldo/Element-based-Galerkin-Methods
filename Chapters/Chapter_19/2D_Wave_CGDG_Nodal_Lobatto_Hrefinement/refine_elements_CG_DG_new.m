%---------------------------------------------------------------------%
%This function refines selected elements and maintains CG connectivity
%given that intma and coord are properly preallocated
%For now the refinement cannot occur at the boundary
%Written by M.A. Kopera on 10/1011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,ft,nft] = refine_elements_CG_DG_new(iref,refel,refpt,coord,...
                                                intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl)

% disp('in refine new');
   [fe,ip,jc,kc,js,je,ks,ke,pi,pj] = get_pointers(ngl);
   
ft = 0; %faces touched
nft = 0; %number of faces touched
                                            
%mark elements for refinement
nref = size(iref,2);

for i=1:nref
    parent = iref(i);
    
    
    % update the element tree
   for j=1:4 
        child(j) = refel+j-1;
        tc(j,parent) = child(j);
        tp(child(j)) = parent;
        tm(child(j)) = 1;
        tl(child(j)) = tl(parent)+1;
   end
   tm(parent) = 0;
   % Create coordinates for new elements
   
   dx = coord(intma(parent,ngl,1),1)-coord(intma(parent,1,1),1);
   dy = coord(intma(parent,1,ngl),2)-coord(intma(parent,1,1),2);
   %a
   xstart(1) = coord(intma(parent,1,1),1); %lower bottom coordinates of each child
   ystart(1) = coord(intma(parent,1,1),2);
   %b
   xstart(2) = coord(intma(parent,1,1),1)+dx/2;
   ystart(2) = coord(intma(parent,1,1),2);
   %c
   xstart(3) = coord(intma(parent,1,1),1)+dx/2;
   ystart(3) = coord(intma(parent,1,1),2)+dy/2;
   %d
   xstart(4) = coord(intma(parent,1,1),1);
   ystart(4) = coord(intma(parent,1,1),2)+dy/2;
   
   
   % create intma and coord for four children elements
   for el=1:4
        for j=js(el):je(el)
            for k=ks(el):ke(el)
                if(j==jc(el))&&(k==kc(el))
                    intma(child(el),j,k) = intma(parent,j,k);
                else
                    coord(refpt,1) = ( xgl(j)+1 )*dx/4 + xstart(el);
                    coord(refpt,2) = ( xgl(k)+1 )*dy/4 + ystart(el);
                    intma(child(el),j,k) = refpt;
                    iperiodic(refpt) = refpt;
                    refpt = refpt + 1;
                end
            end
        end
        refel=refel+1;
   end
   %make sure the points on the internal faces coincide
    intma(child(2),1,:) = intma(child(1),ngl,:);   
    intma(child(3),:,1) = intma(child(2),:,ngl);
    intma(child(4),ngl,:) = intma(child(3),1,:); 
    intma(child(4),:,1) = intma(child(1),:,ngl);
    
   
   %create new internal faces
   [face,iface,facepa,jeside,jesideh] = create_internal_faces(child,ngl,face,facepa,intma,iface,jeside,jesideh);
   
   %external faces
   
   for f=1:4 %go over four external faces
       ifex = jeside(parent,f); %store the number of each parent face
       if face(ifex,5)==parent %determine where parent element is stored in the face structure
           ipar = 5;
           iopp = 6;
       elseif face(ifex,6)==parent
           ipar = 6;
           iopp = 5;
       end
       
       %replace parent element by two children elements
       face(ifex,ipar) = child(fe(f,1));
       face(ifex,ipar+2) = child(fe(f,2));
       bnf = face(ifex,11);       
       
       if face(ifex,9)==0 %if face is conforming
          face(ifex,9)=1; %change status to hanging node
          jeside(child(fe(f,1)),f) = ifex; %update jeside of both new children
          jeside(child(fe(f,2)),f) = ifex;
          
          
          %check whether neighbour elements were previously refined and
          %create new historical faces if necessary
          
          if tc(1,face(ifex,iopp))>0
              for j=1:2
                  ch(j) = child(fe(f,j));
                  nc(j) = tc(fe(ip(f),j),face(ifex,iopp)); %store neighbour children
                  if bnf==0
                        intma(ch(j),:,:) = renumber_face_points(intma(ch(j),:,:),intma(nc(j),:,:),ip(f)); %adjust numbering
                  end
                  
                  %create new face
                  iface = iface+1;
                    face(iface,1) = intma(ch(j),pi(f),pi(5-f));
                    face(iface,2) = intma(ch(j),pj(f),pj(5-f)); 
                    face(iface,3) = f; 
                    face(iface,4) = ip(f);
                    face(iface,5) = ch(j); 
                    face(iface,7) = 0; 
                    face(iface,6) = nc(j); 
                    face(iface,8) = 0; 
                    face(iface,9) = 0;
                    face(iface,10) = face(ifex,10);

                    jesideh(ch(j),f) = iface;
                    jesideh(nc(j),ip(f)) = iface;                
              end
              
          end
          
%         arrange face properly
          face = arrange_face_local(face,ifex);
          %mark for recalculation
          nft=nft+1;
          ft(nft) = ifex;
          
          %deal with boundary conditions

          if(bnf>0)
            if face(bnf,5)==parent %determine where parent element is stored in the face structure
                   iparb = 5;
                   ioppb = 6;
            elseif face(bnf,6)==parent
                   iparb = 6;
                   ioppb = 5;
            end
             face(bnf,iparb) = child(fe(f,1));
             face(bnf,iparb+2) = child(fe(f,2)); 
             face(bnf,9) = 1;
             

             
             %create new historical faces if the boundary neighbour was
             %previously refined
             if tc(1,face(ifex,iopp))>0
                 for j=1:2
                    iface = iface+1;
                    
                    face(iface,1) = intma(nc(j),pi(ip(f)),pi(5-ip(f)));
                    face(iface,2) = intma(nc(j),pj(ip(f)),pj(5-ip(f))); 
                    face(iface,3) = ip(f); 
                    face(iface,4) = f;
                    face(iface,5) = nc(j); 
                    face(iface,7) = 0; 
                    face(iface,6) = ch(j); 
                    face(iface,8) = 0; 
                    face(iface,9) = 0;
                    face(iface,10) = face(bnf,10);
                    
                    face(iface,11) = jesideh(ch(j),f); %?? jesideh or jeside ??
                    face(jesideh(ch(j),f),11) = iface;
                    
                    jesideh(nc(j),ip(f)) = iface; 
                 end
             end
             
             %arrange face properly
             face = arrange_face_local(face,bnf);
             %mark for recalculation
             nft=nft+1;
             ft(nft) = bnf;
          end
       elseif face(ifex,9)==1 %if face has already a hanging node
          face(ifex,9)=2; %update status to double-refined 
          
          if(bnf>0) %if boundary update the boundary neighbour element on the other side (CHECK THAT!!)
               if face(ifex,5)==parent %determine where parent element is stored in the face structure
                   ipar = 5;
                   iopp = 6;
               elseif face(ifex,6)==parent
                   ipar = 6;
                   iopp = 5;
               end
             face(bnf,ipar) = child(fe(f,1)); %here
             face(bnf,ipar+2) = child(fe(f,2)); %here
             face(bnf,9) = 2;
          else %if not boundary then adjust the point numbering along common edges
             intma(child(fe(f,1)),:,:) = renumber_face_points(intma(face(ifex,ipar),:,:),intma(face(ifex,iopp),:,:),ip(f));
             intma(child(fe(f,2)),:,:) = renumber_face_points(intma(face(ifex,ipar+2),:,:),intma(face(ifex,iopp+2),:,:),ip(f));
          end
          
          facepa(ifex) = 0; %disable parent face
            for j=1:2 %create two new faces
                iface = iface+1;
                face(iface,1) = intma(child(fe(f,j)),pi(f),pi(5-f));
                face(iface,2) = intma(child(fe(f,j)),pj(f),pj(5-f)); 
                face(iface,3) = face(ifex,3); 
                face(iface,4) = face(ifex,4);
                face(iface,5) = face(ifex,5+2*(j-1)); 
                face(iface,7) = 0; 
                face(iface,6) = face(ifex,6+2*(j-1)); 
                face(iface,8) = 0; 
                face(iface,9) = 0;
                face(iface,10) = face(ifex,10);
                
                facepa(iface) = 1; %enable the new face
                jeside(child(fe(f,j)),f) = iface; %adjust jeside
                jesideh(child(fe(f,j)),f) = iface;                
            end
          
          
            if (bnf==0) %if not boundary
                for j=1:2 %adjust jeside of opposite elements
                    jeside(face(iface-2+j,iopp),ip(f)) = iface-2+j;  
                    jesideh(face(iface-2+j,iopp),ip(f)) = iface-2+j;    
                end
            else %if boundary 
                             
               facepa(bnf) = 0; %disable boundary neighbour face
               
               ngb(1) = face(ifex,iopp); %store neighbouring elements
               ngb(2) = face(ifex,iopp+2);
               
               for j=1:2 %create two new faces in place of bnf
                    iface = iface+1;
                    face(iface,1) = intma(ngb(j),pi(ip(f)),pi(5-ip(f)));
                    face(iface,2) = intma(ngb(j),pj(ip(f)),pj(5-ip(f))); 
                    face(iface,3) = ip(f); 
                    face(iface,4) = f;
                    face(iface,5) = ngb(j); 
                    face(iface,7) = 0; 
                    face(iface,6) = child(fe(f,j)); 
                    face(iface,8) = 0; 
                    face(iface,9) = 0;
                    face(iface,10) = face(bnf,10);
                    
                    face(iface,11) = iface-2;
                    face(iface-2,11) = iface;
                    
                    facepa(iface) = 1; %enable two new faces
                    jeside(ngb(j),ip(f)) = iface; %update jeside
                    jesideh(ngb(j),ip(f)) = iface; 
               end
            end
        end
          
       
   end
   
end

end