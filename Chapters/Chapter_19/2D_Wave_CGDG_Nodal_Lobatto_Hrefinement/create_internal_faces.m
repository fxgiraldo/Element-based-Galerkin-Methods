%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [face,iface,facepa,jeside,jesideh] = create_internal_faces(child,ngl,face,facepa,intma,iface,jeside,jesideh)

  iface=iface+1;
   face(iface,1) = intma(child(1),1,ngl);
   face(iface,2) = intma(child(1),ngl,ngl);
   face(iface,3) = 2;
   face(iface,4) = 1;
   face(iface,5) = child(1);
   face(iface,7) = 0;
   face(iface,6) = child(4);
   face(iface,8) = 0;
   face(iface,9) = 0;
   
   facepa(iface) = 1;             
   
   jeside(child(1),2)=iface; %local to global face map
   jeside(child(4),1)=iface; %local to global face map
   jesideh(child(1),2)=iface; %local to global face map history
   jesideh(child(4),1)=iface; %local to global face map history
   
   iface=iface+1;
   face(iface,1) = intma(child(1),ngl,1);
   face(iface,2) = intma(child(1),ngl,ngl);
   face(iface,3) = 4;
   face(iface,4) = 3;
   face(iface,5) = child(1);
   face(iface,7) = 0;
   face(iface,6) = child(2);
   face(iface,8) = 0;
   face(iface,9) = 0;
    
   facepa(iface) = 1; 
   
   jeside(child(1),4)=iface; %local to global face map
   jeside(child(2),3)=iface; %local to global face map
   jesideh(child(1),4)=iface; %local to global face map
   jesideh(child(2),3)=iface; %local to global face map
   
   iface=iface+1;
   face(iface,1) = intma(child(2),1,ngl);
   face(iface,2) = intma(child(2),ngl,ngl);
   face(iface,3) = 2;
   face(iface,4) = 1;
   face(iface,5) = child(2);
   face(iface,7) = 0;
   face(iface,6) = child(3);
   face(iface,8) = 0;
   face(iface,9) = 0;
              
   facepa(iface) = 1; 
   
   jeside(child(2),2)=iface; %local to global face map
   jeside(child(3),1)=iface; %local to global face map
   jesideh(child(2),2)=iface; %local to global face map
   jesideh(child(3),1)=iface; %local to global face map
   
   iface=iface+1;
   face(iface,1) = intma(child(4),ngl,1);
   face(iface,2) = intma(child(4),ngl,ngl);
   face(iface,3) = 4;
   face(iface,4) = 3;
   face(iface,5) = child(4);
   face(iface,7) = 0;
   face(iface,6) = child(3);
   face(iface,8) = 0;
   face(iface,9) = 0;
       
   facepa(iface) = 1; 
   
   jeside(child(4),4)=iface; %local to global face map
   jeside(child(3),3)=iface; %local to global face map
   jesideh(child(4),4)=iface; %local to global face map
   jesideh(child(3),3)=iface; %local to global face map
   
end
