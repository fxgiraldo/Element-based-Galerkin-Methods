%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = scatter_to_hanging_nodes_CG_DG(rhs,face,ffc,nffc,ngl)


%go over all faces and find the hanging nodes
for ifa=1:nffc
    iface=ffc(ifa);
    
    plside = 0;
   if face(iface,9)==1 && face(iface,7)==0
       parent = face(iface,5);%identify parent element
       plside = face(iface,3);%parent local face
       child1 = face(iface,6);%identify children elements 
       child2 = face(iface,8);
   elseif face(iface,9)==1 && face(iface,8)==0
        parent = face(iface,6);
        plside = face(iface,4);
        child1 = face(iface,5);
        child2 = face(iface,7);
   else
       continue
   end
   
   parout = zeros(1,ngl); 
   ch1in = zeros(1,ngl);
   ch2in = zeros(1,ngl);
   parloc = zeros(1,ngl);
   
   
   %copy data from parent elements to scatter
   switch plside
       case 0
           continue
       case 1
           for j=1:ngl
               parout(j) = rhs(parent,j,1);
           end
       case 2
           for j=1:ngl
               parout(j) = rhs(parent,ngl,j);
           end
       case 3           
           for j=1:ngl
               parout(j) = rhs(parent,j,ngl);
           end
       case 4
           for j=1:ngl
               parout(j) = rhs(parent,1,j);
           end
   end 
    
   %scatter
    [ch1in,ch2in] = scatter_to_children(parout);
 
%     xgl = legendre_gauss_lobatto(ngl);
%    figure;
%    plot(xgl,parout,'o',(xgl-1)/2,ch1in,'gx',(xgl+1)/2,ch2in,'rx')
%    
   %distribute the data to hanging nodes only (copy, not DSS)
    switch plside
        case 1
            rhs(child1,1:ngl,1)=ch1in(1:ngl);
            rhs(child2,1:ngl,1)=ch2in(1:ngl);
        case 2
            rhs(child1,ngl,1:ngl)=ch1in(1:ngl);
            rhs(child2,ngl,1:ngl)=ch2in(1:ngl);
        case 3
            rhs(child1,1:ngl,ngl)=ch1in(1:ngl);
            rhs(child2,1:ngl,ngl)=ch2in(1:ngl);
        case 4
            rhs(child1,1,1:ngl)=ch1in(1:ngl);
            rhs(child2,1,1:ngl)=ch2in(1:ngl);
    end
    
end

end
