%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Mmatrix = gather_from_hanging_nodes(Mmatrix,intma,face,nface,ngl)

%DSS of hanging nodes
%go over all faces and find the hanging nodes
for iface=1:nface
    
    plside = 0;
   if face(iface,9)==1 && face(iface,7)==0
       parent = face(iface,5);%identify parent element
       plside = face(iface,3);%parent local face
       child1 = face(iface,6);%identify children elements 
       child2 = face(iface,8);
   else if face(iface,9)==1 && face(iface,8)==0
        parent = face(iface,6);
        plside = face(iface,4);
        child1 = face(iface,5);
        child2 = face(iface,7);
       end
   end
   
   %locate parent and children elements in Mmatrix
   switch plside
       case 0
           continue
       case 1
           for j=1:ngl
               parloc(j) = intma(parent,j,1);
               ch1loc(j) = intma(child1,j,ngl);
               ch2loc(j) = intma(child2,j,ngl);
           end
       case 2
           for j=1:ngl
               parloc(j) = intma(parent,ngl,j);
               ch1loc(j) = intma(child1,1,j);
               ch2loc(j) = intma(child2,1,j);
           end
       case 3           
           for j=1:ngl
               parloc(j) = intma(parent,j,ngl);
               ch1loc(j) = intma(child1,j,1);
               ch2loc(j) = intma(child2,j,1);
           end
       case 4
           for j=1:ngl
               parloc(j) = intma(parent,1,j);
               ch1loc(j) = intma(child1,ngl,j);
               ch2loc(j) = intma(child2,ngl,j);
           end
   end
   
   
    %copy the data from children to gather
    parin = zeros(1,ngl);
    ch1out = zeros(1,ngl);
    ch2out = zeros(1,ngl);

    ch1out(1:ngl) = Mmatrix(ch1loc(1:ngl)); 
    ch2out(1:ngl) = Mmatrix(ch2loc(1:ngl));
   
   %gather
    [parin] = gather_from_children(ch1out,ch2out);
   
   %distribute the data and perform DSS
    for j=2:ngl-1
       Mmatrix(parloc(j)) = Mmatrix(parloc(j)) + parin(j);
    end
    
end

end
