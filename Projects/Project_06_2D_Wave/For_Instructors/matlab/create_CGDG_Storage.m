%---------------------------------------------------------------------%
%This function creates the data structures for either CG or DG storage
%Written by F.X. Giraldo on 9/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,intma,iperiodic,DG_to_CG,npoin] = ...
    create_CGDG_Storage(storage_method,npoin_CG,npoin_DG,coord_CG,...
    intma_CG,iperiodic_CG,nelem,ngl)

if strcmp(storage_method,'cg')
    npoin=npoin_CG;
elseif strcmp(storage_method,'dg')
    npoin=npoin_DG;
end
DG_to_CG=zeros(npoin,1);

%Initialize Global Arrays
coord=zeros(2,npoin);
intma=zeros(ngl,ngl,nelem);
iperiodic=zeros(npoin,1);

if strcmp(storage_method,'cg')
    coord=coord_CG;
    intma=intma_CG;
    iperiodic=iperiodic_CG;
    for i=1:npoin_CG
        DG_to_CG(i)=i;
    end
    
elseif strcmp(storage_method,'dg')
    npoin=npoin_DG;
    i_DG=0;

    %Get INTMA
    for e=1:nelem
        for j=1:ngl
            for i=1:ngl
                i_DG=i_DG + 1;
                i_CG=intma_CG(i,j,e);
                intma(i,j,e)=i_DG;
                DG_to_CG(i_DG)=i_CG;
            end
        end
    end
    
    if (i_DG ~= npoin_DG)
        disp(['ERROR in CREATE_CGDG_STORAGE: i_DG /= npoin_DG = ',num2str(i_DG),' ',num2str(npoin_DG)]);
        exit;
    end
    
    %Get COORD
    for i_DG=1:npoin_DG
        i_CG=DG_to_CG(i_DG);
        coord(:,i_DG)=coord_CG(:,i_CG);
    end
    
    %Get IPERIODIC
    for i_DG=1:npoin_DG
        iperiodic(i_DG)=i_DG;
    end
    
end %if