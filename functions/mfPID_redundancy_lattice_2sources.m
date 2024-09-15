function [Re,D,Respec,Dspec] = mfPID_redundancy_lattice_2sources(Ispec,pY)

    %%% Redundancy of each atom

    Respec=zeros(4,size(pY,2)); Re=zeros(4,1);
    Dspec=zeros(4,size(pY,2)); D=zeros(4,1);
    for iY = 1:size(pY,2)
        
        Respec(1,iY)=min(Ispec([1 2],iY));
        Respec(2,iY)=Ispec(1,iY);
        Respec(3,iY)=Ispec(2,iY);
        Respec(4,iY)=Ispec(3,iY);
        
        %%% Partial information of each atom
        Dspec(1,iY) = Respec(1,iY);
        Dspec(2,iY) = Respec(2,iY)-Dspec(1,iY);
        Dspec(3,iY) = Respec(3,iY)-Dspec(1,iY);
        Dspec(4,iY) = Respec(4,iY)-sum(Dspec([1 2 3],iY));
    end

    Re = sum(pY.*Respec,2);
    D = sum(pY.*Dspec,2);
        
end