function [Re,D,Respec,Dspec] = mfPID_redundancy_lattice_3sources(Ispec,pY)

    %%% Redundancy of each atom

    Respec=zeros(18,size(pY,2)); Re=zeros(18,1);
    Dspec=zeros(18,size(pY,2)); D=zeros(18,1);
    for iY = 1:size(pY,2)
        
        Respec(1,iY)=min(Ispec([1 2 3],iY));
        Respec(2,iY)=min(Ispec([1 2],iY));
        Respec(3,iY)=min(Ispec([1 3],iY));
        Respec(4,iY)=min(Ispec([2 3],iY));
        Respec(5,iY)=min(Ispec([1 6],iY));
        Respec(6,iY)=min(Ispec([2 5],iY));
        Respec(7,iY)=min(Ispec([3 4],iY));
        Respec(8,iY)=Ispec(1,iY);
        Respec(9,iY)=Ispec(2,iY);
        Respec(10,iY)=Ispec(3,iY);
        Respec(11,iY)=min(Ispec([4 5 6],iY));
        Respec(12,iY)=min(Ispec([4 5],iY));
        Respec(13,iY)=min(Ispec([4 6],iY));
        Respec(14,iY)=min(Ispec([5 6],iY));
        Respec(15,iY)=Ispec(4,iY);
        Respec(16,iY)=Ispec(5,iY);
        Respec(17,iY)=Ispec(6,iY);
        Respec(18,iY)=Ispec(7,iY);

        %%% Partial information of each atom
        Dspec(1,iY) = Respec(1,iY);
        Dspec(2,iY) = Respec(2,iY)-Dspec(1,iY);
        Dspec(3,iY) = Respec(3,iY)-Dspec(1,iY);
        Dspec(4,iY) = Respec(4,iY)-Dspec(1,iY);
        Dspec(5,iY) = Respec(5,iY)-sum(Dspec([1 2 3],iY));
        Dspec(6,iY) = Respec(6,iY)-sum(Dspec([1 2 4],iY));
        Dspec(7,iY) = Respec(7,iY)-sum(Dspec([1 3 4],iY));
        Dspec(8,iY) = Respec(8,iY)-sum(Dspec([1 2 3 5],iY));
        Dspec(9,iY) = Respec(9,iY)-sum(Dspec([1 2 4 6],iY));
        Dspec(10,iY) = Respec(10,iY)-sum(Dspec([1 3 4 7],iY));
        Dspec(11,iY) = Respec(11,iY)-sum(Dspec(1:7,iY));
        Dspec(12,iY) = Respec(12,iY)-sum(Dspec([1:8,11],iY));
        Dspec(13,iY) = Respec(13,iY)-sum(Dspec([1:7,9,11],iY));
        Dspec(14,iY) = Respec(14,iY)-sum(Dspec([1:7,10,11],iY));
        Dspec(15,iY) = Respec(15,iY)-sum(Dspec([1:9,10:13],iY));
        Dspec(16,iY) = Respec(16,iY)-sum(Dspec([1:8,10:12,14],iY));
        Dspec(17,iY) = Respec(17,iY)-sum(Dspec([1:7,9:11,13,14],iY));
        Dspec(18,iY) = Respec(18,iY)-sum(Dspec(:,iY));
    end

    Re = sum(pY.*Respec,2);
    D = sum(pY.*Dspec,2);
        
end