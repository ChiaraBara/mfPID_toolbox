function out = mfPID_3sources_th(joint_prob,states,base)

    error(nargchk(2,3,nargin));
    if nargin < 3, base=0; end 

    %%% joint ptobability of [Y X1 X2 X3]
    
    nstate = size(states,2);
    statelab = {'Y','X1','X2','X3'};
    
    allcomb = {1,2,3,4,[1 2],[1 3],[1 4],...
        [2 3],[2 4],[3 4],[1 2 3],[1 2 4],[1 3 4],...
        [2 3 4],[1 2 3 4]};
    comb = {1,2,3,[1 2],[1 3],[2 3],[1 2 3]};
   
    for ic = 1:size(allcomb,2)
        
        substat = states(:,allcomb{ic});
        substat_un = unique(substat, 'rows');
        name = [];
        for in = 1:length(allcomb{ic})
            name = [name, statelab{allcomb{ic}(in)}];
        end
        eval(['C',name,' = substat_un;']);
        for i = 1: size(substat_un, 1)
            substat_un_tmp = substat_un(i, :);
            p_ind = find(ismember(substat, substat_un_tmp, 'rows'));
            p_substat(i) = sum(joint_prob(p_ind));
        end
        eval(['p',name,' = p_substat;']);
    end

    Ny=size(CY,2); N1=size(CX1,2); N2=size(CX2,2); N3=size(CX3,2);

    pX = {pX1,pX2,pX3,pX1X2,pX1X3,pX2X3,pX1X2X3};
    cX = {CX1,CX2,CX3,CX1X2,CX1X3,CX2X3,CX1X2X3};
    pYX = {pYX1,pYX2,pYX3,pYX1X2,pYX1X3,pYX2X3,pYX1X2X3};
    cYX = {CYX1,CYX2,CYX3,CYX1X2,CYX1X3,CYX2X3,CYX1X2X3};
    
    Ispec = zeros(size(comb,2),length(CY));
    for iy = 1:length(CY)
        for icomb = 1:size(comb,2)
            I = 0;
            name = [];
            for in = 1:length(comb{icomb})
                name = [name,'X',num2str(comb{icomb}(in))];
            end        
            CYX = eval(['CY',name]);
            CX = eval(['C',name]);
            pYX = eval(['pY',name]);
            pX = eval(['p',name]);
            N = size(CX,2);
            posy = find(ismember(CYX(:,1:Ny),CY(iy,:),'rows'));
            for k = 1:length(posy)
                if pYX(posy(k))~=0
                    [~,posx] = ismember(CYX(posy(k),Ny+1:Ny+N),CX,'rows');
                    if base == 2
                        I = I + (pYX(posy(k))/pY(iy))*log2(pYX(posy(k))/(pY(iy)*pX(posx)));
                    else
                        I = I + (pYX(posy(k))/pY(iy))*log(pYX(posy(k))/(pY(iy)*pX(posx)));
                    end
                end
            end   
            Ispec(icomb,iy) = I;
        end
    end
    
    [Re,D,Respec,Dspec] = mfPID_redundancy_lattice_3sources(Ispec,pY);

    U1spec = sum(Dspec([5,8],:));
    U2spec = sum(Dspec([6,9],:));
    U3spec = sum(Dspec([7,10],:));
    Rspec = sum(Dspec(1:4,:));
    Sspec = sum(Dspec(11:18,:));

    II = D(18)+D(1)-D(5)-D(6)-D(7)-D(12)-D(13)-D(14)-2*D(11);
    I = sum(D);
    I1 = Re(8);
    I2 = Re(9);
    I3 = Re(10);
    U1 = sum(D([5,8]));
    U2 = sum(D([6,9]));
    U3 = sum(D([7,10]));
    R = sum(D(1:4));
    S = sum(D(11:18));

    
    %%% OUTPUI
    out.pY = pY;
    out.I1 = I1;
    out.I2 = I2;
    out.I3 = I3;
    out.I = I;
    out.II = II;
    out.R = R;
    out.U1 = U1;
    out.U2 = U2;
    out.U3 = U3;
    out.S = S;
    out.D = D;
    out.Ispec = Ispec;
    out.Rspec = Rspec;
    out.U1spec = U1spec;
    out.U2spec = U2spec;
    out.U3spec = U3spec;
    out.Sspec = Sspec;
    out.Dspec = Dspec;
    
end
    