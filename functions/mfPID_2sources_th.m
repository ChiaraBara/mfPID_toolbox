function out = mfPID_2sources_th(joint_prob,states,base)

    error(nargchk(2,3,nargin));
    if nargin < 3, base=0; end 

    %%% joint ptobability of [Y X1 X2]
    
    nstate = size(states,2);
    statelab = {'Y','X1','X2','X3'};

    allcomb = {1,2,3,[1 2],[1 3],[2 3],[1 2 3]};
    comb = {1,2,[1 2]};
    
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

    Ny=size(CY,2); N1=size(CX1,2); N2=size(CX2,2); 

    pX = {pX1,pX2,pX1X2};
    cX = {CX1,CX2,CX1X2};
    pYX = {pYX1,pYX2,pYX1X2};
    cYX = {CYX1,CYX2,CYX1X2};
    
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
            posy = find(CYX(:,1:Ny) == CY(iy));
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

    [Re,D,Respec,Dspec] = mfPID_redundancy_lattice_2sources(Ispec,pY);

    Ii12spec = Ispec(3,:);
    Rspec = Dspec(1,:);
    U1spec = Dspec(2,:);
    U2spec = Dspec(3,:);
    Sspec = Dspec(4,:);

    Ii12 = sum(D);
    R = D(1);
    U1 = D(2);
    U2 = D(3);
    S = D(4);
       
    Ii1 = sum(pY.*Ispec(1,:));
    Ii2 = sum(pY.*Ispec(2,:));

    %%% OUTPUT
    out.pY = pY;
    out.I1 = Ii1;
    out.I2 = Ii2;    
    out.I = Ii12;
    out.R = R;
    out.U1 = U1;
    out.U2 = U2;
    out.S = S;
    out.D = D;
    out.Ispec = Ispec;
    out.Rspec = Rspec;
    out.U1spec = U1spec;
    out.U2spec = U2spec;
    out.Sspec = Sspec;
    out.Dspec = Dspec;
    
end