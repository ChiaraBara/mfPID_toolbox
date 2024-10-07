function out = mfPID_2sources_discrete(B,j,i1,i2)

    Y = B(:,j);
    X1 = B(:,i1);
    X2 = B(:,i2);

    Np=size(Y,1); %total number of patterns
    Ny=size(Y,2); N1=size(X1,2); N2=size(X2,2);
    base = 2;
    
    [~,pY,CY] = mfPID_H(Y,base);
    [~,pX1,CX1] = mfPID_H(X1,base);
    [~,pX2,CX2] = mfPID_H(X2,base);
    [~,pYX1,CX1Y] = mfPID_H([X1 Y],base);
    [~,pYX2,CX2Y] = mfPID_H([X2 Y],base);
    [~,pX1X2,CX1X2] = mfPID_H([X1 X2],base);
    [~,pYX1X2,CX1X2Y] = mfPID_H([X1 X2 Y],base);
    
    Ispec = zeros(2,length(CY));
    Ispec12 = zeros(1,length(CY));
    for iy = 1:size(CY,1)
        for ix = 1:2
            I = 0;
            CXY = eval(['CX',num2str(ix),'Y']);
            CX = eval(['CX',num2str(ix)]);
            pYX = eval(['pYX',num2str(ix)]);
            pX = eval(['pX',num2str(ix)]);
            N = eval(['N',num2str(ix)]);
            posy = find(ismember(CXY(:,end-Ny+1:end),CY(iy,:),'rows'));
            for k = 1:length(posy)
                [~,posx] = ismember(CXY(posy(k),1:N),CX,'rows');
                I = I + (pYX(posy(k))/pY(iy))*log2(pYX(posy(k))/(pY(iy)*pX(posx)));
            end   
            Ispec(ix,iy) = I;
        end
        I12 = 0;
        posy = find(ismember(CX1X2Y(:,end-Ny+1:end),CY(iy,:),'rows'));
        for k12 = 1:length(posy)
            [~,posx12] = ismember(CX1X2Y(posy(k12),1:(N1+N2)),CX1X2,'rows');
            I12 = I12 + (pYX1X2(posy(k12))/pY(iy))*log2(pYX1X2(posy(k12))/(pY(iy)*pX1X2(posx12)));
        end
        Ispec12(iy) = I12;
    end
    
    Ispec = [Ispec;Ispec12];

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
    
%     Ii12 = sum(pY.*Ispec(3,:));
%     R = sum(pY.*min(Ispec(1:2,:)));
%     U1 = Ii1 - R;
%     U2 = Ii2 - R;
%     S = Ii12 - Ii1 - Ii2 + R;
    
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