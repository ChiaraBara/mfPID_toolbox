function out = mfPID_3sources_discrete(B,j,i1,i2,i3)

    Y = B(:,j);
    X1 = B(:,i1);
    X2 = B(:,i2);
    X3 = B(:,i3);

    Np=size(Y,1); %total number of patterns
    Ny=size(Y,2); N1=size(X1,2); N2=size(X2,2); N3=size(X3,2);
    base = 2;
    
    comb = {1,2,3,[1 2],[1 3],[2 3],[1 2 3]};
    
    [~,pY,CY] = mfPID_H(Y,base);
    [~,pX1,CX1] = mfPID_H(X1,base);
    [~,pX2,CX2] = mfPID_H(X2,base);
    [~,pX3,CX3] = mfPID_H(X3,base);
    [~,pYX1,CX1Y] = mfPID_H([X1 Y],base);
    [~,pYX2,CX2Y] = mfPID_H([X2 Y],base);
    [~,pYX3,CX3Y] = mfPID_H([X3 Y],base);
    [~,pX1X2,CX1X2] = mfPID_H([X1 X2],base);
    [~,pX1X3,CX1X3] = mfPID_H([X1 X3],base);
    [~,pX2X3,CX2X3] = mfPID_H([X2 X3],base);
    [~,pX1X2X3,CX1X2X3] = mfPID_H([X1 X2 X3],base);
    [~,pYX1X2,CX1X2Y] = mfPID_H([X1 X2 Y],base);
    [~,pYX1X3,CX1X3Y] = mfPID_H([X1 X3 Y],base);
    [~,pYX2X3,CX2X3Y] = mfPID_H([X2 X3 Y],base);
    [~,pYX1X2X3,CX1X2X3Y] = mfPID_H([X1 X2 X3 Y],base);
    
    pX = {pX1,pX2,pX3,pX1X2,pX1X3,pX2X3,pX1X2X3};
    cX = {CX1,CX2,CX3,CX1X2,CX1X3,CX2X3,CX1X2X3};
    pYX = {pYX1,pYX2,pYX3,pYX1X2,pYX1X3,pYX2X3,pYX1X2X3};
    cXY = {CX1Y,CX2Y,CX3Y,CX1X2Y,CX1X3Y,CX2X3Y,CX1X2X3Y};
    
    Ispec = zeros(size(comb,2),length(CY));
    for iy = 1:length(CY)
        for icomb = 1:size(comb,2)
            I = 0;
            name = [];
            for in = 1:length(comb{icomb})
                name = [name,'X',num2str(comb{icomb}(in))];
            end        
            CXY = eval(['C',name,'Y']);
            CX = eval(['C',name]);
            pYX = eval(['pY',name]);
            pX = eval(['p',name]);
            N = size(CX,2);
            posy = find(CXY(:,end-Ny+1:end) == CY(iy));
            for k = 1:length(posy)
                [~,posx] = ismember(CXY(posy(k),1:N),CX,'rows');
                I = I + (pYX(posy(k))/pY(iy))*log2(pYX(posy(k))/(pY(iy)*pX(posx)));
            end   
            Ispec(icomb,iy) = I;
        end
    end
    
    [Re,D,Respec,Dspec] = mfPID_redundancy_lattice_3sources(Ispec,pY);

    %%% PID TERMS
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