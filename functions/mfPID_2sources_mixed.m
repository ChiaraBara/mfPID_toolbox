function out = mfPID_2sources_mixed(B,j,i1,i2,k,warning)

    if nargin < 6, warning=0; end
    
    metric='chebychev';
    
    M_Y = B(:,j);
    M_X1 = B(:,i1);
    M_X2 = B(:,i2);
        
    M_X1X2 = [M_X1 M_X2];
    M_X1X2Y = [M_X1 M_X2 M_Y];
    
    Nall = size(M_Y,1);
    Ny = size(M_Y,2); Ni1=size(M_X1,2); Ni2=size(M_X2,2);   
    
    %%% probability p(y)
    [CY,~,jc] = unique(M_Y,'rows');
    NpY=size(CY,1); 
    for cnt=1:NpY
        Ni(cnt)=sum(jc==cnt);
        pY(cnt)=Ni(cnt)/Nall;
    end
    
    error = zeros(1,size(CY,1));
    for iy = 1:size(CY,1)
        
        posy = ismember(M_X1X2Y(:,end-Ny+1:end),CY(iy,:),'rows'); posy = find(posy==1);
        M_X1X2i = M_X1X2(posy,:);
        
        try
            %%% neighbor search in the lower density space for each yi
            [~, distances] =  knnsearch(M_X1X2i,M_X1X2i,'K',k+1,'Distance',metric);
            dd{iy} = distances(:,end);
        
            %%%%%%%%%%%%% m12, number of samples of [X1 X2] for all possible Y
            [~, distance_X1X2] =  knnsearch(M_X1X2,M_X1X2,'K',Nall,'Distance',metric);
            count_X1X2 = distance_X1X2(posy,2:end) < dd{iy};
            count_X1X2 = max(k-1, sum(count_X1X2,2)); 
        
            for ix = 1:2
    
                M_X_tmp = eval(['M_X',num2str(ix)]);
                M_Xi_tmp = M_X_tmp(posy,:);
            
                %%%%%%%%%%%%% m1 ed m2, number of samples of X1 and X2 for all possible Y
                [~, distance_X] =  knnsearch(M_X_tmp,M_X_tmp,'K',Nall,'Distance',metric);
                count_X = distance_X(posy,2:end) < dd{iy};
                count_X = max(k-1, sum(count_X,2)); 
    
                %%%%%%%%%%%%% m1i ed m2i, number of samples of X1 and X2 for Y=i
                [~, distance_Xi] =  knnsearch(M_Xi_tmp,M_Xi_tmp,'K',Ni(iy),'Distance',metric);
                count_Xi = distance_Xi(:,2:end) < dd{iy};
                count_Xi = max(k-1, sum(count_Xi,2)); 
                
                eval(['count_X',num2str(ix),'{iy} = count_X;']);
                eval(['count_X',num2str(ix),'i{iy} = count_Xi;']);
                
                clear count_X count_Xi count_X_boot
                
            end       
            eval('count_X12{iy} = count_X1X2;');       
        catch
            error(iy) = 1;
            if warning == 1
                disp(['Warning! Not enough samples to counting k=',num2str(k),' neighbors when y=[',num2str(CY(iy,:)),'].'])
            end
        end
    end
    
    %%%%%%%%%%%%% computation of Ispec
    Ispec = zeros(2,length(CY));
    Ispec12 = zeros(1,length(CY));
    for iy = 1:length(CY)
        if error(iy) == 0
            dd2=2*dd{iy}; dd2(dd2==0)=[];
            for ix = 1:2
                eval(['mall = count_X',num2str(ix),'{iy};']);
                eval(['mi = count_X',num2str(ix),'i{iy};']);
                eval(['Nj = Ni',num2str(ix),';']);
                Ispec(ix,iy) = (1/Ni(iy))*sum(psi(mi+1)-psi(mall+1)) - psi(Ni(iy)) + psi(Nall);
            end
            m12 = count_X12{iy};
            Ispec12(iy) = psi(k) - (1/Ni(iy))*sum(psi(m12+1)) - psi(Ni(iy)) + psi(Nall);
        else
            Ispec(:,iy)=nan(2,1);
            Ispec12(iy) = nan;
        end
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