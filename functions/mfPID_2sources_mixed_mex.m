function out = mfPID_2sources_mixed_mex(B,j,i1,i2,k,warning)

    if nargin < 6, warning=0; end
    
    metric='maximum';
    
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
        M_X1X2i_tmp = M_X1X2(posy,:);
        
        %%% neighbor search in the lower density space for each yi
        atria_X1X2i = nn_prepare(M_X1X2i_tmp, metric);
        try
            [~, distances] = nn_search(M_X1X2i_tmp, atria_X1X2i, (1:Ni(iy))', k, 0);
            dd{iy}=distances(:,k);
            
            %%%%%%%%%%%%% m12, number of samples of [X1 X2] for all possible Y
            if ~isempty(M_X1X2)
                 atria_X1X2 = nn_prepare(M_X1X2, metric);
                 [count_X1X2, tmp] = range_search(M_X1X2, atria_X1X2, posy, dd{iy}, 0);
                 tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
                 for n=1:length(tmp)-1
                     count_X1X2(n)=max(k-1,count_X1X2(n)-sum(tmp{n}==dd{iy}(n)));
                 end
             else
                 count_X1X2=(Ni(iy)-1)*ones(Ni(iy),1);
             end 
            
            for ix = 1:2
                
                 M_X_tmp = eval(['M_X',num2str(ix)]);
                 M_Xi_tmp = M_X_tmp(posy,:);
                 
                 %%%%%%%%%%%%% m1 ed m2, number of samples of X1 and X2 for all possible Y
                 if ~isempty(M_X_tmp)
                     atria_X = nn_prepare(M_X_tmp, metric);
                     [count_X, tmp] = range_search(M_X_tmp, atria_X, posy, dd{iy}, 0);
                     tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
                     for n=1:length(tmp)-1
                         count_X(n)=max(k-1,count_X(n)-sum(tmp{n}==dd{iy}(n)));
                     end
                 else
                     count_X=(Ni(iy)-1)*ones(Ni(iy),1);
                 end  
                 
                 %%%%%%%%%%%%% m1i ed m2i, number of samples of X1 and X2 for Y=i
                 if ~isempty(M_Xi_tmp)
                     atria_Xi = nn_prepare(M_Xi_tmp, metric);
                     [count_Xi, tmp] = range_search(M_Xi_tmp, atria_Xi, (1:Ni(iy))', dd{iy}, 0);
                     tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
                     for n=1:length(tmp)-1
                         count_Xi(n)=max(k-1,count_Xi(n)-sum(tmp{n}==dd{iy}(n)));
                     end
                 else
                     count_Xi=(Ni(iy)-1)*ones(Ni(iy),1);
                 end
                          
                 eval(['count_X',num2str(ix),'{iy} = count_X;']);
                 eval(['count_X',num2str(ix),'i{iy} = count_Xi;']);
                 
                 clear count_X count_Xi 
                 
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