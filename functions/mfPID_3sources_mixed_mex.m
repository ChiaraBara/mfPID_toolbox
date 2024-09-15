function out = mfPID_3sources_mixed_mex(B,j,i1,i2,i3,k,warning)

    if nargin < 7, warning=0; end
    
    metric='maximum';
    
    M_Y = B(:,j);
    M_X1 = B(:,i1);
    M_X2 = B(:,i2);
    M_X3 = B(:,i3);
        
    M_X1X2 = [M_X1 M_X2];
    M_X1X3 = [M_X1 M_X3];
    M_X2X3 = [M_X2 M_X3];
    M_X1X2X3 = [M_X1 M_X2 M_X3];
    M_X1X2X3Y = [M_X1 M_X2 M_X3 M_Y];
    
    Nall = size(M_Y,1);
    Ny = size(M_Y,2); Ni1=size(M_X1,2); Ni2=size(M_X2,2); Ni3=size(M_X3,2);   
    
    comb = {1,2,3,[1 2],[1 3],[2 3]};
    Ncomb = [Ni1, Ni2, Ni3, Ni1+Ni2, Ni1+Ni3, Ni2+Ni3];
    
    %%% probability p(y)
    [CY,~,jc] = unique(M_Y,'rows');
    NpY=size(CY,1); 
    for cnt=1:NpY
        Ni(cnt)=sum(jc==cnt);
        pY(cnt)=Ni(cnt)/Nall;
    end

    error = zeros(1,size(CY,1));    
    for iy = 1:size(CY,1)
        
        posy = ismember(M_X1X2X3Y(:,end-Ny+1:end),CY(iy,:),'rows'); posy = find(posy==1);
        M_X1X2X3i_tmp = M_X1X2X3(posy,:);
        
        %%% neighbor search in the lower density space for each yi
        atria_X1X2X3i = nn_prepare(M_X1X2X3i_tmp, metric);

        try 
            [~, distances] = nn_search(M_X1X2X3i_tmp, atria_X1X2X3i, (1:Ni(iy))', k, 0);
            dd{iy}=distances(:,k);
            
            %%%%%%%%%%%%% m123, number of samples of [X1 X2 X3] for all possible Y
            if ~isempty(M_X1X2X3)
                 atria_X1X2X3 = nn_prepare(M_X1X2X3, metric);
                 [count_X1X2X3, tmp] = range_search(M_X1X2X3, atria_X1X2X3, posy, dd{iy}, 0);
                 tmp=tmp(:,2);%subtraction from count of points with distance exactly equal to k-th neighbor
                 for n=1:length(tmp)-1
                     count_X1X2X3(n)=max(k-1,count_X1X2X3(n)-sum(tmp{n}==dd{iy}(n)));
                 end
             else
                 count_X1X2X3=(Ni(iy)-1)*ones(Ni(iy),1);
             end 
            
            for icomb = 1:size(comb,2)
                
                 name = [];
                 for in = 1:length(comb{icomb})
                     name = [name,'X',num2str(comb{icomb}(in))];
                 end 
                 
                 M_X_tmp = eval(['M_',name]);
                 M_Xi_tmp = M_X_tmp(posy,:);
                 
                 %%%%%%%%%%%%% m, number of samples of M for all possible Y
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
                 
                 %%%%%%%%%%%%% mi, number of samples of M for Y=i
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
                          
                 eval(['count_',name,'{iy} = count_X;']);
                 eval(['count_',name,'i{iy} = count_Xi;']);
                 
                 clear count_X count_Xi 
                 
            end    
            
            eval('count_X123{iy} = count_X1X2X3;');

        catch
            error(iy) = 1;
            if warning == 1
                disp(['Warning! Not enough samples to counting k=',num2str(k),' neighbors when y=[',num2str(CY(iy,:)),'].'])
            end
        end
    end

    %%%%%%%%%%%%% computation of Ispec
    Ispec = zeros(size(comb,2),length(CY));
    Ispec123 = zeros(1,length(CY));
    for iy = 1:length(CY)        
        if error(iy) == 0
            for icomb = 1:size(comb,2)          
                name = [];
                for in = 1:length(comb{icomb})
                    name = [name,'X',num2str(comb{icomb}(in))];
                end 
                eval(['mall = count_',name,'{iy};']);
                eval(['mi = count_',name,'i{iy};']);            
                Ispec(icomb,iy) = (1/Ni(iy))*sum(psi(mi+1)-psi(mall+1)) - psi(Ni(iy)) + psi(Nall);
            end
            m123 = count_X123{iy};
            Ispec123(1,iy) = psi(k) - (1/Ni(iy))*sum(psi(m123+1)) - psi(Ni(iy)) + psi(Nall); 
        else
            Ispec(:,iy)=nan(size(comb,2),1);
            Ispec123(1,iy) = nan;  
        end
    end  
    
    Ispec = [Ispec; Ispec123];
 
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