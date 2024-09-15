%%%%%% it works only for consecutive lags

function B = mfPID_ObsMat_lags(Y,lags)

    [N,M]=size(Y);
    Nc=length(lags); % number of candidates

    Lfuture=find(lags>0); %maximum lag in the future 
    Lmax = length(Lfuture);
    Lpast=find(lags<0); %maximum lag in the past
    Lmin = length(Lpast);
    i0 = find(lags==0);

    for n=1:N-Nc+1 %Lmin+1:N-(Lmax+i0)
        for ips=1:Lmin %fill the i-th row of A
            B(n,ips)=Y((n-1)+ips); 
        end
        if i0 == 1
            B(n,Lmin+i0)=Y((n-1)+Lmin+i0);
        end
        for ift=1:Lmax
            B(n,Lmin+i0+ift)=Y((n-1)+Lmin+i0+ift);
        end
    end

end
