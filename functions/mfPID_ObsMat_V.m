%% form observation matrix 
function B=mfPID_ObsMat_V(Y,V)

if isempty(V) %if no conditioning, give back the jth series of Y
    B=Y;
else
    [N,M]=size(Y);
    Nc=size(V,1); % number of candidates
    Lmax=max(V(:,2)); %maximum lag (across all signals)

    B=NaN*ones(N-Lmax,Nc);
    for n=Lmax+1:N
        for i=1:Nc %fill the i-th row of A
            B(n-Lmax,i)=Y(n-V(i,2),V(i,1)); 
        end
    end
end
