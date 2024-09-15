%% quantizes the series y (separately for each column) with c quantization levels
% y: input series, column data
% c: number of quantization levels

function x = mfPID_quantization(y,c)

    [n,nc]=size(y);

    x=zeros(n,nc);
    for ic = 1:nc
        ma=max(y(:,ic)); mi=min(y(:,ic));
        q=(ma-mi)/c; % amplitude of quantization level

        l=zeros(c,1);
        for i=1:c %quantization levels
           l(i)=mi+i*q;
        end

        for i=1:n
           j=1;
           while (y(i,ic)>=l(j))&&(j<c)
              j=j+1;
           end
           x(i,ic)=j;
        end
    end

end

