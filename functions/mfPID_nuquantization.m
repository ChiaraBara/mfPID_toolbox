%% quantizes non uniformly the series y considering bins of size n
% y: input series, column data
% n: bins size

function yq = mfPID_nuquantization(y,n)

N = size(y,1);
yq = zeros(N,1);    
b = N/n;

[~, sort_idx] = sort(y);
yq_sorted = zeros(N,1);
for ib = 1:b
    yq_sorted((ib-1)*n+1:ib*n) = ones(n,1)*ib;
end
for in = 1:N
    yq(in) = yq_sorted(sort_idx==in);
end

end

