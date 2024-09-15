function out = Heaviside(X)

    pos = find(X > 0);
    out = zeros(length(X),1); out(pos)=1;

end