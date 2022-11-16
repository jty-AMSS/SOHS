function g=CZtoS_1(f)
N=f.n;
x=sym('x',[1,length(N)]);
[a,b]=find(f);
g=sym(0);
for i=1:length(b)
    v=sym(b(i));
    t=a(i,:);
    for k=1:length(N)
        if t(k)<N/2
            v=v*x(k)^t(k);
        end
        if t(k)>N/2
            t(k)=N-t(k);
            v=v*conj(x(k))^t(k);
        end
    end
    g=g+v;
end
end