function SOS=SOSPHP(n)
%check the FSOS of PHP
N=ones(1,n)*(n-1);
SOS=CZ(N);
T=eye(n);
for k=1:(N(1)-1)
    h=CZ(N);
    for i=1:length(T)
        h(T(i,:)*k)=1;
    end
    disp(h)
    SOS=SOS+1/(2*N(1))*(h*h');
end
end