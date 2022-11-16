function f=sym2CZ(SoS,X,N)
% X: vector of variables
% N : size of Z_N
% convert f to CZ(N):
% 1/x, conj(x)-> N(x)-1;
X=X(:).';
SoS=expand(SoS);
SoS=subs(SoS, conj(X), X.^(N-1));
T=children(SoS).';%terms
f=CZ(N);
for i=1:length(T)
   c=coeffs(T(i));
   
   f(getdegree(T(i)/c,X))=f(getdegree(T(i)/c,X))+double(c);
end
end


function d=getdegree(T,X)
%T=x^n*y^m sym
n=length(X);
d=zeros(n,1);
for i=1:n
    U=ones(n,1);
    U(i)=2;
    d(i)= double( log2(subs(T,X(:),U)));
end
end