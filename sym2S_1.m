function [A1,A2,B]=sym2S_1(f,x)
%f \in  C[x,conj(x)],i.e. polyomial of x_i and conj(x_i)  , x =sym('x',[1,n])
x=x(:).';
n=length(x);
syms xt
f=expand(f);
y=sym('y',[1,n]);
f=subs(f,conj(x),y);
[B,T]=coeffs(f);
S=char([x,y]);
for i=1:length(T)
    tt=T(i);
    v=mupadmex(['degreevec(' char(tt), ',', S(9:end-2) ' )']);
    A(i,:)=double(v);
end
A1=A(:,1:n);
A2=A(:,n+(1:n));

end