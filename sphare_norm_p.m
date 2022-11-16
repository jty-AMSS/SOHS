% compute the SOHS of 1-real(x)^n-real(y)^n with relaxation order N=2000

clear;clc

N=2000;
et=@(x,y)1-x^10-y^10;
et3=@(x,y,z)1-x^10-y^10-z^10;
syms x1 x2


f=et(1/2*(x1+conj(x1)),(-1i)*1/2*(conj(x1)-x1));
f=expand(f);
g2= restraction(Finite_S_n(f,[x1]),[N]);

% % one can also try: 
% N=200;
% f=et3(1/2*(x1+conj(x1)),(-1i)*1/2*(conj(x1)-x1)*(1/2)*(x2+x2'),(-1i)*1/2*(conj(x1)-x1)*(-i/2)*(x2-x2'));
% f=expand(f);
% g2= restraction(Finite_S_n(f,[x1,x2]),[N,N]);


[Q,I,F,SOS,g]=L1RootSOS(g2);

SOHS=sym(0);
if ~all((I<N/4)+(I>3/4*N))
    error('Cannot Left it')
    return 
end
SF={};
for i=1:length(F)
    SF{i}=Finite_S_n(F{i});
    hi=sym(SF{i});
    SOHS=SOHS+hi*conj(hi);
end
SOHS=expand(SOHS);
%error:
t1=sym('t1','real');
t2=sym('t2','real');
err=subs(expand(SOHS-(f)),[x1,x2],[exp(1i*t1),exp(1i*t2)]);
err=simplify(err);
err=expand(err);
sum(children(err))-err
err=children(err);

err=simplify(err);
vpa(sum(err))

