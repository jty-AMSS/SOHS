function [Q,Index,SF,SOHS,err]=Fun_SOHS(F,N)
%compute the FSOS by the guidence of square root of F
%input: F is an function handle of a polynomial with n inputs (e.g F=@(x) 1-1/2*x(1)-1/2*conj(x(1)) )
% N=[N_1,N_2,...,N_n] is the relaxation order of f 
%output: 
% Q is the Gram matrix indexed by Index
% Index is the index set of Q
% SF is the terms of sum of squares (symbolic object)  with x_i \in S_1 
% SOHS= \sum_|F{i}|^2
% g is the square root of f (in CZ class)
%---------------------------------------------
% example 1: verify 1-x^10-y^10-z^10>=0 for (x^2+y^2+z^2=1) with  relaxation order 50 
%  
% et3=@(x,y,z)1-x^10-y^10-z^10;
% F=@(x)et3(1/2*(x(1)+conj(x(1))),(-1i)*1/2*(conj(x(1))-x(1))*(1/2)*(x(2)+x(2)'),(-1i)*1/2*(conj(x(1))-x(1))*(-i/2)*(x(2)-x(2)'))
% [Q,Index,SF,SOHS,err]=Fun_SOHS(F,[200,200])

% example 2: verify Motzkin polynomial >=0 on [-2,2]*[-2,2]
% M=@(x,y)x^4*y^2+x^2*y^4-3*x^2*y^2+1;
% F=@(x)M(x(1)+conj(x(1)),x(2)+conj(x(2)))
% [Q,Index,SF,SOHS,err]=Fun_SOHS(F,[8,8])

FunF=F;
tol=1e-7;
tolFFt=0;
n=length(N);
N=N(:)';
x=sym('x', [1,n]);
f=F(x);
f=expand(f);
g2= restraction(Finite_S_n(f,x),N);
V=CZifft(g2);
if min(V(:))<-1e-10
    error('~ (f>0)')
end
V=abs(V).^0.5;
g=CZfft(V,tolFFt); %g=sqrt(f)
[a,b]=find(g);
[~, sortb_Index]=sort(abs(b),'descend');
flag=ones(size(sortb_Index,1),1);
for i=1:size(sortb_Index,1)
    t=or( a(sortb_Index(i,:))<=N/4,a(sortb_Index(i,:))>=3*N/4);
    if ~all(t)
        flag(i)=0;
    end
end
sortb_Index=sortb_Index(flag>0);
k=ceil(sqrt(length(f)));
Q=-1;
while max(max(isnan(Q)))||(lambda_min(Q)<-tol)
    Index= a(sortb_Index(1:k),:);
    disp(k)
    if nargout<3
        Q=ComputeSOSByCVX(g2,Index);
    end
    if nargout>=3
        [Q,F]=ComputeSOSByCVX(g2,Index);
    end
    k=min(k+1,length(b)); 
    if k>length(sortb_Index)
        Q=nan;
        Index=nan;
        SF={nan};
        SOHS=nan;
        err=nan;
        warning('order is not enough')
        return 
    end
end

%% SOHS 
if nargout>3
SOHS=sym(0);
for i=1:length(F)
    SF{i}=Finite_S_n(F{i});
    hi=sym(SF{i});
    SOHS=SOHS+hi*conj(hi);
end
SOHS=expand(SOHS);
end
%% error
if nargout>4
t=sym('t',[1,n],'real');
err=expand(SOHS-FunF(x));
err=subs(err,x,exp(1i*[t]));
err=simplify(err);
end
end