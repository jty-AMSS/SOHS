function [Q,Index,F,SOS,g]=L1RootSOS(f,tol,tolFFt)
%compute the FSOS by the guidence of square root of f
%input: f is an element in CZ class (functions on Z_n1*Z_n2*...*Z_nk)
%output: 
% Q is the Gram matrix indexed by Index
% Index is the index set of Q
% F is the terms of sum of squares (in CZ class)
% SOS= \sum_|F{i}|^2
% g is the square root of f (in CZ class)
if nargin==1
    tol=0;
    tolFFt=0;
end
if nargin==2
    tolFFt=0;
end
V=CZifft(f);
if min(V(:))<-1e-10
    error('~ (f>0)')
end
V=abs(V).^0.5;
g=CZfft(V,tolFFt);
[a,b]=find(g);
[~, sortb_Index]=sort(abs(b),'descend');
k=ceil(sqrt(length(f)));
Q=-1;
while max(max(isnan(Q)))||(lambda_min(Q)<-tol)
    Index= a(sortb_Index(1:k),:);
    disp(k)
    if nargout<3
        Q=ComputeSOSByCVX(f,Index);
    end
    if nargout==3
        [Q,F]=ComputeSOSByCVX(f,Index);
    end
    if nargout>3
        [Q,F,SOS]=ComputeSOSByCVX(f,Index);
    end
    k=min(k+1,length(b));
end
end