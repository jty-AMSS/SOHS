function [Q,F,SOS]=ComputeSOSByCVX(f,Index)
% Main function of computation of FSOS:
%input: f is an element in CZ class (functions on Z_n1*Z_n2*...*Z_nk)
%       Index is some elements in dual group of Z_n1*Z_n2*...*Z_nk, which
%is the index set of Gram matrix 
%output: 
% Q is the Gram matrix indexed by Index
% F is the terms of sum of squares (in CZ class)
% SOS= \sum_|F{i}|^2
Eq=containers.Map('KeyType',  'char', 'ValueType', 'any');
m=size(Index,1);
N=f.n;
for i=1:m
    for j=1:m
        t=mod(Index(i,:)-Index(j,:),N);
        if isKey(Eq,char(t))
            Eq(char(t))=[Eq(char(t)),m*(i-1)+j];
        else
            Eq(char(t))=[m*(i-1)+j];
        end
    end
end
G=keys(Eq);
suppf=find(f);
if ~isempty(setdiff(suppf,cell2mat(G(:))+0,'rows'))
    Q=nan;
    F=nan;
    SOS=nan;
    return
end
Ax=[];Ay=[];b=zeros(length(G),1);
for i=1:length(G)
    t=G{i};
    Ax(end+1:end+(length(Eq(t))))=i;
    Ay(end+1:end+(length(Eq(t))))=Eq(t);
    b(i)=f(t);
end
A=sparse(Ax,Ay,1,length(G),m*m);
cvx_begin
% variable  P(m,m) hermitian %hermitian_semidefinite
% maximize(lambda_min(P))
variable  P(m,m) hermitian semidefinite
maximize(-P(1,1))
A*P(:)==b;
cvx_end
Q=P;

if isinf(cvx_optval)
    F={};
    SOS=0;
    return
end
if isnan(cvx_optval)
    Q=nan;
    F={};
    SOS=0;
    return
end
disp('The compute is ended, now check the results')
if nargout>1
    F={};
    H=chol(sym(Q),'nocheck');
    for i=1:length(H)
        h=CZ(N);
        for j=1:length(H)
            h(Index(j,:))=double(H(i,j));
        end
        F{i}=h;
    end
end
if nargout>2
    SOS=CZ(N);
    for i=1:length(H)
        SOS=SOS+F{i}*F{i}';
    end
end
end