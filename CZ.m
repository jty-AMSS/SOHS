classdef CZ %
    % para:
    %   n=[n1,n2,...,nk]  double array means the group is
    % Z_n1*Z_n2*...*Z_nk 
    %   c is the Map<string , double> means the coefficient of element, by
    %  \widehat(f)(t)=c(char(t)) is the  coefficient of term
    %  prod(\chi_{t_i}(x(i))) where t is a double array [t_1,t_2,...,t_k] 
    properties
        n
        c
    end
    methods
        %% Constructor 
        function obj = CZ(n)
            obj.n=n;
            obj.c=containers.Map('KeyType',  'char', 'ValueType', 'any');
        end
        %% Get/Set
        function  x = get(f,t)
            %Get the Coefficient  of  term t
            if ~ischar(t)
                t=char(t);
            end
            if length(t(:))==length(f.n)
                t=t(:)';
                if isKey(f.c,t)
                    x=f.c(t);
                else
                    x=0;
                end
            end
        end
        
        function  f = set(f,t,x)
            t=double(t);
            t=mod(t,f.n);
            t=char(t);
            f.c(t)=x;
        end
        %% subsref
        function x = subsref(this,s)
            if strcmp(s.type,'{}')
                x=cell2mat(s.subs);
                x=PointValue(this,x);
                return
            end
            if strcmp(s.subs,'c')&&strcmp(s.type,'.')
                x=this.c;
                return
            end
            if strcmp(s.subs,'n')&&strcmp(s.type,'.')
                x=this.n;
                return
            end
            s=s.subs;
            t=s{1};
            if ~ischar(t)
                t=mod(t,this.n);
                t=char(t);
            end
            if isKey(this.c,t)
                x=get(this,t);
            else
                x=0;
            end
        end
        
        function this = subsasgn(this,s,b)
            s=s.subs;
            t=s{1};
            if ~ischar(t)
                t=char(t);
            end
            set(this,t,b);
        end
        %% plus
        function h=plus(f,g)
            if isnumeric(f)&&length(f)==1
                h=g;
                N=length(h.n);t=zeros(1,N);
                set(h,t,get(h,t)+f);
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(g,f);
                return
            end
            if ~(all(f.n==g.n))
                error('the addition are not on same group')
            end
            h=CZ(f.n);
            Index=[cell2mat(keys(f.c).');cell2mat(keys(g.c).')];
            Index=unique(Index,'rows');
            for i=1:(size(Index,1))
                t=Index(i,:);
                if get(f,t)+get(g,t)~=0
                    set(h,t,get(f,t)+get(g,t));
                end
            end
        end
        %% Sparsity
        function L=length(this)
            L=length(this.c);
        end
        function L=sparsity(this)
            L=length(this.c);
        end
        %% minus
        function h = minus(f,g)
            if isnumeric(f)&&length(f)==1
                h=minus(CZ(g.n),g)+f;
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(f,-g);
                return
            end
            if f.n~=g.n
                error('the addition are not on same group')
            end
            h=CZ(f.n);
            Index=[cell2mat(keys(f.c).');cell2mat(keys(g.c).')];
            Index=unique(Index,'rows');
            for i=1:(size(Index,1))
                t=Index(i,:);
                if (get(f,t)-get(g,t))~=0
                    set(h,t,get(f,t)-get(g,t));
                end
            end
        end
        %% mtimes
        function h = mtimes(f,g)
            if isnumeric(f)&&length(f)==1% is  scalar
                h=CZ(g.n);
                Index=cell2mat(keys(g.c).');
                for i=1:(size(Index,1))
                    t=Index(i,:);
                    set(h,t,f*get(g,t));
                end
                return
            end
            if isnumeric(g)&&length(g)==1% is  scalar
                h=g*f;
                return
            end
            N=f.n;
            h=CZ(N);
            Index1=cell2mat(keys(f.c).');
            Index2=cell2mat(keys(g.c).');
            for i=1:size(Index1,1)
                x=Index1(i,:);
                for j=1:size(Index2,1)
                    y=Index2(j,:);
                    t=mod(x+y,N);
                    set(h,t,get(h,t)+get(f,x)*get(g,y));
                end
            end
        end
        %%  	ctranspose(a)
        function h = ctranspose(f)
            N=f.n;
            h=CZ(N);
            Index1=cell2mat(keys(f.c).');
            for i=1:size(Index1,1)
                t=Index1(i,:);
                tp=char(N-t);
                set(h,tp,(get(f,t)));
            end
        end
        function h=conj(f)
            h=f';
        end
        %% times .*
        function h=times(f,g)
            h=f*g;
        end
        %% Display
        function disp(f)
            S=sym(f);
            disp(S)
        end
        %% transform to symbolic  object
        function S=sym(f)
            x=sym('x',[length(f.n),1]);
            S=sym(0);
            K=keys(f.c);
            for i=K
                ind=i{1};
                ind=double(ind);
                ind=mod(ind,f.n);
                y=get(f,ind);
                term=y*prod(x(:).^(ind(:)));
                S=S+term;
            end
        end
        
        %% Find
        function [subs,vals] = find(f)
            subs=keys(f.c);
            subs=cell2mat(subs.');
            subs=double(subs);
            vals=values(f.c);
            vals=cell2mat(vals.');
        end
        %% IFFT
        function V=CZifft(f)
            [A,B]=find(f);
            if length(f.n)>1
                F=zeros(f.n);
            else
                F=zeros(f.n,1);
                A=A+1;
                F(A+0)=B;
                V=ifft(F*length(F(:)));
                return 
            end
            N=f.n;
            for i=1:size(A,1)
                t=A(i,:);
                t=mod(t,N);
                t=t+1;
                t(t==0)=N(t==0);
                ndx=mysub2ind(f.n,t);
                F(ndx)=B(i);
            end
            V=ifftn(F*length(F(:)));
        end
        %% compute the value at point x
        function S=PointValue(this,x)
            x=x(:)';
            x=exp(2*pi*1i.*(x(:)')./(this.n));
            [a,b]=find(this);
            a=double(a);
            S=0;
            for i=1:length(b)
                t=prod(x.^[a(i,:)]);
                S=S+b(i)*t;
            end
        end
        %% Simplify
        function f=simplify(f,tol)
            %remove abs(coeffs)<tol
            if nargin==1
                tol=0;
            end
            I=keys(f.c);
            for i=1:length(I)
                if(abs(get(f,I{i}))<tol)
                    remove(f.c,I{i});
                end
            end
        end
        % shows minimal value
        function z=min_value(f)
            V=CZifft(f);
            z=min(V(:));
        end
    end
end



%% rewrite the sub2ind and ind2sub:
function ndx = mysub2ind(siz,v)
v1=v(1);
if length(v)>1
    v2=v(2);
    if length(v)>2
        varargin=cell(length(v)-2,1);
    end
end
for i=3:length(v)
varargin{i-2}=v(i);
end
siz = double(siz);
lensiz = length(siz);
if lensiz < 2
    error(message('MATLAB:sub2ind:InvalidSize'));
end

numOfIndInput = length(v);
if lensiz < numOfIndInput
    %Adjust for trailing singleton dimensions
    siz = [siz, ones(1,numOfIndInput-lensiz)];
elseif lensiz > numOfIndInput
    %Adjust for linear indexing on last element
    siz = [siz(1:numOfIndInput-1), prod(siz(numOfIndInput:end))];
end

if any(min(v1(:)) < 1) || any(max(v1(:)) > siz(1))
    %Verify subscripts are within range
    error(message('MATLAB:sub2ind:IndexOutOfRange'));
end

ndx = double(v1);
s = size(v1);
if numOfIndInput >= 2
    if ~isequal(s,size(v2))
        %Verify sizes of subscripts
        error(message('MATLAB:sub2ind:SubscriptVectorSize'));
    end
    if any(min(v2(:)) < 1) || any(max(v2(:)) > siz(2))
        %Verify subscripts are within range
        error(message('MATLAB:sub2ind:IndexOutOfRange'));
    end
    %Compute linear indices
    ndx = ndx + (double(v2) - 1).*siz(1);
end 
    
if numOfIndInput > 2
    %Compute linear indices
    k = cumprod(siz);
    for i = 3:numOfIndInput
        v = varargin{i-2};
        %%Input checking
        if ~isequal(s,size(v))
            %Verify sizes of subscripts
            error(message('MATLAB:sub2ind:SubscriptVectorSize'));
        end
        if (any(min(v(:)) < 1)) || (any(max(v(:)) > siz(i)))
            %Verify subscripts are within range
            error(message('MATLAB:sub2ind:IndexOutOfRange'));
        end
        ndx = ndx + (double(v)-1)*k(i-1);
    end
end
end



function v = myind2sub(siz,ndx)


nout = max(length(siz),1);
siz = double(siz);
lensiz = length(siz);

if lensiz < nout
    siz = [siz ones(1,nout-lensiz)];
elseif lensiz > nout
    siz = [siz(1:nout-1) prod(siz(nout:end))];
end

if nout > 2
    k = cumprod(siz);
    for i = nout:-1:3
        vi = rem(ndx-1, k(i-1)) + 1;
        vj = (ndx - vi)/k(i-1) + 1;
        varargout{i-2} = double(vj);
        ndx = vi;
    end
end

if nout >= 2
    vi = rem(ndx-1, siz(1)) + 1;
    v2 = double((ndx - vi)/siz(1) + 1);
    v1 = double(vi);
else
    v1 = double(ndx);
end
v=[v1,v2,cell2mat(varargout)];
end


