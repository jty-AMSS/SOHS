classdef Finite_S_n
    % 结构：n是变量个数
    % T:term
    % c: coeffs
    properties
        A1
        A2
        B
    end
    methods
        %% 构造函数
        function obj = Finite_S_n(f,x)
            %% Sym->Sn
            if isa(f,'sym')
                [A1,A2,B]=sym2S_1(f,x);
                obj.A1=A1;
                obj.A2=A2;
                obj.B=B;
                return 
            end
            %% CZ->Sn
            N=f.n;
            [supp,coeff]=find(f);
            for i=1:length(coeff)
                t=supp(i,:);
                tn=N.*(t>N/2)-t.*(t>N/2);
                tp=t.*(t<N/2);
                A1(i,:)=tp;
                A2(i,:)=tn;
                B(i)=coeff(i);
                
            end
            obj.A1=A1;
            obj.A2=A2;
            obj.B=B;
        end
        %% sym disp
        function h=sym(this)
            A1=this.A1;
            A2=this.A2;
            B=this.B;
            x=size(this.A1,2);
            x=sym('x',[1,x]);
            h=sym(0);
            for i=1:length(this.B)
                h=h+prod(x.^A1(i,:))*prod(conj(x).^A2(i,:))*B(i);
            end
        end
        function disp(this)
            disp(sym(this))
        end
        %%  Restraction
        function f=restraction(this,N)
            % this ->C[Z_N]
            f=CZ(N);
            A1=this.A1;
            A2=this.A2;
            B=this.B;
            for i=1:length(B)
                t=mod(A1(i,:)-A2(i,:),N);
                f(t)=f(t)+double(B(i));
            end
        end
    end
end

function i=findterm(A,t)
% like findrow
t=t(:)';
A=(A==t);
i=find(all(A'));
end