% Compute the SOHS of Motzkin polynomial
clear
clc
M=@(x,y)x^4*y^2+x^2*y^4-3*x^2*y^2+1;
syms x1 x2
f=M((x1+conj(x1)),((x2+conj(x2))))
expand(f)
g2= restraction(Finite_S_n(f,[x1,x2]),[8,8])
[Q,I,F,SOS,g]=L1RootSOS(g2);
H=Finite_S_n(F{1});%Since the rank of Gram matrix is 1, otherwise one should compute the compute  \sum_{i=1}^length(F) |Finite_S_n(F{i})|^2, here M=|H|^2

t1=sym('t1','real');
t2=sym('t2','real');

% shows the error of SOHS:
err=subs(expand(sym(H)*sym(H)'-(f)),[x1,x2],[exp(1i*t1),exp(1i*t2)]);

err=simplify(err);
disp('error of SOHS:')
disp(vpa(err))