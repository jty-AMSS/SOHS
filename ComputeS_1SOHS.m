% Generate a random  SOHS on S^1 and compute it by FSOS

%% Generate random  sparse SOHS with support {x^B(1),x^B(2),x^B(3), conj(x)^A(1),conj(x)^A(2),..,conj(x)^A(6)}
syms x x1
A=randi(60,6,1)
B=randi(60,3,1)
SoS=sym(0);
for i=1:10
F{i}=sym(0);
f=sym(0);
for k=1:length(B)
f=f+(randn()+1i*randn())*x^B(k);
end
for k=1:length(A)
f=f+(randn()+1i*randn())*conj(x)^A(k);
end
G{i}=f;
SoS=SoS+f*conj(f);
end
SoS=expand(SoS);
f=sym2CZ(SoS,x,300);

%% Compute the SOHS
[Q,I,F]=L1RootSOS(f);
SoS2=sym(0);
for i=1:length(F)
h=CZtoS_1(F{i});
SoS2=SoS2+h*conj(h);
end
SoS2=expand(SoS2);
SoS2=subs(SoS2,x1,x);
vpa(subs(SoS2-SoS,conj(x),1/x))

