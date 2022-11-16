% MAX_SAT Example
clear
clc
x=sym('x',[10,1]);
for i=1:10
    eval(['x' num2str(i) '=x(i)'])
end
f=234*x3 - 1386*x2 - 1389*x1 + 502*x4 + 3056*x5 - 4692*x6 - 2142*x7 - 1312*x8 - 4645*x9 + 3787*x10 - 3399*x1*x2 - 1140*x1*x3 - 491*x1*x4 - 282*x2*x3 - 2413*x1*x5 - 884*x2*x4 - 2212*x1*x6 + 3457*x2*x5 + 4462*x3*x4 + 4097*x1*x7 + 1707*x2*x6 - 936*x3*x5 + 3419*x1*x8 - 4102*x2*x7 - 976*x3*x6 - 2403*x4*x5 - 1245*x1*x9 - 3786*x2*x8 + 1014*x3*x7 + 3139*x4*x6 + 483*x1*x10 + 4417*x2*x9 - 854*x3*x8 + 6*x4*x7 - 2037*x5*x6 - 1678*x2*x10 + 2057*x3*x9 - 981*x4*x8 + 4848*x5*x7 + 4085*x3*x10 + 1129*x4*x9 - 4936*x5*x8 - 1122*x6*x7 - 2628*x4*x10 + 2787*x5*x9 + 667*x6*x8 - 2002*x5*x10 + 640*x6*x9 + 1874*x7*x8 - 707*x6*x10 + 778*x7*x9 + 3813*x7*x10 - 2764*x8*x9 + 3038*x8*x10 + 2170*x9*x10 + 50450;
h=Finite_S_n(f,x);
N=ones(1,10)*2;
g=restraction(h,N);



% If basis is M_ap
n=10;
I=nchoosek(1:n,2);
subs=[];
for i=1:size(I,1)
subs(i,:)=zeros(1,n);
subs(i,I(i,:))=1;
end
subs=[subs;eye(n);zeros(1,n)];
[Q2]=ComputeSOSByCVX(g,subs);
disp('M_ap cannot compute the FSOS')
pause(1)
% If basis is choosed by suqare root:
[Q,Index,F,SOS]=L1RootSOS(g)
disp('error of FSOS where basis is choosed by suqare root :')
disp(vpa(SOS-g))