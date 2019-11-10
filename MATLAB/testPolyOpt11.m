clear
clc

sdpvar x1 x2;%variable

x=[x1;x2];


f = x1^2*x2^2*(x1^2+x2^2-1);%objective function

%inequalities constraints
g=[]; 


%equalities constraints
h=[];


epsilon=1e-5; %the small parameter

k=2; %the order of relaxation 

PolyOpt11(x,f,g,h,epsilon,k)