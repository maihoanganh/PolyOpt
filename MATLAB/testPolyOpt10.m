clear
clc

sdpvar x1 x2;%variable

x=[x1;x2];


f = x1^3-x2^2;%objective function

%inequalities constraints
g=[x1;x2]; 


%equalities constraints
h=[(x1*x2+1)*(x1-x2)^2];


epsilon=10^(-5); %the small parameter

k=1; %the order of relaxation 

PolyOpt10(x,f,g,h,epsilon,k)