function PolyOpt10(x,f,g,h,epsilon,k)

while 1>0

format long

% Original problem f^* = min{f(x); x in S(g,h)} where S(g,h) is semialgebraic 
% set of inequalities g={g_1,...,g_m} and equalities h={h_1,...,h_l}.
% 
% Solve hierarchy sdp:
%   rho_k(eps) = sup{lambda: theta^k*(f + eps*theta^d - lambda) in   Q(g)_{d+k}}
% where  Q(g,h)_r is truncated quaratic module of index r associated with S(g,h)


tic


n=length(x);

m=length(g);

l=length(h);




%degree
w=[];
for j=1:m
    w=[w;ceil(.5*degree(g(j)))];
end


%degree
u=[];
for j=1:l
    u=[u;ceil(.5*degree(h(j)))];
end

if m==0 && l==0
    d=ceil(.5*degree(f));%unconstraint
else    
    d=floor(.5*degree(f))+1;%constraint
end

lw=[];
for j=1:m
    lw=[lw;nchoosek(n+d+k-w(j),n)];
end

lu=[];
for j=1:l
    lu=[lu;nchoosek(n+2*d+2*k-2*u(j),n)];
end



v = monolist(x,k+d);%vector of monomial

v2=monolist(x,2*d+2*k);

G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

theta=1+x'*x;%quadratic polynomial in (*)

%building the element in quadratic module
quadratic_module=v'*G0*v;
for j=1:m
    G{j}=sdpvar(lw(j),lw(j));
    quadratic_module=quadratic_module+v(1:lw(j))'*G{j}*v(1:lw(j))*g(j);
end

for j=1:l
    c{j}=sdpvar(lu(j),1);
    quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
end

lambda=sdpvar(1);%vaiable of sdp (*)

%identify coefficient of two polynomials
F = [coefficients(theta^k*(f-lambda+epsilon*theta^d)-quadratic_module,x) == 0,G0>=0];
for j=1:m
    F=[F,G{j}>=0];
end

ops = sdpsettings('solver','mosek','verbose',0);

diagnostics=optimize(F,-lambda,ops);

opt_val=double(lambda);%display optimal value

fprintf('optimal value              problem       \n');
fprintf('%15.15f               %0d              \n',opt_val,diagnostics.problem);
fprintf('\n');

if diagnostics.problem~=0 && diagnostics.problem~=4
    fprintf('Increase order k or increase parameter epsilon!\n');
    break
end

toc








% Extraction of global optimizers by using ASEM


%small parameter for omega0_upper
epsilon=1e-2;


%to check zero eigenvalues
tol=1e-2;

% for pivot
TOL=1e-3;

%inequalities 
g=[g;opt_val-f];m=length(g);

w=[w;ceil(.5*degree(g(m)))];%degree

lw=[];
for j=1:m
    lw=[lw;nchoosek(n+2+k-w(j),n)];
end

lu=[];
for j=1:l
    lu=[lu;nchoosek(n+4+2*k-2*u(j),n)];
end

%vector of monomial
v=monolist(x,2+k);

v2=monolist(x,4+2*k);






%finding upper bound of radii r0


%Gram matrix of sos
G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

%building the element in quadratic module
quadratic_module=v'*G0*v;
for j=1:m
  Q{j}=sdpvar(lw(j),lw(j));
  quadratic_module=quadratic_module+v(1:lw(j))'*Q{j}*v(1:lw(j))*g(j);
end

for j=1:l
  c{j}=sdpvar(lu(j),1);
  quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
end

theta=1+x'*x;
lambda=sdpvar(1);

F=[coefficients(theta^k*(x'*x-lambda+epsilon*theta^2)-quadratic_module,x)==0,G0>=0];

for j=1:m
    F=[F,Q{j}>=0];
end



diagnostics=optimize(F,-lambda,ops);

lambda=double(lambda);

omega0_upper=lambda;

  
fprintf('omega0_upper          problem       \n');
fprintf('%0f                 %0d\n',omega0_upper,diagnostics.problem);
fprintf('\n');

if diagnostics.problem~=0 && diagnostics.problem~=4
    fprintf('Increase order k or increase parameter epsilon!\n');
    break
end
    
    
%centers of sequence of balls 
a0=zeros(n,1);


k_min=max([w;u;1]);

lw=[];
for j=1:m
    lw=[lw;nchoosek(n+k_min+k-w(j),n)];
end

lu=[];
for j=1:l
    lu=[lu;nchoosek(n+2*k_min+2*k-2*u(j),n)];
end
lu_sphere=nchoosek(n+2*k_min+2*k-2,n);
%vector of monomial
v=monolist(x,k_min+k);

v2=monolist(x,2*k_min+2*k);

u0=monolist(x,k+k_min-1);

%finding radii r0


%Gram matrix of sos
G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

Q0=sdpvar(length(u0));
 
lambda=sdpvar(1);

%building the element in quadratic module
quadratic_module=v'*G0*v+(omega0_upper-(x-a0)'*(x-a0))*u0'*Q0*u0;
for j=1:m
    Q{j}=sdpvar(lw(j),lw(j));
    quadratic_module=quadratic_module+v(1:lw(j))'*Q{j}*v(1:lw(j))*g(j);
end

for j=1:l
    c{j}=sdpvar(lu(j),1);
    quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
end



F=[coefficients((x-a0)'*(x-a0)-lambda-quadratic_module,x)==0,G0>=0,Q0>=0];

for j=1:m
    F=[F,Q{j}>=0];
end

  
  
diagnostics=optimize(F,-lambda,ops);

lambda=double(lambda);

omega0=lambda;

M = dual(F(2));

dim_subM=nchoosek(n+k,k);

subM=M(1:dim_subM,1:dim_subM);


rank_M=rank(M,tol);

rank_subM=rank(subM,tol);

fprintf('omega0          problem       rank of moment matrix      rank of moment submatrix\n');
fprintf('%0f         %0d              %0d                        %0d\n',omega0,diagnostics.problem,rank_M,rank_subM);
fprintf('\n');

if diagnostics.problem~=0 && diagnostics.problem~=4
    fprintf('Increase order k or increase parameter epsilon!\n');
    break
end

if rank_subM==rank_M
    [V,~] = eig(double(G0));
    r=rank_M;
    V=V(:,1:r);
    V=V';
    [U,pivot] = rref(V,TOL);
    U=U';
    % Figure out multiplying matrices using YALMIP code
    w = v(pivot);
    for i = 1:n
        xw = x(i)*w;
        k = [];
        for j = 1:length(xw)
            k = [k;find(ismember(xw(j),v))];           
        end
        N{i} = U(k,:);
    end
    
    

    % Create random convex combination
    rands = rand(n,1);rands = rands/sum(rands);
    M = 0;
    for i = 1:n
        M = M + rands(i)*N{i};
    end

    [L,T] = schur(M);
    % Extract solution
    for i = 1:r
        solution=zeros(n,1);
        for j = 1:n
            solution(j) =  L(:,i)'*N{j}*L(:,i);
        end
        
        solution=solution
        if m~=0 
            check_ineq=replace(g,x, solution)
        end

        if l~=0
            check_eq=replace(h,x, solution)    
        end
        toc
    end
    break
end
   













%centers of sequence of balls 
a=3*sqrt(omega0)*eye(n);

%find radius
omega=zeros(n,1);
for t=1:n
  
    %Gram matrix of sos
    G0=sdpvar(length(v),length(v));%vaiable of sdp (*)

    Q0=sdpvar(lu_sphere,1);



    lambda=sdpvar(1);

    quadratic_module=v'*G0*v+(omega0-(x-a0)'*(x-a0))*v2(1:lu_sphere)'*Q0;
    for j=1:m
        Q{j}=sdpvar(lw(j),lw(j));
        quadratic_module=quadratic_module+v(1:lw(j))'*Q{j}*v(1:lw(j))*g(j);
    end

    for j=1:l
        c{j}=sdpvar(lu(j),1);
        quadratic_module=quadratic_module+v2(1:lu(j))'*c{j}*h(j);
    end



    if t>1
        for j=1:t-1
           H{j}=sdpvar(lu_sphere,1);
           quadratic_module=quadratic_module+v2(1:lu_sphere)'*H{j}*(omega(j)-(x-a(:,j))'*(x-a(:,j)));
        end
    end

    F=[coefficients((x-a(:,t))'*(x-a(:,t))-lambda-quadratic_module,x)==0,G0>=0];
    for j=1:m
        F=[F,Q{j}>=0];
    end


    diagnostics=optimize(F,-lambda,ops);

    lambda=double(lambda);

    omega(t)=lambda;

    M = dual(F(2));

    dim_subM=nchoosek(n+k,k);

    subM=M(1:dim_subM,1:dim_subM);


    rank_M=rank(M,tol);

    rank_subM=rank(subM,tol);

    fprintf('omega%0d         problem       rank of moment matrix      rank of moment submatrix\n',t);
    fprintf('%0f         %0d              %0d                        %0d\n',omega(t),diagnostics.problem,rank_M,rank_subM);
    fprintf('\n');
    
    if diagnostics.problem~=0 && diagnostics.problem~=4
        fprintf('Increase order k or increase parameter epsilon!\n');
        break
    end

    if rank_subM==rank_M
        [V,~] = eig(double(G0));
        r=rank_M;
        V=V(:,1:r);
        V=V';
        [U,pivot] = rref(V,TOL);
        U=U';
        % Figure out multiplying matrices using YALMIP code
        w = v(pivot);
        for i = 1:n
            xw = x(i)*w;
            k = [];
            for j = 1:length(xw)
                k = [k;find(ismember(xw(j),v))];     
            end
            N{i} = U(k,:);
        end



        % Create random convex combination
        rands = rand(n,1);rands = rands/sum(rands);
        M = 0;
        for i = 1:n
            M = M + rands(i)*N{i};
        end

        [L,T] = schur(M);
        % Extract solution
        solution=zeros(n,1);
        for i = 1:r
            for j = 1:n
                solution(j) =  L(:,i)'*N{j}*L(:,i);
            end
            
            solution=solution
            if m~=0 
                check_ineq=replace(g,x, solution)
            end

            if l~=0
                check_eq=replace(h,x, solution)   
            end
            toc
        end
        break
    end
   

end





 
    
fprintf('Computing approximate root from sphere equations...\n');

%solve linear programming
A=[];
b=[];
for j=1:n
    A=[A;-(a(:,j)-a0)'];
    b=[b;omega(j)-omega0-norm(a(:,j))^2+norm(a0)^2];
end
b=.5*b;

solution=inv(A)*b
   
if m~=0 
    check_ineq=replace(g,x, solution)
end

if l~=0
   check_eq=replace(h,x, solution)    
end
toc
break

end