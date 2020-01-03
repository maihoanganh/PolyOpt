"""
Constrained POP: f*=inf{f(x): x in S(g,h)} where S(g,h) is a semialgebraic set of inequalities g={g1,...,gm} 
and equalities h={h1,...,hl}.

Hierarchy relaxation: rho_k=sup{lambda: theta^k*(f-lambda + eps*theta^d) in quadratic module index k of S(g,h)}.
"""

using DynamicPolynomials

using JuMP

using MosekTools

using CPUTime





while 1>0
    
start = time()





# define polnomial variables
@polyvar x1 x2 x3
x=[x1;x2;x3];n=length(x)

# Polynomial to optimize 
f = (x1+x2)*(x2+x3)
    
#inequalities polynomial
g = x ; m=length(g)


#equalities polynomial
h = [x1*x2*x3*(x1+x2+x3)-1] ; l=length(h)


# small parameter
eps = 1e-5

# index of relaxation
k = 5









# quadraic polynomial
theta=1+x'*x



# Degree of objective polynomial
if m==0 && l==0
    df = ceil(Int,degree(leadingmonomial(f))/2)
else
    df = floor(Int,degree(leadingmonomial(f))/2) + 1
end

# Degree of inequalities polynomials
dg = []
for i = 1:m
    dg = [dg; ceil(Int,degree(leadingmonomial(g[i]))/2)]
end


# Degree of inequalities polynomials
dh = []
for j = 1:l
    dh = [dh; ceil(Int,degree(leadingmonomial(h[j]))/2)]
end


# Define vetor of monomials
v0= monomials(x, 0)
for j in 1:k+df
    v0= [v0;monomials(x, j)]
end

w0= monomials(x, 0)
for j in 1:2*(k+df)
    w0= [w0;monomials(x, j)]
end

length_v=[];
for i in 1:m
    length_v=[length_v;binomial(k+df-dg[i]+n,n)]
end

length_v_max=length(v0)

length_w=[];
for j in 1:l  
    length_w=[length_w;binomial(2*(k+df-dh[j])+n,n)]
end

# Define sum of square cone
model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))

# Weighted SOS matrix
@variable(model, G0[1:length_v_max, 1:length_v_max],PSD)

# Weighted SOS decomposition
wSOS=v0'*G0*v0

for i in 1:m
    
    G = @variable(model, [1:length_v[i], 1:length_v[i]],PSD)
    
    wSOS = wSOS + g[i]*v0[1:length_v[i]]'*G*v0[1:length_v[i]]
    
end


for j in 1:l
    
    q = @variable(model, [1:length_w[j]])
    
    wSOS = wSOS + h[j]*w0[1:length_w[j]]'*q
    
end

@variable(model, lambda)

@constraint(model, coefficients(theta^k*(f-lambda + eps*theta^df) - wSOS) .== 0)

@objective(model, Max, lambda)

optimize!(model)

opt_val = value(lambda)

println("termination status = ", termination_status(model))

println("opt_val = ",opt_val)

elapsed = time() - start
println("elapsed time = ",elapsed)

if termination_status(model) != MOI.SLOW_PROGRESS && termination_status(model)!= MOI.OPTIMAL
    println("Increase order k or increase parameter epsilon!")
    break
end




"""
Find a feasible point in a semialgebraic set S(g,h)  of inequalities g={g1,...,gm} and equalities h={h1,...,hl}.

Adding Sphere Inequalities Method: 

1)  Let a0,a1,...,an in R^n such that a1-a0,...,an-a0 are linear independent. Eg. a0=0 and a1,...,an are natural basis

2)  Get L0 = sup{lambda: theta^k*(|x|^2 - lambda + eps*theta^d) in quadratic module index k of S(g,h)}

    Set g = g cup {L0 - |x|^2} ==> S(g,h) satisfies Archimedian condition.
    
3)  Solve the numerical scheme of SDPs 

          omega_k^j = sup{lambda: |x-aj|^2 - lambda in quadratic module index k of S(g cup {omega_k^i - |x-ai|^2: i=0,...,k-1},h)}, j=1:m

    Check number zero eigenvalue to obtain a feasible solution.
"""




using LinearAlgebra


using RowEchelon





#inequalities polynomial
if m==0
    g =[opt_val-f]; m=length(g)
    dg = [ceil(Int,degree(leadingmonomial(g[m]))/2)]
else
    g =[g; opt_val-f]; m=length(g)
    dg = [dg; ceil(Int,degree(leadingmonomial(g[m]))/2)]
end


# small parameter for L
eps = 1e-2

# for rank of moment matrix
TOL=1e-2

# parameter for pivot of rref
tau=1e-3




# Define centers and square of radius
a0=zeros(Float64,(n, 1)); a = Matrix{Float64}(I, n, n)




dmax=maximum([dg;dh;1])




println("Determine L0:")
 
println("---------------------------------------------------------------------------------------")




# Define vetor of monomials
v0= monomials(x, 0)
for j in 1:k+2
    v0= [v0;monomials(x, j)]
end

w0= monomials(x, 0)
for j in 1:2*(k+2)
    w0= [w0;monomials(x, j)]
end

length_v_max=length(v0)

length_v=[];
for i in 1:m
    length_v=[length_v; binomial(2+(k+dmax-2)-dg[i]+n,n)]
end

length_w=[];
for j in 1:l  
    length_w=[length_w;binomial(2*(2+(k+dmax-2)-dh[j])+n,n)]
end
# Define sum of square cone
model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))

# Weighted SOS matrix
@variable(model, G0[1:length_v_max, 1:length_v_max],PSD)

# Weighted SOS decomposition
wSOS=v0'*G0*v0

for i in 1:m
    
    G = @variable(model, [1:length_v[i], 1:length_v[i]],PSD)
    
    wSOS = wSOS + g[i]*v0[1:length_v[i]]'*G*v0[1:length_v[i]]
    
end


for j in 1:l
    
    q = @variable(model, [1:length_w[j]])
    
    wSOS = wSOS + h[j]*w0[1:length_w[j]]'*q
    
end

@variable(model, lambda)

@constraint(model, coefficients(theta^(k+dmax-2)*(x'*x - lambda + eps*theta^2) - wSOS) .== 0)

@objective(model, Max, lambda)

optimize!(model)

L0 = value(lambda)

status= termination_status(model)

println("termination status = ", status)

println("L0 = ", L0)

println("=======================================================================================")



if termination_status(model) != MOI.SLOW_PROGRESS && termination_status(model)!= MOI.OPTIMAL
    println("Increase order k or increase parameter epsilon!")
    break
end










# Define omegat, t=0,...,n
omega0 = 0; omega = zeros(n)






println("Determine omega",0,":")

println("---------------------------------------------------------------------------------------")

# Define vetor of monomials
v0= monomials(x, 0)
for j in 1:k+dmax
    v0= [v0;monomials(x, j)]
end

w0= monomials(x, 0)
for j in 1:2*(k+dmax)
    w0= [w0;monomials(x, j)]
end

length_v_max=length(v0)

length_v=[];
for i in 1:m
    length_v=[length_v; binomial(dmax+k-dg[i]+n,n)]
end

length_w=[];
for j in 1:l  
    length_w=[length_w;binomial(2*(dmax+k-dh[j])+n,n)]
end

length_v_sphere=binomial(dmax+k-1+n,n)


# Define sum of square cone
model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))


# Weighted SOS matrix
G0=@variable(model, [1:length_v_max, 1:length_v_max])

consG0=@constraint(model, G0 in PSDCone())

@variable(model, lambda)

# Weighted SOS decomposition
wSOS=v0'*G0*v0

for i in 1:m
    
    G = @variable(model, [1:length_v[i], 1:length_v[i]],PSD)
    
    wSOS = wSOS + g[i]*v0[1:length_v[i]]'*G*v0[1:length_v[i]]
    
end




for j in 1:l
    
    q = @variable(model, [1:length_w[j]])
    
    wSOS = wSOS + h[j]*w0[1:length_w[j]]'*q
    
end

G = @variable(model, [1:length_v_sphere, 1:length_v_sphere],PSD)

wSOS=wSOS+(L0-x'*x)*v0[1:length_v_sphere]'*G*v0[1:length_v_sphere]

@constraint(model, coefficients([x-a0]'*[x-a0]-lambda-wSOS).== 0)

@objective(model, Max, lambda)

optimize!(model)

omega0 = value(lambda)

M=dual(consG0)
rM=rank(M, TOL)
dimsubM=binomial(k+n,n)
subM=M[1:dimsubM,1:dimsubM]
rsubM=rank(subM, TOL)

println("termination status = ", termination_status(model))

println("omega",0," = ", value(lambda))

println("rank of moment matrix = ", rM)

println("rank of moment submatrix = ", rsubM)

if termination_status(model) != MOI.SLOW_PROGRESS && termination_status(model)!= MOI.OPTIMAL
    println("Increase order k or increase parameter epsilon!")
    break
end

# extraction of Henrion and Lasserre
if rM==rsubM # check the flat-extension
    F = eigen(value.(G0))
    V = F.vectors
    r=rM
    Ix=sortperm(F.values)

    V=V[:,Ix[1:r]]
    V=Matrix(V')
    V= rref_with_pivots!(V,tau);
    U=V[1]

    U=Matrix(U')
    # Figure out multiplying matrices using YALMIP code
    w=v0[V[2]];
    N=zeros(length(V[2]),r,n)
    for i in 1:n
        xw=x[i]*w
        kk=indexin(xw,v0)
        N[:,:,i]=U[kk,:]
    end



    # Create random convex combination
    rands = rand(n,1);rands = rands/sum(rands);
    M = zeros(length(V[2]),r);
    for i in 1:n
        M=M+rands[i]*N[:,:,i];
    end

    F= schur(M);
    L=F.Z
    # Extract solution
    for i in 1:r
        solution=[]
        for j = 1:n
            solution=[solution;L[:,i]'*N[:,:,j]*L[:,i]];
        end

        println("solution = ",solution)
        #check the feasibility of solution 
        for i in 1:m
            println("check inequality ",i," = ",polynomial(g[i])(x => solution))         
        end

        for i in 1:l
            println("check equality ",i," = ",polynomial(h[i])(x => solution))
        end
        elapsed = time() - start
        println("elapsed time = ",elapsed)
        println("---------------------------------------------------------------------------------------")
    end
    break
 end


 println("=======================================================================================")
    





length_v_sphere=binomial(2*(dmax+k-1)+n, n)

for j in 1:n

    println("Determine omega_k^",j,":")

    println("---------------------------------------------------------------------------------------")

    # Define sum of square cone
    model = Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    
    # Weighted SOS matrix
    G0=@variable(model, [1:length_v_max, 1:length_v_max])

    consG0=@constraint(model, G0 in PSDCone())

    @variable(model, lambda)

    # Weighted SOS decomposition
    wSOS=v0'*G0*v0

    for i in 1:m

        G = @variable(model, [1:length_v[i], 1:length_v[i]],PSD)

        wSOS = wSOS + g[i]*v0[1:length_v[i]]'*G*v0[1:length_v[i]]

    end

    for j in 1:l
        q = @variable(model, [1:length_w[j]])

        wSOS = wSOS + h[j]*w0[1:length_w[j]]'*q

    end

    q = @variable(model, [1:length_v_sphere])

    wSOS = wSOS + (omega0-[x-a0]'*[x-a0])*w0[1:length_v_sphere]'*q


    #G = @variable(model, [1:length_v_sphere, 1:length_v_sphere],PSD)

    #wSOS = wSOS + (L0-x'*x)*v0[1:length_v_sphere]'*G*v0[1:length_v_sphere]
    
    if j>=2
        for i in 1:j-1

            Q= @variable(model, [1:length_v_sphere])

            wSOS = wSOS + (omega[i]-[x-a[:,i]]'*[x-a[:,i]])*w0[1:length_v_sphere]'*Q

        end 
    end
    
    @constraint(model, coefficients([x-a[:,j]]'*[x-a[:,j]] - lambda - wSOS) .== 0)

    @objective(model, Max, lambda)

    optimize!(model)
    
    omega[j] = value(lambda)


    M=dual(consG0)
    rM=rank(M, TOL)
    dimsubM=binomial(k+n,n)
    subM=M[1:dimsubM,1:dimsubM]
    rsubM=rank(subM, TOL)


    println("termination status = ", termination_status(model))
    
    println("omega",j," = ", value(lambda))

    println("rank of moment matrix = ", rM)

    println("rank of moment submatrix = ", rsubM)

    if termination_status(model) != MOI.SLOW_PROGRESS && termination_status(model)!= MOI.OPTIMAL
        println("Increase order k or increase parameter epsilon!")
        break
    end

    if rM==rsubM
        F = eigen(value.(G0))
        V = F.vectors
        r=rM
        Ix=sortperm(F.values)
        
        V=V[:,Ix[1:r]]
        V=Matrix(V')
        V= rref_with_pivots!(V,tau);
        U=V[1]

        U=Matrix(U')
        # Figure out multiplying matrices using YALMIP code
        w=v0[V[2]];
        N=zeros(length(V[2]),r,n)
        for i in 1:n
            xw=x[i]*w
            kk=indexin(xw,v0)
            N[:,:,i]=U[kk,:]
        end



        # Create random convex combination
        rands = rand(n,1);rands = rands/sum(rands);
        M = zeros(length(V[2]),r);
        for i in 1:n
            M=M+rands[i]*N[:,:,i];
        end

        F= schur(M);
        L=F.Z
        # Extract solution
        for i in 1:r
            solution=[]
            for j in 1:n
                solution=[solution;L[:,i]'*N[:,:,j]*L[:,i]];
            end

            println("solution = ",solution)

            for i in 1:m
                println("check inequality ",i," = ",polynomial(g[i])(x => solution))         
            end

            for i in 1:l
                println("check equality ",i," = ",polynomial(h[i])(x => solution))
            end
            elapsed = time() - start
            println("elapsed time = ",elapsed)
            println("---------------------------------------------------------------------------------------")
        end
        break

     end

    
    
     println("=======================================================================================")
    
end
    
break
end
