using IntervalArithmetic

ix = interval(1, 2)
iy = ix^2
@show iy

######################
# Ex.1 Root finding
######################
# Problem f(x) = 0
F(x) = x^2 - 2
DF(x) = 2x
x0 = 1

# Newton iteration
using LinearAlgebra

function newton_F(x)    
    num_itr = 0; p = 1; tol = 5e-15
    Fx = F(x)    
    println("Before iteration: $(norm(Fx,p))")
    while num_itr ≤ 100
        dx = DF(x)\Fx
        x = x - dx
        num_itr += 1
        Fx = F(x)
        println("After $(num_itr) iteration: $(norm(Fx,p))")
        if norm(Fx,p)/ norm(x,p)< tol
            break
        end
    end
    return x
end

x_bar = newton_F(x0)
@show x_bar
@show sqrt(2.);
# ==


# Validating the approximate solution
using IntervalArithmetic
ix = interval(x_bar)
A_dagger = DF(ix); A = 1/A_dagger

# Y0 bound
Y0 = abs(A*F(x_bar))
@show Y0

# Z0, Z1
Z0 = 0.; Z1 = 0.
@show Z0
@show Z1

# Z2 bound
Z2 = A*2
@show Z2

# radii-polynomial
p(r) = Z2*r^2 - (1-Z0-Z1)*r + Y0
# 
a = Z2; b = - (1-Z0-Z1); c = Y0;
if sup(b^2-4*a*c) < 0
  println("error: connot find the root of radii polynomial")
else
  r_min = (-b - sqrt(b^2-4*a*c))/2/a
  r_max = (-b + sqrt(b^2-4*a*c))/2/a
end
@show sup(r_min), inf(r_max)
r0 = sup(r_min);
if sup(p(r0))<0 # validate p(r0) < 0
  sol = interval(x_bar,r0;format=:midpoint)
end
# ==

# Plot radii polynomial
using Plots
midp(r) = mid(p(r))
plot(midp,-0.1,2,xlabel="r",ylabel="p(r)",legend=false)
# ==



######################
# Ex.2 Logistic map
######################
using Plots
# λ = 4.215351654086267
λ = 1.2
P(x) = λ*x*(1-x)

length_x = 49
x = zeros(length_x)
x[1] = 0.01
# x[1] = 0.9520156078132255
 # 0.19256523180939808
 # 0.6554191603773762
for i=1:length_x-1
    x[i+1] = P(x[i])
end
scatter(x,xlabel="n",ylabel="xₙ",legend=false)
# scatter(x,xlabel="n",ylabel="xₙ",legend=false, size=(600,300))
# savefig("Logistic_map.pdf")
# ==


function logisticmap(λ)
    P(x) = λ*x*(1-x)
    length_x = 500
    x = zeros(length_x)
    x[1] = rand()/(1+rand())
    for i=1:length_x-1
        x[i+1] = P(x[i])
    end
    scatter!(plt,λ*ones(Int(length_x/2)),x[Int(length_x/2):end],
    legend=false, mc=:red, ms=2, ma=0.5)
end

plt=plot()
for λ = 2.8:0.01:3.9
    logisticmap(λ)
end
plot(plt,xlabel="λ",ylabel="xₙ",size=(600,300))
# savefig("Logistic_map2.pdf")
# ==


λ = 3.8284
P(x) = λ*x*(1-x)
length_x = 200
x = zeros(length_x)
x[1] = rand()/(1+rand())
for i=1:length_x-1
    x[i+1] = P(x[i])
end
scatter(x,xlabel="n",ylabel="xₙ",legend=false, mc=:red, ms=2, ma=0.5)
# ==

λ = 3.8285
# λ = 3.9
F(x) = [x[1]-λ*x[3]*(1-x[3])
    x[2]-λ*x[1]*(1-x[1])
    x[3]-λ*x[2]*(1-x[2])]
DF(x) = [1 0 -λ*(1-2x[3]);
    -λ*(1-2x[1]) 1 0;
    0 -λ*(1-2x[2]) 1]

x0 = [1., -1., 1.]
# x₀ = [.5, 1.0, 0.1]

x_bar = newton_F(x0)
@show x_bar
# ==


# Define A and A_dagger
using IntervalArithmetic
ix = interval.(x_bar)
A_dagger = DF(ix)
A = interval(inv(mid.(A_dagger)))

# Y0 bound
Y0 = norm(A*F(x_bar),1)
@show Y0

# Z0, Z1
Z0 = opnorm(I - A*A_dagger,1)
Z1 = 0
@show Z0
@show Z1

# Z2 bound
Z2 = 2*abs(interval(λ))*opnorm(A,1)
@show Z2

# radii-polynomial
p(r) = Z2*r^2 - (1-Z0-Z1)*r + Y0
#
a = Z2; b = - (1-Z0-Z1); c = Y0;
if sup(b^2-4*a*c) < 0
  println("error: connot find the root of radii polynomial")
else
  r_min = (-b - sqrt(b^2-4*a*c))/2/a
  r_max = (-b + sqrt(b^2-4*a*c))/2/a
end
@show r0 = sup(r_min);

if sup(p(r0))<0 # validate p(r0) < 0
  sol = interval.(x_bar,r0;format=:midpoint)
end
# ==


# Plot radii polynomial
using Plots
midp(r) = mid(p(r))
plot(midp,-1e-4,2e-3,xlabel="r",ylabel="p(r)",legend=false)
# ==

######################
# Ex.3 Periodic orbits with unknown period
######################
eta = 1.6311
# eta = 1.8
# eta = 2.5
# eta = 1.2

F((lam, x1, x2, x3)) = [x1 + x2 + x3 - eta
    x1 - lam*x3*(1-x3)
    x2 - lam*x1*(1-x1)
    x3 - lam*x2*(1-x2)]

DF((lam, x1, x2, x3)) = [0 1 1 1;
    -x3*(1-x3) 1 0 -lam*(1-2*x3);
    -x1*(1-x1) -lam*(1-2*x1) 1 0;
    -x2*(1-x2) 0 -lam*(1-2*x2) 1]

x0  = [3.0,1.0,-1.0,1.0]
x_bar = newton_F(x0)
@show x_bar
# ==


# Define A and A_dagger
using IntervalArithmetic
ix = interval.(x_bar)
A_dagger = DF(ix)
A = interval(inv(mid.(A_dagger)))

# Y0 bound
Y0 = norm(A*F(x_bar),1)
@show Y0
# Z0, Z1
Z0 = opnorm(I - A*A_dagger,1)
Z1 = 0
@show Z0
@show Z1
# Z2 bound
λ, x1, x2, x3 = abs.(x_bar)
B(r) = [
    0 0 0 0
    1+2*(x3)+r 0 0 (1+(x3)+r)+λ+r
    1+2*(x1)+r (1+(x1)+r)+λ+r 0 0
    1+2*(x2)+r 0 (1+(x2)+r)+λ+r 0]
Z2_func(r) = opnorm(A*B(r),1)

# radii-polynomial
p(r) = Z2_func(r)*r^2 - (1-Z0-Z1)*r + Y0
# ==

using ForwardDiff
A = inv(DF(x_bar)) # for ForwardDiff
Z2_funcmid(r) = opnorm(A*B(r),1)
pmid(r) = Z2_funcmid(r)*r^2 - (1-sup(Z0)-sup(Z1))*r + sup(Y0)
F(x) = pmid(x)
DF(x)= ForwardDiff.derivative(pmid,x)

r_min = 1.1*newton_F(0)
r_max = 0.9*newton_F(0.03)
@show r_min, r_max

# ==

if sup(p(r_min))<0 # validate p(r0) < 0
  sol = interval.(x_bar,r_min;format=:midpoint)
end
# ==

# Plot radii polynomial
using Plots
midp(r) = mid(p(r))
plot(midp,-1e-4,0.05,xlabel="r",ylabel="p(r)",legend=false,size=(600,300))
# savefig("radii-poly.pdf")