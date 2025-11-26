#######################
# Aliasing in Fourier Series
#######################
using Plots

f1 = x->sinpi(x)
f2 = x->sinpi(9x)

plot(f1,-1,1,
    line       = 2,
    size       = (600,300),
    legend     = false,
)

plot!(f2,-1,1,line=1.6)
# ==

x = -1:0.125:1
scatter!(x,f1.(x))
# savefig("aliasing.pdf")

# ==
include("FourierChebyshev.jl")
#f(x)の概形
f(x) = exp(sin(5x))/(1+sin(cos(x)))
plot(f,0,2π,size=(800,400),legend=false,xlabel="\$x\$", ylabel="\$f(x)\$")