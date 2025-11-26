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

#######################
# Fourier Series Convolution
#######################
include("FourierChebyshev.jl")
#f(x)の概形
f(x) = exp(sin(5x))/(1+sin(cos(x)))
plot(f,0,2π,size=(800,400),legend=false,xlabel="\$x\$", ylabel="\$f(x)\$")

g(x) = f(x)^2
plot(g,0,2π,size=(800,400),xlabel="\$x\$", ylabel="\$f(x)^2\$",label = "Julia")


#######################
# naive FFT Algorithm
#######################
M = 70
a = fouriercoeffs(f,M)
p = 2
N = (p-1)*M
ta = [zeros(N);a;zeros(N)] # 1. Padding zeros
tb = ifft(ifftshift(ta)) # 2. IFFT of ta
tbp = tb.^p # 3. tb*^tb
tc = fftshift(fft(tbp))*(2.0*p*M-1)^(p-1) # 4. FFT of tb2 
plot_fourier!(tc[p:end-(p-1)],label="FFT algorithm")


# ==
using SpecialFunctions
f(x) = erf(sin(3x)+cos(2x))^4
plot(f, 0, 2π, legend=false, size=(800,400))

M = 150
p = 5
# f(x) = erf(sin(3x)+cos(2x))^4
g(x) = f(x)^p
plot(g,0,2π,size=(800,400),legend=false,xlabel="\$x\$", ylabel="\$f(x)^{$(p)}\$")

a = fouriercoeffs(f,M) # size(a) = 2M-1
# plot(abs.(a),yscale=:log10,)
function powerconvfourier(a::Vector{Complex{T}},p) where T
    M = Int((length(a)+1)/2)
    N = (p-1)*M
    ta = [zeros(N);a;zeros(N)] # 1. Padding zeros: size(ta) = 2pM-1
    tb = ifft(ifftshift(ta)) # 2. IFFT of ta
    tbᵖ = tb.^p # 3. tb*^tb
    cᵖ = fftshift(fft(tbᵖ))*(2.0*p*M-1)^(p-1)
    return cᵖ[N+1:end-N], cᵖ[p:end-(p-1)]# return (truncated, full) version
end

ap, ap_full = powerconvfourier(a,p)
plot_fourier!(ap_full)

plot_fouriercoeffs(ap_full)

# ==
f(x) = exp(-(100*(x-0.5)^2))
a = cheb(f) # Two-sided

using RadiiPolynomial
a_cheb = Sequence(Chebyshev(length(a)-1),[a[1];0.5*a[2:end]])
# space(a_cheb)
# plot(x -> a_cheb(x),-1,1)
a_cheb_ext = project(a_cheb,Chebyshev(500))

ap = (a_cheb_ext)^(0.4)
# plot(x->ap(x),-1,1)
plot(abs.((ap[:])),yscale=:log10)