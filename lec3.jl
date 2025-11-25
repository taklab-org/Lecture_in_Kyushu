######################
# Fourier Series
######################

using FFTW, Plots
f(x) = exp(sin(5x))/(1+sin(cos(x)))
# f(x) = sin(x^3*(x-2π)^2/10)
# f(x) = exp(erf(x^2)+x^5)*sinpi(3x) + x
plot(f,0,2π,lw=1.6,legend=false,
    xlabel = "\$x\$",
    ylabel = "\$f(x)\$",
    # size=(600,300)
)

# ==
# setprecision(104)
N = 100

function fouriercoeffs(f,N)
    h = (2.0)*π/(2N-1)
    xj = h*(0:2N-2)
    fj = f.(xj);
    return fftshift(fft(fj))/(2N-1)
end

ck = fouriercoeffs(f,N)

plot(-N+1:N-1,abs.(ck),yscale=:log10,
    legend=false,
    xlabel = "\$k\$",
    ylabel = "\$|\\bar{c}_k\\,|\$",
    lw = 1.6,
    # yticks = [1e-30,1e-16,1e-8,1],
    title = "Absolute values of Fourier coefficients"
)
# ==

function plot_fourier(ck)
    N = (length(ck)+1)/2# 2N-1
    n_pad = 200
    ck_pad = [zeros(n_pad);ck;zeros(n_pad)]
    N_pad = N + n_pad
    h_pad = 2.0π/(2N_pad-1)
    xj_pad = h_pad*(0:(2N_pad-2))
    
    fNj_pad = real((2N_pad-1)*ifft(ifftshift(ck_pad)))
    plot(xj_pad, fNj_pad, legend=false, 
        title = "Plot f(x) with padding",
        xlabel = "\$x\$",
        ylabel = "\$f(x)\$",lw=1.6)
end
plot_fourier(ck)
plot!(f,0,2π,lw=1.6,size=(600,300))
# savefig("plot_fourier.pdf")
# ==

function fouriercoeffs(f, N, I=[0,2π])
  a = I[1]; b = I[2]
  h = (b-a)/(2N-1)
  j = 0:2N-2
  xj = a .+ j*h
  fj = f.(xj);
  return fftshift(fft(fj))/(2N-1)
end
function plot_fourier(ck, I=[0,2π])
  a = I[1]; b = I[2]
  N = (length(ck)+1)/2 # 2N-1
  n_pad = 200
  ck_pad = [zeros(n_pad);ck;zeros(n_pad)]
  N_pad = N + n_pad
  h_pad = (b-a)/(2N_pad-1)
  xj_pad = a .+ h_pad*(0:2N_pad-2)

  fNj_pad = real((2N_pad-1)*ifft(ifftshift(ck_pad)))
  plot(xj_pad, fNj_pad, legend=false,
  xlabel = "\$x\$", ylabel = "\$f(x)\$")
end

N = 30
g(x) = exp(sin(5x))/(1+sin(cos(x)))
ω = 1.3
a = 1
f(t) = g(ω*(t-a))
c = fouriercoeffs(g,N)
d = fouriercoeffs(f,N,[a;a+2π/ω])
k = (-N+1):(N-1)
reshape([c;d],2N-1,2)

function plot_fourier!(ck, I=[0,2π])
  a = I[1]; b = I[2]
  N = (length(ck)+1)/2 # 2N-1
  n_pad = 200
  ck_pad = [zeros(n_pad);ck;zeros(n_pad)]
  N_pad = N + n_pad
  h_pad = (b-a)/(2N_pad-1)
  xj_pad = a .+ h_pad*(0:2N_pad-2)

  fNj_pad = real((2N_pad-1)*ifft(ifftshift(ck_pad)))
  plot!(xj_pad, fNj_pad, legend=false,
    xlabel = "\$x\$", ylabel = "\$f(x)\$",size=(600,300))
end

plot_fourier(c) # plot f(t)
plot_fourier!(d,[a,a+2π/ω]) # plot g(x)
# savefig("plot_profiles.pdf")
# ==


######################
# Chebyshev Series
######################
using Polynomials

function T(n)
    return ChebyshevT([zeros(n);1.])
end

plot(
    plot(T(0), ylabel = "\$T_0(x)\$"),
    plot(T(1), ylabel = "\$T_1(x)\$"),
    plot(T(2), ylabel = "\$T_2(x)\$"),
    plot(T(3), ylabel = "\$T_3(x)\$"),
    plot(T(4), ylabel = "\$T_4(x)\$"),
    plot(T(5), ylabel = "\$T_5(x)\$"),
    layout    = (3, 2),
    line        = 2,
    color      = 2, 
    xlabel = "\$x\$",
    xlims      = (-1, 1),
    ylims      = (-1.1, 1.1),
    legend   = false,
    lw = 3,
)
# ==

n = 10
tt= range(0,stop=π,length=n+1)
zz =exp.(im*tt)

plot(zz,
    title = "Equispaced  points on the unit circle",
    legend = false,
    size = (800,400),
)
scatter!(zz)

# ==
function chebpts(n, a=-1, b=1) # n: order of Chebyshev polynomials
    m = -n:2:n
    x = sinpi.(m/(2*n))
    return (1.0 .- x).*a/2 + (1.0 .+ x).*b/2
end
x = chebpts(n)
for j = 2:n
    plot!([x[j];(zz[end-j+1])],
        title = "Chebyshev points",
        legend = false,
        color=1,
    )
end
scatter!(x.+0*im)
# ==

plot!(T(n), color=2,line= 1.6,)
# ==

using SpecialFunctions
M = 50 # number of Chebyshev polynimals
# fc = Fun(x->exp(erf(x^2)+x.^5).*sin(3*pi*x) + x, Chebyshev(),M)
# n = ncoefficients(fc) - 1 # maximum order of Chebyshev polynomials (n = M-1)
f(x) = exp(erf(x^2)+x^5)*sinpi(3x) + x
n = M-1
plot(f,-1,1,
    xlabel     = "\$x\$",
    ylabel     = "\$f(x)\$",
    line       = 3,
    title      = "\$f(x) = \\exp\\,(\\mathrm{erf}\\,(x^2)+x^5)\\sin\\,(3\\pi x) + x\$", 
    size       = (800,400),
    legend     = false,
)

xj = chebpts(n)
fvals = f.(xj)
scatter!(xj,fvals)

# ==

valsUnitDisc = [reverse(fvals);fvals[2:end-1]]
FourierCoeffs = real(fft(valsUnitDisc))
ChebCoeffs = FourierCoeffs[1:n+1]/n
ChebCoeffs[1] = ChebCoeffs[1]/2
ChebCoeffs[end] = ChebCoeffs[end]/2
# reshape([ChebCoeffs; Cheb_twosided],n+1,2)
ChebCoeffs

# ==
function chebcoeffs(f,M,I=[-1;1])
    a = I[1]; b = I[2]
    n = M-1
    cpts  = chebpts(n, a, b)
    fvals = f.(cpts)
    FourierCoeffs = real(fft([reverse(fvals);fvals[2:end-1]]))
    ChebCoeffs = FourierCoeffs[1:n+1]/n
    ChebCoeffs[1] = ChebCoeffs[1]/2
    ChebCoeffs[end] = ChebCoeffs[end]/2
    return ChebCoeffs # return Two-sided Chebyshev
end

# ==
f(x) = exp(x)
chebcoeffs(f,2^3+1) # M=9として計算

# ==
chebcoeffs(f,2^4+1)

# ==
function cheb(f,I=[-1;1];tol = 5e-15,Nmax = 10000)
    a = I[1]; b = I[2]; m = 0.5*(a+b); r = 0.5*(b-a); x = rand(5)
    x1 = m .+ x*r; x2 = m .- x*r
    if f.(x1) ≈ f.(x2)
        odd_even = 1 # even function: 1
    elseif f.(x1) ≈ -f.(x2)
        odd_even = -1 #  odd function: -1
    else
        odd_even = 0 # otherwise: 0
    end
    i = 3
    schbc = 0 # sampling chebyshev coefficients
    while true
        schbc = chebcoeffs(f,2^i+1,I)
        if all(abs.(schbc[end-2:end]) .< tol) || (2^i+1 > Nmax) 
            break
        end
        i += 1
    end    
    M = findlast(abs.(schbc) .> tol)
    cc = schbc[1:M]
    # cc = chebcoeffs(f,M,I)
    if odd_even == 1 # even function
        cc[2:2:end] .= 0
    elseif odd_even == -1 # odd function
        cc[1:2:end] .= 0
    end
    return cc # return Two-sided Chebyshev
end

# ==
f(x) = sin(x^3*(x-2π)^2/10)
Cheb_twosided = cheb(f,[0,2π])

# ==
f(x) = 1/(1+1000*(x+.5)^2)+1/sqrt(1+1000*(x-.5)^2)
Cheb_twosided = cheb(f)

# ==
plot(abs.(cheb(f)),
    yscale=:log10,
    legend = false,
    xlabel="\$n\$",
    ylabel="\$|a_n\\,|\$",
    title="Chebyshev coefficients",
    size = (600,300)
)
# savefig("chebcoeffs.pdf")
# ==
function plot_cheb(ChebCoeffs_twosided;I=[-1,1],title="") # Input: Two-sided Chebyshev
    M = length(ChebCoeffs_twosided) # M: size of chebyshev
    a = I[1]; b = I[2];
    k = 0:M-1
    ξⱼ = range(-1, stop=1, length=5000)
    xc = (1.0 .- ξⱼ)*a/2 + (1.0 .+ ξⱼ)*b/2 # points in [a,b]
    fxc = cos.(Vector(k)' .* acos.(ξⱼ)) * ChebCoeffs_twosided
    plot(xc, fxc, label="plot_cheb", title=title, size=(600,300),
        xlabel="\$x\$", ylabel="\$f(x)\$")
end
plot_cheb(Cheb_twosided)
plot!(f,label="original")
# savefig("two_peaks.pdf") 
# ==