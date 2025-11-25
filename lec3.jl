using FFTW, Plots
f(x) = exp(sin(5x))/(1+sin(cos(x)))
# f(x) = sin(x^3*(x-2π)^2/10)
# f(x) = exp(erf(x^2)+x^5)*sinpi(3x) + x
plot(f,0,2π,lw=3)

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
    lw = 3,
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
    plot(xj_pad, fNj_pad, legend=false, title = "Plot f(x) with padding",
        xlabel = "\$x\$",
        ylabel = "\$f(x)\$",lw=3)
end
plot_fourier(ck)
plot!(f,0,2π,lw=3,size=(600,300))
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
    xlabel = "\$x\$", ylabel = "\$f(x)\$")
end

plot_fourier(c) # plot f(t)
plot_fourier!(d,[a,a+2π/ω]) # plot g(x)
# ==


