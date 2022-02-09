using ApproxFun, QuantumTimeSteppers, Plots, MTFun

F = Fourier(-4..4)
u0F = Fun(x -> exp(10im*x)*exp(-30*(x+2.0)^2),F,1024)
Δt = 0.01
ε = 0.01
V = x -> .5*x^2
T = 20
tF, uF = Strang_evolve(u0F,V,ε,Δt,T,0);
@gif for k=1:5:length(tF)
    plot(uF[k],ylims=(-2,2))
end

MT = MalmquistTakenaka(10.0im)
u0MT = Fun(x -> exp(10im*x)*exp(-30*(x+2.0)^2),MT,1024)
Δt = 0.01
ε = 0.01
V = x -> .5*x^2
T = 20
tMT, uMT = Strang_evolve(u0MT,V,ε,Δt,T,0);
@gif for k=1:5:length(tMT)
    plot(uMT[k],xlims=(-4,4),ylims=(-2,2))
end


H = HermiteFSE()
u0H = Fun(x -> exp(10im*x)*exp(-30*(x+2.0)^2),H,1024)
Δt = 0.01
ε = 0.01
V = x -> .5*x^2
T = 20
tH, uH = Strang_evolve(u0H,V,ε,Δt,T,0);
@gif for k=1:5:length(tH)
    plot(uH[k],xlims=(-4,4),ylims=(-2,2))
end