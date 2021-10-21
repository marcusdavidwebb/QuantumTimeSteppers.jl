
function plan_R0_step(u::Fun,V,ε,Δt)
    itransplan = plan_itransform(space(u),coefficients(u))
    Vmultiplier = exp.(-1im*(Δt/ε)*V.(points(u)))
    transplan = plan_transform(space(u),values(u))
    return transplan, itransplan, Vmultiplier
end

function R0_step(u::Fun, transplan, itransplan, Vmultiplier)
    vals = itransplan*coefficients(u)
    newcoeffs = transplan*(Vmultiplier.*vals)
    return Fun(space(u),newcoeffs)
end

function plan_R1_step(u::Fun,ε,Δt,order = 0)
    n = length(coefficients(u))
    D2 = Derivative(space(u),2)[1:n,1:n]
    if order == 2
        R1D2top = 2I + 1im*ε*Δt*D2
        R1D2bottom = 2I - 1im*ε*Δt*D2
    elseif order == 4
        R1D2top = 12I + 6im*ε*Δt*D2 - (ε*Δt*D2)^2
        R1D2bottom = 12I - 6im*ε*Δt*D2 - (ε*Δt*D2)^2
    elseif order == 6
        R1D2top = 120I + 60im*ε*Δt*D2 - 12*(ε*Δt*D2)^2 - 1im*(ε*Δt*D2)^3
        R1D2bottom = 120I - 60im*ε*Δt*D2 - 12*(ε*Δt*D2)^2 + 1im*(ε*Δt*D2)^3
    elseif order == 0
        if bandwidths(D2) == (0,0) # for example, Fourier
            R1D2top = complex(D2)
            R1D2top.data = exp.(1im*ε*Δt*R1D2top.data)
        else
            R1D2top = exp(Matrix(1im*ε*Δt*D2))
        end
        R1D2bottom = I
    end
    return (R1D2bottom, R1D2top)
end

function R1_step(u::Fun,R1D2)
    newcoeffs = R1D2[1] \ (R1D2[2] * coefficients(u))
    return Fun(space(u),newcoeffs)
end

function plan_R2_step(u0::Fun,V1,V2,ε,Δt)
    n = length(coefficients(u0))
    R2D2 = (Δt^3 * ε / 6 ) * Derivative(space(u0),2)[1:n,1:n]
    # it is assumed that the transforms were produced in plan_R1_step
    #itransplan = plan_itransform(space(u0),coefficients(u0))
    #transplan = plan_transform(space(u0),values(u0))
    R2V1 = (Δt^3 / 12ε) * V1.(points(u0)).^2
    R2V2 = V2.(points(u0))
    return R2D2, R2V1, R2V2
end

function R2_step(u::Fun, transplan, itransplan, R2D2, R2V1, R2V2, Krylov_dim=10)
    q1 = complex(coefficients(u))
    T = eltype(q1)
    normq1 = norm(q1)
    n = length(q1)
    m = Krylov_dim
    Q = zeros(T,n,m) 
    Q[:,1] = q1/normq1
    H = zeros(T,m,m)
        
    for k=1:m
        z = one(T)im * R2D2 * (transplan * (R2V2 .* (itransplan * Q[:,k])))
        z = z + one(T)im * transplan * (R2V2 .* (itransplan * (R2D2 * Q[:,k])))
        z = z + one(T)im * transplan * (R2V1 .* (itransplan * Q[:,k]))
        
        for i=1:k
            H[i,k] = Q[:,i]'*z
            z = z - H[i,k]*Q[:,i]
        end
        if k < m
            H[k+1,k] = norm(z)
            Q[:,k+1] = z / H[k+1,k]
        end
    end
    newcoeffs = normq1*(Q*exp(H))[:,1]
    return Fun(space(u), newcoeffs)
end

function plan_Lie_Trotter_step(u::Fun,V,ε,Δt,order=0)
    transplan, itransplan, Vmultiplier = plan_R0_step(u,V,ε,Δt)
    R1D2 = plan_R1_step(u,ε,Δt,order)
    return transplan, itransplan, Vmultiplier, R1D2
end

function Lie_Trotter_step(u::Fun, transplan, itransplan, Vmultiplier, R1D2)
    newu = R0_step(u,transplan, itransplan, Vmultiplier)
    newu = R1_step(newu, R1D2)
    return newu
end

function Lie_Trotter_evolve(u0::Fun,V,ε,Δt,T,order=0)
    M = Integer(floor(T/Δt))
    u0 = Fun(space(u0),complex(coefficients(u0)))
    retFuns = Vector{typeof(u0)}(undef,M+1)
    retFuns[1] = u0
    transplan, itransplan, Vmultiplier, R1D2 = plan_Lie_Trotter_step(retFuns[1],V,ε,Δt,order)
    for k = 1:M
        retFuns[k+1] = Lie_Trotter_step(retFuns[k], transplan, itransplan, Vmultiplier, R1D2)
    end
    return 0:Δt:T, retFuns
end


function plan_Strang_step(u::Fun,V,ε,Δt,order=0)
    transplan, itransplan, Vmultiplier = plan_R0_step(u,V,ε,Δt/2)
    R1D2 = plan_R1_step(u,ε,Δt,order)
    return transplan, itransplan, Vmultiplier, R1D2
end

function Strang_step(u::Fun, transplan, itransplan, Vmultiplier, R1D2)
    newu = R0_step(u, transplan, itransplan, Vmultiplier)
    newu = R1_step(newu, R1D2)
    newu = R0_step(newu, transplan, itransplan, Vmultiplier)
    return newu
end

function Strang_evolve(u0::Fun,V,ε,Δt,T,order=0)
    M = Integer(floor(T/Δt))
    u0 = Fun(space(u0),complex(coefficients(u0)))
    retFuns = Vector{typeof(u0)}(undef,M+1)
    retFuns[1] = u0
    transplan, itransplan, Vmultiplier, R1D2 = plan_Strang_step(retFuns[1],V,ε,Δt,order)
    for k = 1:M
        retFuns[k+1] = Strang_step(retFuns[k], transplan, itransplan, Vmultiplier, R1D2)
    end
    return 0:Δt:T, retFuns
end

function plan_Zassenhaus2_step(u::Fun,V,V1,V2,ε,Δt,order=0)
    transplan, itransplan, Vmultiplier = plan_R0_step(u,V,ε,Δt/2)
    R1D2 = plan_R1_step(u,ε,Δt/2,order)
    R2D2, R2V1, R2V2 = plan_R2_step(u,V1,V2,ε,Δt)
    return transplan, itransplan, Vmultiplier, R1D2, R2D2, R2V1, R2V2
end

function Zassenhaus2_step(u::Fun, transplan, itransplan, Vmultiplier, R1D2, R2D2, R2V1, R2V2, Krylov_dim = 10)
    newu = R0_step(u, transplan, itransplan, Vmultiplier)
    newu = R1_step(newu,R1D2)
    newu = R2_step(newu, transplan, itransplan, R2D2, R2V1, R2V2, Krylov_dim)
    newu = R1_step(newu, R1D2)
    newu = R0_step(newu, transplan, itransplan, Vmultiplier)
    return newu
end

function Zassenhaus2_evolve(u0::Fun,V,V1,V2,ε,Δt,T,order=0, Krylov_dim = 10)
    M = Integer(floor(T/Δt))
    u0 = Fun(space(u0), complex(coefficients(u0)))
    retFuns = Vector{typeof(u0)}(undef,M+1)
    retFuns[1] = u0
    transplan, itransplan, Vmultiplier, R1D2, R2D2, R2V1, R2V2 = plan_Zassenhaus2_step(retFuns[1],V,V1,V2,ε,Δt,order)
    for k = 1:M
        retFuns[k+1] = Zassenhaus2_step(retFuns[k], transplan, itransplan, Vmultiplier, R1D2, R2D2, R2V1, R2V2, Krylov_dim)
    end
    return 0:Δt:T, retFuns
end
