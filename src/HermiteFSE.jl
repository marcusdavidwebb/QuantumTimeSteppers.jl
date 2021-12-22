struct HermiteFSE{T<:Real} <: Space{Line{false,T},T}
    t :: T
end

HermiteFSE{T}() where T = HermiteFSE{T}(zero(T))
HermiteFSE() = HermiteFSE{Float64}()

domain(::HermiteFSE{T}) where T = Line{false,T}()
canonicaldomain(::HermiteFSE{T}) where T = Line{false,T}()
normalization(::Type{T}, ::HermiteFSE, ::Int) where T = one(T)

spacescompatible(a::HermiteFSE,b::HermiteFSE) = a.t == b.t
canonicalspace(::HermiteFSE{T}) where T = HermiteFSE(zero(T))

function Derivative(sp::HermiteFSE,order::Integer)
    ConcreteDerivative(sp,order)
end

space(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = D.space
rangespace(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = domainspace(D)
bandwidths(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = D.order,D.order
Base.stride(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = D.order

function getindex(D::ConcreteDerivative{HermiteFSE{T},K,KK},k::Integer,j::Integer) where {T,K,KK}
    # differentiation matrix does not depend on HermiteFSE.t (see Kropielnicka, Iserles, Schratz, Webb 2021)
    
    m = D.order
    bandind = j-k

    if m == 1
        if bandind == 1
            convert(T,sqrt(k/2))
        elseif bandind == -1
            convert(T,-sqrt(j/2))
        else
            zero(T)
        end
    elseif m == 2
        if bandind == 0
            convert(T,.5-k)
        elseif bandind == 2
            convert(T,sqrt(k*(k+1))/2)
        elseif bandind == -2
            convert(T,sqrt(j*(j+1))/2)
        else
            zero(T)
        end
    else
        error("Higher order HermiteFSE derivatives are not supported yet")
    end
end

function evaluate(cfs::AbstractVector,Hspace::HermiteFSE,x)
    # Clenshaw's algorithm (working implementation given here) can have overflow/underflow problems:
    # T = eltype(x)
    # n = length(cfs)   
    # bk2 = zero(T)
    # bk1 = convert(T,cfs[n])
    # for k = n-1:-1:1
    #     bk1,bk2 = convert(T,cfs[k]) + x*sqrt(2/k)*bk1 - sqrt(k/(k+1)) * bk2 , bk1
    # end
    # return bk1 * exp(-x^2 / 2) * convert(T,π)^(-1/4)

    # Rescaled forward recurrence to avoid overflow/underflow
    T = eltype(x)
    t = Hspace.t
    x = x/sqrt(1+4*t^2)
    cfs = cfs .* ((1+2*one(T)im*t)/(1-2*one(T)im*t)).^(zero(T) : one(T)/2 : n*one(T)/2)
    n = length(cfs)
    ret = zero(T)
    if (n > 0) ret += cfs[1] * exp(-x^2 / 2) * convert(T,π)^(-1/4) end
    if (n > 1) ret += cfs[2] * x * exp(-x^2 / 2) * (4/convert(T,π))^(1/4) end
    
    if (n > 2)
        hkm1 = convert(T,π)^(-1/4)       #h_0(x)
        hk = sqrt(one(T)*2) * x * hkm1  #h_1(x)
        sum_log_scale = zero(T)
        for k = 2:n-1
            # hk = h_k(x), hkm1 = h_{k-1}(x) (actually, recaled versions thereof)
            hkm1, hk = hk, sqrt(one(T)*2/k)*x*hk - sqrt((k-one(T))/k)*hkm1
            # rescale values if necessary
            scale = (x->(x<one(T)*100) ? one(T) : inv(x))(abs(hk))
            hk *= scale
            hkm1 *= scale
            # keep track of final rescaling factor
            sum_log_scale += log(scale)
            ret += cfs[k+1] * hk * exp(-x^2 / 2 - sum_log_scale)
        end
    end
    ret = ret * exp(-one(T)im*t*x^2)/sqrt(1-2*one(T)im*t)
    return ret
end

eval_Hermite_function(n::Integer,x::T, t=0.0) where T <: Number = evaluate([zeros(T,n);one(T)],HermiteFSE(t),x)

## Transforms for HermiteFSE space
hasfasttransform(::HermiteFSE) = false
# Gauss-Hermite points
points(Hspace::HermiteFSE,n::Int) = gausshermite(n)[1] * sqrt(1+4*(Hspace.t)^2)

function plan_Hermite_transform(n::Integer)
   # builds the orthogonal matrix Q and weight vector 'valweights' such that
   # the val2coeffs transform is Q * (valweights .* vals) and
   # the coeffs2vals transform is (Q' * cfs) ./ valweights.

    Q = zeros(n+1,n+1)
    x = points(HermiteFSE(),n+1)

    hkm1 = ones(n+1)*π^(-1/4)       #h_0(x)
    hk = sqrt(2) * x .* hkm1  #h_1(x)
    cum_log_scale = zeros(n+1,n+1)

    Q[1,:] = hkm1
    if (n ≥ 1) Q[2,:] = hk end

    for k = 2:n 
        # hk = h_k(x), hkm1 = h_{k-1}(x) (actually, recaled versions thereof)
        hkm1, hk = hk, sqrt(2/k)* x .* hk - sqrt((k-1)/k)*hkm1
        # rescale values if necessary
        scale = (x->(x<100) ? 1 : inv(x)).(abs.(hk))
        hk .*= scale
        hkm1 .*=  scale
        Q[k+1,:] = hk
        # keep track of rescaling factors
        cum_log_scale[k+1,:] = cum_log_scale[k,:] + log.(scale)
    end
    
    # the function values should be weighted by these
    valweights = exp.(cum_log_scale[n+1,:] + x.^2 / 2) ./ (sqrt(n+1)*abs.(Q[n+1,:]))
    
    # Q contains rescaled Hermite polynomials. Need to convert to unscaled ψ_k(x_j)/(sqrt(n+1)ψ_n(x_j)),
    # where ψ_n is the degree n Hermite function
    for k = 1:n+1
        Q[k,:] = Q[k,:]./abs.(Q[n+1,:]) .* exp.(cum_log_scale[n+1,:]-cum_log_scale[k,:] .- .5*log(n+1))
    end

    return valweights, Q
end

function plan_Hemite_transform_coeff_scaling(n::Integer, t = 0.0)
    # the coefficients should be weighted by these
    if t == 0.0
        coeffweights = ones(n+1)
    else
        x = points(HermiteFSE(t),n+1)
        coeffweights = ((1+2.0im*t)/(1-2.0im*t)).^(0:.5:n/2) * exp(-1.0im*t*x.^2)/sqrt(1-2.0im*t)
    end
    return coeffweights
end

function plan_transform(S::HermiteFSE,vals::AbstractVector)
    valweights, Q = plan_Hermite_transform(length(vals)-1)
    coeffweights = plan_Hemite_transform_coeff_scaling(length(vals)-1, S.t)
    TransformPlan(S,(valweights,Q,coeffweights),Val{false})
end
function plan_itransform(S::HermiteFSE,cfs::AbstractVector)
    valweights, Q, coeffweights = plan_Hermite_transform(length(cfs)-1)
    coeffweights = plan_Hemite_transform_coeff_scaling(length(vals)-1, S.t)
    ITransformPlan(S,(valweights,Q,coeffweights),Val{false})
end

function *(P::TransformPlan{T,S,false},vals::AbstractVector) where {T,S<:HermiteFSE}
    valweights, Q, coeffweights = P.plan
    (Q*(valweights.*vals)) ./ coeffweights
end
function *(P::ITransformPlan{T,S,false},cfs::AbstractVector) where {T,S<:HermiteFSE}
    valweights, Q, coeffweights = P.plan
    (Q' * (coeffweights.*cfs)) ./ valweights
end