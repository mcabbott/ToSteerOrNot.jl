export zero_eigvec
export zero_eigvec_iter, zero_eigvec_iter2, zero_eigvec_iter3
export second, third

someeltype(xyz...) = promote_type(map(_eltype, xyz)...)
_eltype(::Nothing) = Bool
_eltype(x::AbstractArray) = eltype(x)

second(xs) = xs[2]
third(xs) = xs[3]

#####
##### Zero-eigval eigenvectors
#####

"""
    zero_eigvec(M, nrm=sum)

Assumes matrix `M` has at least one zero eigenvalue,
and returns the corresponding eigenvector, assumed to be real.

Result is normalised `nrm(v) == 1`, by default `sum`.

Several implementations via dispatch!

Keyword `realsort=true` replaces default `sortby=abs2` with `sortby=(-)∘real`.
Folk wisdom is that a small positive eigenvalue is likely the 0 + numerical noise,
while all the rest from F-P will typically be -ve
"""
zero_eigvec(mat::AbstractMatrix{<:Real}, nrm::Function=sum; kw...) = zero_eigvec_LA(mat, nrm; kw...)
zero_eigvec(mat::AbstractMatrix{<:ForwardDiff.Dual}, nrm::Function=sum; kw...) = zero_eigvec_DE(mat, nrm; kw...)

function zero_eigvec_LA(mat::AbstractMatrix, nrm::Function=sum; verbose::Bool=true, realsort::Bool=false)
    if realsort
        eig = LinearAlgebra.eigen(mat; sortby=(-)∘real)
    else
        eig = LinearAlgebra.eigen(mat; sortby=abs2) # doesn't allow ForwardDiff
    end
    # @assert abs(eig.values[1]) < 1e-10
    val = eig.values[1]
    if verbose && abs(val) > 1e-10
        @info "zero eigenvalue (LA) is:" val = ForwardDiff.value(val)
    end
    prob = real(eig.vectors[:,1])
    @assert sum(abs2, imag(eig.vectors[:,1])) < 1e-10

    prob ./ nrm(prob)
end

function zero_eigvec_DE(mat::AbstractMatrix, nrm::Function=sum; verbose::Bool=true, realsort::Bool=false)
    n, _ = size(mat)
    _vals, _vecs = DifferentiableEigen.eigen(mat)
    eig_values = _vals[1:2:end]
    # @assert all(iszero, _vals[2:2:end])
    eig_vectors = reshape(_vecs[1:2:end], n, n)
    if realsort
        val, i = findmax(real, eig_values)
    else
        val, i = findmin(abs2, eig_values)
    end
    # @assert abs(val) < 1e-10
    if verbose && abs(val) > 1e-10
        @info "zero eigenvalue (DE) is:" val = ForwardDiff.value(val)
    end
    prob = eig_vectors[:,i]
    prob ./ nrm(prob)
end

function zero_eigvec_SVD(mat::AbstractMatrix, nrm::Function=sum; verbose::Bool=true, realsort::Bool=false)  # faster, but no ForwardDiff and wrong Zygote?
    S = svd(mat)
    if realsort
        val, i = findmax(real, S.S)
    else
        val, i = findmin(abs2, S.S)
    end
    if verbose && abs(val) > 1e-10
        @info "zero eigenvalue (SVD) is:" val = ForwardDiff.value(val)
    end
    prob = S.Vt[i, :]
    prob ./ nrm(prob)
end

"""
    zero_eigvec_iter(M, nrm=sum)

Same as `zero_eigvec(M)` but works by an iterative algorithm.
This is much faster, but of dubious accuracy?
Accepts `ForwardDiff.Dual`, again accuracy unclear.
"""
function zero_eigvec_iter(mat::AbstractMatrix{T}, nrm::Function=sum; tol::Real=1e-10, maxiter::Int=50) where T
    shifted_lu = lu(mat - T(tol)*I)  # smallest eig now slightly negative?
    v_guess = T.(rand(size(mat, 2))) |> normalize!
    err = T(0)
    for _ in 1:maxiter
        v_new = shifted_lu \ v_guess |> normalize!
        err = 1 - abs(dot(v_guess, v_new))
        if err < tol
            return v_new ./ nrm(v_new)
        end
        v_guess = v_new
    end
    @warn "`zero_eigvec_iter` failed to converge" maxiter err = ForwardDiff.value(err)
    return v_guess ./ nrm(v_guess)
end


function zero_eigvec_iter2(mat::AbstractMatrix, nrm::Function=sum; tol::Real=1e-10, maxiter::Int=50, scale::Real=2)
    T = float(eltype(mat))
    T <: ForwardDiff.Dual && @warn "zero_eigvec_iter2 gives wrong answers with ForwardDiff!"
    v = T.(rand(size(mat, 2))) |> normalize!
    err = T(0)
    shift = T(tol)

    for _ in 1:maxiter
        shifted = mat - I * shift * T(scale)
        v_new = lu(shifted) \ v |> normalize!

        err = 1 - abs(dot(v, v_new))
        if err < tol
            return v_new ./ nrm(v_new)
        end

        shift = dot(v, mat, v) / dot(v, v)  # this is an eigenvalue estimate?
        v = v_new
    end

    @warn "`zero_eigvec_iter2` failed to converge" maxiter err = ForwardDiff.value(err)
    return v ./ nrm(v)
end

"""
    zero_eigvec_iter3(M, nrm=sum; [tol, maxiter, shift])

Same as `zero_eigvec(M)` but works by by calling `IterativeSolvers.invpowm!`.

Has a special method for `ForwardDiff.Dual`.
"""
function zero_eigvec_iter3(mat::AbstractMatrix, nrm::Function=sum; verbose::Bool=true, kw...)
    val, eig = _zero_eigvec_iter3(mat; kw...)
    if verbose && abs(val) > 1e-10
        @info "zero eigenvalue (iter 3) is:" val = ForwardDiff.value(val)
    end
    eig ./ nrm(eig)
end

function _zero_eigvec_iter3(A::AbstractMatrix; tol::Real=1e-10, maxiter::Int=50, shift::Real=1e-12)
    n = size(A, 1)
    T = float(eltype(A))
    # B = LinearMaps.LinearMap{T}(x -> A\x, n)
    F = lu(A - shift * I)  # consider allowsingular = true here
    B = LinearMaps.LinearMap{T}((y, x) -> ldiv!(y, F, x), n; ismutating=true)
    λ, v = IterativeSolvers.invpowm!(B, T.(rand(n)); shift, tol, maxiter)
    v ./= sum(v)  # this is needed for ForwardDiff.Dual method to know what to expect?
    (; λ, v)
end

function _zero_eigvec_iter3(A::AbstractMatrix{<:ForwardDiff.Dual}; kw...)
    A0 = ForwardDiff.value.(A)
    λ, v0 = _zero_eigvec_iter3(A0; kw...)  # really λ0, don't bother to work out perturbation of λ
    b = one.(v0)'
    @assert b * v0 ≈ 1  # this is the normalisation we assume for the eigenvector
    dvdλ = [A0 -v0; b 0] \ [A * v0; 0]  # here we only really want partials(A)
    v = @views v0 .- dvdλ[1:end-1] .+ ForwardDiff.value.(dvdλ[1:end-1])
    (; λ, v)
end


# https://julialinearalgebra.github.io/ArnoldiMethod.jl/stable/
#=
using ArnoldiMethod, LinearMaps

function zero_eigvec4(A::AbstractMatrix, nrm::Function=sum; kw...)
    # F = factorize(A)
    F = lu(A + 1e-16 * I; allowsingular = true)
    B = LinearMaps.LinearMap{eltype(A)}((y, x) -> ldiv!(y, F, x), size(A,1), ismutating=true)

    decomp, history = ArnoldiMethod.partialschur(B; nev=1, tol=1e-10, which=:LM, restarts=100);
    λs_inv, X = ArnoldiMethod.partialeigen(decomp)

    λs = 1 ./ λs_inv
    # @show λs
    eig = @view X[:,1]
    eig ./ nrm(eig)
end

=#
