export flatten, maximise, grad

#=

Because I was writing this function too many times, I wrapped it up nicely here.

=#

"""
    maximise(fun, params::NamedTuple; method=LBFGS(), options...)
    maximise(fun, vector; kw...) = maximise(fun, (; vector); kw...)

Wrapper for Optim.optimise which:
- Maximises instead of minimising the objective function.
- Calls `flatten(params)` to handle NamedTuples of parameters. This gets keywords `limit`, `log` passed through.
- If the initial parameter is positive (or all positive numbers) then this is imposed as a constraint.
- Allows `y = fun(params)` to return a NamedTuple with aux info, provided the `first` part is ``-loss``.
- Passes most keywords into Optim.Options, except `autodiff=:forward` which is passed separately, and `method=LBFGS()` which becomes positional.
- Keeps an honest count of how many times `fun` has been called during optimisation, and records the time taken.

Returns `(; y..., params..., positive, length, calls, time)`.
Best you avoid these 4 words as parameter names!

## Examples:

```
julia> maximise(sum, [1,2,3]; limit=100)  # :objective, :vector are default names
(objective = 300.0, vector = [100.0, 100.0, 100.0], positive = (vector = +,), length = 3, calls = 197, time = 1.8)

julia> maximise(x -> (; minus = -sum(x), also="aux"), [1,2,3])  # with NamedTuple output, notice x .≥ 0 constraint
(minus = -0.0, also = "aux", vector = [0.0, 0.0, 0.0], positive = (vector = +,), length = 3, calls = 7, time = 1.0)

julia> maximise(nt -> -sum(nt.x .+ nt.y), (x=[1,2], y=[3,-4]); limit=100)  # NamedTuple input, and y allows negative entries
(objective = 200.0, x = [0.0, 0.0], y = [-100.0, -100.0], positive = (x = +, y = nothing), length = 4, calls = 153, time = 1.2)
```
"""
function maximise(obj::Function, param::Union{NamedTuple,AbstractVector};
    # KW for flatten:
        limit=10^6, log=false,
    # KW for Optim:
        # method=Optim.LBFGS(),
        method=Optim.LBFGS(linesearch=Optim.LineSearches.BackTracking()),  # https://github.com/JuliaNLSolvers/Optim.jl/issues/713
        autodiff=:forward,
        g_abstol=1e-8, # 1e-8 default https://julianlsolvers.github.io/Optim.jl/stable/user/config/
        iterations=1000, # 1000 default
        )
    v, re, positive = flatten(param; limit, log)
    calls = Ref(0)
    # cnt(x) = begin calls[]+=1; x end
    function loss(x)
        calls[] += 1
        -(_first(obj(re(x))))
    end
    t1 = time()
    # This line was missing `method` for a while, very confused by why nothing converged...
    # res = Optim.optimize((-)∘_first∘cnt∘obj∘re, v, method,
    res = Optim.optimize(loss, v, method,
        Optim.Options(; g_abstol, outer_g_abstol=g_abstol, iterations); autodiff = optim_autosymbol(autodiff))
    t2 = time()
    maxparam = re(res.minimizer)
    maxval = _ensure_nt(obj(maxparam), :objective)
    (; maxval..., _ensure_nt(maxparam, :vector)..., positive, length=length(v), calls=calls[], time=round(t2-t1, digits=1))
end

_first(nt::NamedTuple) = _first(first(nt))
_first(y::Real) = y
_first(y) = error(string("expected function to return a real number or a NamedTuple, but got this: ", y))

_ensure_nt(nt::NamedTuple, _::Symbol) = nt
_ensure_nt(x::Union{Number,AbstractVector}, sy::Symbol) = NamedTuple{(sy,)}((x,))
_ensure_nt(tup::Tuple, sy::Symbol) = NamedTuple{(sy, _auxnames(Val(length(tup)-1))...)}(tup)
_auxnames(n) = ntuple(i -> Symbol(:aux, i), n)  # ToSteerOrNot._ensure_nt((1,2,3), :sy) == (sy = 1, aux1 = 2, aux2 = 3)

"""
    flatten(param::Namedtuple; limit=10^6, log=false)

Returns `vector, rebuilder, positive`. The rebuild function does this:
- If the initial parameter was all positiver, then it imposes `clamp.(xs, 0, limit)`
- If not, i.e. any zeros or negative, then it imposes `clamp.(xs, -limit, limit)`.

With `log=true`, the flattened vector contains `log.(xs)` for positive parameters.

## Examples:

```
julia> v, re, sig = flatten((a=0.0, b=[2 3; 4 5.]); limit=10);

julia> println(v)  # integers become floats
[0.0, 2.0, 4.0, 3.0, 5.0]

julia> re([-1, 2, -3, 40, -50])  # constrains 0 .≤ b .≤ 10
(a = -1, b = [2 10; 0 0])

julia> sig
(a = nothing, b = +)
```
This example should be type-stable, `Test.@inferred re(randn(Float32, 5))`.
But construction is not, e.g. `@code_warntype flatten(rand(3))`.

Example with `log=true`:
```
julia> v, re, sig = flatten(1:3f0; log=true);

julia> println(v)
Float32[0.0, 0.6931472, 1.0986123]

julia> println(re([-1, 0, 1f0]))
Float32[0.36787945, 1.0, 2.7182817]

julia> sig
(vector = log,)
```
"""
function flatten(param::NamedTuple; limit::Real=10^6, log::Bool=false, signs=map(_describe(log), param))
    function rebuild(v::AbstractVector, sig::NamedTuple=signs)
        ind = Ref(0)
        map(param, sig) do p, s
            p isa Number && return _limit(v[ind[]+=1], s, limit)
            i = ind[]
            l = length(p)
            ind[] += l
            return _limit(reshape(view(v, i+1:i+l), size(p)), s, limit)
        end
    end
    w = vcat(map(param, signs) do p, s
        s === (Base.log) ? Base.log.(_vec(p)) : float.(_vec(p))
    end...)
    w, rebuild, signs
end
function flatten(param::AbstractVector; kw...)
    w, re, signs = flatten((; vector=param); kw...)
    w, first∘re, signs
end

function _describe(param::Union{AbstractArray, Number}, log::Bool=false)
    sig = all(>(0), param)
    (sig && log) ? Base.log : sig ? (+) : nothing
end
_describe(log::Bool) = Base.Fix2(_describe, log)

_limit(x, ::Nothing, limit::Real) = clamp.(x, -limit, limit)
_limit(x, ::typeof(+), limit::Real) = clamp.(x, 0, limit)
_limit(x, ::typeof(log), limit::Real) = clamp.(exp.(x), 0, limit)

_vec(x::Number) = x
_vec(x::AbstractArray) = vec(x)


"""
    grad(obj, param; method=LBFGS(), options...)

Wrapper for ForwardDiff.gradient which behaves like maximise... calling
`flatten(params)` to handle NamedTuples of parameters.

Forms a NamedTuple of the result, without imposing clamp / sign / log transformations of the forward reconstruction.

# Examples:
```
julia> grad(nt -> sum(nt.x ./ nt.y),  (;x = [1,2], y=[4,5]))
(objective = 0.65, x = [0.25, 0.2], y = [-0.0625, -0.08], positive = (x = +, y = +), length = 4, calls = 1)

julia> grad(nt -> sum(nt.x ./ nt.y),  (;x = [1,2], y=[4,5]); log=true)
(objective = 0.65, x = [0.25, 0.4000000000000001], y = [-0.25, -0.4000000000000001], positive = (x = log, y = log), length = 4, calls = 1)
```
"""
function grad(obj, param::NamedTuple; kw...)
    v, re, positive = flatten(param; kw...)
    calls = Ref(0)
    cnt(x) = begin calls[]+=1; x end
    y = _ensure_nt(obj(param), :objective)
    dv = ForwardDiff.gradient(first∘cnt∘obj∘re, v)
    (; y..., re(dv, map(_ -> nothing, positive))..., positive, length=length(v), calls=calls[])
end


function optim_autosymbol(sy::Symbol)
    sy == :forward && return :forward  # Optim.jl understands this
    sy in (:reverse, :zygote) && return ADTypes.AutoZygote()
    error("not sure what to do with this: $sy")
end
function optim_autosymbol(m::Module)
    m === ForwardDiff && return :forward
    m === Zygote && return ADTypes.AutoZygote()
    error("not sure what to do with this: $m")
end
