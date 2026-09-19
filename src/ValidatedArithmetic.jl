module ValidatedArithmetic

using IntervalArithmetic
using LinearAlgebra

export IV, iv, checked, upper, magnitude, norminf, symbox, eye, pointmatrix
export ValidationFailure, require, TJet, constantjet, timejet, coefficient
export pown, polynomial_value, jacobian, verified_solve, verified_inverse, positive_definite
export matmul, checked_hull, checked_intersection

const IV = Interval{BigFloat}
iv(x::AbstractString) = parse(IV, x)
iv(x::Real) = interval(BigFloat, x)
iv(x::IV) = x
iv(a::Real, b::Real) = interval(BigFloat, a, b)

struct ValidationFailure <: Exception
    message::String
end
Base.showerror(io::IO, e::ValidationFailure) = print(io, e.message)
require(ok, message) = ok || throw(ValidationFailure(message))

function checked(x::IV)
    require(isguaranteed(x) && decoration(x) == com &&
            isfinite(inf(x)) && isfinite(sup(x)) && inf(x) <= sup(x),
            "Invalid, non-common, unbounded, or non-guaranteed interval: $x")
    x
end
checked(xs::AbstractArray) = (foreach(checked, xs); xs)
# IEEE set operations default to trv even for valid inputs. :auto inherits
# the input decorations; invalid inputs are rejected BEFORE the operation.
# This never upgrades a failed computation's decoration or guarantee flag.
checked_hull(x::IV,y::IV) = checked(hull(checked(x),checked(y);dec=:auto))
checked_intersection(x::IV,y::IV) =
    checked(intersect_interval(checked(x),checked(y);dec=:auto))
upper(x::IV) = sup(checked(x))
magnitude(x::IV) = (checked(x); max(abs(inf(x)), abs(sup(x))))
symbox(r::Real) = iv(-r, r)
eye(n) = [iv(i == j ? 1 : 0) for i in 1:n, j in 1:n]
pointmatrix(A) = [iv(A[i,j]) for i in axes(A,1), j in axes(A,2)]
norminf(A::AbstractMatrix) = maximum(upper(sum(iv(magnitude(A[i,j]))
    for j in axes(A,2))) for i in axes(A,1))
norminf(x::AbstractVector) = maximum(magnitude, x)

# Also works for jets, without requiring Julia's Number promotion machinery.
function matmul(A::AbstractMatrix,B::AbstractMatrix)
    size(A,2)==size(B,1) || throw(DimensionMismatch("Matrix product"))
    [sum(A[i,k]*B[k,j] for k in axes(A,2)) for i in axes(A,1),j in axes(B,2)]
end

# Finite time derivative jets. These are NOT Taylor models and carry no
# remainder. They are used only to enclose finitely many derivatives;
# every use as a function approximation adds a separate Taylor remainder.
struct TJet
    c::Vector{IV}                 # c[k+1] encloses the kth derivative / k!
end
constantjet(x::IV, n) = TJet(vcat([x], fill(iv(0), n)))
timejet(x::IV, n) = TJet(vcat([x, iv(1)], fill(iv(0), n-1)))
coefficient(x::TJet, k) = x.c[k+1]
Base.zero(x::TJet) = constantjet(iv(0), length(x.c)-1)
Base.one(x::TJet) = constantjet(iv(1), length(x.c)-1)
function same_order(a::TJet, b::TJet)
    length(a.c) == length(b.c) || throw(DimensionMismatch("Time-jet orders differ"))
end
Base.:+(a::TJet,b::TJet) = (same_order(a,b); TJet(a.c+b.c))
Base.:-(a::TJet) = TJet(-a.c)
Base.:-(a::TJet,b::TJet) = a+(-b)
Base.:+(a::TJet,b::Real) = a+constantjet(iv(b),length(a.c)-1)
Base.:+(a::Real,b::TJet) = b+a
Base.:-(a::TJet,b::Real) = a+(-iv(b))
Base.:-(a::Real,b::TJet) = iv(a)+(-b)
Base.:*(a::TJet,b::TJet) = (same_order(a,b);
    TJet([sum(a.c[j+1]*b.c[k-j+1] for j in 0:k) for k in 0:length(a.c)-1]))
Base.:*(a::TJet,b::Real) = TJet(a.c .* iv(b))
Base.:*(a::Real,b::TJet) = b*a
function Base.inv(a::TJet)
    checked(a.c[1])
    require(inf(a.c[1]) > 0 || sup(a.c[1]) < 0, "Zero in reciprocal jet")
    b=fill(iv(0),length(a.c)); b[1]=inv(a.c[1])
    for k in 1:length(b)-1
        b[k+1]=-sum(a.c[j+1]*b[k-j+1] for j in 1:k)/a.c[1]
    end
    TJet(checked(b))
end
Base.:/(a::TJet,b::TJet) = a*inv(b)
Base.:/(a::TJet,b::Real) = a*inv(iv(b))
Base.:/(a::Real,b::TJet) = iv(a)*inv(b)
function Base.sqrt(a::TJet)
    require(inf(checked(a.c[1])) > 0, "Square-root jet is not strictly positive")
    b=fill(iv(0),length(a.c)); b[1]=sqrt(a.c[1])
    for k in 1:length(b)-1
        s=iv(0)
        for j in 1:k-1
            s+=b[j+1]*b[k-j+1]
        end
        b[k+1]=(a.c[k+1]-s)/(iv(2)*b[1])
    end
    TJet(checked(b))
end
function Base.:^(a::TJet, n::Integer)
    n<0 && return inv(a)^(-n)
    b=one(a)
    for _ in 1:n
        b=b*a
    end
    b
end
function polynomial_value(p::TJet,x)
    # Lift the leading coefficient explicitly, without a mixed-type zero.
    # constantlike has methods for intervals and for time/derivative jets.
    y=constantlike(x,p.c[end])
    for k in length(p.c)-1:-1:1
        y=y*x+p.c[k]
    end
    y
end

# First derivative jets, nestable to obtain Hessians. Their scalar values
# may themselves be interval time jets. No finite differences are used.
struct DJet{T}
    v::T
    d::Vector{T}
end
depth(x) = 0
depth(x::DJet) = 1+depth(x.v)
constantlike(x::IV, b::Real) = iv(b)
constantlike(x::TJet, b::Real) = constantjet(iv(b),length(x.c)-1)
constantlike(x::TJet, b::TJet) = b
constantlike(x::DJet, b) = DJet(constantlike(x.v,b), [zero(t) for t in x.d])
Base.zero(x::DJet) = DJet(zero(x.v),[zero(t) for t in x.d])
Base.one(x::DJet) = DJet(one(x.v),[zero(t) for t in x.d])
function align(a::DJet,b::DJet)
    if depth(a)>depth(b)
        return a,constantlike(a,b)
    elseif depth(a)<depth(b)
        return constantlike(b,a),b
    end
    length(a.d)==length(b.d) || throw(DimensionMismatch("Derivative-jet dimensions differ"))
    a,b
end
constantlike(x::DJet{T}, b::DJet{T}) where T = b
function Base.:+(a::DJet,b::DJet)
    a,b=align(a,b); DJet(a.v+b.v,[a.d[j]+b.d[j] for j in eachindex(a.d)])
end
Base.:-(a::DJet) = DJet(-a.v,[-x for x in a.d])
Base.:-(a::DJet,b::DJet) = a+(-b)
function Base.:*(a::DJet,b::DJet)
    a,b=align(a,b)
    DJet(a.v*b.v,[a.d[j]*b.v+a.v*b.d[j] for j in eachindex(a.d)])
end
Base.inv(a::DJet) = DJet(inv(a.v),[-x/(a.v*a.v) for x in a.d])
Base.:/(a::DJet,b::DJet) = a*inv(b)
function Base.sqrt(a::DJet)
    b=sqrt(a.v); DJet(b,[x/(b+b) for x in a.d])
end
for op in (:+,:-,:*,:/)
    @eval Base.$op(a::DJet,b::Union{Real,TJet}) = Base.$op(a,constantlike(a,b))
    @eval Base.$op(a::Union{Real,TJet},b::DJet) = Base.$op(constantlike(b,a),b)
end
function Base.:^(a::DJet,n::Integer)
    n<0 && return inv(a)^(-n)
    b=one(a)
    for _ in 1:n
        b=b*a
    end
    b
end
function jacobian(f,x)
    n=length(x)
    z=[DJet(x[i],[i==j ? one(x[i]) : zero(x[i]) for j in 1:n]) for i in 1:n]
    y=f(z)
    [y[i] isa DJet ? y[i].d[j] : zero(x[1]) for i in eachindex(y),j in 1:n]
end

# Floating matrices below only propose preconditioners. All acceptance
# inequalities and returned bounds use outward-rounded intervals.
function verified_solve(A::AbstractMatrix,b::AbstractVector)
    checked(A); checked(b); n=length(b)
    C=pointmatrix(inv(Float64.(mid.(A))))
    E=eye(n)-C*A; q=norminf(E)
    require(q<1,"Interval linear system: Neumann test failed")
    x=iv.(mid.(C*b))
    r=C*(b-A*x)
    e=upper(iv(norminf(r))/(iv(1)-iv(q)))
    X=x .+ symbox(e)
    for _ in 1:4
        Z=x+r+E*(X-x)
        X=checked_intersection.(X,Z)
    end
    checked(X)
end
function verified_inverse(A)
    n=size(A,1); E=eye(n)
    hcat([verified_solve(A,E[:,j]) for j in 1:n]...)
end
function positive_definite(A)
    checked(A); n=size(A,1); L=eye(n); d=fill(iv(0),n)
    for j in 1:n
        d[j]=A[j,j]
        for k in 1:j-1
            d[j]-=L[j,k]^2*d[k]
        end
        require(inf(checked(d[j]))>0,"Interval LDL positivity test failed at pivot $j")
        for i in j+1:n
            a=A[i,j]
            for k in 1:j-1
                a-=L[i,k]*L[j,k]*d[k]
            end
            L[i,j]=a/d[j]
        end
    end
    minimum(inf.(d))
end

end
