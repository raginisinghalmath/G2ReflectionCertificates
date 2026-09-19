module G2CapdCertificates

using IntervalArithmetic
using LinearAlgebra
using SHA

IntervalArithmetic.configure(;rounding=:correct,matmul=:slow)

include("ValidatedArithmetic.jl")
include("G2Equations.jl")
include("SingularStart.jl")
include("IntervalKrawczyk.jl")
include("CapdFlow.jl")

using .ValidatedArithmetic
using .G2Equations
using .SingularStart
using .IntervalKrawczyk
using .CapdFlow

export certify_capd, build_capd_driver

function print_environment(kind)
    println("G2 CAPD C1 certificate: ",kind)
    println("Julia version: ",VERSION)
    println("IntervalArithmetic version: ",pkgversion(IntervalArithmetic))
    println("BigFloat/CAPD precision: ",precision(BigFloat)," bits")
    for name in ("Project.toml","Manifest.toml")
        path=normpath(joinpath(@__DIR__,"..",name))
        println("SHA256 ",name," ",bytes2hex(sha256(read(path))))
    end
    for directory in ("src","deps")
        root=normpath(joinpath(@__DIR__,"..",directory))
        for name in sort(readdir(root))
            path=joinpath(root,name)
            isfile(path) || continue
            println("SHA256 ",directory,"/",name," ",
                    bytes2hex(sha256(read(path))))
        end
    end
    println("Flow engine: CAPD MpIOdeSolver with MpC1Rect2Set.")
    flush(stdout)
end

"""
Certify the scaled shooting map using the Julia Fuchsian start, one CAPD C1
flow from delta to the uncertain terminal time, and strict interval Krawczyk.
"""
function certify_capd(kind; precision_bits=256,start_order=8,
                      start_delta=BigInt(1)//BigInt(1024),radius_scale=1,
                      capd_order=24,centre_tolerance="1e-32",
                      box_tolerance="1e-24")
    setprecision(BigFloat,precision_bits) do
        print_environment(kind)
        lam=lambda_value()
        d,k=dimensions(kind)
        Pc=parameters(kind;center=true,radius_scale=radius_scale)
        PX=parameters(kind;radius_scale=radius_scale)
        require(length(matching(fill(iv(0),d),kind))==k+1,
                "Shooting map and parameter dimensions differ")
        println("Parameter centre = ",Pc.center)
        println("Parameter box = ",PX.box)
        println("Singular start: delta=",start_delta," order=",start_order)

        println("Computing centre Fuchsian enclosure")
        sc=singular_start(Pc.box[1:k],Pc.radii[1:k],lam,kind;
                          N=start_order,delta=start_delta)
        println("Computing centre CAPD C1 flow")
        fc=capd_flow(sc,Pc.box[end],lam,kind;order=capd_order,
                     tolerance=centre_tolerance,bits=precision_bits)
        F0=matching(fc.q,kind)

        println("Computing full-box Fuchsian enclosure")
        sx=singular_start(PX.box[1:k],PX.radii[1:k],lam,kind;
                          N=start_order,delta=start_delta)
        println("Computing full-box CAPD C1 flow")
        fx=capd_flow(sx,PX.box[end],lam,kind;order=capd_order,
                     tolerance=box_tolerance,bits=precision_bits)

        initial_derivative=fill(iv(0),d+1,k+1)
        initial_derivative[1:d,1:k]=sx.V
        initial_derivative[d+1,k+1]=PX.radii[end]
        endpoint_derivative=fx.M[1:d,:]*initial_derivative
        checked(endpoint_derivative)
        J=hcat([matching(endpoint_derivative[:,j],kind)
                for j in 1:k+1]...)
        certificate=certify_krawczyk(F0,J)
        println("SUCCESS: strict Krawczyk inclusion certified with CAPD C1 flow.")
        flush(stdout)
        (certificate=certificate,centre_flow=fc,box_flow=fx,
         endpoint_derivative=endpoint_derivative)
    end
end

end
