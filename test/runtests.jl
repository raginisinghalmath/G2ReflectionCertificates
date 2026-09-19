using Test
using IntervalArithmetic
using LinearAlgebra

include(joinpath(@__DIR__,"..","src","G2CapdCertificates.jl"))
include(joinpath(@__DIR__,"..","src","ApproximateShooting.jl"))
include(joinpath(@__DIR__,"..","src","ApproximateS7Gluing.jl"))
using .G2CapdCertificates.ValidatedArithmetic
using .G2CapdCertificates.G2Equations
using .G2CapdCertificates.CapdFlow
using .G2CapdCertificates.IntervalKrawczyk

containszero(x)=inf(x)<=0<=sup(x)
encloses(x,r)=issubset_interval(iv(r),checked(x))

setprecision(BigFloat,256) do
    @testset "Interval arithmetic" begin
        @test encloses(iv("0.1"),1//10)
        @test_throws ValidationFailure checked(sqrt(iv(-1,1)))
        A=iv.([2 1;1 3]); b=iv.([1,2])
        x=verified_solve(A,b)
        @test encloses(x[1],1//5) && encloses(x[2],3//5)
        @test positive_definite(A)>0
    end

    @testset "G2 equations and U(1) reduction" begin
        lam=lambda_value()
        for kind in (:tau4,:tau23_u1)
            p=parameters(kind;center=true).center[1:end-1]
            b,v=initial_pair(p,lam,kind)
            @test all(containszero,regularized(b,p[1],p[2],iv(0),lam,kind))
            L=jacobian(y->regularized(y,p[1],p[2],iv(0),lam,kind),b)
            @test all(containszero,L*v-v)
        end
        u=iv.([1//5,3//100,-1//10,-3//10,2//25,2//5])
        f=vectorfield(lift(u,:tau23_u1),lam,:tau4)
        @test containszero(f[2]-f[3])
        @test containszero(f[6]-f[7])
        @test all(containszero,
                  select(f,:tau23_u1)-vectorfield(u,lam,:tau23_u1))
    end

    @testset "CAPD field construction" begin
        f4=capd_field(:tau4)
        f23=capd_field(:tau23_u1)
        @test startswith(f4,"time:s;par:lam,delta;var:")
        @test startswith(f23,"time:s;par:lam,delta;var:")
        @test occursin("q1,q2,q3,q4,q5,q6,q7,q8,tt",f4)
        @test occursin("q1,q2,q4,q5,q6,q8,tt",f23)
        @test count(==(','),split(f4,";fun:")[2])==8
        @test count(==(','),split(f23,";fun:")[2])==6
        @test endswith(f4,",0;") && endswith(f23,",0;")
    end

    @testset "Strict Krawczyk acceptance and rejection" begin
        J=iv.([2 1;1 3])
        @test certify_krawczyk(iv.([1//100,-1//100]),J).margin>0
        @test_throws ValidationFailure certify_krawczyk(iv.([100,100]),J)
    end

    @testset "Approximate shooting-map smoke tests" begin
        z4=[0.042329982765393195,-0.06469940434457604,
            0.006759263945977987,0.7656914313714904]
        z23=[0.025880416406250467,-0.07764124921873236,
             0.0,0.8578017401911906]
        shot4=ApproximateShooting.approximate_shooting(:tau4,z4)
        shot23=ApproximateShooting.approximate_shooting(:tau23,z23)
        @test norm(shot4.residual,Inf)<1e-8
        @test norm(shot23.residual,Inf)<1e-8
        @test all(>(0),shot4.chamber)
        @test all(>(0),shot23.chamber)
    end

    @testset "Approximate S7 gluing smoke tests" begin
        seed=ApproximateS7Gluing.squashed_s7_seed()
        transformed=ApproximateS7Gluing.pack_s7_parameters(
            seed,1/64,3e-5,
        )
        recovered=ApproximateS7Gluing.unpack_s7_parameters(
            transformed,1/64,3e-5,
        )
        @test isapprox(recovered,seed;rtol=2e-15,atol=2e-15)

        xi_start=ApproximateS7Gluing.xi_singular_start(
            seed[5:8],3e-5,
        )
        @test all(>(0),ApproximateShooting.chamber_margins(xi_start,:tau4))

        match=ApproximateS7Gluing.approximate_s7_matching(seed)
        @test norm(match.residual,Inf)<1e-5
        @test all(>(0),match.eta_chamber)
        @test all(>(0),match.xi_chamber)
    end
end

println("Julia tests passed. Run scripts/test_capd_setup.jl next.")
