module SingularStart

using IntervalArithmetic
using LinearAlgebra
using ..ValidatedArithmetic
using ..G2Equations

export singular_start, augmented_G, initial_augmented, formal_coefficients

# The sensitivity columns are derivatives with respect to scaled parameters
# xi_j, p_j=p_j^0+r_j xi_j. The terminal time is not a singular-start parameter.
function augmented_G(z,p,r,t,lam,kind)
    d,k=dimensions(kind); y=z[1:d]; V=reshape(z[d+1:end],d,k)
    x=vcat(y,[zero(y[1]) for _ in 1:k])
    fun=x->regularized(x[1:d],p[1]+r[1]*x[d+1],p[2]+r[2]*x[d+2],t,lam,kind)
    J=jacobian(fun,x)
    g=regularized(y,p[1],p[2],t,lam,kind)
    dv=matmul(J[:,1:d],V)+J[:,d+1:end]
    vcat(g,vec(dv))
end
function initial_augmented(p,r,lam,kind)
    d,k=dimensions(kind)
    b,v=initial_pair(p,lam,kind)
    J=jacobian(x->vcat(initial_pair(x,lam,kind)...),p)
    Z0=vcat(b,vec([J[i,j]*r[j] for i in 1:d,j in 1:k]))
    Z1=vcat(v,vec([J[d+i,j]*r[j] for i in 1:d,j in 1:k]))
    checked(Z0); checked(Z1)
    Z0,Z1
end
function formal_coefficients(p,r,lam,kind,N)
    Z0,Z1=initial_augmented(p,r,lam,kind); n=length(Z0)
    A0=jacobian(z->augmented_G(z,p,r,iv(0),lam,kind),Z0)
    checked(A0)
    C=hcat(Z0,Z1,fill(iv(0),n,N-1))
    # Sanity checks, not proofs by interval overlap. The exact identities
    # are algebraic; proof.tex gives their radical reduction and explains
    # how the verified stability test then supplies analytic existence.
    g0=augmented_G(Z0,p,r,iv(0),lam,kind)
    require(all(x->inf(x)<=0<=sup(x),g0),"Initial equilibrium identity excluded")
    require(all(x->inf(x)<=0<=sup(x),A0*Z1-Z1),"Resonant initial identity excluded")
    for m in 2:N
        t=timejet(iv(0),m)
        z=[TJet(C[i,1:m+1]) for i in 1:n]
        F=augmented_G(z,p,r,t,lam,kind)
        b=[coefficient(F[i],m) for i in 1:n]
        C[:,m+1]=verified_solve(iv(m)*eye(n)-A0,b)
    end
    checked(C); C,A0
end

function evaluate_polynomial(C,t)
    [polynomial_value(TJet(collect(C[i,:])),t) for i in axes(C,1)]
end

"""
Prove a singular-start enclosure, INCLUDING parameter derivatives.
The exact initial identities and the verified matrix inequality imply local
analytic existence by the regular-singular argument in proof.tex. Smooth
geometric extension still requires the draft's orbit-extension argument.
No start-error, amplification or derivative-error constant is imported.
"""
function singular_start(p,r,lam,kind; N=8, delta=BigInt(1)//BigInt(1024), m=4)
    require(N>=2 && 1<m<N+1,"Need 1 < m < N+1")
    require(inf(p[1])>0 && sup(p[2])<0,"Singular-orbit parameter chamber failed")
    deltaI=iv(delta); require(inf(deltaI)>0 && sup(deltaI)<1,"Bad start time")
    println("Start: computing order ",N+1," coefficient enclosures, including sensitivities")
    flush(stdout)
    C,A0=formal_coefficients(p,r,lam,kind,N+1)
    Cbar=C[:,1:N+1]; nz=size(C,1)
    # Taylor's theorem applied to the FULL analytic residual, not to a
    # truncated residual polynomial: R(t)/t^(N+1) is enclosed by R^(N+1)/(N+1)!.
    t=timejet(iv(0,sup(deltaI)),N+1)
    z=evaluate_polynomial(Cbar,t)
    F=augmented_G(z,p,r,t,lam,kind)
    Cp=copy(Cbar)
    for j in 0:N
        Cp[:,j+1].*=iv(j)
    end
    R=evaluate_polynomial(Cp,t)-F
    Rbound=norminf([coefficient(x,N+1) for x in R])

    # A floating Lyapunov solve proposes one fixed, exactly binary matrix P.
    # Both its positivity and the full-tube matrix inequality are verified below.
    A=Float64.(mid.(A0))-m*Matrix{Float64}(I,nz,nz)
    Pf=lyap(transpose(A),Matrix{Float64}(I,nz,nz))
    Pf=(Pf+transpose(Pf))/2
    P=pointmatrix(Pf)
    ppivot=positive_definite(P)
    Pinv=verified_inverse(P)
    CL=upper(sqrt(iv(nz)*iv(norminf(P))*iv(norminf(Pinv))))
    H0=upper(iv(CL)*iv(Rbound)/iv(N+1-m))
    # The next exact coefficient bounds the limiting error/t^(N+1).
    # This makes the first-exit bootstrap valid from some positive initial time.
    H=upper(iv(2)*(iv(max(H0,norminf(C[:,N+2])))+iv("1e-50")))
    T=iv(0,sup(deltaI))
    zbox=evaluate_polynomial(Cbar,T)
    # N+1 is an exact integer; generic mixed ^ may lose the guarantee flag.
    delta_power=checked(pown(deltaI,N+1))
    tube=zbox .+ symbox(upper(iv(H)*delta_power))
    J=jacobian(z->augmented_G(z,p,r,T,lam,kind),tube)
    B=J-iv(m)*eye(nz)
    energy=-(transpose(B)*P+P*B)
    epivot=positive_definite(energy)
    require(H0<H,"Singular-start first-exit inequality failed")
    endpoint=evaluate_polynomial(Cbar,deltaI) .+
        symbox(upper(iv(H0)*delta_power))
    checked(endpoint)
    d,k=dimensions(kind); y=endpoint[1:d]; V=reshape(endpoint[d+1:end],d,k)
    q=q_from_y(y,p[1],p[2],deltaI,kind)
    x=vcat(y,fill(iv(0),k))
    QJ=jacobian(x->q_from_y(x[1:d],p[1]+r[1]*x[d+1],
                         p[2]+r[2]*x[d+2],deltaI,kind),x)
    Vq=QJ[:,1:d]*V+QJ[:,d+1:end]
    checked(q); checked(Vq); chamber(q,kind)
    println("Start PASS: Rbound=",Rbound," CL=",CL," H0=",H0," H=",H)
    println("Interval LDL lower pivots: P=",ppivot," energy=",epivot)
    flush(stdout)
    (q=q,V=Vq,y=y,Vy=V,delta=delta,coefficients=C,Rbound=Rbound,CL=CL,
     H0=H0,H=H,P=P,energy=energy)
end

end
