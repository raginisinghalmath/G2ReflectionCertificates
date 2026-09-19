module IntervalKrawczyk

using IntervalArithmetic
using LinearAlgebra
using ..ValidatedArithmetic

export certify_krawczyk

"""F0 encloses F(0); J encloses DF([-1,1]^n), in SCALED coordinates."""
function certify_krawczyk(F0,J)
    checked(F0); checked(J); n=length(F0)
    require(size(J)==(n,n),"Shooting Jacobian has the wrong dimensions")
    A=pointmatrix(inv(Float64.(mid.(J))))
    E=eye(n)-A*J
    q=norminf(E)
    # This also proves the fixed preconditioner A is nonsingular.
    require(q<1,"Krawczyk contraction/nonsingularity test failed")
    K=-A*F0+E*fill(iv(-1,1),n)
    checked(K)
    margin=minimum(inf(iv(1)-iv(magnitude(x))) for x in K)
    println("F(0) enclosure = ",F0)
    println("Scaled Jacobian enclosure = ")
    show(stdout,"text/plain",J); println()
    println("Fixed binary preconditioner = ")
    show(stdout,"text/plain",A); println()
    println("||I-AJ||_infinity <= ",q)
    println("Krawczyk image = ",K)
    println("Strict inclusion margin >= ",margin)
    require(margin>0,"Krawczyk image is not strictly inside the unit cube")
    (K=K,A=A,J=J,F0=F0,contraction=q,margin=margin)
end

end
