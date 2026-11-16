# Check of Lemma lem:projection / Corollary cor:ckpbms in ordered.tex:
# the projected functional  P_i[F] = <∂̄^3 F, φ> = -∫ F ∂̄^3 φ d^2x
# depends on the residue numerator g of F = g/(z-z_i) only modulo (z-z_i),
# and the CKP term S B^2/2 and the Sahoo/BMS form give the same value.
# Run: julia check_projection.jl

include(joinpath(@__DIR__, "quad.jl"))
using Printf

failures = 0
function check(name, ok)
    global failures
    println(ok ? "PASS  " : "FAIL  ", name)
    ok || (failures += 1)
end

# punctures, charges (sum zero), anchor i = 1
zs = ComplexF64[0.1 + 0.2im, 1.3 - 0.4im, -0.9 + 1.1im, 0.4 - 1.5im]
Qs = [1.0, -0.6, 1.7, -2.1]
i  = 1
zi = zs[i]
ε  = 0.35                       # support radius of the test function (< distances)

# test function φ = f(r^2) * m(Z, Zb), Z = z - z_i, f(s) = exp(-1/(1-s/ε^2))
# represent terms  coef * f^{(n)}(r^2) * Z^a * Zb^b   as Dict((n,a,b) => coef)
const Term = Tuple{Int,Int,Int}
function dbar(T::Dict{Term,ComplexF64})
    R = Dict{Term,ComplexF64}()
    for ((n, a, b), c) in T
        # ∂̄ f^{(n)}(Z Zb) = f^{(n+1)} Z ;  ∂̄ Zb^b = b Zb^{b-1}
        R[(n + 1, a + 1, b)] = get(R, (n + 1, a + 1, b), 0) + c
        b > 0 && (R[(n, a, b - 1)] = get(R, (n, a, b - 1), 0) + b * c)
    end
    R
end
# derivatives of g(x) = exp(-1/(1-x)): g^{(n)} = g * Qn(y), y = 1/(1-x), Q_{n+1} = y^2 (Qn' - Qn)
function fderivs(nmax)
    polys = [[1.0]]                     # coefficients in y (ascending)
    for n in 1:nmax
        Pn = polys[end]
        dP = [k * Pn[k+1] for k in 1:length(Pn)-1]
        L = max(length(dP), length(Pn)) + 2
        R = zeros(L)
        for k in 1:length(dP); R[k+2] += dP[k]; end
        for k in 1:length(Pn); R[k+2] -= Pn[k]; end
        push!(polys, R)
    end
    polys
end
const FP = fderivs(6)
function fder(n, s)     # d^n/ds^n exp(-1/(1-s/ε^2))
    x = s / ε^2
    x >= 1 && return 0.0
    y = 1 / (1 - x)
    val = exp(-y) * sum(FP[n+1][k] * y^(k - 1) for k in 1:length(FP[n+1]))
    val / ε^(2n)
end
function evalterm(T, Z)
    s = abs2(Z); Zb = conj(Z)
    sum(c * fder(n, s) * Z^a * Zb^b for ((n, a, b), c) in T)
end

# projection: -∫ F ∂̄^3 φ d^2x in polar coordinates around z_i
function project(F, mono::Dict{Term,ComplexF64})
    D3 = dbar(dbar(dbar(mono)))
    rb = collect(range(0, ε, length=25)); θb = collect(range(0, 2pi, length=33))
    integrate(r -> r * integrate(θ -> begin
        Z = r * cis(θ); -F(zi + Z) * evalterm(D3, Z) end, θb), rb)
end

# CKP:  F = Q_i P1(z̄)^2/2 /(z-z_i),  P1 = Σ_{j≠i} Q_j (z̄-z̄_j)/(z̄_i-z̄_j)
P1(z) = sum(Qs[j] * (conj(z) - conj(zs[j])) / (conj(zi) - conj(zs[j])) for j in eachindex(zs) if j != i)
Rb(z) = sum(Qs[c] * (z - zs[c]) * (conj(z) - conj(zs[c])) / ((zi - zs[c]) * (conj(zi) - conj(zs[c])))
            for c in eachindex(zs) if c != i)
F_ckp(z) = Qs[i] * P1(z)^2 / 2 / (z - zi)
F_bms(z) = Qs[i] * P1(z) * Rb(z) / 2 / (z - zi)

tests = [Dict{Term,ComplexF64}((0, 0, 0) => 1.0),
         Dict{Term,ComplexF64}((0, 0, 0) => 1.0, (0, 1, 0) => 0.7 - 0.2im, (0, 0, 1) => -0.4im,
                               (0, 0, 2) => 1.3, (0, 1, 1) => 0.5, (0, 2, 3) => -0.8 + 0.1im)]
for (k, φ) in enumerate(tests)
    a = project(F_ckp, φ); b = project(F_bms, φ)
    @printf("   test %d: P[CKP] = % .10f %+.10fi   P[BMS] = % .10f %+.10fi\n", k, real(a), imag(a), real(b), imag(b))
    check("projection CKP == BMS (test $k)", isapprox(a, b; rtol=1e-8, atol=1e-12))
    # difference numerator is divisible by (z - z_i): projection of the difference vanishes
    check("projection of (z-z_i)h/(z-z_i) with ∂̄^3 h = 0 vanishes (test $k)",
          abs(project(z -> Qs[i] * P1(z) * (Rb(z) - P1(z)) / 2 / (z - zi), φ)) < 1e-10)
end
# radial test: value π φ(z_i) ∂̄^2 g(z_i) = π e^{-1} Q_i P1'^2
P1p = sum(Qs[j] / (conj(zi) - conj(zs[j])) for j in eachindex(zs) if j != i)
pred = pi * exp(-1) * Qs[i] * P1p^2
a = project(F_ckp, tests[1])
check("radial test: P[F] = π φ(z_i) ∂̄^2 g(z_i)", isapprox(a, pred; rtol=1e-8))
# a single-B term (degree-1 numerator) has a nonzero projection
Fb(z) = Qs[i] * P1(z) / (z - zi)
check("single-B term S B has nonzero projection (non-radial test)", abs(project(Fb, tests[2])) > 1e-3)

println(failures == 0 ? "ALL CHECKS PASSED" : "$failures CHECK(S) FAILED")
exit(failures == 0 ? 0 : 1)
