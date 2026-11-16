# Checks for the Euclidean kinked-Wilson-line model of ordered.tex
# (Definition def:kinkmodel, Lemmas lem:density, lem:taildecay,
#  Theorems thm:kinkexact, thm:kinkhalf, thm:running).
# Run: julia check_kinked_line.jl     (exits nonzero on failure)

include(joinpath(@__DIR__, "quad.jl"))
using LinearAlgebra
using Printf

failures = 0
function check(name, ok)
    global failures
    println(ok ? "PASS  " : "FAIL  ", name)
    ok || (failures += 1)
end

# ---------------- model data (R^4, Euclidean) ----------------
const p  = [1.0, 0.0, 0.0, 0.0]
const us = [normalize([cos(1.1), sin(1.1), 0.0, 0.0]),
            normalize([cos(2.0), 0.0, sin(2.0), 0.0]),
            normalize([cos(2.6), 0.3, 0.0, sin(2.6)])]
const cs = [0.7, -1.3, 0.6]          # charge products e_a e_b (arbitrary reals)
const ws = cs ./ (4pi^2)              # weights including the 1/(4π^2) of the propagator
const q  = [0.35, 0.5, -0.2, 0.15]
const c0 = 1.0                       # emission rate: weight exp(-c0 ω s)

# inner sigma-integral, closed form:  ∫_0^∞ dσ (ẋ·u)/|x-σu|^2
function Isig(x, xd, u)
    a = dot(x, u); ρ = sqrt(max(dot(x, x) - a^2, 0.0))
    return dot(xd, u) / ρ * (pi / 2 + atan(a / ρ))
end
# same with log kernel: ∫_0^∞ dσ (ẋ·u) ln|x-σu| / |x-σu|^2
function IsigD(x, xd, u)
    a = dot(x, u); ρ = sqrt(max(dot(x, x) - a^2, 0.0))
    φ0 = atan(-a / ρ); b = pi / 2 - φ0        # ∫_{φ0}^{π/2} ln cos φ dφ = ∫_0^b ln sin ψ dψ
    lnsin_int = (b * log(b) - b) + integrate(ψ -> ψ == 0 ? 0.0 : log(sin(ψ) / ψ), [0.0, b / 2, b])
    return dot(xd, u) / ρ * (b * log(ρ) - lnsin_int)
end
Cf(v, u) = Isig(v, v, u)
Df(v, u) = IsigD(v, v, u)

# check closed forms against direct numerical σ-integration
function Isig_num(x, xd, u; logk=false)
    f(σ) = begin
        r = x - σ * u; d2 = dot(r, r)
        dot(xd, u) * (logk ? 0.5 * log(d2) : 1.0) / d2
    end
    integrate(f, vcat(collect(0:0.5:20), [50.0, 200.0, 1000.0])) + integrate_tail(f, 1000.0)
end
let x = [0.8, -0.3, 0.5, 0.2], xd = [0.4, 0.9, -0.1, 0.3]
    for u in us
        check(@sprintf("closed-form inner integral (u·p=%.2f)", dot(u, p)),
              isapprox(Isig(x, xd, u), Isig_num(x, xd, u); rtol=1e-9))
        check("closed-form inner log integral",
              isapprox(IsigD(x, xd, u), Isig_num(x, xd, u; logk=true); rtol=1e-8, atol=1e-10))
    end
end
# Euclidean cusp function: C(v,u) = cotθ (π-θ) for unit vectors
for u in us
    θ = acos(dot(p, u))
    check("cusp formula C = cot(θ)(π-θ)", isapprox(Cf(p, u), cot(θ) * (pi - θ); rtol=1e-12))
    check("C scale invariant in v", isapprox(Cf(3.7 * p, u), Cf(p, u); rtol=1e-12))
    check("D(λv)=C(v)lnλ+D(v)", isapprox(Df(3.7 * p, u), Cf(p, u) * log(3.7) + Df(p, u); rtol=1e-9))
end

# ---------------- kinked line ----------------
Pv(ω) = p .+ ω .* q
Ctot(v) = sum(ws[b] * Cf(v, us[b]) for b in eachindex(us))
α(ω) = Ctot(Pv(ω)) - Ctot(p)
# tail functions
hC(λ, ω, u) = Isig(Pv(ω) .+ (λ - 1) .* p, p, u)
hD(λ, ω, u) = IsigD(Pv(ω) .+ (λ - 1) .* p, p, u)
τC(ω, u) = integrate_tail(λ -> hC(λ, ω, u) - Cf(p, u) / λ, 1.0)
τD(ω, u) = integrate_tail(λ -> hD(λ, ω, u) - (Cf(p, u) * log(λ) + Df(p, u)) / λ, 1.0)
tt(ω) = sum(ws[b] * τC(ω, us[b]) for b in eachindex(us))

# Lemma lem:taildecay: λ^2 |h(λ)-C(p)/λ| ≤ K ω
let ω = 0.05
    vals = [λ^2 * abs(hC(λ, ω, us[1]) - Cf(p, us[1]) / λ) / ω for λ in exp.(range(0, log(1e6), length=60))]
    check(@sprintf("tail decay O(ω/λ^2): sup λ^2|h-C/λ|/ω = %.3f bounded", maximum(vals)),
          maximum(vals) < 10 && vals[end] < 2 * vals[end-10])
end

# direct evaluation of ΔΦ(s) from the τ-integral (Λ = 1), for testing scale covariance
function ΔΦ_direct(s, ω; Λ=1.0)
    tot = 0.0
    for b in eachindex(us)
        u = us[b]; P = Pv(ω)
        dens(τ) = τ < s ? Isig(τ .* P, P, u) : Isig(s .* P .+ (τ - s) .* p, p, u)
        f(τ) = dens(τ) - Cf(p, u) / τ
        lo = 1 / Λ
        if s > lo
            v = integrate(f, exp.(range(log(lo), log(s), length=30)))
            v += integrate_tail(f, s)
        else
            v = integrate_tail(f, lo)
        end
        tot += ws[b] * v
    end
    tot
end
ΔΦ_formula(s, ω; Λ=1.0) = s ≥ 1 / Λ ? α(ω) * log(Λ * s) + tt(ω) :
    sum(ws[b] * integrate_tail(λ -> hC(λ, ω, us[b]) - Cf(p, us[b]) / λ, 1 / (Λ * s)) for b in eachindex(us))
let ω = 0.07
    for s in (0.3, 1.0, 4.0, 50.0, 3000.0)
        check(@sprintf("ΔΦ(s) = α ln(Λs)+τ (scale covariance), s=%g", s),
              isapprox(ΔΦ_direct(s, ω), ΔΦ_formula(s, ω); rtol=1e-8, atol=1e-12))
    end
end

# Theorem thm:kinkexact: F(ω) = e^τ Γ(1+α) (Λ/β)^α + G,  |G| ≤ K ω^2
function Fnum(ω; Λ=1.0, g=1.0)
    β = c0 * ω
    αω = α(ω); tω = tt(ω)
    Φ(s) = s ≥ 1 / Λ ? αω * log(Λ * s) + tω : ΔΦ_formula(s, ω; Λ=Λ)
    y0 = β / Λ
    f(y) = exp(-y) * exp(g * Φ(y / β))
    integrate(f, vcat(range(0, y0, length=8), exp.(range(log(y0), log(60.0), length=60))[2:end]))
end
Fexact_main(ω; Λ=1.0, g=1.0) = exp(g * tt(ω)) * gammafn(1 + g * α(ω)) * (Λ / (c0 * ω))^(g * α(ω))
Gs = Float64[]
for ω in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
    G = Fnum(ω) - Fexact_main(ω)
    push!(Gs, abs(G) / ω^2)
    @printf("   ω=%-6g  F=%.12f  G/ω^2=%.5f\n", ω, Fnum(ω), G / ω^2)
end
check("G(ω) = O(ω^2)", maximum(Gs) < 5 * Gs[end] + 1 && maximum(Gs) < 10)

# Expansion coefficients A, A2 (α = ωA + ω^2 A2 + O(ω^3)) by Richardson differences
hstep = 1e-3
A  = (α(hstep) - α(-hstep)) / (2hstep)
A2 = (α(hstep) + α(-hstep)) / (2hstep^2)
@printf("   A = %.10f   A2 = %.10f\n", A, A2)
check("A2 ≠ 0 for this configuration (running counterexample is nontrivial)", abs(A2) > 1e-4)

# Theorem thm:kinkhalf: at fixed ω the s≥1/Λ part of [g^2]F is an exact quadratic in
# T' = ln(Λ/β); its T'^2 coefficient is α(ω)^2/2, and α(ω)^2/(2ω^2) → A^2/2.
function g2F_num(ω; Λ=1.0, β0=0.0)
    β = c0 * ω
    αω = α(ω); tω = tt(ω)
    ΔD = sum(ws[b] * (Df(Pv(ω), us[b]) - Df(p, us[b])) for b in eachindex(us))
    dω = ΔD + tω
    tDω = sum(ws[b] * τD(ω, us[b]) for b in eachindex(us))
    function integrand(s)
        if s ≥ 1 / Λ
            L = log(Λ * s)
            Ψ1 = αω * L + tω
            Ψ2 = -(αω * L^2 / 2 + dω * L + tDω)
        else
            λ0 = 1 / (Λ * s); Ls = log(Λ * s)
            Ψ1 = ΔΦ_formula(s, ω; Λ=Λ)
            Ψ2 = -sum(ws[b] * integrate_tail(λ -> Ls * (hC(λ, ω, us[b]) - Cf(p, us[b]) / λ) +
                       hD(λ, ω, us[b]) - (Cf(p, us[b]) * log(λ) + Df(p, us[b])) / λ, λ0) for b in eachindex(us))
        end
        Ψ1^2 / 2 + β0 * Ψ2
    end
    y0 = β / Λ
    f(y) = exp(-y) * integrand(y / β)
    integrate(f, vcat(range(0, y0, length=6), exp.(range(log(y0), log(60.0), length=50))[2:end]))
end

# Direct check of the log-kernel tail being affine in ln(Λs) (used in thm:running)
function Ψ2_direct(s, ω)
    tot = 0.0
    for b in eachindex(us)
        u = us[b]; P = Pv(ω)
        densD(τ) = τ < s ? IsigD(τ .* P, P, u) : IsigD(s .* P .+ (τ - s) .* p, p, u)
        dens(τ)  = τ < s ? Isig(τ .* P, P, u)  : Isig(s .* P .+ (τ - s) .* p, p, u)
        # kernel ln(Λ|x-σu|) with Λ=1
        f(τ) = (densD(τ)) - (Cf(p, u) * log(τ) + Df(p, u)) / τ
        v = integrate(f, exp.(range(0.0, log(s), length=30))) + integrate_tail(f, s)
        tot += ws[b] * v
    end
    -tot
end
let ω = 0.06
    αω = α(ω); tω = tt(ω)
    dω = sum(ws[b] * (Df(Pv(ω), us[b]) - Df(p, us[b])) for b in eachindex(us)) + tω
    tDω = sum(ws[b] * τD(ω, us[b]) for b in eachindex(us))
    for s in (1.0, 7.0, 400.0)
        L = log(s)
        check(@sprintf("Ψ2(s) = -(αL^2/2 + dL + t_D), s=%g", s),
              isapprox(Ψ2_direct(s, ω), -(αω * L^2 / 2 + dω * L + tDω); rtol=1e-7, atol=1e-11))
    end
end

# Fit the T'^2 coefficient at fixed ω by varying Λ (large Λ makes the s<1/Λ part O(ω/Λ)).
function T2coef(ω; β0=0.0)
    Λs = [1e3, 1e5, 1e7, 1e9]
    Tp = [log(Λ / (c0 * ω)) for Λ in Λs]
    y  = [g2F_num(ω; Λ=Λ, β0=β0) for Λ in Λs]
    # least squares quadratic
    M = hcat(ones(length(Tp)), Tp, Tp .^ 2)
    coef = M \ y
    coef[3]
end
println("   [g^2]F: T'^2 coefficient κ2(ω) vs α^2/2 - β0 α/2")
for β0 in (0.0, 0.8)
    rat = Float64[]
    for ω in (0.04, 0.02, 0.01)
        κ2 = T2coef(ω; β0=β0)
        pred = α(ω)^2 / 2 - β0 * α(ω) / 2
        @printf("   β0=%.1f ω=%-5g κ2=% .10e  pred=% .10e\n", β0, ω, κ2, pred)
        check(@sprintf("T'^2 coefficient matches α^2/2-β0α/2 (β0=%.1f, ω=%g)", β0, ω),
              isapprox(κ2, pred; rtol=1e-6, atol=1e-12))
        push!(rat, (κ2 + β0 * A * ω / 2) / ω^2)
    end
    lim = A^2 / 2 - β0 * A2 / 2
    @printf("   (κ2+β0Aω/2)/ω^2 → %.6f ... %.6f ; predicted limit A^2/2-β0 A2/2 = %.6f\n", rat[1], rat[end], lim)
    check(@sprintf("ω→0 limit of ω^2T^2 coefficient (β0=%.1f)", β0),
          abs(rat[end] - lim) < 0.5 * abs(rat[1] - lim) + 1e-9 || abs(rat[end] - lim) < 2e-2 * abs(lim))
end

# Moments of the exponential law used in the proofs
let T = 7.3
    m1 = integrate(y -> exp(-y) * (T + log(y)), vcat([0.0], exp.(range(log(1e-12), log(80.0), length=80))))
    m2 = integrate(y -> exp(-y) * (T + log(y))^2, vcat([0.0], exp.(range(log(1e-12), log(80.0), length=80))))
    check("E[T+lnY] = T-γ", isapprox(m1, T - EULERGAMMA; rtol=1e-9))
    check("E[(T+lnY)^2] = (T-γ)^2+π^2/6", isapprox(m2, (T - EULERGAMMA)^2 + pi^2 / 6; rtol=1e-9))
end

println(failures == 0 ? "ALL CHECKS PASSED" : "$failures CHECK(S) FAILED")
exit(failures == 0 ? 0 : 1)
