# verify_sampling.jl -- numerical checks for Section "Certified discrete sampling".
# Standard library only.  Run:  julia verify_sampling.jl
using Printf, LinearAlgebra

const NPASS = Ref(0); const NFAIL = Ref(0)
function check(name, ok; info = "")
    if ok
        NPASS[] += 1; println("PASS  ", name, "   ", info)
    else
        NFAIL[] += 1; println("FAIL  ", name, "   ", info)
    end
end
const TOL = 1e-12            # absolute floating-point allowance in inequality checks

# ---------------------------------------------------------------- special functions
const BERN = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330]
function lgammac(z::ComplexF64)               # log Gamma up to 2*pi*i, Re z > 0
    s = zero(z)
    while real(z) < 20
        s += log(z); z += 1
    end
    r = (z - 0.5) * log(z) - z + 0.5 * log(2pi)
    zi = 1 / z; z2 = zi * zi; p = zi
    for k in 1:length(BERN)
        r += BERN[k] / (2k * (2k - 1)) * p
        p *= z2
    end
    r - s
end
gammac(z) = exp(lgammac(complex(float(z))))
gammar(x::Real) = real(gammac(x))

function simpson(g, lo, hi, n)                # composite Simpson, n even
    hstep = (hi - lo) / n
    s = g(lo) + g(hi)
    for i in 1:n-1
        s += (isodd(i) ? 4 : 2) * g(lo + i * hstep)
    end
    s * hstep / 3
end
erfc_(x) = x < 0 ? 2 - erfc_(-x) : 2 / sqrt(pi) * simpson(s -> exp(-(x + s)^2), 0.0, 12.0, 24000)

function gausslegendre(n)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cos(pi * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
            end
            dp = n * (z * p1 - p2) / (z^2 - 1)
            dz = p1 / dp; z -= dz
            abs(dz) < 1e-16 && break
        end
        p1, p2 = 1.0, 0.0
        for j in 1:n
            p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
        end
        dp = n * (z * p1 - p2) / (z^2 - 1)
        x[i] = z; w[i] = 2 / ((1 - z^2) * dp^2)
    end
    x, w
end
const GLX, GLW = gausslegendre(10)
function glint(g, lo, hi, panels)            # composite Gauss-Legendre
    s = zero(g(lo))
    hp = (hi - lo) / panels
    for p in 0:panels-1
        c = lo + (p + 0.5) * hp
        for i in eachindex(GLX)
            s += GLW[i] * hp / 2 * g(c + hp / 2 * GLX[i])
        end
    end
    s
end

# ---------------------------------------------------------------- constants of the section
Aconst(al, sg) = exp(al^2 * sg^2 / 2)
Bconst(c, sg) = exp(c^2 * sg^2 / 2) / (sg * sqrt(2pi))
theta_exact(s) = sum(exp(-2pi^2 * n^2 / s^2) for n in -50:50)
theta_bound(s) = 1 + 2 / (exp(2pi^2 / s^2) - 1)
const THETA0 = 1 + 2 / (exp(2) - 1)

kappa(a, sg, h, nu) = exp(-(a + im * nu) * log(h) + sg^2 * (a + im * nu)^2 / 2)

# discrete Gaussian functional:  (Delta/2pi) sum_{|k|<=K} Y_k kappa_h(k Delta)
function Edisc(Yfun, a, sg, h, Dl, K)
    s = zero(Yfun(0))
    for k in -K:K
        s += Yfun(k) * kappa(a, sg, h, k * Dl)
    end
    s * Dl / (2pi)
end
# continuum Gaussian functional for exact data:  E f(h e^{sigma Z}) = int phi_sigma(w) f(h e^w) dw
function Econt(f, sg, h; shift = 0.0, jump = nothing)
    phi(w) = exp(-w^2 / (2sg^2)) / (sg * sqrt(2pi))
    g(w) = phi(w) * f(h * exp(w + shift))
    W = 11sg
    if jump === nothing
        return simpson(g, -W, W, 40000)
    else
        wj = log(jump / h) - shift
        if wj <= -W || wj >= W
            return simpson(g, -W, W, 40000)
        end
        return simpson(g, -W, wj - 1e-13, 40000) + simpson(g, wj + 1e-13, W, 40000)
    end
end
noise_sum(a, sg, h, Dl, K) = Dl / (2pi) * sum(abs(kappa(a, sg, h, k * Dl)) for k in -K:K)

# four-term defect bound for one Gaussian functional (Lemma lem:onefunctional)
function defect_terms(; a, sg, dl, h, Dl, K, eta, Qm, Qp, Qa)
    noise = theta_exact(sg * Dl) * Bconst(a, sg) * eta * h^(-a)
    trunc = Bconst(a, sg) * Qa * h^(-a) * erfc_(sg * K * Dl / sqrt(2))
    alias = (Bconst(a - dl, sg) * Qm * h^(-(a - dl)) + Bconst(a + dl, sg) * Qp * h^(-(a + dl))) /
            (exp(2pi * dl / Dl) - 1)
    (noise, trunc, alias)
end

println("== 0. special-function sanity ==")
check("Gamma(1/2)=sqrt(pi)", abs(gammac(0.5) - sqrt(pi)) < 1e-13)
check("|Gamma(1+i)|^2 = pi/sinh(pi)", abs(abs2(gammac(1 + 1im)) - pi / sinh(pi)) < 1e-13)
check("|Gamma(1/2+5i)|^2 = pi/cosh(5pi)", abs(abs2(gammac(0.5 + 5im)) / (pi / cosh(5pi)) - 1) < 1e-11)
check("erfc(1)", abs(erfc_(1.0) - 0.15729920705028513) < 1e-14)
check("erfc(y) <= exp(-y^2) on grid", all(erfc_(y) <= exp(-y^2) + 1e-16 for y in 0:0.05:6))

println("\n== 1. Gaussian Poisson identities (Lemma lem:poisson) ==")
for (a, sg, h, Dl, u) in [(1.0, 1.0, 0.3, 1.7, 0.4), (0.7, 0.8, 2.0, 2.5, -1.1), (1.5, 1.3, 0.05, 0.9, 2.3)]
    lhs = Dl / (2pi) * sum(exp(sg^2 * (a + im * k * Dl)^2 / 2 + im * k * Dl * u) for k in -400:400)
    p(w) = exp(-a * w - w^2 / (2sg^2)) / (sg * sqrt(2pi))
    T = 2pi / Dl
    rhs = sum(p(u - n * T) for n in -40:40)
    check(@sprintf("Poisson kernel identity a=%.1f s=%.1f D=%.1f", a, sg, Dl), abs(lhs - rhs) < 1e-12 * max(1, abs(rhs)),
          info = @sprintf("|lhs-rhs|=%.2e", abs(lhs - rhs)))
    lhs2 = Dl / (2pi) * sum(exp(-sg^2 * (k * Dl)^2 / 2) for k in -400:400)
    rhs2 = theta_exact(sg * Dl) / (sg * sqrt(2pi))
    check("theta identity", abs(lhs2 - rhs2) < 1e-13, info = @sprintf("%.15f vs %.15f", lhs2, rhs2))
    check("theta <= 1+2/(e^g-1)", theta_exact(sg * Dl) <= theta_bound(sg * Dl) + 1e-15)
end
check("theta(s) <= THETA0 for s <= pi", all(theta_exact(s) <= THETA0 + 1e-15 for s in 0.01:0.01:pi))

# ---------------------------------------------------------------- scalar profiles
# each: f, G, residue, alpha, M, Q(s) = Q_s(f) (exact or upper bound), jump location
profiles = [
    (name = "c*exp(-x)", f = x -> 0.7exp(-x), G = z -> 0.7gammac(z), r = 0.7, al = 1.0, M = 0.7,
     Q = s -> 0.7gammar(s), jump = nothing),
    (name = "(1+x)^-6", f = x -> (1 + x)^(-6), G = z -> gammac(z) * gammac(6 - z) / gammac(6.0), r = 1.0,
     al = 1.0, M = 6.0, Q = s -> gammar(s) * gammar(6 - s) / gammar(6.0), jump = nothing),
    (name = "exp(-sqrt x)", f = x -> exp(-sqrt(x)), G = z -> 2gammac(2z), r = 1.0, al = 0.5, M = 1.0,
     Q = s -> 2gammar(2s), jump = nothing),
]

println("\n== 2. exact aliasing formula and the three defect terms (Lemmas lem:alias, lem:onefunctional) ==")
for P in profiles, (a, sg, dl) in [(1.0, 1.0, 0.5), (0.8, 0.7, 0.4)], h in (0.05, 0.4), Dl in (0.6, 1.5, 2.5)
    T = 2pi / Dl
    Kbig = ceil(Int, 12 / (sg * Dl)) + 5
    Yex = k -> P.G(a + im * k * Dl)
    Dfull = Edisc(Yex, a, sg, h, Dl, Kbig)
    Ec = Econt(P.f, sg, h)
    # aliasing series  sum_n e^{a n T} E f(h e^{nT + sigma Z})
    ser = sum(exp(a * n * T) * Econt(P.f, sg, h; shift = n * T) for n in -30:30)
    tag = @sprintf("%s a=%.1f s=%.1f d=%.1f h=%.2f D=%.1f", P.name, a, sg, dl, h, Dl)
    check("alias formula " * tag, abs(Dfull - ser) < 1e-9 * max(1, abs(ser)), info = @sprintf("|diff|=%.1e", abs(Dfull - ser)))
    Qm, Qp = P.Q(a - dl), P.Q(a + dl); Qa = P.Q(a)
    check("Q_a <= sqrt(Q_- Q_+) " * tag, Qa <= sqrt(Qm * Qp) * (1 + 1e-13))
    for K in (1, 3, 8)
        nz, tr, al = defect_terms(; a, sg, dl, h, Dl, K, eta = 0.0, Qm, Qp, Qa)
        DK = Edisc(Yex, a, sg, h, Dl, K)
        okA = abs(Dfull - Ec) <= al + TOL
        okT = abs(DK - Dfull) <= tr + TOL
        K == 1 && check("aliasing bound " * tag, okA, info = @sprintf("|D-E|=%.3e <= %.3e", abs(Dfull - Ec), al))
        check("truncation bound K=$K " * tag, okT, info = @sprintf("%.3e <= %.3e", abs(DK - Dfull), tr))
    end
    for K in (2, 6)
        ns = noise_sum(a, sg, h, Dl, K)
        check("noise kernel sum K=$K " * tag, ns <= theta_exact(sg * Dl) * Bconst(a, sg) * h^(-a) * (1 + 1e-13))
    end
end

println("\n== 3. scalar residue estimator, adversarial noise, full four-term bound (Theorem thm:discrete (i)) ==")
for P in profiles, (a, sg, dl) in [(1.0, 1.0, 0.5), (1.0, 0.7071067811865476, 0.5)], h in (0.03, 0.2),
    (Dl, K) in [(0.4, 10), (1.0, 3), (2.0, 2), (0.25, 30)], eta in (0.0, 1e-4, 1e-2)
    Yex = k -> P.G(a + im * k * Dl)
    det = Edisc(Yex, a, sg, h, Dl, K) - P.r
    ph = abs(det) > 0 ? det / abs(det) : 1.0 + 0im
    # adversarial noise: e_k = eta * ph * conj(kappa_k)/|kappa_k|  (aligns with the deterministic error)
    Y = k -> Yex(k) + eta * ph * conj(kappa(a, sg, h, k * Dl)) / abs(kappa(a, sg, h, k * Dl))
    err = abs(Edisc(Y, a, sg, h, Dl, K) - P.r)
    Qm, Qp, Qa = P.Q(a - dl), P.Q(a + dl), P.Q(a)
    Q = max(Qm, Qp)
    bias = Aconst(P.al, sg) * P.M * h^P.al
    nz, tr, al = defect_terms(; a, sg, dl, h, Dl, K, eta, Qm = Q, Qp = Q, Qa = Q)
    bound = bias + nz + tr + al
    check(@sprintf("residue %s s=%.2f h=%.2f D=%.2f K=%d eta=%.0e", P.name, sg, h, Dl, K, eta), err <= bound + TOL,
          info = @sprintf("err=%.3e bound=%.3e", err, bound))
end

println("\n== 4. vector-valued (C^2) log data, m=2: top-log and full-jet estimators (Theorem thm:discrete (ii),(iii)) ==")
u = [ComplexF64[1.0, 0.5im], ComplexF64[-0.3, 0.2], ComplexF64[0.25 + 0.1im, -0.4]]
wv = ComplexF64[0.6, -0.2im]
mlog = 2
function vecprof(al)
    Ck(k) = k == 0 ? 1.0 : (k / (al * exp(1)))^k
    M = norm(wv) + sum(norm(u[k+1]) * Ck(k) for k in 0:mlog)
    f = x -> (x <= 1 ? sum(u[k+1] * log(1 / x)^k for k in 0:mlog) : zeros(ComplexF64, 2)) + wv * x^al * exp(-x)
    G = z -> sum(factorial(k) * u[k+1] / z^(k + 1) for k in 0:mlog) + wv * gammac(z + al)
    Q = s -> sum(norm(u[k+1]) * factorial(k) / s^(k + 1) for k in 0:mlog) + norm(wv) * gammar(s + al)
    (f = f, G = G, M = M, Q = Q, al = al)
end
function lagrange_coeffs(j, m, L, d)        # coefficients (ascending) of ell_j(t)
    c = [1.0]
    for r in 0:m
        r == j && continue
        root = L + r * d; den = d * (j - r)
        nc = zeros(length(c) + 1)
        for i in eachindex(c)
            nc[i+1] += c[i] / den; nc[i] -= root * c[i] / den
        end
        c = nc
    end
    c
end
function fulljet_weights(m, L, d, sg)
    W = zeros(m + 1, m + 1)                    # W[k+1, j+1]
    for j in 0:m
        c = lagrange_coeffs(j, m, L, d)
        for k in 0:m, r in 0:div(m - k, 2)
            W[k+1, j+1] += (-sg^2 / 2)^r / factorial(r) * factorial(k + 2r) * c[k+2r+1]
        end
    end
    W
end
Hpoly(n, A, sg) = sum((sg^2 / 2)^r * A^(n - 2r) / (factorial(r) * factorial(n - 2r)) for r in 0:div(n, 2))

for al in (1.0, 0.5), (a, sg, dl, d) in [(1.0, 1.0, 0.5, 1.0), (1.2, 0.8, 0.6, 0.7)], h in (0.05, 0.3),
    (Dl, K) in [(0.5, 8), (1.2, 4), (0.3, 25)], eta in (0.0, 1e-3)
    V = vecprof(al)
    Qm, Qp = V.Q(a - dl), V.Q(a + dl); Q = max(Qm, Qp)
    hs = [h * exp(-j * d) for j in 0:mlog]
    tag = @sprintf("al=%.1f a=%.1f s=%.1f h=%.2f D=%.1f K=%d eta=%.0e", al, a, sg, h, Dl, K, eta)
    # --- top-log estimator: weights d^-m (-1)^{m-j} binom(m,j)
    wt = [d^(-mlog) * (-1)^(mlog - j) * binomial(mlog, j) for j in 0:mlog]
    # combined coefficient c_k (vector-free) for adversarial noise
    for (target, w) in [("top", wt)]
        Yex = k -> V.G(a + im * k * Dl)
        det = sum(w[j+1] * Edisc(Yex, a, sg, hs[j+1], Dl, K) for j in 0:mlog) - factorial(mlog) * u[mlog+1]
        uvec = norm(det) > 0 ? det / norm(det) : ComplexF64[1, 0]
        ck = [Dl / (2pi) * sum(w[j+1] * kappa(a, sg, hs[j+1], k * Dl) for j in 0:mlog) for k in -K:K]
        err = norm(det) + eta * sum(abs.(ck))     # = norm of error under e_k = eta*uvec*conj(c_k)/|c_k|
        # direct evaluation with that noise, to confirm
        Yn = k -> Yex(k) + eta * uvec * conj(ck[k+K+1]) / max(abs(ck[k+K+1]), 1e-300)
        err2 = norm(sum(w[j+1] * Edisc(Yn, a, sg, hs[j+1], Dl, K) for j in 0:mlog) - factorial(mlog) * u[mlog+1])
        U = Aconst(V.al, sg) * (1 + exp(-V.al * d))^mlog / d^mlog
        Vc = Bconst(a, sg) * (1 + exp(a * d))^mlog / d^mlog
        bound = U * V.M * h^V.al + theta_exact(sg * Dl) * Vc * eta * h^(-a) +
                Vc * Q * h^(-a) * erfc_(sg * K * Dl / sqrt(2)) +
                Q / ((exp(2pi * dl / Dl) - 1) * d^mlog) *
                (Bconst(a - dl, sg) * h^(-(a - dl)) * (1 + exp((a - dl) * d))^mlog +
                 Bconst(a + dl, sg) * h^(-(a + dl)) * (1 + exp((a + dl) * d))^mlog)
        check("top-log " * tag, err <= bound + TOL && abs(err - err2) < 1e-9 * max(1, err),
              info = @sprintf("err=%.3e bound=%.3e", err, bound))
    end
    # --- full jet
    L = -log(h)
    W = fulljet_weights(mlog, L, d, sg)
    Aj = abs(L) + mlog * d
    for k in 0:mlog
        wbound = all(abs(W[k+1, j+1]) <= binomial(mlog, j) * d^(-mlog) * Hpoly(mlog - k, Aj, sg) * (1 + 1e-12) for j in 0:mlog)
        w = W[k+1, :]
        Yex = kk -> V.G(a + im * kk * Dl)
        det = sum(w[j+1] * Edisc(Yex, a, sg, hs[j+1], Dl, K) for j in 0:mlog) - factorial(k) * u[k+1]
        ck = [Dl / (2pi) * sum(w[j+1] * kappa(a, sg, hs[j+1], kk * Dl) for j in 0:mlog) for kk in -K:K]
        err = norm(det) + eta * sum(abs.(ck))
        Hn = Hpoly(mlog - k, Aj, sg) / d^mlog
        bound = Hn * (Aconst(V.al, sg) * V.M * h^V.al * (1 + exp(-V.al * d))^mlog +
                      (theta_exact(sg * Dl) * eta + Q * erfc_(sg * K * Dl / sqrt(2))) * Bconst(a, sg) * h^(-a) * (1 + exp(a * d))^mlog +
                      Q / (exp(2pi * dl / Dl) - 1) *
                      (Bconst(a - dl, sg) * h^(-(a - dl)) * (1 + exp((a - dl) * d))^mlog +
                       Bconst(a + dl, sg) * h^(-(a + dl)) * (1 + exp((a + dl) * d))^mlog))
        check("full-jet k=$k " * tag, wbound && err <= bound + TOL, info = @sprintf("err=%.3e bound=%.3e", err, bound))
    end
end

println("\n== 5. normalized branch estimator (Theorem thm:discrete (iv)) ==")
function Htmpl(z, lam, gm, a)   # int_0^infty e^{-zt}(1+lam t)^{-gamma} dt, Re z = a > 0
    Tm = 45 / a
    glint(t -> exp(-z * t) * (1 + lam * t)^(-gm), 0.0, Tm, 3000)
end
function Phi(h, sg, lam, gm)
    L = -log(h)
    lo = max(0.0, L - 12sg); hi = L + 12sg
    glint(t -> exp(-(t - L)^2 / (2sg^2)) * (1 + lam * t)^(-gm), lo, hi, 400) / (sg * sqrt(2pi))
end
lam, gm = 1.0, 0.5
for (bb, Mr, al) in [(0.8 + 0.3im, 0.5, 1.0), (1.5, 2.0, 0.7)], (a, sg, dl) in [(1.0, 1.0, 0.5)], h in (0.02, 0.2),
    (Dl, K) in [(0.5, 8), (1.1, 4), (0.35, 20)], eta in (0.0, 1e-3)
    Hcache = Dict{Int,ComplexF64}()
    G = k -> get!(Hcache, k) do
        z = a + im * k * Dl
        bb * Htmpl(z, lam, gm, a) + Mr * gammac(z + al)
    end
    Qs(s) = abs(bb) * real(Htmpl(s + 0im, lam, gm, s)) + Mr * gammar(s + al)
    Q = max(Qs(a - dl), Qs(a + dl))
    Ph = Phi(h, sg, lam, gm)
    det = Edisc(G, a, sg, h, Dl, K) / Ph - bb
    ns = noise_sum(a, sg, h, Dl, K)
    err = abs(det) + eta * ns / Ph
    nz, tr, alv = defect_terms(; a, sg, dl, h, Dl, K, eta, Qm = Q, Qp = Q, Qa = Q)
    bound = (Aconst(al, sg) * Mr * h^al + nz + tr + alv) / Ph
    check(@sprintf("branch b=%s M=%.1f h=%.2f D=%.2f K=%d eta=%.0e", string(bb), Mr, h, Dl, K, eta), err <= bound + TOL,
          info = @sprintf("err=%.3e bound=%.3e", err, bound))
    # template response identity  Phi_h = E phi(h e^{sigma Z})
    phif = x -> x <= 1 ? (1 + lam * log(1 / x))^(-gm) : 0.0
    check(@sprintf("Phi_h identity h=%.2f", h), abs(Ph - Econt(phif, sg, h; jump = 1.0)) < 1e-8)
end

# ---------------------------------------------------------------- calibrated node sets and rates
function calibrate(; a, sg, dl, Q, eta, eps, hlo, hhi)
    Om = Q > eps * eta ? sqrt(2) / sg * sqrt(log(Q / (eps * eta))) : 0.0
    X = Q * (Bconst(a - dl, sg) * hhi^dl + Bconst(a + dl, sg) * hlo^(-dl)) / (eps * Bconst(a, sg) * eta)
    Dl = min(pi / sg, 2pi * dl / log(1 + X))
    K = ceil(Int, Om / Dl)
    (Dl = Dl, K = K, Om = Om, X = X)
end

println("\n== 6. calibrated nodes: scalar residue, rate and sample count (Corollary cor:calibrated) ==")
let a = 1.0, al = 1.0, sg = 1.0, dl = 0.5, M = 1.0, Q = 2.0, eps = 0.5
    th = al / (a + al)
    U = Aconst(al, sg); V = Bconst(a, sg)
    C = U^(a / (a + al)) * V^th * ((a / al)^th + (al / a)^(a / (a + al)))
    for eta in 10.0 .^ (-2:-1:-10)
        et = (THETA0 + 2eps) * eta
        h = (a * V * et / (al * U * M))^(1 / (a + al))
        cal = calibrate(; a, sg, dl, Q, eta, eps, hlo = h, hhi = h)
        Dl, K = cal.Dl, cal.K
        rate = C * M^(a / (a + al)) * et^th
        # the three defect terms at this node set are each at most eps*B*eta*h^-a (except the noise term theta)
        nz, tr, alv = defect_terms(; a, sg, dl, h, Dl, K, eta, Qm = Q, Qp = Q, Qa = Q)
        Nbound = 3 + 2 * cal.Om * (sg / pi + log(1 + cal.X) / (2pi * dl))
        check(@sprintf("sample-count bound eta=%.0e", eta), 2K + 1 <= Nbound,
              info = @sprintf("N=%d <= %.1f ; asymptotic constant sqrt2(a+al+dl)/(pi s dl (a+al)) = %.3f", 2K + 1, Nbound,
                              sqrt(2) * (a + al + dl) / (pi * sg * dl * (a + al))))
        check(@sprintf("calibration eta=%.0e", eta),
              tr <= eps * V * eta * h^(-a) * (1 + 1e-12) && alv <= eps * V * eta * h^(-a) * (1 + 1e-12) &&
              theta_exact(sg * Dl) <= THETA0, info = @sprintf("Delta=%.4f K=%d N=%d N/log(1/eta)^1.5=%.3f", Dl, K, 2K + 1, (2K + 1) / log(1 / eta)^1.5))
        # profiles: exp(-x) (Q_s = Gamma(s) <= 2) and the lower-bound profile of thm:direct
        c = (2al * eta * M^(a / al) / gammar(a / al))^th
        prof = [("exp(-x)", z -> gammac(z), 1.0),
                ("lower f_c", z -> (c / al) * (c / M)^(z / al) * gammac(z / al), c)]
        Qc = s -> (c / al) * (c / M)^(s / al) * gammar(s / al)
        check(@sprintf("lower-bound profile satisfies Q prior eta=%.0e", eta), max(Qc(a - dl), Qc(a + dl)) <= Q,
              info = @sprintf("Q_- = %.2e, Q_+ = %.2e", Qc(a - dl), Qc(a + dl)))
        for (nm, G, r) in prof
            Yex = k -> G(a + im * k * Dl)
            det = Edisc(Yex, a, sg, h, Dl, K) - r
            err = abs(det) + eta * noise_sum(a, sg, h, Dl, K)
            check(@sprintf("rate %s eta=%.0e", nm, eta), err <= rate + TOL, info = @sprintf("err=%.3e rate bound=%.3e  lower=%.3e", err, rate, c / 2))
        end
    end
end

println("\n== 7. calibrated nodes: full jet m=2 (vector) and branch, rates with log factors ==")
let a = 1.0, al = 1.0, sg = 1.0, dl = 0.5, d = 1.0, eps = 0.5
    V = vecprof(al)
    M = V.M; Q = max(V.Q(a - dl), V.Q(a + dl))
    for eta in 10.0 .^ (-3:-1:-9)
        h = (eta / M)^(1 / (a + al)); L = -log(h)
        hs = [h * exp(-j * d) for j in 0:mlog]
        cal = calibrate(; a, sg, dl, Q, eta, eps, hlo = hs[end], hhi = hs[1])
        Dl, K = cal.Dl, cal.K
        S = M^(a / (a + al)) * eta^(al / (a + al))
        W = fulljet_weights(mlog, L, d, sg)
        for k in 0:mlog
            w = W[k+1, :]
            Yex = kk -> V.G(a + im * kk * Dl)
            det = sum(w[j+1] * Edisc(Yex, a, sg, hs[j+1], Dl, K) for j in 0:mlog) - factorial(k) * u[k+1]
            ck = [Dl / (2pi) * sum(w[j+1] * kappa(a, sg, hs[j+1], kk * Dl) for j in 0:mlog) for kk in -K:K]
            err = norm(det) + eta * sum(abs.(ck))
            bound = Hpoly(mlog - k, abs(L) + mlog * d, sg) / d^mlog * S *
                    (Aconst(al, sg) * (1 + exp(-al * d))^mlog + (THETA0 + 2eps) * Bconst(a, sg) * (1 + exp(a * d))^mlog)
            check(@sprintf("full-jet rate k=%d eta=%.0e", k, eta), err <= bound + TOL,
                  info = @sprintf("N=%d err=%.3e bound=%.3e bound/(S log^(m-k))=%.2f", 2K + 1, err, bound, bound / (S * log(M / eta)^(mlog - k))))
        end
        # lower-bound log profile (thm:loglower, thm:fulljet) satisfies the extra prior
        Cm = (mlog / (al * exp(1)))^mlog
        h0 = (2eta * Cm * a^(mlog + 1) / (factorial(mlog) * M))^(1 / (a + al))
        c0 = M * h0^al / Cm
        Qlog = s -> c0 * factorial(mlog) * h0^s / s^(mlog + 1)
        check(@sprintf("log lower profile Q prior eta=%.0e", eta), max(Qlog(a - dl), Qlog(a + dl)) <= Q,
              info = @sprintf("Q_-=%.2e", Qlog(a - dl)))
    end
end
let a = 1.0, al = 1.0, sg = 1.0, dl = 0.5, eps = 0.5, bb = 0.8 + 0.3im, Mr = 0.5
    Qs(s) = abs(bb) * real(Htmpl(s + 0im, lam, gm, s)) + Mr * gammar(s + al)
    Q = max(Qs(a - dl), Qs(a + dl)); M = Mr
    for eta in 10.0 .^ (-3:-1:-9)
        h = (eta / M)^(1 / (a + al)); L = -log(h)
        L >= sg || continue
        cal = calibrate(; a, sg, dl, Q, eta, eps, hlo = h, hhi = h)
        Dl, K = cal.Dl, cal.K
        S = M^(a / (a + al)) * eta^(al / (a + al))
        G = k -> (z = a + im * k * Dl; bb * Htmpl(z, lam, gm, a) + Mr * gammac(z + al))
        Ph = Phi(h, sg, lam, gm)
        p0 = 0.6826894921370859
        err = abs(Edisc(G, a, sg, h, Dl, K) / Ph - bb) + eta * noise_sum(a, sg, h, Dl, K) / Ph
        bound = (Aconst(al, sg) + (THETA0 + 2eps) * Bconst(a, sg)) / p0 * S * (1 + lam * (L + sg))^gm
        check(@sprintf("branch rate eta=%.0e", eta), Ph >= p0 * (1 + lam * (L + sg))^(-gm) && err <= bound + TOL,
              info = @sprintf("N=%d err=%.3e bound=%.3e", 2K + 1, err, bound))
        # branch lower-bound profile: Q_s(f_L) <= (M/s)(2 a eta/M)^{(al+s)/(a+al)}
        Le = log(M / (2a * eta)) / (a + al)
        if Le >= 0 && exp(-al * Le) * (1 + lam * Le)^gm <= 1
            bL = M * exp(-al * Le) * (1 + lam * Le)^gm
            Qexact = s -> bL * glint(t -> exp(-s * t) * (1 + lam * t)^(-gm), Le, Le + 60 / s, 600)
            Qbd = s -> (M / s) * (2a * eta / M)^((al + s) / (a + al))
            check(@sprintf("branch lower profile Q prior eta=%.0e", eta),
                  Qexact(a - dl) <= Qbd(a - dl) * (1 + 1e-10) && Qexact(a + dl) <= Qbd(a + dl) * (1 + 1e-10) &&
                  max(Qbd(a - dl), Qbd(a + dl)) <= Q, info = @sprintf("Q_-<=%.2e Q_+<=%.2e Q=%.2f", Qbd(a - dl), Qbd(a + dl), Q))
        end
    end
end

println("\n== 7b. calibrated nodes: top log coefficient m=2 at h_* (Theorem thm:discreterates (a)) ==")
let a = 1.0, al = 1.0, sg = 1.0, dl = 0.5, d = 1.0, eps = 0.5
    V = vecprof(al)
    M = V.M; Q = max(V.Q(a - dl), V.Q(a + dl))
    th = al / (a + al)
    U = Aconst(al, sg) * (1 + exp(-al * d))^mlog / d^mlog
    Vc = Bconst(a, sg) * (1 + exp(a * d))^mlog / d^mlog
    Cm = U^(a / (a + al)) * Vc^th * ((a / al)^th + (al / a)^(a / (a + al)))
    for eta in 10.0 .^ (-3:-1:-9)
        et = (THETA0 + 2eps) * eta
        h = (a * Vc * et / (al * U * M))^(1 / (a + al))
        hs = [h * exp(-j * d) for j in 0:mlog]
        cal = calibrate(; a, sg, dl, Q, eta, eps, hlo = hs[end], hhi = hs[1])
        Dl, K = cal.Dl, cal.K
        w = [d^(-mlog) * (-1)^(mlog - j) * binomial(mlog, j) for j in 0:mlog]
        Yex = kk -> V.G(a + im * kk * Dl)
        det = sum(w[j+1] * Edisc(Yex, a, sg, hs[j+1], Dl, K) for j in 0:mlog) - factorial(mlog) * u[mlog+1]
        ck = [Dl / (2pi) * sum(w[j+1] * kappa(a, sg, hs[j+1], kk * Dl) for j in 0:mlog) for kk in -K:K]
        err = norm(det) + eta * sum(abs.(ck))
        bound = Cm * (THETA0 + 2eps)^th * M^(a / (a + al)) * eta^th
        check(@sprintf("top-log rate eta=%.0e", eta), err <= bound + TOL,
              info = @sprintf("N=%d err=%.3e bound=%.3e", 2K + 1, err, bound))
    end
end

println("\n== 8. negative control: the bound must fail if the aliasing term is dropped at coarse spacing ==")
let a = 1.0, sg = 1.0, h = 0.2, Dl = 3.0, K = 3
    G = z -> gammac(z)
    det = abs(Edisc(k -> G(a + im * k * Dl), a, sg, h, Dl, K) - 1.0)
    bias = Aconst(1.0, sg) * 1.0 * h
    tr = Bconst(a, sg) * sqrt(pi) * h^(-a) * erfc_(sg * K * Dl / sqrt(2))
    check("coarse grid error exceeds bias+truncation (aliasing is real)", det > bias + tr,
          info = @sprintf("err=%.3e > %.3e", det, bias + tr))
end

println()
@printf("SUMMARY: %d PASS, %d FAIL\n", NPASS[], NFAIL[])
NFAIL[] == 0 || exit(1)
