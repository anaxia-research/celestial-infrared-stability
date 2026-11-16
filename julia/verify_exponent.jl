# verify_exponent.jl -- numerical checks for Section "Recovery of an unknown branch exponent".
# Standard library only. Own Gauss--Legendre quadrature (composite panels).
# Run: julia verify_exponent.jl

using Printf
using LinearAlgebra

const NFAIL = Ref(0)
const NPASS = Ref(0)
function check(name, ok)
    if ok
        NPASS[] += 1
        println("PASS  ", name)
    else
        NFAIL[] += 1
        println("FAIL  ", name)
    end
end

# ---------------- Gauss--Legendre nodes (Newton on Legendre recursion) -------------
function gauleg(n)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cos(pi * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p1 = 1.0; p2 = 0.0
            for j in 1:n
                p3 = p2; p2 = p1
                p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            z1 = z; z = z1 - p1 / pp
            abs(z - z1) < 1e-15 && break
        end
        p1 = 1.0; p2 = 0.0
        for j in 1:n
            p3 = p2; p2 = p1
            p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
        end
        pp = n * (z * p1 - p2) / (z^2 - 1)
        x[i] = z; w[i] = 2 / ((1 - z^2) * pp^2)
    end
    x, w
end
const GX, GW = gauleg(30)

# composite GL on [lo,hi] with np panels; f may return real or complex
function quad(f, lo, hi; np = 40)
    hstep = (hi - lo) / np
    s = zero(f(lo + 0.5 * hstep))
    for k in 0:np-1
        c = lo + (k + 0.5) * hstep; r = 0.5 * hstep
        for i in eachindex(GX)
            s += GW[i] * r * f(c + r * GX[i])
        end
    end
    s
end

gdens(z) = exp(-z^2 / 2) / sqrt(2pi)
const p0 = quad(gdens, -1.0, 1.0; np = 4)                 # erf(1/sqrt2)
const m2 = quad(z -> z^2 * gdens(z), -1.0, 1.0; np = 4)    # E Z^2 1{|Z|<=1}
Qtail(x) = x > 40 ? 0.0 : quad(gdens, x, x + 40.0; np = 80)  # P(Z > x), x >= 0

println("p0 = ", p0, "   m2 = ", m2, "   p0-2phi(1) = ", p0 - 2 * gdens(1.0))
check("p0 = erf(1/sqrt2) value", abs(p0 - 0.6826894921370859) < 1e-13)
check("m2 = p0 - 2 phi(1)", abs(m2 - (p0 - 2 * gdens(1.0))) < 1e-13)
check("Gaussian tail bound Q(x) <= exp(-x^2/2)/2 on grid",
      all(Qtail(x) <= 0.5 * exp(-x^2 / 2) * (1 + 1e-12) for x in 0:0.25:8))

# ---------------- template quantities --------------------------------------------
struct Par
    a::Float64; al::Float64; lam::Float64; sig::Float64; M::Float64
    g0::Float64; g1::Float64; b0::Float64; b1::Float64; s::Float64
end
Kg(P, g, t) = t < 0 ? 0.0 : (1 + P.lam * t)^(-g)
u(P, t) = log1p(P.lam * t)
Acon(P) = exp(P.al^2 * P.sig^2 / 2)
Bcon(P) = exp(P.a^2 * P.sig^2 / 2) / (sqrt(2pi) * P.sig)

# integrate a function of t against the N(L, sig^2) density restricted to t >= 0
function gint(P, L, F)
    lo = max(0.0, L - 14P.sig); hi = L + 14P.sig
    quad(t -> exp(-(t - L)^2 / (2P.sig^2)) / (sqrt(2pi) * P.sig) * F(t), lo, hi; np = 60)
end
Phi(P, L, g) = gint(P, L, t -> Kg(P, g, t))
Emu(P, L, g, F) = gint(P, L, t -> Kg(P, g, t) * F(t)) / Phi(P, L, g)
covTu(P, L, g) = Emu(P, L, g, t -> t * u(P, t)) - Emu(P, L, g, t -> t) * Emu(P, L, g, t -> u(P, t))
kappa(P, L1) = 3.0^(-P.g1) / (1 + 0.5 * exp(-L1^2 / (8P.sig^2)) * (1 + P.lam * L1 / 2)^P.g1)

# ---------------- Lemma checks ---------------------------------------------------
PA = Par(1.0, 1.0, 1.0, 1.0, 1.0, 0.2, 0.8, 1.0, 3.0, 1.0)
PB = Par(0.5, 1.5, 2.0, 0.7, 2.0, 0.3, 0.9, 0.5, 2.0, 0.6)
PC = Par(2.0, 0.5, 0.3, 1.3, 0.5, 0.1, 0.6, 2.0, 2.5, 1.5)

for (nm, P) in (("A", PA), ("B", PB), ("C", PC))
    gs = range(0.0, P.g1, length = 7)
    Ls = [P.sig, 2P.sig, 3P.sig, 5.0, 10.0, 25.0, 60.0]
    ok1 = true; ok2 = true; ok3 = true; ok4 = true; ok5 = true; ok6 = true
    for L in Ls, g in gs
        ph = Phi(P, L, g)
        L >= P.sig && (ok1 &= ph >= p0 * (1 + P.lam * (L + P.sig))^(-g))
        # d/dgamma log Phi = -E_mu u  (central finite difference)
        dg = 1e-5
        fd = (log(Phi(P, L, g + dg)) - log(Phi(P, L, g - dg))) / (2dg)   # K_g is defined for g<0 too
        ok2 &= abs(fd + Emu(P, L, g, t -> u(P, t))) < 1e-7 * (1 + abs(fd))
        # E_mu u <= log(1+lam(L+sig))
        ok3 &= Emu(P, L, g, t -> u(P, t)) <= log1p(P.lam * (L + P.sig))
        # Phi upper bound K(L/2)+Q(L/(2 sig))
        ok4 &= ph <= Kg(P, g, L / 2) + Qtail(L / (2P.sig)) + 1e-14
        if L >= 2P.sig
            cL = Kg(P, g, L + P.sig) / ph
            ok5 &= cL >= kappa(P, L)
            # covariance lower bound
            ok6 &= covTu(P, L, g) / P.sig^2 >= P.lam * kappa(P, L)^2 * p0 * m2 / (1 + P.lam * (L + P.sig))
        end
    end
    check("[$nm] Phi_L >= p0 (1+lam(L+sig))^(-g)  (L>=sig)", ok1)
    check("[$nm] d/dg log Phi_L(g) = -E_mu[log(1+lam T)]", ok2)
    check("[$nm] E_mu[log(1+lam T)] <= log(1+lam(L+sig))", ok3)
    check("[$nm] Phi_L <= K(L/2) + Q(L/(2sig))", ok4)
    check("[$nm] c_L = K(L+sig)/Phi_L >= kappa(L)  (L>=2sig)", ok5)
    check("[$nm] Cov(T,u)/sig^2 >= lam kappa^2 p0 m2/(1+lam(L+sig))", ok6)
    # d/dL E_mu[u] = Cov/sig^2
    okL = true
    for L in (2P.sig, 4.0, 12.0), g in (0.1, P.g1)
        dL = 1e-4
        fd = (Emu(P, L + dL, g, t -> u(P, t)) - Emu(P, L - dL, g, t -> u(P, t))) / (2dL)
        okL &= abs(fd - covTu(P, L, g) / P.sig^2) < 1e-6 * (1 + abs(fd))
    end
    check("[$nm] d/dL E_mu[u] = Cov(T,u)/sig^2", okL)
    # monotonicity of tau term on L>=2sig
    taus = [0.5 * exp(-L^2 / (8P.sig^2)) * (1 + P.lam * L / 2)^P.g1 for L in range(2P.sig, 30.0, length = 400)]
    check("[$nm] tau(L) decreasing on L>=2sig", all(diff(taus) .<= 0))
    # rho'(g) >= D for several (L1,L2)
    okD = true
    for L1 in (2P.sig, 5.0, 15.0)
        L2 = L1 + 2P.s
        D = p0 * m2 * kappa(P, L1)^2 * log((1 + P.lam * (L2 + P.sig)) / (1 + P.lam * (L1 + P.sig)))
        okD &= D >= 2 * p0 * m2 * kappa(P, L1)^2 * P.lam * P.s / (1 + P.lam * (L2 + P.sig))
        for g in range(0.0, P.g1, length = 9)
            rp = Emu(P, L2, g, t -> u(P, t)) - Emu(P, L1, g, t -> u(P, t))
            okD &= rp >= D
        end
    end
    check("[$nm] rho'(g) >= D = p0 m2 kappa^2 log(...) >= 2 p0 m2 kappa^2 lam s/(1+lam(L2+sig))", okD)
end

# ---------------- Gaussian identity and kernel norm (full frequency pipeline) -------
# f = b phi_g + R,  R(x) = r1 M x^al 1{x<=1} + r2 M x^(-c) 1{x>1},  c > a
Hg(P, g, z) = quad(t -> exp(-z * t) * (1 + P.lam * t)^(-g), 0.0, 60.0 / P.a; np = 600)
GR(P, r1, r2, c, z) = r1 * P.M / (z + P.al) + r2 * P.M / (c - z)
function Ehs_line(P, Yfun, h)      # (1/2pi) int Y(nu) h^{-a-i nu} exp(sig^2 (a+i nu)^2/2) dnu
    W = 13 / P.sig
    quad(nu -> Yfun(nu) * h^(-P.a - im * nu) * exp(P.sig^2 * (P.a + im * nu)^2 / 2), -W, W; np = 200) / (2pi)
end
kern_abs(P, h, nu) = abs(h^(-P.a - im * nu) * exp(P.sig^2 * (P.a + im * nu)^2 / 2)) / (2pi)
function Ef_direct(P, b, g, r1, r2, c, h)
    Rf(x) = x <= 1 ? r1 * P.M * x^P.al : r2 * P.M * x^(-c)
    L = -log(h); z0 = L / P.sig        # R jumps at x = 1, i.e. z = z0: split there
    b * Phi(P, L, g) + quad(z -> gdens(z) * Rf(h * exp(P.sig * z)), -14.0, z0; np = 400) +
        quad(z -> gdens(z) * Rf(h * exp(P.sig * z)), z0, 14.0; np = 400)
end
let P = PA, b = 1.7, g = 0.55, r1 = -0.8, r2 = 0.6, c = 1.6
    okI = true; okK = true
    for L in (1.0, 2.5, 4.0)
        h = exp(-L)
        Gline(nu) = b * Hg(P, g, P.a + im * nu) + GR(P, r1, r2, c, P.a + im * nu)
        lhs = Ehs_line(P, Gline, h)
        rhs = Ef_direct(P, b, g, r1, r2, c, h)
        @printf("   Gaussian identity L=%.1f : line=%.12f%+.2ei  direct=%.12f\n", L, real(lhs), imag(lhs), rhs)
        okI &= abs(lhs - rhs) < 1e-8
        W = 13 / P.sig
        kn = quad(nu -> kern_abs(P, h, nu), -W, W; np = 200)
        okK &= abs(kn - Bcon(P) * h^(-P.a)) < 1e-10 * Bcon(P) * h^(-P.a)
    end
    check("Gaussian identity E_{h,sig}(G) = E f(h e^{sig Z}) (frequency quadrature)", okI)
    check("L1 norm of frequency kernel = B_{a,sig} h^{-a}", okK)
end

# ---------------- estimator ------------------------------------------------------
# bias of remainder R(x) = sgn*M*min(x^al, x^(-a-1)) (|R| <= M x^al, Q_a finite)
function biasR(P, h, sgn)
    Rf(x) = sgn * P.M * min(x^P.al, x^(-P.a - 1))
    z0 = -log(h) / P.sig              # kink of R at x = 1
    quad(z -> gdens(z) * Rf(h * exp(P.sig * z)), -14.0, z0; np = 400) +
        quad(z -> gdens(z) * Rf(h * exp(P.sig * z)), z0, 14.0; np = 400)
end
rho(P, L1, L2, g) = log(Phi(P, L1, g)) - log(Phi(P, L2, g))
function gamma_hat(P, L1, L2, E1, E2)
    lr = log(norm(E1) / norm(E2))
    lo, hi = P.g0, P.g1
    rlo = rho(P, L1, L2, lo); rhi = rho(P, L1, L2, hi)
    lr <= rlo && return lo
    lr >= rhi && return hi
    for _ in 1:70
        mid = (lo + hi) / 2
        rho(P, L1, L2, mid) < lr ? (lo = mid) : (hi = mid)
    end
    (lo + hi) / 2
end
Sfun(P, eta) = P.M^(P.a / (P.a + P.al)) * eta^(P.al / (P.a + P.al))
function constants(P, eta)
    Lc = log(P.M / eta) / (P.a + P.al); L1 = Lc - P.s; L2 = Lc + P.s
    A = Acon(P); B = Bcon(P); m = max(P.a, P.al); S = Sfun(P, eta)
    k = kappa(P, L1)
    Cg = 2 * (A + B) * exp(m * P.s) / (p0^2 * m2 * k^2 * P.lam * P.s)
    Lam = 1 + P.lam * (L2 + P.sig)
    Uc = log1p(P.lam * (Lc + P.sig))
    R1 = L1 >= 2P.sig
    R2 = 2 * (A + B) * exp(m * P.s) * S * Lam^P.g1 <= p0 * P.b0
    R3 = Cg * S * Lam^(1 + P.g1) * Uc / P.b0 <= 1
    (; Lc, L1, L2, A, B, S, Cg, Lam, Uc, R1, R2, R3)
end
function gamma_bound(P, eta, g, bn)
    c = constants(P, eta); c.Cg * c.S * c.Lam^(1 + g) / bn
end
function b_bound(P, eta, g)
    c = constants(P, eta)
    2 * c.Uc * c.Cg * c.S * c.Lam^(1 + g) + exp(1) * (c.A + c.B) * c.S * (1 + P.lam * (c.Lc + P.sig))^g / p0
end

# adversarial data at scale L: E = b Phi + bias + noise, worst-case sign/phase patterns
function run_estimator(P, eta, b::Vector{ComplexF64}, g, sgnR, nz)
    c = constants(P, eta)
    Ls = (c.L1, c.L2, c.Lc)
    Es = Vector{Vector{ComplexF64}}(undef, 3)
    for (j, L) in enumerate(Ls)
        h = exp(-L)
        noise = nz[j] * Bcon(P) * eta * h^(-P.a)       # nz[j]: vector of norm <= 1
        Es[j] = b .* Phi(P, L, g) .+ biasR(P, h, sgnR[j]) .* (b ./ norm(b)) .+ noise
    end
    gh = gamma_hat(P, c.L1, c.L2, Es[1], Es[2])
    bh = Es[3] ./ Phi(P, c.Lc, gh)
    gh, bh
end

println("\nEstimator under adversarial perturbations (scalar and C^2-valued amplitudes):")
# Float64 arithmetic: noise levels are chosen so that the proved bounds stay far above
# the ~1e-14 relative accuracy of the quadrature of Phi (otherwise roundoff dominates).
for (nm, P, etas) in (("A", PA, (1e-14, 1e-20, 1e-30)),
                      ("B", PB, (1e-12, 1e-16, 1e-22)),
                      ("C", PC, (1e-30, 1e-50, 1e-80)))
    for eta in etas
        c = constants(P, eta)
        regime = c.R1 && c.R2 && c.R3
        check("[$nm] regime (R1)-(R3) holds at eta=$eta", regime)
        regime || continue
        okg = true; okb = true; worstg = 0.0; worstb = 0.0
        for g in (P.g0, (P.g0 + P.g1) / 2, P.g1), bmag in (P.b0, P.b1)
            for sR in Iterators.product((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)),
                sN in Iterators.product((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0))
                # scalar (embedded as first coordinate) with real signs
                b = ComplexF64[bmag, 0]
                nz = [ComplexF64[sN[j], 0] for j in 1:3]
                gh, bh = run_estimator(P, eta, b, g, sR, nz)
                eg = abs(gh - g); eb = norm(bh - b)
                okg &= eg <= gamma_bound(P, eta, g, bmag)
                okb &= eb <= b_bound(P, eta, g)
                worstg = max(worstg, eg / gamma_bound(P, eta, g, bmag))
                worstb = max(worstb, eb / b_bound(P, eta, g))
            end
            # Hilbert-valued: b in C^2, noise orthogonal or complex-phased
            b = ComplexF64[bmag * cis(0.3), 0]
            for nzset in ([ComplexF64[0, 1], ComplexF64[0, -1], ComplexF64[0, 1]],
                          [ComplexF64[cis(1.1), 0], ComplexF64[-cis(0.3), 0], ComplexF64[cis(2.0), 0]])
                gh, bh = run_estimator(P, eta, b, g, (1.0, -1.0, 1.0), nzset)
                okg &= abs(gh - g) <= gamma_bound(P, eta, g, bmag)
                okb &= norm(bh - b) <= b_bound(P, eta, g)
            end
        end
        @printf("   [%s] eta=%.0e  worst |g^-g|/bound=%.3e  worst |b^-b|/bound=%.3e\n", nm, eta, worstg, worstb)
        check("[$nm] |gamma_hat-gamma| <= C_g S (1+lam(L2+sig))^{1+g}/|b|  at eta=$eta", okg)
        check("[$nm] |b_hat-b| <= 2U_c C_g S Lam^{1+g} + e(A+B)S(1+lam(Lc+sig))^g/p0  at eta=$eta", okb)
    end
end

# one end-to-end run with genuine frequency-line data and adversarial noise function
let P = Par(1.0, 1.0, 1.0, 0.6, 1.0, 0.2, 0.8, 1.0, 3.0, 1.0), eta = 1e-14
    c = constants(P, eta)
    check("[E2E] regime holds", c.R1 && c.R2 && c.R3)
    b = 2.2; g = 0.7; r1 = 1.0; r2 = -1.0; cc = 1.5
    Es = ComplexF64[]
    for (j, L) in enumerate((c.L1, c.L2, c.Lc))
        h = exp(-L)
        sgn = j == 1 ? 1.0 : -1.0
        Gline(nu) = b * Hg(P, g, P.a + im * nu) + GR(P, r1, r2, cc, P.a + im * nu)
        # adversarial noise: e(nu) = sgn*eta*conj(k)/|k|, k the frequency kernel
        knl(nu) = h^(-P.a - im * nu) * exp(P.sig^2 * (P.a + im * nu)^2 / 2)
        Yfun(nu) = Gline(nu) + sgn * eta * conj(knl(nu)) / abs(knl(nu))
        push!(Es, Ehs_line(P, Yfun, h))
    end
    gh = gamma_hat(P, c.L1, c.L2, [Es[1]], [Es[2]])
    bh = Es[3] / Phi(P, c.Lc, gh)
    @printf("   [E2E] gamma=%.4f gamma_hat=%.8f (bound %.3e) ; b=%.4f b_hat=%.8f%+.1ei (bound %.3e)\n",
            g, gh, gamma_bound(P, eta, g, b), b, real(bh), imag(bh), b_bound(P, eta, g))
    check("[E2E] gamma error within bound", abs(gh - g) <= gamma_bound(P, eta, g, b))
    check("[E2E] b error within bound", abs(bh - b) <= b_bound(P, eta, g))
end

# ---------------- lower bound two-point construction -----------------------------
function Calpha(P)
    # sup_{s>=0} lam s (1+lam s)^{1+g1} e^{-al s}: grid + golden refinement
    f(s) = P.lam * s * (1 + P.lam * s)^(1 + P.g1) * exp(-P.al * s)
    smax = 50 * (1 + P.g1 + 2) / P.al
    grid = range(0, smax, length = 200001)
    vals = f.(grid); i = argmax(vals)
    lo = grid[max(i - 1, 1)]; hi = grid[min(i + 1, length(grid))]
    for _ in 1:200
        m1 = lo + (hi - lo) * 0.381966; m2_ = lo + (hi - lo) * 0.618034
        f(m1) < f(m2_) ? (lo = m1) : (hi = m2_)
    end
    f((lo + hi) / 2) * (1 + 1e-12)
end
Calpha_closed(P) = max(1.0, (3P.lam / (exp(1) * P.al))^3 * exp(P.al / P.lam))
function lower_pair(P, eta, gs)
    # pair: f1 = beta*phi_{gs} (R1=0), beta = b0 (1+lam L)^Del ; f2 = b0*phi_{gs-Del} + R2
    Ca = Calpha(P); th = P.al / (P.a + P.al)
    cflat = (2P.a^2 / P.lam)^th * Ca^(-(1 - th))
    Lf = log(P.lam * P.M / (2P.a^2 * Ca * eta)) / (P.a + P.al)
    S = Sfun(P, eta)
    Del = cflat * S * (1 + P.lam * Lf)^(1 + gs) / (exp(1) * P.b0)
    beta = P.b0 * (1 + P.lam * Lf)^Del
    bdiff = P.b0 * expm1(Del * log1p(P.lam * Lf))     # beta - b0 without cancellation
    (; Ca, cflat, Lf, Del, beta, bdiff, S)
end
println("\nTwo-point lower-bound constructions:")
for (nm, P, etas) in (("A", PA, (1e-14, 1e-20, 1e-30, 1e-45)),
                      ("B", PB, (1e-16, 1e-24, 1e-36)),
                      ("C", PC, (1e-30, 1e-50, 1e-80)))
    Ca = Calpha(P)
    check("[$nm] C_alpha <= closed-form bound max(1,(3lam/(e al))^3 e^{al/lam})", Ca <= Calpha_closed(P))
    rows = []
    for eta in etas
        gs = P.g1
        lp = lower_pair(P, eta, gs)
        cond = lp.Lf >= 0 && lp.Del <= gs - P.g0 && lp.beta <= P.b1 && lp.Del * log1p(P.lam * lp.Lf) <= 1
        check("[$nm] lower-bound conditions at eta=$eta", cond)
        cond || continue
        L = lp.Lf; Del = lp.Del; gp = gs - Del; beta = lp.beta
        w(t) = log((1 + P.lam * t) / (1 + P.lam * L))
        D(t) = -beta * Kg(P, gs, t) * expm1(Del * w(t))   # = beta K_g - b0 K_{g-Del}
        okid = all(abs(D(t) - (beta * Kg(P, gs, t) - P.b0 * Kg(P, gp, t))) <= 1e-13 * beta * Kg(P, gs, t) for t in (0.0, L / 3, L, 2L))
        check("[$nm] D(t) = beta K_g(t) - b0 K_{g-Del}(t) at eta=$eta", okid)
        check("[$nm] beta Del (1+lam L)^{-1-g} <= c_flat S_eta", beta * Del * (1 + P.lam * L)^(-1 - gs) <= lp.cflat * lp.S * (1 + 1e-12))
        # remainder constraint on [0,L]
        ts = range(0, L, length = 200001)
        supR = maximum(abs(D(t)) * exp(P.al * t) for t in ts)
        check("[$nm] sup_{0<=t<=L} |R_2| e^{al t} <= M at eta=$eta  (ratio $(round(supR/P.M,digits=4)))", supR <= P.M)
        # transform constraint: int_L^inf e^{-a t}|D| dt <= 2 eta
        T = L + 80 / P.a
        I = quad(t -> exp(-P.a * t) * abs(D(t)), L, T; np = 2000)
        check("[$nm] int_L^inf e^{-at}|D| <= 2 eta at eta=$eta (ratio $(round(I/(2eta),digits=4)))", I <= 2eta)
        # rate comparison
        c = constants(P, eta)
        lowg = Del / 2
        lowb = lp.bdiff / 2
        lowb2 = lp.cflat / (2 * exp(1)) * lp.S * (1 + P.lam * L)^(1 + gs) * log1p(P.lam * L)
        check("[$nm] (beta-b0)/2 >= (c_flat/(2e)) S (1+lam L)^{1+g} log(1+lam L)", lowb >= lowb2 * (1 - 1e-12))
        upg = gamma_bound(P, eta, gs, P.b0); upb = b_bound(P, eta, gs)
        check("[$nm] lower <= upper (gamma, b) at eta=$eta", lowg <= upg && lowb <= upb)
        Lr = log(P.M / eta)
        push!(rows, (eta, lowg / (lp.S * Lr^(1 + gs) / P.b0), upg / (lp.S * Lr^(1 + gs) / P.b0),
                     lowb / (lp.S * Lr^(1 + gs) * log(Lr)), upb / (lp.S * Lr^(1 + gs) * log(Lr))))
    end
    println("   eta      low_g/rate   up_g/rate   low_b/rate   up_b/rate   (rate_g=S L^{1+g1}/b0, rate_b=S L^{1+g1} log L, L=log(M/eta))")
    for r in rows
        @printf("   %.0e   %.4e  %.4e  %.4e  %.4e\n", r...)
    end
end

# asymptotic ratio upper/lower tends to a finite constant: check stabilisation for very small eta
let P = PA
    rat(eta) = begin
        lp = lower_pair(P, eta, P.g1)
        (gamma_bound(P, eta, P.g1, P.b0) / (lp.Del / 2), b_bound(P, eta, P.g1) / (lp.bdiff / 2))
    end
    r1 = rat(1e-100); r2 = rat(1e-200); r3 = rat(1e-300)
    @printf("   upper/lower ratios: eta=1e-100: %.3f, %.3f ; 1e-200: %.3f, %.3f ; 1e-300: %.3f, %.3f\n", r1..., r2..., r3...)
    check("[A] upper/lower ratio for gamma bounded and converging", abs(r3[1] / r2[1] - 1) < abs(r2[1] / r1[1] - 1) + 1e-3 && r3[1] < 2 * r1[1])
    check("[A] upper/lower ratio for b bounded and converging", r3[2] < 2 * r1[2])
end

println()
@printf("Summary: %d PASS, %d FAIL\n", NPASS[], NFAIL[])
println(NFAIL[] == 0 ? "ALL PASS" : "SOME CHECKS FAILED")
