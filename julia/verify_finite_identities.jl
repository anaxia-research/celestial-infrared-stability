#!/usr/bin/env julia
# Exact replays (Rational{BigInt}, standard library only) of the statements
# formalized in CelestialIR.lean, over many instances, plus an exact check that
# the full-jet weights W_{k,j} of the paper (mdpi-preprints/CelestialIRPaper.pdf) reproduce beta_k.
# This is a finite computational check, not a proof.

using Random

const Q = Rational{BigInt}
q(a, b=1) = BigInt(a) // BigInt(b)
fact(n) = factorial(BigInt(n))
C(n, k) = (k < 0 || k > n) ? BigInt(0) : binomial(BigInt(n), BigInt(k))

rng = MersenneTwister(20261116)
randq(lo=-9, hi=9, den=7) = q(rand(rng, lo:hi), rand(rng, 1:den))

counts = Dict{String,Int}()
function check(name, cond)
    cond || error("FAIL: $name")
    counts[name] = get(counts, name, 0) + 1
end
pass(name) = println("PASS  $(rpad(name, 52)) ($(counts[name]) exact checks)")

# ---------------------------------------------------------------------------
# 1. Pole-jet shear of arbitrary length r (u-basis: beta_k = k! u_k)
# ---------------------------------------------------------------------------
shear(sig, u) = [sum(C(j, k) * (-sig)^(j - k) * u[j+1] for j in k:length(u)-1) for k in 0:length(u)-1]
evalpoly(u, t) = sum(u[k+1] * t^k for k in 0:length(u)-1)
# eq:jet verbatim in the beta basis (neutral prefactor w+s0=0, and a general prefactor c)
jet(sig, beta) = [sum((-sig)^(j - k) // fact(j - k) * beta[j+1] for j in k:length(beta)-1) for k in 0:length(beta)-1]

for trial in 1:400
    r = rand(rng, 1:9)
    u = [randq() for _ in 1:r]
    s, t, ell = randq(), randq(), randq()
    check("1 group law S_t S_s = S_(s+t)", shear(t, shear(s, u)) == shear(s + t, u))
    check("1 identity S_0 = id", shear(q(0), u) == u)
    check("1 inverse S_(-s) S_s = id", shear(-s, shear(s, u)) == u)
    check("1 compensator S_(l-s) S_s = S_l", shear(ell - s, shear(s, u)) == shear(ell, u))
    check("1 top coefficient invariant", shear(s, u)[end] == u[end])
    x = randq()
    check("1 translation semantics P(t) -> P(t - s)", evalpoly(shear(s, u), x) == evalpoly(u, x - s))
    beta = [fact(k) * u[k+1] for k in 0:r-1]
    c = randq(1, 9)
    check("1 eq:jet (beta basis, prefactor c)", c .* [fact(k) * shear(s, u)[k+1] for k in 0:r-1] == c .* jet(s, beta))
    # pole order preserved: lower actual order d
    d = rand(rng, 0:r-1)
    ud = [k <= d ? u[k+1] : q(0) for k in 0:r-1]
    sd = shear(s, ud)
    check("1 pole order and top coefficient preserved", all(iszero, sd[d+2:end]) && sd[d+1] == ud[d+1])
    # integer powers of the unit shear
    n = rand(rng, 0:6)
    v = u; for _ in 1:n; v = shear(q(1), v); end
    check("1 integer powers S_1^n = S_n", v == shear(q(n), u))
end
# length-2 reproduction of ScaleJets.shift and residue negative control
check("1 r=2: (a1,a2) -> (a1 - s a2, a2)", shear(q(5), [q(3), q(2)]) == [q(3) - 5 * 2, q(2)])
check("1 residue not invariant (control)", shear(q(5), [q(3), q(2), q(1)])[1] != 3)
for nm in ["1 group law S_t S_s = S_(s+t)", "1 identity S_0 = id", "1 inverse S_(-s) S_s = id",
           "1 compensator S_(l-s) S_s = S_l", "1 top coefficient invariant",
           "1 translation semantics P(t) -> P(t - s)", "1 eq:jet (beta basis, prefactor c)",
           "1 pole order and top coefficient preserved", "1 integer powers S_1^n = S_n",
           "1 r=2: (a1,a2) -> (a1 - s a2, a2)", "1 residue not invariant (control)"]
    pass(nm)
end

# ---------------------------------------------------------------------------
# 2. Forward-difference identity and estimator exactness
# ---------------------------------------------------------------------------
fd(m, f) = sum((-1)^(m - j) * C(m, j) * f(j) for j in 0:m)
for m in 0:30, k in 0:m
    val = fd(m, j -> BigInt(j)^k)
    check("2 sum (-1)^(m-j) C(m,j) j^k = 0 (k<m), m! (k=m)", val == (k == m ? fact(m) : 0))
end
for m in 1:12
    # control: k = m+1 gives m! m(m+1)/2, not 0
    check("2 control k=m+1 gives m!*m(m+1)/2", fd(m, j -> BigInt(j)^(m + 1)) == fact(m) * m * (m + 1) ÷ 2)
end
for trial in 1:300
    m = rand(rng, 0:10)
    a = [randq() for _ in 0:m]
    L, d = randq(), randq(1, 9)
    val = fd(m, j -> sum(a[k+1] * (L + j * d)^k for k in 0:m))
    check("2 sum (-1)^(m-j)C(m,j)P(L+jd) = m! d^m a_m", val == fact(m) * d^m * a[m+1])
end
pass("2 sum (-1)^(m-j) C(m,j) j^k = 0 (k<m), m! (k=m)")
pass("2 control k=m+1 gives m!*m(m+1)/2")
pass("2 sum (-1)^(m-j)C(m,j)P(L+jd) = m! d^m a_m")

# ---------------------------------------------------------------------------
# 3. Bondi norm and H^2 domination
# ---------------------------------------------------------------------------
for l in -60:300
    L = BigInt(l); lam = L * (L + 1)
    check("3 l(l+1)(l(l+1)-2) = (l-1)l(l+1)(l+2)", lam * (lam - 2) == (L - 1) * L * (L + 1) * (L + 2))
    if l >= 0
        check("3 H_l >= 0, = 0 iff l <= 1", ((L - 1) * L * (L + 1) * (L + 2) >= 0) &&
              (((L - 1) * L * (L + 1) * (L + 2) == 0) == (l <= 1)))
    end
end
for trial in 1:1000
    lam = trial <= 500 ? q(trial - 1) : abs(randq(0, 50, 11))
    check("3 lam(lam-2) <= (1+lam)^2, gap 1+4lam", lam * (lam - 2) <= (1 + lam)^2 && (1 + lam)^2 - lam * (lam - 2) == 1 + 4lam)
end
for trial in 1:100
    n = rand(rng, 1:40)
    c = [abs(randq()) for _ in 0:n-1]
    lhs = sum(BigInt(l) * (l + 1) * (BigInt(l) * (l + 1) - 2) * c[l+1] for l in 0:n-1)
    rhs = sum((1 + BigInt(l) * (l + 1))^2 * c[l+1] for l in 0:n-1)
    check("3 spectral sum bound (harmonic eigenvalues)", lhs <= rhs)
end
for nm in ["3 l(l+1)(l(l+1)-2) = (l-1)l(l+1)(l+2)", "3 H_l >= 0, = 0 iff l <= 1",
           "3 lam(lam-2) <= (1+lam)^2, gap 1+4lam", "3 spectral sum bound (harmonic eigenvalues)"]
    pass(nm)
end

# ---------------------------------------------------------------------------
# 4a. Ordered-limit residues: two expansions of 1/(w(z+w))
# ---------------------------------------------------------------------------
eA(i, j) = (i >= 0 && j == -i - 2) ? BigInt(-1)^i : BigInt(0)
eB(i, j) = (j + 1 >= 0 && i == -j - 2) ? BigInt(-1)^(j + 1) : BigInt(0)
for (nm, e) in (("eA", eA), ("eB", eB)), i in -25:25, j in -25:25
    check("4 (zw+w^2) e = 1 for both expansions", e(i - 1, j - 1) + e(i, j - 2) == ((i == 0 && j == 0) ? 1 : 0))
end
for trial in 1:200
    N = rand(rng, 1:8)
    p = [randq() for a in 1:N, b in 1:N]
    p[1, 1] = trial % 2 == 0 ? q(1) : p[1, 1]
    rA = sum(p[a+1, b+1] * eA(-1 - a, -1 - b) for a in 0:N-1, b in 0:N-1)
    rB = sum(p[a+1, b+1] * eB(-1 - a, -1 - b) for a in 0:N-1, b in 0:N-1)
    check("4 Res_w Res_z = 0, Res_z Res_w = p(0,0)", rA == 0 && rB == p[1, 1])
end
for m in 0:5, T in 0:9
    tuples = Vector{Vector{Int}}()
    function rec(prefix, left, k)
        if k == 0; push!(tuples, copy(prefix)); return; end
        for u in 0:left; push!(prefix, u); rec(prefix, left - u, k - 1); pop!(prefix); end
    end
    rec(Int[], T, m)
    brute = [t for t in Iterators.product(fill(0:T, m)...) if sum(t; init=0) <= T]
    check("4 #gap tuples = C(T+m,m)", length(tuples) == C(T + m, m) && length(brute) == C(T + m, m))
    # ordered times 0 <= t_1 <= ... <= t_m <= T via partial sums: a bijection
    times = Set([cumsum(t) for t in tuples])
    incr = [t for t in Iterators.product(fill(0:T, m)...) if all(t[i] <= t[i+1] for i in 1:m-1)]
    check("4 #ordered tuples = C(T+m,m) (partial sums)", length(times) == length(tuples) == length(incr))
    if m >= 1
        w = [randq() for _ in 0:T]
        lhs = sum(w[t[1]+1] for t in tuples)
        rhs = sum(w[u+1] * C(T - u + m - 1, m - 1) for u in 0:T)
        check("4 gap identity sum w(d_0) = sum_u w(u) C(T-u+m-1,m-1)", lhs == rhs)
        j = rand(rng, 1:m)   # any gap may be singled out
        check("4 gap identity for an arbitrary gap index", sum(w[t[j]+1] for t in tuples) == rhs)
    end
end
for nm in ["4 (zw+w^2) e = 1 for both expansions", "4 Res_w Res_z = 0, Res_z Res_w = p(0,0)",
           "4 #gap tuples = C(T+m,m)",
           "4 #ordered tuples = C(T+m,m) (partial sums)",
           "4 gap identity sum w(d_0) = sum_u w(u) C(T-u+m-1,m-1)", "4 gap identity for an arbitrary gap index"]
    pass(nm)
end

# ---------------------------------------------------------------------------
# 5. Chebyshev: leading coefficients and exact weights
# ---------------------------------------------------------------------------
function chebcoeffs(M)
    T = [BigInt[1], BigInt[0, 1]]
    for n in 2:M
        a = vcat(0, 2 .* T[n]); b = vcat(T[n-1], zeros(BigInt, length(a) - length(T[n-1])))
        push!(T, a .- b)
    end
    return T
end
Tc = chebcoeffs(60)
compose_2u_minus_1(p) = begin   # coefficients of p(2u-1)
    res = zeros(BigInt, length(p)); pw = BigInt[1]
    for k in 1:length(p)
        res[1:length(pw)] .+= p[k] .* pw
        pw = vcat(0, 2 .* pw) .- vcat(pw, 0)
    end
    res
end
for m in 1:60
    check("5 lead T_m = 2^(m-1), deg T_m = m", Tc[m+1][end] == BigInt(2)^(m - 1) && length(Tc[m+1]) == m + 1)
    check("5 lead T_m(2u-1) = 2^(2m-1)", compose_2u_minus_1(Tc[m+1])[end] == BigInt(2)^(2m - 1))
    check("5 T_m(1) = 1, T_m(-1) = (-1)^m", evalpoly(Tc[m+1], 1) == 1 && evalpoly(Tc[m+1], -1) == (-1)^m)
end
pass("5 lead T_m = 2^(m-1), deg T_m = m"); pass("5 lead T_m(2u-1) = 2^(2m-1)"); pass("5 T_m(1) = 1, T_m(-1) = (-1)^m")

# Exact quadratic field Q(sqrt D): x = a + b sqrt(D)
struct QF; a::Q; b::Q; D::Int; end
QF(a, D) = QF(Q(a), q(0), D)
import Base: +, -, *, /, ==, zero, one, ^
+(x::QF, y::QF) = QF(x.a + y.a, x.b + y.b, x.D)
-(x::QF, y::QF) = QF(x.a - y.a, x.b - y.b, x.D)
-(x::QF) = QF(-x.a, -x.b, x.D)
*(x::QF, y::QF) = QF(x.a * y.a + x.D * x.b * y.b, x.a * y.b + x.b * y.a, x.D)
/(x::QF, y::QF) = (n = y.a^2 - y.D * y.b^2; x * QF(y.a / n, -y.b / n, y.D))
==(x::QF, y::QF) = x.a == y.a && x.b == y.b
^(x::QF, n::Int) = n == 0 ? QF(1, x.D) : x * x^(n - 1)
function sgn(x::QF)    # exact sign of a + b sqrt(D), D > 0
    sa, sb = sign(x.a), sign(x.b)
    sa == 0 && return sb; sb == 0 && return sa; sa == sb && return sa
    return x.a^2 > x.D * x.b^2 ? sa : sb
end
# exact Chebyshev nodes cos(j pi/m), m = 1..6
r2 = QF(q(0), q(1, 2), 2)    # sqrt(2)/2
s3 = QF(q(0), q(1, 2), 3)    # sqrt(3)/2
c5a = QF(q(1, 4), q(1, 4), 5) # cos(pi/5)  = (1+sqrt5)/4
c5b = QF(q(-1, 4), q(1, 4), 5) # cos(2pi/5) = (sqrt5-1)/4
nodes = Dict(
    1 => [QF(1, 2), QF(-1, 2)],
    2 => [QF(1, 2), QF(0, 2), QF(-1, 2)],
    3 => [QF(1, 2), QF(q(1, 2), 2), QF(q(-1, 2), 2), QF(-1, 2)],
    4 => [QF(1, 2), r2, QF(0, 2), -r2, QF(-1, 2)],
    5 => [QF(1, 5), c5a, c5b, -c5b, -c5a, QF(-1, 5)],
    6 => [QF(1, 3), s3, QF(q(1, 2), 3), QF(0, 3), QF(q(-1, 2), 3), -s3, QF(-1, 3)])
evalT(m, x) = (acc = QF(0, x.D); for k in 0:m; acc = acc + QF(Q(Tc[m+1][k+1]), x.D) * x^k; end; acc)
for m in 1:6
    xs = nodes[m]; D = xs[1].D
    # check nodes really are cos(j pi/m): T_m(x_j) = (-1)^j and the nodes are distinct
    for j in 0:m
        check("5 T_m(x_j) = (-1)^j at exact nodes (m<=6)", evalT(m, xs[j+1]) == QF((-1)^j, D))
    end
    for (L, T) in [(q(0), q(1)), (q(-1), q(2)), (randq(), randq(1, 9)), (randq(), randq(1, 9)), (q(3, 7), q(5, 3))]
        ts = [QF(L, D) + QF(T, D) * (QF(1, D) + x) / QF(2, D) for x in xs]
        w = [QF(1, D) / prod(ts[j] - ts[k] for k in 1:m+1 if k != j) for j in 1:m+1]
        check("5 sign w_j = (-1)^j", all(sgn(w[j]) == (-1)^(j - 1) for j in 1:m+1))
        absum = sum(sgn(w[j]) == 1 ? w[j] : -w[j] for j in 1:m+1)
        check("5 sum |w_j| = K_m = 2^(2m-1)/T^m", absum == QF(q(2)^(2m - 1) / T^m, D))
        for trial in 1:20
            c = [randq() for _ in 0:m]
            P(t) = (acc = QF(0, D); for k in 0:m; acc = acc + QF(c[k+1], D) * t^k; end; acc)
            check("5 sum w_j P(t_j) = a_m (deg P <= m)", sum(w[j] * P(ts[j]) for j in 1:m+1) == QF(c[m+1], D))
        end
    end
end
# the Lean-stated weights for m = 1..4 on [-1,1]
leanw = Dict(1 => [q(1, 2), q(-1, 2)], 2 => [q(1, 2), q(-1), q(1, 2)],
             3 => [q(2, 3), q(-4, 3), q(4, 3), q(-2, 3)], 4 => [q(1), q(-2), q(2), q(-2), q(1)])
for m in 1:4
    xs = nodes[m]; D = xs[1].D
    w = [QF(1, D) / prod(xs[j] - xs[k] for k in 1:m+1 if k != j) for j in 1:m+1]
    check("5 Lean weights for m=1..4 equal Lagrange weights", all(w[j] == QF(leanw[m][j], D) for j in 1:m+1))
end
for nm in ["5 T_m(x_j) = (-1)^j at exact nodes (m<=6)", "5 sign w_j = (-1)^j", "5 sum |w_j| = K_m = 2^(2m-1)/T^m",
           "5 sum w_j P(t_j) = a_m (deg P <= m)", "5 Lean weights for m=1..4 equal Lagrange weights"]
    pass(nm)
end

# ---------------------------------------------------------------------------
# 6. Full-jet weights W_{k,j} (the paper (mdpi-preprints/CelestialIRPaper.pdf), eq:fullweights) reproduce beta_k
# ---------------------------------------------------------------------------
polymul(p, r) = (res = zeros(Q, length(p) + length(r) - 1);
                 for i in eachindex(p), j in eachindex(r); res[i+j-1] += p[i] * r[j]; end; res)
function lagrange_l(j, m, L, d)   # l_j(t) = prod_{r != j} (t - L - r d)/(d (j - r))
    p = Q[1]
    for r in 0:m
        r == j && continue
        p = polymul(p, Q[-(L + r * d), q(1)]) ./ (d * (j - r))
    end
    p
end
coef(p, n) = n + 1 <= length(p) ? p[n+1] : q(0)
Wkj(k, j, m, L, d, s2) = sum((-s2 / 2)^r / fact(r) * fact(k + 2r) * coef(lagrange_l(j, m, L, d), k + 2r)
                             for r in 0:div(m - k, 2))
# Gaussian average E[(t - sigma Z)^n] via normal moments E Z^(2r) = (2r)!/(2^r r!)
gauss_moment(n, t, s2) = sum(C(n, 2r) * t^(n - 2r) * s2^r * fact(2r) // (BigInt(2)^r * fact(r)) for r in 0:div(n, 2))
for trial in 1:250
    m = rand(rng, 1:5)
    u = [randq() for _ in 0:m]
    L, d, s2 = randq(), randq(1, 9), abs(randq(1, 9))
    E = [sum(u[n+1] * gauss_moment(n, L + j * d, s2) for n in 0:m) for j in 0:m]
    for k in 0:m
        bhat = sum(Wkj(k, j, m, L, d, s2) * E[j+1] for j in 0:m)
        check("6 full-jet weights: sum_j W_kj E_j = beta_k = k! u_k", bhat == fact(k) * u[k+1])
    end
    if m >= 2   # wrong-sign heat control must fail
        k = m - 2
        bad = sum(sum((s2 / 2)^r / fact(r) * fact(k + 2r) * coef(lagrange_l(j, m, L, d), k + 2r)
                      for r in 0:div(m - k, 2)) * E[j+1] for j in 0:m)
        check("6 control: wrong heat sign gives error m! s2 u_m", bad - fact(k) * u[k+1] == fact(m) * s2 * u[m+1])
    end
    # top coefficient: W_{m,j} = m! * (forward-difference weight)/d^m
    check("6 top weights W_mj = (-1)^(m-j) C(m,j)/d^m", all(Wkj(m, j, m, L, d, s2) == (-1)^(m - j) * C(m, j) * fact(m) // (fact(m) * d^m) for j in 0:m))
end
pass("6 full-jet weights: sum_j W_kj E_j = beta_k = k! u_k")
pass("6 control: wrong heat sign gives error m! s2 u_m")
pass("6 top weights W_mj = (-1)^(m-j) C(m,j)/d^m")

println("ALL PASS  total exact checks: ", sum(values(counts)))
