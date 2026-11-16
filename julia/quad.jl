# Self-contained quadrature and special-function helpers (no packages).

"Gauss-Legendre nodes and weights on [-1,1] by Newton iteration."
function gausslegendre(n::Int)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cos(pi * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p3 = p2; p2 = p1
                p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            dz = p1 / pp
            z -= dz
            abs(dz) < 1e-15 && break
        end
        p1, p2 = 1.0, 0.0
        for j in 1:n
            p3 = p2; p2 = p1
            p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
        end
        pp = n * (z * p1 - p2) / (z^2 - 1)
        x[i] = z; w[i] = 2 / ((1 - z^2) * pp^2)
    end
    return x, w
end

const GLX, GLW = gausslegendre(48)

"Composite Gauss-Legendre on the breakpoints `bs` (sorted)."
function integrate(f, bs::AbstractVector{<:Real})
    s = 0.0
    for k in 1:length(bs)-1
        a, b = bs[k], bs[k+1]
        h = (b - a) / 2; m = (a + b) / 2
        for i in eachindex(GLX)
            s += GLW[i] * h * f(m + h * GLX[i])
        end
    end
    return s
end

"Integral over (0,1] with geometric breakpoints toward 0 (handles log endpoint singularities)."
integrate01(f; levels=40) = integrate(f, vcat([0.0], [2.0^(-k) for k in levels:-1:0]))

"∫_a^∞ f, via λ = a/u, dλ = a du/u^2, u ∈ (0,1]."
integrate_tail(f, a) = integrate01(u -> u == 0 ? 0.0 : f(a / u) * a / u^2)

"Lanczos approximation to Γ(x) for real x > 0."
function gammafn(x::Float64)
    g = 7.0
    c = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
    if x < 0.5
        return pi / (sin(pi * x) * gammafn(1 - x))
    end
    x -= 1
    a = c[1]; t = x + g + 0.5
    for i in 1:8
        a += c[i+1] / (x + i)
    end
    return sqrt(2pi) * t^(x + 0.5) * exp(-t) * a
end

const EULERGAMMA = 0.57721566490153286061
