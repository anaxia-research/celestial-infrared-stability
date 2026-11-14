#!/usr/bin/env julia
# Independent finite exact certificate for full Laurent-jet recovery.
# Requires only Julia's standard library. This is not a Lean formalization
# and does not verify the analytical Mellin or minimax statements.

using Dates
using TOML

const Q = Rational{BigInt}
rat(a::Integer, b::Integer=1) = BigInt(a) // BigInt(b)
fact(n::Int) = factorial(BigInt(n))

# Polynomials are coefficient vectors in ascending powers.
function polynomial_product(p::Vector{Q}, q::Vector{Q})
    result = fill(rat(0), length(p) + length(q) - 1)
    for i in eachindex(p), j in eachindex(q)
        result[i+j-1] += p[i] * q[j]
    end
    return result
end

function polynomial_value(p::Vector{Q}, t::Q)
    value = rat(0)
    for coefficient in reverse(p)
        value = value * t + coefficient
    end
    return value
end

function derivative(p::Vector{Q}, order::Int)
    result = copy(p)
    for iteration in 1:order
        if length(result) == 1
            return Q[rat(0)]
        end
        result = Q[j * result[j+1] for j in 1:(length(result)-1)]
    end
    return result
end

function lagrange_basis(nodes::Vector{Q}, index::Int)
    result = Q[rat(1)]
    for other in eachindex(nodes)
        if other != index
            result = polynomial_product(result, Q[-nodes[other], rat(1)])
            result ./= nodes[index] - nodes[other]
        end
    end
    return result
end

# Independent route 1: average each monomial by exact normal moments.
# E[Z^(2r)] = (2r)!/(2^r*r!), and all odd moments vanish.
# Thus this computes E[P(t-sigma*Z)] without applying a heat operator.
function averaged_monomial(power::Int, t::Q, sigma_squared::Q)
    value = rat(0)
    for r in 0:div(power, 2)
        normal_moment = rat(fact(2*r), BigInt(2)^r * fact(r))
        value += rat(binomial(BigInt(power), 2*r)) * t^(power-2*r) *
                 sigma_squared^r * normal_moment
    end
    return value
end

function averaged_polynomial(p::Vector{Q}, t::Q, sigma_squared::Q)
    return sum(p[n+1] * averaged_monomial(n, t, sigma_squared)
               for n in 0:(length(p)-1))
end

# Independent route 2: undo averaging on the Lagrange polynomials using
# actual polynomial derivatives. sign=+1 is a deliberate wrong-sign control.
function heat_polynomial(p::Vector{Q}, sigma_squared::Q, sign::Int)
    result = fill(rat(0), length(p))
    for r in 0:div(length(p)-1, 2)
        term = derivative(p, 2*r)
        factor = (rat(sign) * sigma_squared / 2)^r / fact(r)
        for index in eachindex(term)
            result[index] += factor * term[index]
        end
    end
    return result
end

function hermite_bound(n::Int, A::Q, sigma_squared::Q)
    return sum((sigma_squared/2)^r * A^(n-2*r) /
               (fact(r) * fact(n-2*r)) for r in 0:div(n, 2))
end

function main()
    designs = [
        (name="positive_origin", L=rat(3,7), d=rat(5,3), sigma_squared=rat(2,5)),
        (name="negative_origin", L=rat(-4,5), d=rat(7,4), sigma_squared=rat(3,2)),
        (name="zero_origin", L=rat(0), d=rat(2,3), sigma_squared=rat(1,7)),
    ]
    rows = Dict{String,Any}[]
    basis_identity_count = 0
    interpolation_identity_count = 0
    mixed_coefficient_count = 0
    weight_bound_count = 0
    wrong_sign_count = 0
    missing_factorial_count = 0

    for design in designs, m in 1:6
        L, d, sigma_squared = design.L, design.d, design.sigma_squared
        nodes = Q[L + j*d for j in 0:m]
        bases = [lagrange_basis(nodes, j) for j in eachindex(nodes)]
        inverse_bases = [heat_polynomial(basis, sigma_squared, -1) for basis in bases]
        wrong_sign_bases = [heat_polynomial(basis, sigma_squared, 1) for basis in bases]

        # Check the interpolation construction separately before recovery.
        for j in eachindex(nodes), other in eachindex(nodes)
            expected = j == other ? rat(1) : rat(0)
            @assert polynomial_value(bases[j], nodes[other]) == expected
            interpolation_identity_count += 1
        end

        # Check the full linear recovery map on every monomial, not merely
        # on a single chosen polynomial. beta_k = k!*u_k.
        for power in 0:m
            observations = Q[averaged_monomial(power, t, sigma_squared) for t in nodes]
            for k in 0:m
                recovered = fact(k) * sum(inverse_bases[j][k+1] * observations[j]
                                          for j in eachindex(nodes))
                expected = power == k ? rat(fact(k)) : rat(0)
                @assert recovered == expected
                basis_identity_count += 1
            end
        end

        # A mixed polynomial exercises all recovered coefficients together.
        coefficients = Q[rat((-1)^k * (k+2), k+3) for k in 0:m]
        observations = Q[averaged_polynomial(coefficients, t, sigma_squared) for t in nodes]
        recovered_beta = Q[]
        target_beta = Q[]
        A = abs(L) + m*d
        for k in 0:m
            recovered = fact(k) * sum(inverse_bases[j][k+1] * observations[j]
                                      for j in eachindex(nodes))
            target = fact(k) * coefficients[k+1]
            @assert recovered == target
            push!(recovered_beta, recovered)
            push!(target_beta, target)
            mixed_coefficient_count += 1

            # Explicit H_(m-k) coefficient bound from the manuscript.
            H = hermite_bound(m-k, A, sigma_squared)
            for j in eachindex(nodes)
                beta_weight = fact(k) * inverse_bases[j][k+1]
                upper = rat(binomial(BigInt(m), j-1)) * H / d^m
                @assert abs(beta_weight) <= upper
                weight_bound_count += 1
            end
        end

        if m >= 2
            # Smoothing followed by +heat instead of -heat doubles the
            # heat time. At k=m-2 its exact error is m!*sigma^2*u_m.
            k = m-2
            wrong = fact(k) * sum(wrong_sign_bases[j][k+1] * observations[j]
                                  for j in eachindex(nodes))
            expected_error = fact(m) * sigma_squared * coefficients[m+1]
            @assert wrong - target_beta[k+1] == expected_error
            @assert wrong != target_beta[k+1]
            wrong_sign_count += 1

            # Omitting the Laurent factorial recovers u_k instead of beta_k.
            k = m
            missing_factorial = sum(inverse_bases[j][k+1] * observations[j]
                                     for j in eachindex(nodes))
            @assert missing_factorial == coefficients[k+1]
            @assert missing_factorial != target_beta[k+1]
            missing_factorial_count += 1
        end

        push!(rows, Dict(
            "design"=>design.name, "degree"=>m,
            "L"=>string(L), "d"=>string(d), "sigma_squared"=>string(sigma_squared),
            "polynomial_coefficients"=>string.(coefficients),
            "gaussian_averaged_nodes"=>string.(observations),
            "target_beta"=>string.(target_beta), "recovered_beta"=>string.(recovered_beta),
        ))
    end

    @assert basis_identity_count == 417
    @assert interpolation_identity_count == 417
    @assert mixed_coefficient_count == 81
    @assert weight_bound_count == 417
    @assert wrong_sign_count == 15
    @assert missing_factorial_count == 15
    report = Dict{String,Any}(
        "schema"=>"celestial-julia-log-jet-v1",
        "generated_utc"=>string(now(UTC)), "julia_version"=>string(VERSION),
        "passed"=>true, "formal_proof"=>false,
        "arithmetic"=>"Rational{BigInt}; standard library only",
        "scope"=>"Finite exact full-jet Gaussian interpolation and inverse-heat identities, degrees 1 through 6",
        "gaussian_observation_method"=>"Direct normal-moment sum, independent of the inverse-heat derivative implementation",
        "interpolation_identity_count"=>interpolation_identity_count,
        "monomial_recovery_identity_count"=>basis_identity_count,
        "mixed_coefficient_recovery_count"=>mixed_coefficient_count,
        "explicit_weight_bound_count"=>weight_bound_count,
        "negative_controls"=>Dict("wrong_heat_sign_rejected_count"=>wrong_sign_count,
                                   "missing_Laurent_factorial_rejected_count"=>missing_factorial_count),
        "limits"=>["Finite exact identities do not establish the analytical Bochner, Mellin, or minimax theorems.",
                    "This certificate does not establish a physical soft-charge Ward identity."],
        "cases"=>rows,
    )
    output = length(ARGS)>0 ? abspath(ARGS[1]) :
             normpath(joinpath(@__DIR__, "..", "reports", "julia-log-jet.toml"))
    mkpath(dirname(output))
    open(output, "w") do io
        TOML.print(io, report; sorted=true)
    end
    println("PASS: $basis_identity_count exact monomial identities, $mixed_coefficient_count mixed coefficients, inverse-heat negative controls; $output")
end

main()
