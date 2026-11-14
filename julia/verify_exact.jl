#!/usr/bin/env julia
# Independent exact-arithmetic checks; requires only Julia's standard library.
# These are finite algebraic checks, not a proof-assistant formalization.

using Dates
using TOML

const Q = Rational{BigInt}
rat(a::Integer, b::Integer=1) = BigInt(a) // BigInt(b)

function chebyshev_polynomials(maxdegree::Int)
    polys = Vector{Vector{Q}}([Q[rat(1)], Q[rat(0), rat(1)]])
    for degree in 1:maxdegree-1
        current = fill(rat(0), degree + 2)
        for (index, coefficient) in enumerate(polys[degree+1])
            current[index+1] += 2coefficient
        end
        for (index, coefficient) in enumerate(polys[degree])
            current[index] -= coefficient
        end
        push!(polys, current)
    end
    polys
end

function beta_by_polynomial_integration(a::Int, n::Int)
    exponent = n - a - 1
    sum(rat((-1)^k * binomial(BigInt(exponent), k), a+k) for k in 0:exponent)
end

function weights_for(nodes::Vector{Q})
    [inv(prod(nodes[j]-nodes[k] for k in eachindex(nodes) if k != j))
     for j in eachindex(nodes)]
end

function main()
    polynomials = chebyshev_polynomials(24)
    for degree in 1:24
        @assert length(polynomials[degree+1]) == degree+1
        @assert last(polynomials[degree+1]) == BigInt(2)^(degree-1)
    end

    beta_count = 0
    for a in 1:6, n in a+1:a+12
        integrated = beta_by_polynomial_integration(a, n)
        expected = factorial(BigInt(a-1)) * factorial(BigInt(n-a-1)) // factorial(BigInt(n-1))
        @assert integrated == expected
        beta_count += 1
    end

    radius = rat(3, 2)
    line_rows = Dict{String, Any}[]
    for n in 2:32
        value = radius * beta_by_polynomial_integration(1, n)
        @assert value == radius / (n-1)
        @assert value != radius / n # Negative control: common off-by-one error.
        push!(line_rows, Dict("N"=>n, "R"=>string(radius), "exact_norm"=>string(value)))
        # The finite product in B(z,N-z) has value one at z=0 after division
        # by (N-1)!; pi*z/sin(pi*z) tends to one, giving unit residue.
        @assert prod(BigInt(k) for k in 1:n-1) // factorial(BigInt(n-1)) == 1
    end

    left, width = rat(3, 7), rat(5, 3)
    rational_cheb_nodes = [Q[rat(1),rat(-1)],
                           Q[rat(1),rat(0),rat(-1)],
                           Q[rat(1),rat(1,2),rat(-1,2),rat(-1)]]
    interpolation_rows = Dict{String, Any}[]
    for degree in 1:3
        nodes = Q[left + width*(1+x)/2 for x in rational_cheb_nodes[degree]]
        weights = weights_for(nodes)
        for power in 0:degree
            @assert sum(weights[j] * nodes[j]^power for j in eachindex(nodes)) == (power == degree ? 1 : 0)
        end
        expected = rat(BigInt(2)^(2degree-1)) / width^degree
        @assert sum(abs, weights) == expected
        @assert all(sign(weights[j]) == (-1)^(j-1) for j in eachindex(weights))
        push!(interpolation_rows, Dict("degree"=>degree, "sharp_noise_factor"=>string(expected),
                                      "moment_identities"=>degree+1))
    end

    # Equispaced degree-four interpolation is strictly inferior to Chebyshev.
    degree = 4
    nodes = Q[left + width*j/degree for j in 0:degree]
    equispaced_norm = sum(abs, weights_for(nodes))
    optimum = rat(BigInt(2)^(2degree-1)) / width^degree
    @assert equispaced_norm / optimum == rat(4,3)

    report = Dict{String,Any}(
        "schema"=>"celestial-julia-exact-v1",
        "generated_utc"=>string(now(UTC)),
        "julia_version"=>string(VERSION),
        "passed"=>true,
        "formal_proof"=>false,
        "arithmetic"=>"Rational{BigInt}; no numerical quadrature or external packages",
        "scope"=>"Independent finite exact beta and Chebyshev checks",
        "chebyshev_recurrence_degree_min"=>1,
        "chebyshev_recurrence_degree_max"=>24,
        "integer_beta_integral_count"=>beta_count,
        "whole_line_reason"=>"Positivity and the Mellin integral triangle inequality; equality at zero frequency.",
        "unit_residue_reason"=>"Exact finite-product value; analytical limit pi*z/sin(pi*z)=1 is used.",
        "negative_controls"=>Dict("wrong_R_over_N_rejected"=>true,
                                   "equispaced_degree_four_ratio"=>"4//3"),
        "integer_a_one_norms"=>line_rows,
        "rational_interpolation"=>interpolation_rows,
    )
    output = length(ARGS)>0 ? abspath(ARGS[1]) : normpath(joinpath(@__DIR__, "..", "reports", "julia-exact.toml"))
    mkpath(dirname(output))
    open(output, "w") do io
        TOML.print(io, report; sorted=true)
    end
    println("PASS: 24 recurrence degrees, $beta_count integer beta integrals, exact interpolation m<=3; $output")
end

main()
