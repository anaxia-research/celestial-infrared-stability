# Proof status and research targets

Updated 15 November 2026. This record distinguishes mathematical proofs from computations, physical assumptions, and open targets. It is not a claim to have completed celestial holography.

## Statements proved in the manuscript

| Statement | Data and assumptions | Proof and verification |
| --- | --- | --- |
| Fixed nonzero soft residue with vanishing line data | Scalar beta-profile family; a vertical line strictly inside its convergence strip | Exact beta integral, positivity, gamma asymptotics; exact Python/Julia controls |
| Instability on the leading-soft boundary contour | Abel boundary distributions; the standard frequency-variable norm `H^(-q)`, `q > 1/2` | Fourier representation, exact weighted norm, dominated convergence, and a two-tail asymptotic argument; numerical quadrature controls |
| Direct soft-residue recovery with optimal exponent | Complete displaced-line data with bounded deterministic noise; global endpoint modulus `|f(x)-f(0)| <= M x^alpha`; weighted integrability | Explicit Gaussian estimator, rigorous bias/noise bound, and matching two-profile minimax lower exponent; high-precision identity and negative controls |
| Finite frequency-band recovery bound | Continuum-band observations and an additional known bound on the weighted energy integral | Gaussian tail bound; numerical truncation controls. No discrete sampling guarantee |
| Sharp highest soft-log coefficient recovery | Energy-space observations; a fixed-degree polynomial on the measured interval | Chebyshev alternation and an indistinguishable pair; exact algebraic controls. Pole interpretation also needs the stated infrared tail condition |
| Numerical-energy pole-jet transformation and reference cancellation | Energy-independent local changes of units; stated integration and continuation assumptions | Written analytical derivation, exact symbolic controls, and limited integer double-pole Lean formalization |

The new analytical theorems are **not formalized in Lean**. The existing `lean/ScaleJets.lean` proves only the documented symbolic double-pole identities. A successful computation is not a substitute for the written infinite-dimensional proof.

## What changed in this revision

The previous draft left its leading-soft boundary contour untreated and obtained its constructive estimator from energy-space observations. The new sections give a boundary result in a specified distribution norm and a separate estimator using Mellin observations directly. They use different data models, which the paper does not identify with one another.

At displaced-line distance `a` and endpoint exponent `alpha`, the direct-data minimax error is bounded above and below by explicit multiples of

`M^(a/(a+alpha)) eta^(alpha/(a+alpha))`.

The error exponent matches. The upper and lower constants do not match, and no exact minimax-constant claim is made.

## Open targets with concrete completion criteria

| Target | What a proof would have to supply | Current status |
| --- | --- | --- |
| A physical scattering realization | A family satisfying the relevant unitarity, crossing, gauge, and kinematic constraints, plus the required Mellin bounds | Open. Scalar positivity and complete monotonicity do not establish these conditions |
| Bounds for observable soft charges | Specify the angular/kinematic smearing and the physical measurement norm, then prove its relation to this paper's norms | Open. The residue is a scalar model of soft data |
| Loop-resummed data | Replace a finite-order pole ansatz by the actual branch-cut/spectral structure and establish remainder bounds | Open. Finite soft-log proofs do not establish resummation |
| Certified discrete recovery | A sampling scheme with a rigorously bounded quadrature and truncation error in addition to the continuum-band bound | Open. Numerical examples do not supply a uniform sampling theorem |
| Full formalization | Mellin/Fourier integration, distribution boundary values, the Sobolev limit, and the minimax inequalities in a proof assistant | Open. Existing Lean scope is deliberately much smaller |
| Originality and publication readiness | Systematic comparison with the nearest theorems, specialist review, and approval by every listed author | Pending. No priority or breakthrough certification |

The [principal-series completeness result](https://arxiv.org/abs/1705.01027), [contour treatment of conformally soft modes](https://arxiv.org/abs/2005.08990), and [distributional celestial transform](https://arxiv.org/abs/2401.08877) are existing work. Our stability statements do not refute them. [Brown and Grabovsky](https://arxiv.org/abs/2401.15178v2) provide close prior work on extrapolation of completely monotone functions, with different observations and norms.

No external submission, release, DOI registration, or contact with cited researchers is part of this revision.
