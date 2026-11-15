# Proof status and research targets

Updated 15 November 2026. This record distinguishes mathematical proofs from computations, physical assumptions, and open targets. It is not a claim to have completed celestial holography.

## Statements proved in the manuscript

| Statement | Data and assumptions | Proof and verification |
| --- | --- | --- |
| Fixed nonzero soft residue with vanishing line data | Scalar beta-profile family; a vertical line strictly inside its convergence strip | Exact beta integral, positivity, gamma asymptotics; exact Python/Julia controls |
| Instability on the leading-soft boundary contour | Abel boundary distributions; the standard frequency-variable norm `H^(-q)`, `q > 1/2` | Fourier representation, exact weighted norm, dominated convergence, and a two-tail asymptotic argument; numerical quadrature controls |
| Direct soft-residue recovery with optimal exponent | Complete displaced-line data with bounded deterministic noise; global endpoint modulus `|f(x)-f(0)| <= M x^alpha`; weighted integrability | Explicit Gaussian estimator, rigorous bias/noise bound, and matching two-profile minimax lower exponent; high-precision identity and negative controls |
| Finite frequency-band recovery bound | Continuum-band observations and an additional known bound on the weighted energy integral | Gaussian tail bound; numerical truncation controls. No discrete sampling guarantee |
| Angularly smeared coefficient recovery | Banach-valued stripped kernels, global endpoint prior, bounded dual test, and Mellin-line noise | Bochner Gaussian identity and duality; exact spherical-mode and numerical controls |
| Direct highest soft-log pole recovery | Fixed log degree and global remainder bound; arbitrary lower polynomial coefficients | Finite differences cancel lower degrees; explicit cutoff profiles give the matching minimax exponent. Constants depend on degree |
| Full principal-part recovery | Same fixed-degree Banach class, unconstrained lower polynomial coefficients, fixed convention | Exact interpolation/inverse-Gaussian weights and matching marginal minimax rates, including logarithmic losses for lower poles; independent Python and rational Julia finite-identity replays |
| Normalized gravitational soft charge and time-window error | Explicit smooth classical Bondi news, fixed angular test | Sphere integrations, exact Lorentzian/Fourier/Mellin formulas and minimax exponent under the uniform modulus |
| Recovery of a prescribed fractional branch coefficient | Known template parameters, global power remainder bound, displaced-line noise | Normalized Gaussian recovery and matching cutoff-profile lower rate |
| Noncommuting ordered residues | Explicit scalar two-energy kernel | Exact double Mellin formula, mixed divisor, and iterated residues |
| Bounded two-channel transform change | Positive-energy compact smooth profiles, both signs of retarded time, real Mellin weight 1/2 | Gamma reflection formula and sharp linear norm. One-channel unboundedness by modulation |
| Shifted-pole and branch-point model transforms | The two explicit scalar kernels in the paper | Complete integral derivations and numerical checks. They illustrate failure of the finite-polynomial resummation inference |
| Sharp highest soft-log coefficient recovery | Energy-space observations; a fixed-degree polynomial on the measured interval | Chebyshev alternation and an indistinguishable pair; exact algebraic controls. Pole interpretation also needs the stated infrared tail condition |
| Numerical-energy pole-jet transformation and reference cancellation | Energy-independent local changes of units; stated integration and continuation assumptions | Written analytical derivation, exact symbolic controls, and limited integer double-pole Lean formalization |

The new analytical theorems, including the angular and logarithmic extensions, are **not formalized in Lean**. The existing `lean/ScaleJets.lean` proves only the documented symbolic double-pole identities. A successful computation is not a substitute for the written infinite-dimensional proof.

## What changed in this revision

The previous draft left its leading-soft boundary contour untreated and obtained its constructive estimator from energy-space observations. The new sections give a boundary result in a specified distribution norm and a separate estimator using Mellin observations directly. They use different data models, which the paper does not identify with one another.

At displaced-line distance `a` and endpoint exponent `alpha`, the direct-data minimax error is bounded above and below by explicit multiples of

`M^(a/(a+alpha)) eta^(alpha/(a+alpha))`.

The error exponent matches. The upper and lower constants do not match, and no exact minimax-constant claim is made.

## Progress on the two selected targets

1. **Stable recovery of angular soft data:** proved for bounded angular pairings in the specified Banach-space observation model. A normalized physical soft-charge interpretation is conditional on an independently established charge map and a controlled hard contribution.
2. **Loop corrections and logarithmic modes:** proved a direct-data recovery theorem and matching minimax exponent for every fixed log degree under the global remainder prior. The entire principal part is also recoverable with matching componentwise logarithmic losses for lower poles. Existing loop operators and Ward identities are credited to the primary sources. Two explicit resummation models prove why finite polynomial truncations cannot certify all-orders singularity structure.

The new estimates do not close all problems of celestial holography, Sabrina Pasterski, or her coauthors. The [bounded problem map](OPEN-PROBLEMS.md) and [calculation ledger](RESEARCH-ATTEMPTS.md) give the status of the eight research topics.

## Open targets with concrete completion criteria

The [massive angular theorem](paper/massive-angular.tex) proves two-sided spectral bounds and a matching minimax recovery rate for the leading massive QED angular kernel at a known speed. The [translation theorem](paper/mellin-translation.tex) proves sharp band leakage for physical time translations on the massless scalar principal series. Both are written analytical proofs. The massive Ward tower, interacting higher-loop amplitudes, and a complete holographic duality are still open.

| Target | What a proof would have to supply | Current status |
| --- | --- | --- |
| A physical scattering realization | A family satisfying the relevant unitarity, crossing, gauge, and kinematic constraints, plus the required Mellin bounds | Open. Scalar positivity and complete monotonicity do not establish these conditions |
| Bounds for physical soft charges | Identify the normalized charge map and hard contribution in a physical theory; verify its angular and Mellin priors | A conventional Bondi soft functional is proved for specified classical radiative data. Global scattering, hard-charge balance and a physical instrument noise model remain open |
| Loop-resummed data | Recover the actual branch-cut/spectral structure with a justified QFT remainder and renormalization prescription | Known-template branch recovery is proved. Unknown physical templates and an all-orders theorem remain open |
| Certified discrete recovery | A sampling scheme with a rigorously bounded quadrature and truncation error in addition to the continuum-band bound | Open. Numerical examples do not supply a uniform sampling theorem |
| Full formalization | Mellin/Fourier integration, distribution boundary values, the Sobolev limit, and the minimax inequalities in a proof assistant | Open. Existing Lean scope is deliberately much smaller |
| Originality and publication readiness | Systematic comparison with the nearest theorems, specialist review, and approval by every listed author | Pending. No priority or breakthrough certification |

The [principal-series completeness result](https://arxiv.org/abs/1705.01027), [contour treatment of conformally soft modes](https://arxiv.org/abs/2005.08990), and [distributional celestial transform](https://arxiv.org/abs/2401.08877) are existing work. Our stability statements do not refute them. [Brown and Grabovsky](https://arxiv.org/abs/2401.15178v2) provide close prior work on extrapolation of completely monotone functions, with different observations and norms.

No external submission, release, DOI registration, or contact with cited researchers is part of this revision.
