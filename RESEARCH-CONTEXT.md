# Research context and claim boundaries

Bounded primary-source audit, 15 November 2026. This document supports the scalar and Banach-valued recovery statements in the manuscript. It establishes relevant antecedents, not a priority certificate or an exhaustive literature review.

## Celestial dictionary and loop corrections

| Primary source | Relevance and boundary |
| --- | --- |
| Sabrina Pasterski, [*A Chapter on Celestial Holography*, 2310.04932v2](https://arxiv.org/html/2310.04932v2), “Holographic Dictionary / Soft Charges as Currents,” equations 29–31 | The soft expansion/Mellin-pole dictionary motivates the target. The chapter's loop qualification does not provide a uniform remainder bound or a stable noisy-data reconstruction theorem. |
| Laura Donnay, Sabrina Pasterski and Andrea Puhm, [*Asymptotic Symmetries and Celestial CFT*, 2005.08990](https://arxiv.org/pdf/2005.08990), §§4.1–4.5 and §5.3.3 | Regulated principal-series contour superpositions already explain analytically continued conformal primaries; soft charges require symplectic and contact-term analysis. Exact representability does not imply stability in the noise norm used here. |
| Biswajit Sahoo and Ashoke Sen, [*Classical and Quantum Results on Logarithmic Terms in the Soft Theorem in Four Dimensions*, 1808.03288v3](https://arxiv.org/abs/1808.03288v3), JHEP 02 (2019) 086, §§2, 5–7 | Supplies physical logarithmic soft terms from classical scattering and one-loop amplitudes. Much of the quantum calculation concerns charged scalars with electromagnetic/gravitational interactions; extending it to another amplitude requires its own hypotheses. |
| Rishabh Bhardwaj and Akshay Yelleshpur Srikant, [*Celestial soft currents at one-loop and their OPEs*, 2403.10443](https://arxiv.org/pdf/2403.10443), JHEP 07 (2024) 034, §§3–5 | One-loop Yang–Mills has double and triple conformal-dimension poles and corresponding soft operators. Multiple-soft limits require care. Our integer **m is log degree, not loop order**; the recovery theorem does not discover these operators or their OPEs. |
| Shreyansh Agrawal, Laura Donnay, Kevin Nguyen and Romain Ruzziconi, [*Logarithmic soft graviton theorems from superrotation Ward identities*, 2309.11220](https://arxiv.org/abs/2309.11220), JHEP 02 (2024) 120 | Their Ward-identity construction uses radiative phase space, dressed fields and fluxes at null/timelike infinity. Those physical inputs are stronger than the functional-analytic coefficient class in this repository. |
| Sangmin Choi, Ameya Kadhe and Andrea Puhm, [*Long-Range Interactions in Celestial CFT*, 2601.05951v1](https://arxiv.org/html/2601.05951v1), §4.1.2, equations 4.8–4.11; §§5–6 | Defines loop soft operators through higher Mellin poles for the universal QED/gravity logarithmic tower. Symmetry-governed terms exponentiate, and ordered limits matter. It does not establish that arbitrary resummed amplitudes possess finite Laurent jets or our global remainder prior. |
| Lorenzo Magnea and Enrico Zunino, [*Non-abelian soft radiation data for a celestial theory*, 2512.22104v1](https://arxiv.org/html/2512.22104v1), §4.1, equations 40–43, §4.2 and §7 | Examines non-abelian single/multiple-emission data. Single-current logarithms can be absorbed into dipole-dependent choices of the dimensionally regulated running-coupling scale. Multiple-emission limits retain energy fractions and ordering issues. This qualifies claims that an observed log necessarily requires new logarithmic CCFT degrees of freedom; it is not a general no-go theorem for logarithmic CCFTs. |
| Jiayin Gu, Ying-Ying Li and Lian-Tao Wang, [*A Note on the Analytic Structure of Celestial Amplitudes*, 2206.01757](https://arxiv.org/pdf/2206.01757), JHEP 08 (2022) 257, §§2–5 | Leading-log resummation in massless EFTs generically yields branch cuts, with special shifted-pole or terminating cases. Fixed-order multipoles cannot be promoted to an all-orders pole theorem. The elementary shifted-pole and incomplete-gamma examples here are scope controls, not newly computed physical loop amplitudes. |

The existing logarithmic-doublet literature, including [Bissi–Donnay–Valsesia, 2407.17123](https://arxiv.org/abs/2407.17123), must also be credited. A triangular coefficient transformation alone does not establish a new logarithmic celestial theory.

## Inverse-problem antecedents

* Flajolet–Gourdon–Dumas, [*Mellin transforms and asymptotics: harmonic sums*](https://algo.inria.fr/flajolet/Publications/mellin-harm.pdf), Theoretical Computer Science 144 (1995), 3–58: direct/converse mapping principles are classical antecedents of power-log/pole relations.
* Trefethen, [*Quantifying the ill-conditioning of analytic continuation*](https://people.maths.ox.ac.uk/trefethen/trefethen_ac.pdf), BIT 60 (2020), 901–915, Theorems 5.1–5.3: unrestricted instability and conditional continuity under a global analytic bound are established. Fixed analyticity radius without a uniform magnitude bound is materially weaker.
* Grabovsky–Hovsepyan, [*Optimal error estimates for analytic continuation in the upper half-plane*, 1812.05715](https://arxiv.org/abs/1812.05715): sharp conditional point recovery with a Hardy-class prior and noisy curve data. Compare data norms and target functionals before claiming a new optimal continuation principle.
* Brown–Grabovsky, [*On Feasibility of Extrapolation of Completely Monotone Functions*, 2401.15178](https://arxiv.org/abs/2401.15178), SIAM J. Math. Anal. 56 (2024), 7713–7747: impossible left extrapolation and conditional right extrapolation already occur within completely monotone classes. Their relative interval-data norm differs from our full Mellin-line norm and fixed-radius family.
* Gerth–Hofmann–Hofmann–Kindermann, [*The Hausdorff Moment Problem in the light of ill-posedness of type I*, 2104.06029v2](https://arxiv.org/abs/2104.06029v2): full-function and endpoint-functional recovery have different conditional stability. Merely imposing “moments” does not select a stable class; state the actual norm and bounds.

The finite-window leading-polynomial-coefficient constant follows Chebyshev extremality. The direct-Mellin logarithmic result adds a specified observation operator, global tail class, angular Banach norm, estimator and matching fixed-degree minimax exponent. That combination is a candidate mathematical contribution; originality remains subject to specialist comparison. A matching exponent does not mean sharp constants, a physical amplitude construction or a Ward-identity proof.

## When a smeared coefficient becomes a physical charge

The current target is a Laurent coefficient beta, valued in a stated angular space X. A bounded dual pairing ell(beta) is a smeared coefficient. For spin, the test belongs to the dual spin bundle; a scalar L2 pairing need not accommodate a singular collinear kernel.

To transfer a recovery estimate to a physical soft charge, independently establish:

1. The physical angular map D, normalization, gauge/helicity/color convention and reference-energy convention.
2. A dual test section D*epsilon in X*, including integration-by-parts boundary/contact terms and any patch or distributional prescriptions.
3. The hard contribution and its error, plus the admissible radiative phase space, dressings and matching conditions.
4. The relevant charge conservation or Ward identity for that theory and perturbative/resummed setting.

Only then does the coefficient error transfer with factor ||D*epsilon||, plus the independently controlled hard-sector error. The manuscript's global remainder bound is assumed, not inferred from a few energy samples. Neither the numerical replays nor formal certificates supply the missing physical inputs.
