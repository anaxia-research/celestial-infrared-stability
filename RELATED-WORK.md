# Sources and originality review

Reviewed on 15 November 2026. This is a bounded search of primary sources and the supplied chapter.

| Source | Relationship to this draft |
| --- | --- |
| Sabrina Pasterski, *A Chapter on Celestial Holography*, arXiv:2310.04932v2, especially equations 29-31 and the loop footnote | Supplies the celestial Mellin and soft-pole motivation. The local copy supplied by the user was read in full. [Paper](https://arxiv.org/abs/2310.04932v2) |
| Perimeter Institute, *Celestial Holography Initiative* | Identifies the research programme linking scattering, asymptotic symmetries, and celestial correlators. [Initiative](https://perimeterinstitute.ca/research/research-centres-and-initiatives/celestial-holography-initiative) |
| Benjamin Muntz, *Dimensional Analysis is a Gauge Theory*, arXiv:2607.21695v1, section 3 | Motivates the local-unit connection. Our calculation explicitly keeps track of a dimensionless reference energy. [Paper](https://arxiv.org/html/2607.21695v1) · [Author's publications](https://benjaminmuntz.com/publications/) |
| Jiayin Gu, Ying-Ying Li, Lian-Tao Wang, *A Note on the Analytic Structure of Celestial Amplitudes*, JHEP 08 (2022) 257 | Establishes prior work on fixed-order multipoles and on branch cuts from resummation. This prevents extending the finite-log theorem to a resummed correlator without further estimates. [Paper](https://arxiv.org/abs/2206.01757) |
| Agnese Bissi, Laura Donnay, Beniamino Valsesia, *Logarithmic doublets in CCFT*, arXiv:2407.17123 | Establishes prior celestial logarithmic multiplets. A triangular pole-jet law by itself cannot establish discovery of logarithmic celestial CFT. [Paper](https://arxiv.org/abs/2407.17123) |
| NIST DLMF, sections 5.12 and 18.5 | Standard beta integral and Chebyshev identities used in the proofs. [Beta integral](https://dlmf.nist.gov/5.12) · [Chebyshev polynomials](https://dlmf.nist.gov/18.5) |

The review found extensive prior work on instability of analytic continuation, Mellin asymptotics, and Chebyshev optimal recovery. The explicit fixed-radius completely monotone example and its exact vertical-line norm were derived here. Establishing whether their celestial formulation is new requires a broader literature review and specialist feedback. The manuscript therefore describes verified mathematical results and a candidate research synthesis.

The profiles used in the counterexample form a specified mathematical class. Positivity of these scalar functions supplies a mathematical restriction. Unitarity, crossing symmetry, gauge constraints, and a complete celestial CFT are additional requirements for a physical scattering model.

Author-name provenance is the six-person `CITATION.cff` on the `paper` branch of `creelie/4Y-M`, read during preparation. This supplies spellings, while approval of the new draft remains pending. [Citation record](https://github.com/creelie/4Y-M/blob/paper/CITATION.cff)

## Additional primary sources in the stability revision

| Source | What is established already; relationship to this revision |
| --- | --- |
| S. Pasterski and S.-H. Shao, *A Conformal Basis for Flat Space Amplitudes*, Phys. Rev. D 96 (2017) 065022. [Primary paper](https://arxiv.org/abs/1705.01027) | Principal-series completeness and the massless Mellin basis are prior work. A discontinuous extraction functional in a specified topology does not contradict completeness. |
| L. Donnay, S. Pasterski, A. Puhm, *Asymptotic Symmetries and Celestial CFT*, JHEP 09 (2020) 176. [Primary paper](https://arxiv.org/abs/2005.08990) | Analytically continued primaries and conformally soft modes have principal-series contour representations. This is the relevant existing extraction setting; our boundary theorem addresses stability for a scalar family in specified negative Sobolev norms. |
| Y. Pano, M. Borji, *Distributional Celestial Amplitudes*, arXiv:2401.08877 (2024). [Primary paper](https://arxiv.org/abs/2401.08877) | Distributional Mellin transforms are an established subject. Our use of log-energy tempered distributions does not introduce distributional celestial amplitudes. |
| H. J. Brown, Y. Grabovsky, *On feasibility of extrapolation of completely monotone functions*, arXiv:2401.15178v2 (2024). [Primary paper](https://arxiv.org/abs/2401.15178v2) | Close prior art: impossibility of left extrapolation despite complete monotonicity and power-law right extrapolation. Their observation interval and norms differ from the Mellin-line operator used here. Positivity alone is not a new instability principle. |

The new direct-data result combines a classical Gaussian smoothing identity with a standard bias/noise balance and an explicit two-profile lower bound. The candidate contribution is the precise residue-recovery observation model and its connection to the celestial problem; no claim of discovering regularization, minimax methods, or the general instability of analytic continuation is made. Priority of the exact combined theorems remains unresolved.
