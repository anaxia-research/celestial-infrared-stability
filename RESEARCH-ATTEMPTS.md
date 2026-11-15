# Attempts across the eight mapped research categories

Bounded attempt ledger, 15 November 2026. All eight categories in [OPEN-PROBLEMS.md](OPEN-PROBLEMS.md) were considered. They are broad research categories, not eight precisely posed conjectures. The outcomes below distinguish proved mathematics, replays of established results, conditional consequences, and unresolved physical obligations. They do not establish all of celestial holography or claim to solve every problem of the cited researchers.

## 1. Recover soft coefficients from Mellin observations

**Attempt and outcome — conditional mathematical results.** The manuscript constructs direct Gaussian recovery from noisy Mellin-line observations, a matching minimax exponent for a fixed global infrared modulus, and angular/fixed-log-degree extensions. In its notation the highest-coefficient rate is proportional to M^(a/(a+alpha)) eta^(alpha/(a+alpha)); lower Laurent coefficients incur explicit logarithmic losses in the full-jet problem. These are results for the stated observation spaces and priors, not unconditional continuation. The fixed-radius completely monotone family and boundary Sobolev family supply negative controls.

**Evidence:** [stability-theorems.tex](paper/stability-theorems.tex), [angular-logs.tex](paper/angular-logs.tex), [full-jet.tex](paper/full-jet.tex), and their verification reports. Exact principal-series superpositions already exist in [Donnay–Pasterski–Puhm, §§4.1–4.5](https://arxiv.org/pdf/2005.08990). Originality of the present combination relative to inverse-problem literature remains under review.

## 2. Identify one conventional gravitational soft functional

**Attempt and outcome — explicit radiative-data realization, limited scope.** With Bondi news N_AB=partial_u C_AB and memory b_AB=integral N_AB du, two sphere integrations by parts give

Q_S[F] = -<E[F],b>/(16 pi G),

where E_AB[F]=D_A D_B F - q_AB Delta F/2. For a unit normalized l-harmonic Y, ||E[Y]||²=(l-1)l(l+1)(l+2)/2. This supplies an actual dual test and exact charge-error factor for the conventional soft contribution.

The constructed news is b_AB tau/[pi(u²+tau²)]. Its memory is b, Fourier news is b exp(-tau|omega|), and gravitational energy is ||b||²/(64 pi² G tau). The pulse duration can grow with fixed memory while both flux and the displaced Mellin-line norm vanish. For unknown duration under a uniform derivative bound, a pair of opposite harmonic pulses gives the same optimal exponent as the general recovery theorem. A known exact duration gives a linear noise bound with an exact constant. The source of the charge convention is [Pasterski, equation 12 and the adjacent Bondi constraint](https://arxiv.org/html/2310.04932v2); the integrations and family are worked out in [radiative-example.tex](paper/radiative-example.tex).

**Not completed:** Global Einstein evolution, endpoint constraints for independently chosen data, antipodal matching, a dressed quantum state space, and the full soft-plus-hard Ward identity. Smooth finite-flux null radiative data are not a complete scattering solution.

## 3. Construct further logarithmic symmetry charges

**Already established:** The leading logarithmic triangle has important solved cases; it must not be presented as a newly discovered open target. See [massless sQED, §4.2](https://arxiv.org/html/2512.21239v1), and the November massless-gravity derivation [§4.1.1, equations 4.2–4.10](https://arxiv.org/html/2611.15921v1). Further subleading logarithmic charges are a different question.

**Actual attempted algebra — conditional contact terms.** Taking the assumed higher-loop sQED factor in [Choi–Kadhe–Puhm, equation 4.24](https://arxiv.org/html/2601.05951v1), each angular term has the form P_l(bar z)/(z-z_i), with polynomial degree at most l. Away from punctures its (l+1)st bar derivative vanishes. Distributionally, in the convention bar-partial(1/z)=pi delta²(z), our calculation gives

bar-partial^(l+1)[P_l(bar z)/(z-z_i)]
= pi sum_(k=0)^l P_l^(k)(bar z_i) bar-partial^(l-k) delta²(z-z_i).

For l=2 and P_2=(A(bar z-bar z_i)+B)²/2, the contact term is pi[(B²/2)bar-partial² delta + AB bar-partial delta + A² delta]. Pairing against a smooth test gives pi[(B²/2)bar-partial² epsilon - AB bar-partial epsilon + A² epsilon] at z_i. For a monomial in the shifted coordinate, the first k derivatives of `(bar z-bar z_i)^k/(z-z_i)` reduce it to `k!/(z-z_i)`. The next derivative is `k! pi delta²`. This follows from the product rule and the vanishing of positive powers of `bar z-bar z_i` against the delta. The [exact arithmetic check](scripts/verify_ward_controls.py) independently verifies the Leibniz coefficients for degrees 1–8.

**Obstruction:** This follows from an assumed soft-factor ansatz. It does not compute an unverified quantum coefficient, produce the canonical charge, or prove its hard-sector action and matching. It is not a new physical symmetry theorem.

## 4. Test higher-loop and massive extensions

**Scope correction:** First-log massive Ward identities are established, for example [Campiglia–Laddha, 1903.09133](https://arxiv.org/abs/1903.09133). The massive extension of the particular celestial tower and the quantum-beyond-one-loop check requested in [2601.05951v1, §6](https://arxiv.org/html/2601.05951v1) are more specific obligations. The bounded latest-paper search did not certify their general closure.

**Actual consistency replay — a known counterexample.** [Karan et al., October 2026, §4.3, equations 4.21–4.27](https://arxiv.org/html/2610.07134v1) compares a Newtonian three-particle tail coefficient with a previous full multiparticle conjecture. Write K=G^4 m^4/(R alpha^7 v^7). Their coefficients are (125/16)K e_x e_x and (225/16)K e_x e_x. Subtracting and projecting transverse-traceless to detector direction e_z gives

Delta A_TT = -(25/8)K diag(1,-1,0),

The [exact projection check](reports/ward-controls.json) reproduces this nonzero discrepancy from the published coefficients. The source's relativistic replacement in §5 remains conjectural. A failure of a full waveform ansatz does not automatically refute a symmetry-projected formula: terms annihilated by the angular Ward operator must be compared explicitly.

**Next check:** The [ordered-scale criterion](QUANTUM-LOOP-CRITERION.md) proves when integrable corrections between successive virtual scales can contribute at most one logarithm. The physical kernel representation and bound must be derived from the regulated amplitude. The [massive angular inverse](paper/massive-angular.tex) quantifies recovery from the established leading kernel. The higher-log massive Ward construction and the two-loop quantum calculation remain open.

## 5. Account for a finite observation time

**Attempt and outcome — exact idealized window bound.** For the radiative family above, a news measurement on [-U,U] captures b_U=(2/pi)arctan(U/tau)b. With L² time-noise norm epsilon, the soft-charge error is at most

||E[F]||/(16 pi G) [sqrt(2U)epsilon + (1-(2/pi)arctan(U/tau))||b||].

The second term is the exact missing memory. At fixed U it approaches the whole memory as tau grows. [radiative-example.tex](paper/radiative-example.tex)

**Obstruction:** This is a specified news-measurement model, not a derivation of detector strain response, quantum dressing, or a regulator prescription. The newer gravitational paper itself flags a regulated-slice problem before [equation 3.17](https://arxiv.org/html/2611.15921v1).

## 6. Make ordered soft targets well defined

**Attempt and outcome — exact negative diagnostic.** The bounded stripped kernel b(x,y)=x exp(-x-y)/(x+y) has ordered soft values 0 and 1. Its double Mellin transform is Gamma(z+1)Gamma(w)/(z+w). The mixed divisor z+w=0 obstructs a common axes-only Laurent jet. Any single number approximating both ordered values has worst error at least 1/2, even for exact data. [ordered-limits.tex](paper/ordered-limits.tex)

This scalar example illustrates a target-definition issue; it is not a scattering-amplitude counterexample. Actual source ordering issues include [Bhardwaj–Srikant, §4.2](https://arxiv.org/pdf/2403.10443) and [Magnea–Zunino, equations 99–102](https://arxiv.org/html/2512.22104v1). Choosing an order identifies one target; it does not prove OPE associativity or the appropriate physical bracket.

## 7. Replace a finite-log model when resummation changes the target

**Attempt and outcome — exact templates and rejection controls.** The kernels exp(-lambda t) and (1+lambda t)^(-gamma), 0<gamma<1, have transforms 1/(z+lambda) and lambda^(-gamma)z^(gamma-1)exp(z/lambda)Gamma(1-gamma,z/lambda). Finite Taylor/binomial truncations have higher poles at zero; their global error does not decay exponentially in t. [resummation verification](scripts/verify_resummation.py) checks the formulas and negative controls.

These illustrate the analytic structures of [Gu–Li–Wang, §§2–5](https://arxiv.org/pdf/2206.01757). For the known running-log template, [the new recovery proof](paper/branch-recovery.tex) divides a Gaussian Mellin average by the exact template response and uses a cutoff profile as the lower-bound example. It gives matching small-noise error order M^(a/(a+alpha)) eta^(alpha/(a+alpha)) [log(M/eta)]^gamma under a global power bound on the remainder. It also gives a finite-band bound and an obstruction to estimating an unknown exponent uniformly when the amplitude can approach zero. Deriving the template and a uniform remainder for a concrete quantum amplitude remains open.

## 8. Connect transform presentations without claiming a full duality

**Recent restricted physics:** [Caminiti–Pasterski, November 2026, §1, equations 1.1–1.2](https://arxiv.org/html/2611.15854v1) gives a special-state swing-length/mutual-information relation and an inner-horizon reformulation. It does not resolve every flat-space holographic framework.

**Actual derived calculation — classical kinematic bridge.** For compact smooth positive-energy f supported away from zero, let C_±(u)=integral exp(∓i omega u)f(omega)d omega, u>0. An Abel-regulated Fubini argument gives

U_±(s)=Mellin_u C_±(s)=Gamma(s)exp(∓i pi s/2)F(1-s).

On s=1/2+i nu the squared response norms sum to exactly 2 pi. The conjugate-weighted two-channel inverse therefore has sharp linear noise norm 1/sqrt(2 pi), pointwise and on joint L² frequency errors. A single channel has an exponentially bad frequency direction; compact smooth source modulations prove that it has no uniform L² inverse bound. [framework-bridge.tex](paper/framework-bridge.tex) · [verification](reports/framework-bridge.json)

This independently derived application of classical Fourier–Mellin identities is not a new duality. Both time channels are assumed, and Re(q)=1/2 is not automatically the physical Re(Delta)=1 convention. No boundary dynamics, OPE, bulk reconstruction, or entanglement observable follows.

## Result of the bounded campaign

The campaign yields conditional recovery mathematics, one explicit conventional soft functional on radiative data, finite-window and ordered-target diagnostics, a minimax theorem for one prescribed branch template, and a stable two-channel transform bridge. Higher symmetries, massive/higher-quantum-loop celestial identities and a complete holographic theory retain concrete unresolved inputs. Numerical checks replay finite cases; written proofs and source-established physics are identified separately.
