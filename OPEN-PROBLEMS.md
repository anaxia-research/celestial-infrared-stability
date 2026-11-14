# Bounded celestial research map

Prepared 15 November 2026 from the specific primary-source versions linked below. This is not a list of every open problem of Sabrina Pasterski, her colleagues, or celestial holography. Here “colleagues” means the coauthors of the cited papers, without implying that these researchers collaborate on or endorse this repository.

**Status vocabulary:** “established in source” records a result with that paper's scope. “Open in source” records an explicit question or future direction at the cited version, not a certification that no later paper has answered it. “Project question” identifies our proposed next step. None of the categories below is marked fully solved by this project.

## 1. Exact continuation versus stable soft-data recovery

**Established in source:** Donnay–Pasterski–Puhm express continued primaries through regulated principal-series superpositions (§4.2, equations 4.12–4.13); §§4.3–4.5 explain pairings and modes. These are not missing basis states. [2005.08990, §§4.1–4.5](https://arxiv.org/pdf/2005.08990)

**Project question:** Under which physically justified uniform priors can noisy Mellin-line data recover a prescribed soft coefficient? Our manuscript proves conditional scalar/angular coefficient estimates and fixed-log-degree minimax exponents in its stated classes. It does not convert exact contour representability into unconditional stable continuation.

## 2. From a coefficient pairing to a physical soft charge

**Established in sources:** The chapter's “Soft Charges as Currents,” equations 30–32, maps soft coefficients to specified current constructions. Donnay–Pasterski–Puhm §5.3.3 and Appendix C address contact terms and symplectic renormalization. [Pasterski chapter](https://arxiv.org/html/2310.04932v2) · [Donnay–Pasterski–Puhm](https://arxiv.org/pdf/2005.08990)

**Project question:** Embed the present recovery class in a concrete physical scattering model. Specify an angular operator D, normalization, dual spin test D*epsilon, boundary/contact terms, hard contribution, phase space and Ward identity. Bounded pairing estimates alone do not complete this step.

## 3. Further subleading modes, memory and logarithmic symmetries

**Open/future directions in sources:** Donnay–Pasterski–Puhm §1, printed page 4, leaves further subleading non-pure-gauge modes and their memories for future work. Choi–Kadhe–Puhm §1 identifies the asymptotic-symmetry derivation of subleading logarithmic soft theorems as open. [2005.08990, §1](https://arxiv.org/pdf/2005.08990) · [2601.05951v1, §1](https://arxiv.org/html/2601.05951v1)

**Project progress:** Fixed-order log-coefficient reconstruction, assuming the global remainder prior. No new asymptotic symmetry or memory law is established here.

## 4. Higher quantum loops and massive matter

**Open in source:** Choi–Kadhe–Puhm §6 seeks checks of the quantum logarithmic structure beyond one loop and describes massive-matter celestial Ward identities as work in progress. This qualifies the scope of its exponentiation and massless celestial calculations. [2601.05951v1, §6](https://arxiv.org/html/2601.05951v1)

**Project progress:** The reconstruction theorem treats any fixed finite log degree satisfying its assumptions. It neither computes higher loops nor supplies massive-particle Ward identities. Log degree must be obtained from the actual calculation.

## 5. Detector scales, dressing and the soft limit

**Future work in source:** The detector/system-scale interpretation and dressed-state argument in Choi–Kadhe–Puhm §3.4 is revisited in §6, which requests a fuller derivation. [2601.05951v1, §§3.4, 6](https://arxiv.org/html/2601.05951v1)

**Project question:** Derive an observation norm, infrared cutoff and uniform remainder bound for a specified detector/dressing prescription. Our fixed reference-energy convention and finite-band error term are mathematical inputs, not this derivation.

## 6. Loop OPEs, ordering and non-abelian energy dependence

**Established complications, with completion still open in the cited analyses:** Bhardwaj–Yelleshpur Srikant §4.2 studies order-dependent multiple-soft operations. Magnea–Zunino §§5–6 and §7 examine energy-fraction dependence and limitations of associativity in multi-emission limits. These results prevent assuming a universal holomorphic local algebra without specifying limits. [2403.10443, §§4.2, 6](https://arxiv.org/pdf/2403.10443) · [2512.22104v1, §§5–7, especially equations 99–102](https://arxiv.org/html/2512.22104v1)

**Project progress:** None toward an OPE associativity proof. Recovering individual Laurent coefficients does not determine their product algebra or justify interchanging soft limits.

## 7. Fixed-order multipoles versus resummed analytic structure

**Established in source:** Gu–Li–Wang §§2–5 gives massless-EFT resummation examples with branch cuts and special shifted/remaining poles. It leaves broader higher-order analytic structures unresolved. [2206.01757, §5](https://arxiv.org/pdf/2206.01757)

**Project progress:** The shifted-pole and incomplete-gamma scalar models in `paper/angular-logs.tex` and `scripts/verify_resummation.py` show why finite log truncations can violate the global decaying-remainder prior. They are analytic scope controls, not newly solved physical loops. A resummed charge-recovery theorem needs its own target and function class.

## 8. Completing and connecting holographic descriptions

**Research directions in source:** The chapter's “Merging Frameworks” discusses twistor/twisted descriptions, Carrollian and light-ray viewpoints, memory observables and horizon applications. These are broad programme directions, not a single named conjecture with a demonstrated closure. [2310.04932v2, “Merging Frameworks” and its three subsections](https://arxiv.org/html/2310.04932v2)

**Project progress:** No full bulk/boundary duality, black-hole information solution or general celestial theory is constructed. The current mathematical results concern noisy soft-coefficient recovery and convention dependence.

## Next reviewable target

Choose one regulated amplitude and one normalized physical charge. Prove membership in the theorem's angular space and a global remainder bound uniform over a stated kinematic set; derive the corresponding D and admissible dual test; then propagate the data and hard-sector errors. This is a concrete research target. Its hypotheses and Ward identity remain obligations, rather than outcomes implied by this map.

See [RESEARCH-CONTEXT.md](RESEARCH-CONTEXT.md) for inverse-problem antecedents and further loop-source qualifications.
