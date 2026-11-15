# An ordered-scale criterion for the quantum leading logarithm

Choi, Kadhe, and Puhm ask whether the classical structure of their logarithmic soft factors survives beyond one loop in quantum scattering. For massless scalar QED, the first test is the coefficient of $\omega\log^2(\omega/\Lambda)$ in the symmetry-governed single-soft-photon insertion. The projection relevant to this test pairs the angular Ward descendant against a smooth test at a hard puncture, keeping its contact distributions. Equations (3.25), (4.24), and (4.25) specify the predicted term and projection. Equation (4.24) amounts to the half-square of the one-loop anchored operator, with the same anchored tree insertion as a prefactor. Section 6 gives the status of the quantum assumption. [Choi–Kadhe–Puhm](https://arxiv.org/html/2601.05951v1)

A calculation based on ordered virtual energy or momentum scales can use the following sufficient test. Its hypotheses concern the **sum** of the regulated contributions after the prescribed projection and subtractions. The scale reduction and the kernel bound need independent proofs for the QED amplitude.

## Criterion

Let $X$ be a Banach space, $m\ge2$, $T=\log(\Lambda/\omega)>0$, and let
$$
 S_m(T)=\{(t_1,\ldots,t_m):0<t_1<\cdots<t_m<T\}.
$$
If the scales are $\kappa_j=\omega e^{t_j}$, each $dt_j$ is $d\kappa_j/\kappa_j$. Suppose a coefficient extracted from a regulated calculation, after dividing by its specified overall soft-energy power, admits the representation
$$
 F_m(T)=\int_{S_m(T)}K_T(t)\,dt+E_T,\qquad
 \|E_T\|\le C T^{m-1},
$$
with a Bochner-integrable kernel $K_T$. Any multiplicities from the original ordered regions must be included in the definition of $K_T$.

Define $m+1$ gaps by
$$
 d_0=t_1,\quad d_j=t_{j+1}-t_j\ (1\le j<m),\quad d_m=T-t_m.
$$
Suppose there is $A_m\in X$, independent of $T$, and nonnegative
$w_j\in L^1([0,\infty))$, independent of $T$, such that almost everywhere
$$
 \|K_T(t)-A_m\|\le \sum_{j=0}^m w_j(d_j).
$$
Then
$$
 \left\|F_m(T)-\frac{T^m}{m!}A_m\right\|
 \le T^{m-1}\left(C+\frac{1}{(m-1)!}\sum_{j=0}^m\|w_j\|_1\right).
$$
In particular, if the well-separated region yields $A_m=A_1^m$ as a composition on a specified hard insertion, the leading coefficient is $A_1^m/m!$ acting on that insertion. An omitted integration region of size $T^m$ or a connected contribution to $A_m$ would change this conclusion.

## Proof

Changing from the ordered $t_j$ to the independent gaps $d_0,\ldots,d_{m-1}$ has unit Jacobian. The remaining gap is $d_m=T-\sum_{j=0}^{m-1}d_j$. In these coordinates the domain is the simplex $\{d_j\ge0,\ \sum_{j=0}^{m-1}d_j\le T\}$. Any gap may serve as the dependent coordinate, since swapping it with any independent gap has Jacobian of absolute value one.

Hold the gap $d_j=u$ fixed. The remaining independent coordinates form a simplex of dimension $m-1$ and volume $(T-u)^{m-1}/(m-1)!$. Thus
$$
 \int_{S_m(T)}w_j(d_j)\,dt=
 \frac1{(m-1)!}\int_0^T(T-u)^{m-1}w_j(u)\,du
 \le\frac{T^{m-1}}{(m-1)!}\|w_j\|_1 .
$$
The volume of $S_m(T)$ is $T^m/m!$, giving the constant term. Apply the norm inequality for the integral and add the bound on $E_T$.

## Sharpness and a two-loop diagnostic

- For a fixed unit vector $v\in X$, the perturbation $K_T-A_m=e^{-d_j}v$ has integrated norm asymptotic to $T^{m-1}/(m-1)!$. Divide the one-gap formula by $T^{m-1}$ and use dominated convergence. The power $m-1$ in the bound is sharp.
- The profile $(1+d_j)^{-1}v$ tends to zero at large separation, but its integrated norm is $T^{m-1}(\log(1+T)+O(1))/(m-1)!$. Indeed $0\le1-(1-u/T)^{m-1}\le(m-1)u/T$ for $0\le u\le T$. Large-gap pointwise decay alone therefore gives no $O(T^{m-1})$ remainder.
- The constant perturbation $K_T=A_1^2+v$ at order two gives a correction $v/2$ to the coefficient of $T^2$ at that order, with the order-one operator held fixed.

For the massless scalar QED target, one would have to derive the projected and infrared-subtracted two-loop kernel in a fixed regulator and soft-limit prescription. This fixes the hard normalization and includes collinear and overlap subtractions and both orderings of the two virtual scales. The remaining checks are to match its separated-scale limit against the ordered composition in equation (4.24), and to bound or calculate the integral of its connected remainder, including comparable virtual scales and the endpoints. A bound of the form above would show that comparable-scale and endpoint corrections contribute at most one logarithm in the chosen projection. An independent contribution there of order $T^2$ would require a different coefficient test.

Exponentiation of the ordinary virtual infrared divergence is established in related settings. For example, Naculich and Schnitzer derive it for the graviton scattering infrared factor in dimensional regularization using gravitational Wilson lines. The quantity tested here is the subleading soft-energy coefficient of a projected radiative insertion. Its required ordered representation and kernel estimates remain open in this project. [Naculich–Schnitzer](https://arxiv.org/abs/1101.1524)
