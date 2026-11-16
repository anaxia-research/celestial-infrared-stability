import Std

/-!
# CelestialIR: machine-checked finite/algebraic parts of the manuscript

Lean 4 core + `Std` only (no Mathlib), no `sorry`, no new axioms.

1. General pole-jet shear of Theorem `thm:units`, eq. `(eq:jet)`, for jets of any
   length `r` over `Int` (and the verbatim factorial form over `Rat`): group law,
   identity, inverse, compensator, top-coefficient and pole-order invariance,
   integer powers of one unipotent matrix, translation semantics, `Fin r` version.
2. Forward-difference identity of Theorem `thm:vectorlogs` for general `m`, and the
   exactness `Σ (-1)^{m-j} C(m,j) P(L+jd) = m! d^m a_m` on polynomials of degree `≤ m`.
3. Bondi norm `λ(λ-2) = (l-1)l(l+1)(l+2)`, `H^2` domination, spectral bound.
4. Ordered-limit residues as coefficient extraction in the two expansions of
   `1/(w(z+w))`; mixed divisor; discrete simplex/gap counts.
5. Chebyshev: degree and leading coefficient `2^{n-1}` of `T_n` and `2^{2n-1}` of
   `T_n(2u-1)` for all `n`; exact weights of Theorem `thm:cheb` for `m = 1,2,3,4`.

Not formalized: Mellin integrals, meromorphic continuation, Gaussian averages,
Bochner integrals, minimax lower bounds, real-analytic `cos(jπ/m)` for general `m`.
-/

set_option linter.deprecated false

namespace CelestialIR

/-! ## 0. A small finite-sum / binomial library (Std only) -/

def fact : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * fact n

def choose : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

/-- `sumR f n = f 0 + ... + f (n-1)`. -/
def sumR (f : Nat → Int) : Nat → Int
  | 0 => 0
  | n + 1 => sumR f n + f n

@[simp] theorem fact_zero : fact 0 = 1 := rfl
@[simp] theorem fact_succ (n : Nat) : fact (n + 1) = (n + 1) * fact n := rfl
theorem fact_pos (n : Nat) : 0 < fact n := by
  induction n with
  | zero => decide
  | succ n ih => show 0 < (n + 1) * fact n; exact Nat.mul_pos (Nat.succ_pos n) ih

@[simp] theorem choose_zero_right (n : Nat) : choose n 0 = 1 := by cases n <;> rfl
@[simp] theorem choose_zero_succ (k : Nat) : choose 0 (k + 1) = 0 := rfl
theorem choose_succ_succ (n k : Nat) : choose (n + 1) (k + 1) = choose n k + choose n (k + 1) := rfl

theorem choose_eq_zero_of_lt : ∀ {n k : Nat}, n < k → choose n k = 0
  | 0, 0, h => absurd h (by decide)
  | 0, _ + 1, _ => rfl
  | n + 1, 0, h => absurd h (by omega)
  | n + 1, k + 1, h => by
      rw [choose_succ_succ, choose_eq_zero_of_lt (by omega : n < k),
        choose_eq_zero_of_lt (by omega : n < k + 1)]

@[simp] theorem choose_self : ∀ n : Nat, choose n n = 1
  | 0 => rfl
  | n + 1 => by rw [choose_succ_succ, choose_self n, choose_eq_zero_of_lt (by omega)]

/-- Factorial formula, multiplicative form: `C(k+a,k) k! a! = (k+a)!`. -/
theorem choose_mul_fact' : ∀ (k a : Nat), choose (k + a) k * fact k * fact a = fact (k + a)
  | 0, a => by simp
  | k + 1, 0 => by simp
  | k + 1, a + 1 => by
      have h1 := choose_mul_fact' k (a + 1)
      have h2 := choose_mul_fact' (k + 1) a
      have e : k + 1 + (a + 1) = (k + a + 1) + 1 := by omega
      have e1 : k + (a + 1) = k + a + 1 := by omega
      have e2 : k + 1 + a = k + a + 1 := by omega
      rw [e, choose_succ_succ]
      rw [e1] at h1; rw [e2] at h2
      simp only [fact_succ] at h1 h2 ⊢
      generalize choose (k + a + 1) k = A at *
      generalize choose (k + a + 1) (k + 1) = B at *
      generalize fact k = F at *
      generalize fact a = G at *
      generalize fact (k + a) = H at *
      have h2' : B * ((k + 1) * F) * G = A * F * ((a + 1) * G) := by rw [h2, h1]
      rw [← h1]
      grind

theorem choose_mul_fact (n k : Nat) (h : k ≤ n) : choose n k * fact k * fact (n - k) = fact n := by
  have := choose_mul_fact' k (n - k)
  rwa [Nat.add_sub_cancel' h] at this

/-! ### sums -/

@[simp] theorem sumR_zero (f : Nat → Int) : sumR f 0 = 0 := rfl
@[simp] theorem sumR_succ (f : Nat → Int) (n : Nat) : sumR f (n + 1) = sumR f n + f n := rfl

theorem sumR_congr {f g : Nat → Int} : ∀ {n : Nat}, (∀ i, i < n → f i = g i) → sumR f n = sumR g n
  | 0, _ => rfl
  | n + 1, h => by
      simp only [sumR_succ]
      rw [sumR_congr (fun i hi => h i (by omega)), h n (by omega)]

theorem sumR_eq_zero {f : Nat → Int} {n : Nat} (h : ∀ i, i < n → f i = 0) : sumR f n = 0 := by
  have := sumR_congr (g := fun _ => 0) h
  rw [this]; clear this h
  induction n with
  | zero => rfl
  | succ n ih => simp [ih]

theorem sumR_add (f g : Nat → Int) : ∀ n, sumR (fun i => f i + g i) n = sumR f n + sumR g n
  | 0 => rfl
  | n + 1 => by simp only [sumR_succ, sumR_add f g n]; omega

theorem sumR_sub (f g : Nat → Int) : ∀ n, sumR (fun i => f i - g i) n = sumR f n - sumR g n
  | 0 => rfl
  | n + 1 => by simp only [sumR_succ, sumR_sub f g n]; omega

theorem sumR_mul_left (c : Int) (f : Nat → Int) : ∀ n, sumR (fun i => c * f i) n = c * sumR f n
  | 0 => by simp
  | n + 1 => by simp only [sumR_succ, sumR_mul_left c f n, Int.mul_add]

theorem sumR_mul_right (c : Int) (f : Nat → Int) (n : Nat) :
    sumR (fun i => f i * c) n = sumR f n * c := by
  rw [Int.mul_comm, ← sumR_mul_left]; exact sumR_congr (fun i _ => Int.mul_comm _ _)

/-- Split off the first term. -/
theorem sumR_succ' (f : Nat → Int) : ∀ n, sumR f (n + 1) = f 0 + sumR (fun i => f (i + 1)) n
  | 0 => by simp
  | n + 1 => by
      rw [sumR_succ, sumR_succ' f n, sumR_succ]; omega

theorem sumR_split (f : Nat → Int) (a : Nat) : ∀ b, sumR f (a + b) = sumR f a + sumR (fun i => f (a + i)) b
  | 0 => by simp
  | b + 1 => by
      rw [← Nat.add_assoc, sumR_succ, sumR_split f a b, sumR_succ]; omega

/-- Zeros beyond `n` do not change the sum. -/
theorem sumR_extend {f : Nat → Int} {n m : Nat} (hnm : n ≤ m) (h : ∀ i, n ≤ i → i < m → f i = 0) :
    sumR f m = sumR f n := by
  obtain ⟨b, rfl⟩ := Nat.exists_eq_add_of_le hnm
  rw [sumR_split, sumR_eq_zero (fun i hi => h (n + i) (by omega) (by omega))]; omega

/-- A sum with a single possibly nonzero term. -/
theorem sumR_single {f : Nat → Int} {n d : Nat} (hd : d < n) (h : ∀ i, i < n → i ≠ d → f i = 0) :
    sumR f n = f d := by
  rw [sumR_extend (n := d + 1) (by omega) (fun i h1 h2 => h i h2 (by omega)), sumR_succ,
    sumR_eq_zero (fun i hi => h i (by omega) (by omega))]; omega

/-- Interchange of a rectangular double sum. -/
theorem sumR_comm (f : Nat → Nat → Int) : ∀ n m : Nat,
    sumR (fun i => sumR (fun j => f i j) m) n = sumR (fun j => sumR (fun i => f i j) n) m
  | 0, m => by simp; exact (sumR_eq_zero (fun _ _ => rfl)).symm
  | n + 1, m => by
      rw [sumR_succ, sumR_comm f n m, ← sumR_add]; rfl

/-- **Binomial theorem** over `Int`. -/
theorem binomial (x y : Int) : ∀ n : Nat,
    (x + y) ^ n = sumR (fun i => (choose n i : Int) * x ^ i * y ^ (n - i)) (n + 1)
  | 0 => by simp
  | n + 1 => by
      have ih := binomial x y n
      have hA : x * sumR (fun i => (choose n i : Int) * x ^ i * y ^ (n - i)) (n + 1)
          = sumR (fun i => (choose n i : Int) * x ^ (i + 1) * y ^ (n - i)) (n + 1) := by
        rw [← sumR_mul_left]; apply sumR_congr; intro i _; rw [Int.pow_succ]; grind
      have hB : y * sumR (fun i => (choose n i : Int) * x ^ i * y ^ (n - i)) (n + 1)
          = y ^ (n + 1) + sumR (fun i => (choose n (i + 1) : Int) * x ^ (i + 1) * y ^ (n - i)) (n + 1) := by
        have hc : ∀ i, i < n → y * ((choose n (i + 1) : Int) * x ^ (i + 1) * y ^ (n - (i + 1)))
            = (choose n (i + 1) : Int) * x ^ (i + 1) * y ^ (n - i) := by
          intro i hi
          have : n - i = (n - (i + 1)) + 1 := by omega
          rw [this, Int.pow_succ]; grind
        rw [← sumR_mul_left, sumR_succ', sumR_succ (fun i => (choose n (i + 1) : Int) * _ * _),
          choose_eq_zero_of_lt (by omega : n < n + 1), sumR_congr hc]
        simp only [choose_zero_right, Int.natCast_one, Int.pow_zero, Nat.sub_zero,
          Int.pow_succ, Int.natCast_zero]
        grind
      have hR : sumR (fun i => (choose (n + 1) i : Int) * x ^ i * y ^ (n + 1 - i)) (n + 1 + 1)
          = y ^ (n + 1) + (sumR (fun i => (choose n i : Int) * x ^ (i + 1) * y ^ (n - i)) (n + 1)
            + sumR (fun i => (choose n (i + 1) : Int) * x ^ (i + 1) * y ^ (n - i)) (n + 1)) := by
        rw [sumR_succ', ← sumR_add]
        congr 1
        · simp
        · apply sumR_congr; intro i _
          have : n + 1 - (i + 1) = n - i := by omega
          rw [this, choose_succ_succ]; push_cast; grind
      rw [hR, Int.pow_succ, Int.mul_add, ih, Int.mul_comm _ x, Int.mul_comm _ y, hA, hB]
      grind

/-- `C(k+n,k+l) C(k+l,k) = C(k+n,k) C(n,l)` for `l ≤ n`. -/
theorem choose_mul_choose (k n l : Nat) (hl : l ≤ n) :
    choose (k + n) (k + l) * choose (k + l) k = choose (k + n) k * choose n l := by
  have c1 := choose_mul_fact' k l
  have c2 := choose_mul_fact (k + n) (k + l) (by omega)
  have c3 := choose_mul_fact n l hl
  have c4 := choose_mul_fact' k n
  have e : k + n - (k + l) = n - l := by omega
  rw [e] at c2
  have hpos : 0 < fact k * fact l * fact (n - l) :=
    Nat.mul_pos (Nat.mul_pos (fact_pos k) (fact_pos l)) (fact_pos (n - l))
  apply Nat.eq_of_mul_eq_mul_right hpos
  generalize fact (n - l) = G at *
  calc choose (k + n) (k + l) * choose (k + l) k * (fact k * fact l * G)
      = choose (k + n) (k + l) * (choose (k + l) k * fact k * fact l) * G := by ac_rfl
    _ = fact (k + n) := by rw [c1, c2]
    _ = choose (k + n) k * fact k * (choose n l * fact l * G) := by rw [c3, c4]
    _ = choose (k + n) k * choose n l * (fact k * fact l * G) := by ac_rfl

theorem choose_mul_fact_int (j k : Nat) (h : k ≤ j) :
    (choose j k : Int) * fact k * fact (j - k) = fact j := by
  have := choose_mul_fact j k h
  exact_mod_cast this

/-! ## 1. The numerical-energy pole-jet shear for jets of arbitrary length `r`

Paper, Theorem `thm:units`, eq. `(eq:jet)`:
`a'_{-k} = e^{-(w+s_0)σ} Σ_{j=k}^r (-σ)^{j-k}/(j-k)! a_{-j}`.

Exact integer model.  Write the principal part through the log-polynomial
coefficients `u` (paper `(eq:vectorclass)`): `β_k = a_{-(k+1)} = k! u_k`,
`0 ≤ k < r`.  In the divided basis `u` the shear by an integer `σ` is the
binomial (Taylor) matrix `C(j,k) (-σ)^{j-k}`, i.e. `exp(-σ D)` for the
nilpotent derivation `D`; it has integer entries for every integer `σ`,
and in the `β` basis it is exactly `exp(-σ N)` with the nilpotent shift
`(Nβ)_k = β_{k+1}` (theorem `jet_weight` / `jet_formula`).  The scalar
prefactor `e^{-(w+s_0)σ}` commutes with everything (`shear_smul`) and equals
`1` in the neutral case `w+s_0 = 0`. -/

/-- The shear on a length-`r` jet (`u : Nat → Int`, only `u 0, …, u (r-1)` are read). -/
def shear (r : Nat) (σ : Int) (u : Nat → Int) (k : Nat) : Int :=
  sumR (fun j => (choose j k : Int) * (-σ) ^ (j - k) * u j) r

/-- A jet of length `r`: coefficients vanish from index `r` on. -/
def IsJet (r : Nat) (u : Nat → Int) : Prop := ∀ j, r ≤ j → u j = 0

theorem shear_vanish (r : Nat) (σ : Int) (u : Nat → Int) {k : Nat} (hk : r ≤ k) :
    shear r σ u k = 0 :=
  sumR_eq_zero (fun j hj => by rw [choose_eq_zero_of_lt (by omega)]; simp)

theorem shear_isJet (r : Nat) (σ : Int) (u : Nat → Int) : IsJet r (shear r σ u) :=
  fun _ hk => shear_vanish r σ u hk

/-- The convolution kernel identity behind the group law. -/
theorem shear_kernel (a b : Int) (j k r : Nat) (hj : j < r) :
    sumR (fun i => (choose i k : Int) * (choose j i : Int) * a ^ (i - k) * b ^ (j - i)) r
      = (choose j k : Int) * (a + b) ^ (j - k) := by
  by_cases hkj : k ≤ j
  · obtain ⟨n, rfl⟩ := Nat.exists_eq_add_of_le hkj
    rw [sumR_extend (n := k + n + 1) (by omega) (fun i h1 _ => by
          rw [choose_eq_zero_of_lt (n := k + n) (k := i) (by omega)]; simp),
      show k + n + 1 = k + (n + 1) from rfl, sumR_split,
      sumR_eq_zero (fun i hi => by rw [choose_eq_zero_of_lt (n := i) (k := k) hi]; simp),
      Int.zero_add, binomial, ← sumR_mul_left]
    have e0 : k + n - k = n := by omega
    rw [e0]
    apply sumR_congr; intro l hl
    have e1 : k + l - k = l := by omega
    have e2 : k + n - (k + l) = n - l := by omega
    rw [e1, e2]
    have hc : (choose (k + n) (k + l) : Int) * choose (k + l) k = choose (k + n) k * choose n l := by
      exact_mod_cast choose_mul_choose k n l (by omega)
    grind
  · rw [choose_eq_zero_of_lt (by omega : j < k)]
    simp only [Int.natCast_zero, Int.zero_mul]
    apply sumR_eq_zero; intro i _
    by_cases hik : i < k
    · rw [choose_eq_zero_of_lt hik]; simp
    · rw [choose_eq_zero_of_lt (n := j) (k := i) (by omega)]; simp

/-- **Group law** of the jet shear, for all integers `σ, τ` and every length `r`. -/
theorem shear_add (r : Nat) (σ τ : Int) (u : Nat → Int) :
    shear r τ (shear r σ u) = shear r (σ + τ) u := by
  funext k
  unfold shear
  have h1 : ∀ i, i < r → (choose i k : Int) * (-τ) ^ (i - k) *
        sumR (fun j => (choose j i : Int) * (-σ) ^ (j - i) * u j) r
      = sumR (fun j => (choose i k : Int) * (choose j i : Int) * (-τ) ^ (i - k) * (-σ) ^ (j - i) * u j) r := by
    intro i _; rw [← sumR_mul_left]; apply sumR_congr; intro j _; grind
  rw [sumR_congr h1, sumR_comm]
  apply sumR_congr; intro j hj
  rw [sumR_mul_right (u j) (fun i => (choose i k : Int) * (choose j i : Int) * (-τ) ^ (i - k) * (-σ) ^ (j - i)),
    shear_kernel _ _ j k r hj, show -(σ + τ) = -τ + -σ by omega]

theorem shear_zero_apply (r : Nat) (u : Nat → Int) {k : Nat} (hk : k < r) : shear r 0 u k = u k := by
  unfold shear
  rw [sumR_single hk]
  · simp
  · intro j _ hjk
    by_cases h : j < k
    · rw [choose_eq_zero_of_lt h]; simp
    · obtain ⟨e, he⟩ : ∃ e, j - k = e + 1 := ⟨j - k - 1, by omega⟩
      rw [he, Int.pow_succ]; simp

/-- The identity element acts trivially on jets. -/
theorem shear_zero (r : Nat) (u : Nat → Int) (hu : IsJet r u) : shear r 0 u = u := by
  funext k
  by_cases hk : k < r
  · exact shear_zero_apply r u hk
  · rw [shear_vanish r 0 u (by omega), hu k (by omega)]

/-- **Inverse**: the shear by `-σ` undoes the shear by `σ`. -/
theorem shear_inverse (r : Nat) (σ : Int) (u : Nat → Int) (hu : IsJet r u) :
    shear r (-σ) (shear r σ u) = u := by
  rw [shear_add, show σ + -σ = 0 by omega, shear_zero r u hu]

theorem shear_inverse' (r : Nat) (σ : Int) (u : Nat → Int) (hu : IsJet r u) :
    shear r σ (shear r (-σ) u) = u := by
  rw [shear_add, show -σ + σ = 0 by omega, shear_zero r u hu]

/-- **Reference-scale compensator**: a log-reference shifted by `-s` cancels the
whole nilpotent jet action, `S_{ℓ-s} ∘ S_s = S_ℓ`. -/
theorem shear_compensator (r : Nat) (ell s : Int) (u : Nat → Int) :
    shear r (ell - s) (shear r s u) = shear r ell u := by
  rw [shear_add, show s + (ell - s) = ell by omega]

/-- **Invariance of the top coefficient** (highest pole coefficient `a_{-r}`). -/
theorem shear_top (n : Nat) (σ : Int) (u : Nat → Int) : shear (n + 1) σ u n = u n := by
  unfold shear
  rw [sumR_single (Nat.lt_succ_self n)]
  · simp
  · intro j hj hjn; rw [choose_eq_zero_of_lt (by omega)]; simp

/-- **Pole order is preserved**: if the jet has actual top index `d`, the sheared
jet has the same top index and the same top coefficient. -/
theorem shear_preserves_order (r : Nat) (σ : Int) (u : Nat → Int) (d : Nat) (hd : d < r)
    (hu : ∀ j, d < j → u j = 0) :
    (∀ k, d < k → shear r σ u k = 0) ∧ shear r σ u d = u d := by
  constructor
  · intro k hk
    apply sumR_eq_zero; intro j _
    by_cases hj : j ≤ d
    · rw [choose_eq_zero_of_lt (by omega : j < k)]; simp
    · rw [hu j (by omega)]; simp
  · unfold shear
    rw [sumR_single hd]
    · simp
    · intro j _ hjd
      by_cases hj : j < d
      · rw [choose_eq_zero_of_lt hj]; simp
      · rw [hu j (by omega)]; simp

/-- Linearity (the prefactor `e^{-(w+s_0)σ}` is a scalar and commutes with the shear). -/
theorem shear_smul (r : Nat) (σ c : Int) (u : Nat → Int) :
    shear r σ (fun j => c * u j) = fun k => c * shear r σ u k := by
  funext k; unfold shear; rw [← sumR_mul_left]; apply sumR_congr; intro j _; grind

theorem shear_add_jets (r : Nat) (σ : Int) (u v : Nat → Int) :
    shear r σ (fun j => u j + v j) = fun k => shear r σ u k + shear r σ v k := by
  funext k; unfold shear; rw [← sumR_add]; apply sumR_congr; intro j _; grind

/-- Integer powers of the unit shear `S_1` (and of its inverse `S_{-1}`). -/
def shearPow (r : Nat) (σ : Int) : Nat → (Nat → Int) → (Nat → Int)
  | 0, u => u
  | n + 1, u => shear r σ (shearPow r σ n u)

theorem shearPow_eq (r : Nat) (σ : Int) (u : Nat → Int) (hu : IsJet r u) :
    ∀ n : Nat, shearPow r σ n u = shear r (n * σ) u
  | 0 => by simp [shearPow, shear_zero r u hu]
  | n + 1 => by
      rw [shearPow, shearPow_eq r σ u hu n, shear_add]
      congr 1; push_cast; grind

/-- `S_1^n = S_n` and `S_{-1}^n = S_{-n}`: the shear is the integer-power family of one
unipotent matrix. -/
theorem shear_integer_powers (r : Nat) (u : Nat → Int) (hu : IsJet r u) (n : Nat) :
    shearPow r 1 n u = shear r n u ∧ shearPow r (-1) n u = shear r (-(n : Int)) u := by
  constructor
  · rw [shearPow_eq r 1 u hu]; simp
  · rw [shearPow_eq r (-1) u hu]; congr 1; omega

/-- Evaluation semantics: the shear is the translation `t ↦ t - σ` of the logarithmic
polynomial `Σ_k u_k t^k` (`t = log(1/x)`; a unit change `x ↦ e^{σ}x`). -/
def evalJ (r : Nat) (u : Nat → Int) (t : Int) : Int := sumR (fun k => u k * t ^ k) r

theorem shear_eval (r : Nat) (σ : Int) (u : Nat → Int) (t : Int) :
    evalJ r (shear r σ u) t = evalJ r u (t - σ) := by
  unfold evalJ shear
  rw [sumR_congr (fun k _ => (sumR_mul_right (t ^ k) _ r).symm), sumR_comm]
  apply sumR_congr; intro j hj
  rw [sumR_extend (n := j + 1) (by omega) (fun k h1 _ => by
        rw [choose_eq_zero_of_lt (n := j) (k := k) (by omega)]; simp),
    show t - σ = t + -σ by omega, binomial t (-σ) j, ← sumR_mul_left]
  apply sumR_congr; intro i _; grind

/-- Bridge to eq. `(eq:jet)`: in the `β` basis (`β_j = j! u_j`) the coefficient of `β_j`
in `β'_k` is `(-σ)^{j-k}/(j-k)!`; integer form, `(j-k)! · W_{kj} = (-σ)^{j-k} · j!`. -/
def jetW (σ : Int) (k j : Nat) : Int := (fact k : Int) * choose j k * (-σ) ^ (j - k)

theorem jet_formula (r : Nat) (σ : Int) (u : Nat → Int) (k : Nat) :
    (fact k : Int) * shear r σ u k = sumR (fun j => jetW σ k j * u j) r := by
  unfold shear jetW; rw [← sumR_mul_left]; apply sumR_congr; intro j _; grind

theorem jet_weight (σ : Int) (k j : Nat) (h : k ≤ j) :
    (fact (j - k) : Int) * jetW σ k j = (-σ) ^ (j - k) * fact j := by
  unfold jetW
  have := choose_mul_fact_int j k h
  rw [← this]; grind

theorem jet_weight_zero (σ : Int) (k j : Nat) (h : j < k) : jetW σ k j = 0 := by
  unfold jetW; rw [choose_eq_zero_of_lt h]; simp

/-- Eq. `(eq:jet)` verbatim over `Rat`: `β'_k = Σ_{j=k}^{r-1} (-σ)^{j-k}/(j-k)! β_j`. -/
def sumQ (f : Nat → Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => sumQ f n + f n

theorem sumQ_cast (f : Nat → Int) : ∀ n, ((sumR f n : Int) : Rat) = sumQ (fun i => (f i : Rat)) n
  | 0 => by simp [sumR, sumQ]
  | n + 1 => by simp [sumR, sumQ, sumQ_cast f n]

theorem jet_formula_rat (r : Nat) (σ : Int) (u : Nat → Int) (k : Nat) :
    (((fact k : Int) * shear r σ u k : Int) : Rat)
      = sumQ (fun j => if k ≤ j then ((-σ : Int) : Rat) ^ (j - k) / (fact (j - k) : Rat)
                        * (((fact j : Int) * u j : Int) : Rat) else 0) r := by
  rw [jet_formula, sumQ_cast]
  congr 1; funext j
  split
  · rename_i h
    have hw := jet_weight σ k j h
    have hF : (fact (j - k) : Rat) ≠ 0 := by
      have := fact_pos (j - k); intro h0; simp at h0; omega
    have hw' : ((fact (j - k) : Nat) : Rat) * (jetW σ k j : Rat) = ((-σ : Int) : Rat) ^ (j - k) * ((fact j : Nat) : Rat) := by
      exact_mod_cast hw
    generalize ((-σ : Int) : Rat) ^ (j - k) = S at hw' ⊢
    push_cast
    generalize ((fact (j - k) : Nat) : Rat) = F at hF hw' ⊢
    grind
  · rename_i h; rw [jet_weight_zero σ k j (by omega)]; simp

/-! ### The same action on `Fin r → Int` -/

def extJ {r : Nat} (v : Fin r → Int) : Nat → Int := fun j => if h : j < r then v ⟨j, h⟩ else 0

def shearFin {r : Nat} (σ : Int) (v : Fin r → Int) : Fin r → Int := fun i => shear r σ (extJ v) i

theorem extJ_shearFin {r : Nat} (σ : Int) (v : Fin r → Int) : extJ (shearFin σ v) = shear r σ (extJ v) := by
  funext j
  by_cases h : j < r
  · simp [extJ, shearFin, h]
  · simp only [extJ, h, dite_false]; rw [shear_vanish r σ _ (by omega)]

theorem extJ_isJet {r : Nat} (v : Fin r → Int) : IsJet r (extJ v) := by
  intro j hj; simp [extJ]; omega

theorem shearFin_add {r : Nat} (σ τ : Int) (v : Fin r → Int) :
    shearFin τ (shearFin σ v) = shearFin (σ + τ) v := by
  funext i; simp only [shearFin]; rw [extJ_shearFin, shear_add]

theorem shearFin_zero {r : Nat} (v : Fin r → Int) : shearFin 0 v = v := by
  funext i; simp only [shearFin]; rw [shear_zero_apply r _ i.isLt]; simp [extJ]

theorem shearFin_inverse {r : Nat} (σ : Int) (v : Fin r → Int) : shearFin (-σ) (shearFin σ v) = v := by
  rw [shearFin_add, show σ + -σ = 0 by omega, shearFin_zero]

theorem shearFin_compensator {r : Nat} (ell s : Int) (v : Fin r → Int) :
    shearFin (ell - s) (shearFin s v) = shearFin ell v := by
  rw [shearFin_add, show s + (ell - s) = ell by omega]

theorem shearFin_top {n : Nat} (σ : Int) (v : Fin (n + 1) → Int) :
    shearFin σ v (Fin.last n) = v (Fin.last n) := by
  simp only [shearFin, Fin.val_last]; rw [shear_top]; simp [extJ]; rfl

/-- Length 2 recovers `ScaleJets.shift`: `(a_{-1}, a_{-2}) ↦ (a_{-1} - σ a_{-2}, a_{-2})`. -/
theorem shearFin_two (σ : Int) (v : Fin 2 → Int) :
    shearFin σ v 0 = v 0 - σ * v 1 ∧ shearFin σ v 1 = v 1 := by
  simp [shearFin, shear, extJ, sumR, choose]; grind

/-- The residue is not invariant (negative control, length 3). -/
theorem shear_residue_moves : shear 3 5 (fun j => if j = 0 then 3 else if j = 1 then 2 else 1) 0 ≠ 3 := by
  decide

/-! ## 2. The forward-difference identity of Theorem `thm:vectorlogs`

`Σ_{j=0}^m (-1)^{m-j} C(m,j) j^k = 0` for `k < m` and `= m!` for `k = m`,
for every `m`; and its consequence for the estimator `(eq:logestimator)`:
on nodes `L + j d` the `m`-th difference of any polynomial of degree `≤ m`
is `m! d^m a_m`. -/

/-- `m`-th forward difference at `0`: `Σ_{j=0}^m (-1)^{m-j} C(m,j) f(j)`. -/
def fdiff (m : Nat) (f : Nat → Int) : Int :=
  sumR (fun j => (-1 : Int) ^ (m - j) * choose m j * f j) (m + 1)

/-- `F m k = Σ_j (-1)^{m-j} C(m,j) j^k`. -/
def fdPow (m k : Nat) : Int := fdiff m (fun j => (j : Int) ^ k)

theorem absorb (m j : Nat) : (j + 1) * choose (m + 1) (j + 1) = (m + 1) * choose m j := by
  by_cases hj : j ≤ m
  · have c1 := choose_mul_fact (m + 1) (j + 1) (by omega)
    have c2 := choose_mul_fact m j hj
    have e : m + 1 - (j + 1) = m - j := by omega
    rw [e] at c1
    simp only [fact_succ] at c1
    have hpos : 0 < fact j * fact (m - j) := Nat.mul_pos (fact_pos j) (fact_pos (m - j))
    apply Nat.eq_of_mul_eq_mul_right hpos
    generalize fact j = A at *
    generalize fact (m - j) = B at *
    calc (j + 1) * choose (m + 1) (j + 1) * (A * B)
        = choose (m + 1) (j + 1) * ((j + 1) * A) * B := by ac_rfl
      _ = (m + 1) * fact m := c1
      _ = (m + 1) * (choose m j * A * B) := by rw [c2]
      _ = (m + 1) * choose m j * (A * B) := by ac_rfl
  · rw [choose_eq_zero_of_lt (by omega : m + 1 < j + 1), choose_eq_zero_of_lt (by omega : m < j)]
    simp

theorem fdiff_mul_left (m : Nat) (c : Int) (f : Nat → Int) :
    fdiff m (fun j => c * f j) = c * fdiff m f := by
  unfold fdiff; rw [← sumR_mul_left]; apply sumR_congr; intro j _; grind

theorem fdiff_sum (m n : Nat) (g : Nat → Nat → Int) :
    fdiff m (fun j => sumR (fun i => g i j) n) = sumR (fun i => fdiff m (g i)) n := by
  unfold fdiff
  rw [sumR_congr (fun j _ => (sumR_mul_left _ (fun i => g i j) n).symm), sumR_comm]

theorem fdPow_succ_succ (m k : Nat) :
    fdPow (m + 1) (k + 1) = (m + 1) * sumR (fun l => (choose k l : Int) * fdPow m l) (k + 1) := by
  unfold fdPow fdiff
  dsimp only
  rw [sumR_succ']
  have h0 : ((0 : Nat) : Int) ^ (k + 1) = 0 := by simp [Int.pow_succ]
  rw [h0, Int.mul_zero, Int.zero_add]
  -- termwise rewrite of the left side
  have hterm : ∀ i, i < m + 1 →
      (-1 : Int) ^ (m + 1 - (i + 1)) * (choose (m + 1) (i + 1) : Int) * ((i + 1 : Nat) : Int) ^ (k + 1)
        = (m + 1) * sumR (fun l => (-1 : Int) ^ (m - i) * choose m i * ((choose k l : Int) * (i : Int) ^ l)) (k + 1) := by
    intro i _
    rw [sumR_mul_left, show m + 1 - (i + 1) = m - i by omega]
    have hb := binomial (i : Int) 1 k
    simp only [Int.one_pow, Int.mul_one] at hb
    rw [← hb]
    have ha : ((i : Int) + 1) * choose (m + 1) (i + 1) = (m + 1) * choose m i := by
      exact_mod_cast absorb m i
    push_cast
    rw [Int.pow_succ]
    generalize ((i : Int) + 1) ^ k = P
    grind
  rw [sumR_congr hterm, sumR_mul_left, sumR_comm]
  congr 1
  apply sumR_congr; intro l _
  rw [← sumR_mul_left]; apply sumR_congr; intro i _; grind

theorem fdPow_zero_left (k : Nat) : fdPow 0 k = if k = 0 then 1 else 0 := by
  cases k <;> simp [fdPow, fdiff, Int.pow_succ]

/-- `(1 + (-1))^n = 0` for `n ≥ 1`, i.e. `fdPow (m+1) 0 = 0`. -/
theorem fdPow_succ_zero (m : Nat) : fdPow (m + 1) 0 = 0 := by
  have hb := binomial 1 (-1) (m + 1)
  rw [show (1 : Int) + -1 = 0 by rfl, Int.pow_succ, Int.mul_zero] at hb
  unfold fdPow fdiff
  rw [hb]
  apply sumR_congr; intro j _
  simp only [Int.pow_zero, Int.mul_one, Int.one_pow]
  grind

/-- **Forward-difference identity**, general `m`. -/
theorem fdPow_eq : ∀ m k : Nat, k ≤ m → fdPow m k = if k = m then (fact m : Int) else 0
  | 0, k, hk => by
      have : k = 0 := by omega
      subst this; simp [fdPow_zero_left]
  | m + 1, 0, _ => by rw [fdPow_succ_zero]; simp
  | m + 1, k + 1, hk => by
      rw [fdPow_succ_succ]
      have ih := fdPow_eq m
      by_cases hkm : k < m
      · rw [sumR_eq_zero (fun l hl => by rw [ih l (by omega), if_neg (by omega)]; simp)]
        rw [if_neg (by omega)]; simp
      · have hk' : k = m := by omega
        subst hk'
        rw [sumR_single (d := k) (by omega) (fun l hl hlk => by rw [ih l (by omega), if_neg hlk]; simp),
          ih k (Nat.le_refl k), if_pos rfl, if_pos rfl, choose_self]
        simp [fact_succ]

theorem forward_difference_lt (m k : Nat) (h : k < m) :
    sumR (fun j => (-1 : Int) ^ (m - j) * choose m j * (j : Int) ^ k) (m + 1) = 0 := by
  have := fdPow_eq m k (by omega); rw [if_neg (by omega)] at this; exact this

theorem forward_difference_eq (m : Nat) :
    sumR (fun j => (-1 : Int) ^ (m - j) * choose m j * (j : Int) ^ m) (m + 1) = fact m := by
  have := fdPow_eq m m (Nat.le_refl m); rw [if_pos rfl] at this; exact this

/-- On the geometric scale grid `t_j = L + j d` (log variable) the `m`-th difference of
`t^k`, `k ≤ m`, is `0` for `k < m` and `m! d^m` for `k = m`. -/
theorem fdiff_affine_pow (m k : Nat) (hk : k ≤ m) (L d : Int) :
    fdiff m (fun j => (L + (j : Int) * d) ^ k) = if k = m then (fact m : Int) * d ^ m else 0 := by
  have hb : ∀ j : Nat, (L + (j : Int) * d) ^ k
      = sumR (fun i => ((choose k i : Int) * d ^ i * L ^ (k - i)) * (j : Int) ^ i) (k + 1) := by
    intro j
    rw [Int.add_comm, binomial]
    apply sumR_congr; intro i _
    rw [Int.mul_pow]; grind
  rw [show (fun j : Nat => (L + (j : Int) * d) ^ k) = (fun j : Nat =>
      sumR (fun i => ((choose k i : Int) * d ^ i * L ^ (k - i)) * (j : Int) ^ i) (k + 1)) from funext hb,
    fdiff_sum]
  simp only [fdiff_mul_left]
  by_cases hkm : k < m
  · rw [if_neg (by omega)]
    apply sumR_eq_zero; intro i hi
    show _ * fdPow m i = 0
    rw [fdPow_eq m i (by omega), if_neg (by omega)]; simp
  · have : k = m := by omega
    subst this
    rw [if_pos rfl, sumR_single (d := k) (by omega)]
    · show _ * fdPow k k = _
      rw [fdPow_eq k k (Nat.le_refl k), if_pos rfl]; simp; grind
    · intro i hi hik
      show _ * fdPow k i = 0
      rw [fdPow_eq k i (by omega), if_neg hik]; simp

/-- **Exactness of the estimator `(eq:logestimator)` on polynomials**: for
`P(t) = Σ_{k≤m} a_k t^k`, `Σ_j (-1)^{m-j} C(m,j) P(L + j d) = m! d^m a_m`. -/
theorem fdiff_polynomial (m : Nat) (a : Nat → Int) (L d : Int) :
    fdiff m (fun j => sumR (fun k => a k * (L + (j : Int) * d) ^ k) (m + 1))
      = (fact m : Int) * d ^ m * a m := by
  rw [fdiff_sum]
  simp only [fdiff_mul_left]
  rw [sumR_single (d := m) (by omega)]
  · rw [fdiff_affine_pow m m (Nat.le_refl m), if_pos rfl]; grind
  · intro k hk hkm
    rw [fdiff_affine_pow m k (by omega), if_neg hkm]; simp

/-! ## 3. Bondi harmonic norm and `H^2` domination (Prop. `prop:harmoniccharge`) -/

/-- `λ_l(λ_l - 2) = (l-1) l (l+1) (l+2)` with `λ_l = l(l+1)`. -/
theorem bondi_norm (l : Int) : l * (l + 1) * (l * (l + 1) - 2) = (l - 1) * l * (l + 1) * (l + 2) := by
  grind

/-- `λ(λ-2) ≤ (1+λ)^2` for `λ ≥ 0`, with exact gap `1 + 4λ`. -/
theorem h2_domination (lam : Int) (h : 0 ≤ lam) : lam * (lam - 2) ≤ (1 + lam) ^ 2 := by
  grind

theorem h2_gap (lam : Int) : (1 + lam) ^ 2 - lam * (lam - 2) = 1 + 4 * lam := by grind

/-- For integer `l ≥ 0` the norm `H_l` is nonnegative and vanishes exactly for `l = 0, 1`. -/
theorem bondi_norm_nonneg (l : Nat) : 0 ≤ ((l : Int) - 1) * l * (l + 1) * (l + 2) := by
  rcases l with _ | _ | l
  · decide
  · decide
  · have h1 : (0 : Int) ≤ ((l + 1 + 1 : Nat) : Int) - 1 := by omega
    have h2 : (0 : Int) ≤ ((l + 1 + 1 : Nat) : Int) := by omega
    exact Int.mul_nonneg (Int.mul_nonneg (Int.mul_nonneg h1 h2) (by omega)) (by omega)

theorem bondi_norm_eq_zero_iff (l : Nat) : ((l : Int) - 1) * l * (l + 1) * (l + 2) = 0 ↔ l ≤ 1 := by
  rcases l with _ | _ | l
  · decide
  · decide
  · constructor
    · intro h
      have h1 : (0 : Int) < ((l + 1 + 1 : Nat) : Int) - 1 := by omega
      have h2 : (0 : Int) < ((l + 1 + 1 : Nat) : Int) := by omega
      have := Int.mul_pos (Int.mul_pos (Int.mul_pos h1 h2)
        (by omega : (0 : Int) < ((l + 1 + 1 : Nat) : Int) + 1))
        (by omega : (0 : Int) < ((l + 1 + 1 : Nat) : Int) + 2)
      omega
    · intro h; omega

theorem sumR_le {f g : Nat → Int} : ∀ {n : Nat}, (∀ i, i < n → f i ≤ g i) → sumR f n ≤ sumR g n
  | 0, _ => Int.le_refl 0
  | n + 1, h => by
      simp only [sumR_succ]
      have := sumR_le (n := n) (fun i hi => h i (by omega))
      have := h n (by omega)
      omega

/-- Spectral `H^2` bound: `Σ λ_l(λ_l-2)|F_l|^2 ≤ Σ (1+λ_l)^2 |F_l|^2` for any finite
spectral data with `λ_l ≥ 0` and weights `c_l = |F_l|^2 ≥ 0`. -/
theorem h2_spectral_bound (n : Nat) (lam c : Nat → Int) (hl : ∀ i, 0 ≤ lam i) (hc : ∀ i, 0 ≤ c i) :
    sumR (fun i => lam i * (lam i - 2) * c i) n ≤ sumR (fun i => (1 + lam i) ^ 2 * c i) n := by
  apply sumR_le; intro i _
  have h := Int.mul_nonneg (show 0 ≤ 1 + 4 * lam i by have := hl i; omega) (hc i)
  have e : (1 + lam i) ^ 2 * c i - lam i * (lam i - 2) * c i = (1 + 4 * lam i) * c i := by grind
  omega

/-- The same bound for the harmonic eigenvalues `λ_l = l(l+1)`, `l ≥ 0`. -/
theorem h2_harmonic_bound (n : Nat) (c : Nat → Int) (hc : ∀ i, 0 ≤ c i) :
    sumR (fun l => ((l : Int) - 1) * l * (l + 1) * (l + 2) * c l) n
      ≤ sumR (fun l => (1 + (l : Int) * (l + 1)) ^ 2 * c l) n := by
  have := h2_spectral_bound n (fun l => (l : Int) * (l + 1)) c
    (fun l => Int.mul_nonneg (by omega) (by omega)) hc
  rw [sumR_congr (g := fun l => (l : Int) * (l + 1) * ((l : Int) * (l + 1) - 2) * c l)
    (fun l _ => by rw [← bondi_norm])]
  exact this

/-! ## 4a. Ordered-limit residues (paper `ordered-limits.tex`), algebraic form

`G(z,w) = Γ(z+1)Γ(w)/(z+w) = p(z,w) / (w (z+w))` with `p = Γ(z+1)Γ(w+1)`
holomorphic, `p(0,0) = 1`.  Taking `Res_{z=0}` first (at fixed `w ≠ 0`, so
`|z| < |w|`) uses the expansion `e_A` of `1/(w(z+w))` in powers of `z/w`;
taking `Res_{w=0}` first (`|w| < |z|`) uses the expansion `e_B` in powers of
`w/z`.  Both are genuine two-sided inverses of the Laurent polynomial
`zw + w^2` (coefficientwise), they have different supports, and the iterated
residue — the coefficient of `z^{-1} w^{-1}` of `p · e` — is `0` for `e_A` and
`p(0,0)` for `e_B`, for every polynomial truncation `p`. -/

/-- `1/(w(z+w)) = Σ_{n≥0} (-1)^n z^n w^{-n-2}` (valid for `|z| < |w|`); coefficient of `z^i w^j`. -/
def eA (i j : Int) : Int := if 0 ≤ i ∧ j = -i - 2 then (-1) ^ i.toNat else 0

/-- `1/(w(z+w)) = Σ_{n≥0} (-1)^n w^{n-1} z^{-n-1}` (valid for `|w| < |z|`). -/
def eB (i j : Int) : Int := if 0 ≤ j + 1 ∧ i = -j - 2 then (-1) ^ (j + 1).toNat else 0

/-- Coefficient of `z^i w^j` in `(z w + w^2) · e`. -/
def mulDen (e : Int → Int → Int) (i j : Int) : Int := e (i - 1) (j - 1) + e i (j - 2)

def delta (i j : Int) : Int := if i = 0 ∧ j = 0 then 1 else 0

theorem eA_inverse (i j : Int) : mulDen eA i j = delta i j := by
  unfold mulDen eA delta
  by_cases h1 : 1 ≤ i ∧ j = -i
  · rw [if_pos (by omega), if_pos (by omega), if_neg (by omega)]
    have : i.toNat = (i - 1).toNat + 1 := by omega
    rw [this, Int.pow_succ]; omega
  · by_cases h2 : i = 0 ∧ j = 0
    · rw [if_neg (by omega), if_pos (by omega), if_pos h2]; obtain ⟨rfl, rfl⟩ := h2; rfl
    · rw [if_neg (by omega), if_neg (by omega), if_neg h2]; rfl

theorem eB_inverse (i j : Int) : mulDen eB i j = delta i j := by
  unfold mulDen eB delta
  by_cases h1 : 1 ≤ j ∧ i = -j
  · rw [if_pos (by omega), if_pos (by omega), if_neg (by omega)]
    have : (j - 1 + 1).toNat = (j - 2 + 1).toNat + 1 := by omega
    rw [this, Int.pow_succ]; omega
  · by_cases h2 : i = 0 ∧ j = 0
    · rw [if_pos (by omega), if_neg (by omega), if_pos h2]; obtain ⟨rfl, rfl⟩ := h2; rfl
    · rw [if_neg (by omega), if_neg (by omega), if_neg h2]; rfl

/-- The two expansions live in different completions: `e_A` has only
nonnegative powers of `z`, `e_B` has `w`-exponents `≥ -1`. -/
theorem eA_support (i j : Int) (h : eA i j ≠ 0) : 0 ≤ i := by
  unfold eA at h; split at h <;> omega

theorem eB_support (i j : Int) (h : eB i j ≠ 0) : -1 ≤ j := by
  unfold eB at h; split at h <;> omega

/-- Iterated residue: coefficient of `z^{-1} w^{-1}` in `p · e` for a polynomial `p`
with coefficients `p a b`, `a, b < N`. -/
def iterRes (p : Nat → Nat → Int) (N : Nat) (e : Int → Int → Int) : Int :=
  sumR (fun a => sumR (fun b => p a b * e (-1 - a) (-1 - b)) N) N

/-- `Res_{w=0} Res_{z=0} G = 0`. -/
theorem res_w_after_z (p : Nat → Nat → Int) (N : Nat) : iterRes p N eA = 0 := by
  apply sumR_eq_zero; intro a _
  apply sumR_eq_zero; intro b _
  unfold eA; rw [if_neg (by omega)]; simp

/-- `Res_{z=0} Res_{w=0} G = p(0,0)`. -/
theorem res_z_after_w (p : Nat → Nat → Int) (N : Nat) (hN : 0 < N) : iterRes p N eB = p 0 0 := by
  unfold iterRes
  rw [sumR_single (d := 0) hN]
  · rw [sumR_single (d := 0) hN]
    · unfold eB; rw [if_pos (by omega)]; simp
    · intro b _ hb; unfold eB; rw [if_neg (by omega)]; simp
  · intro a _ ha
    apply sumR_eq_zero; intro b _
    unfold eB; rw [if_neg (by omega)]; simp

/-- For the numerator `Γ(z+1)Γ(w+1)` (constant term `1`) the ordered residues are `0` and `1`. -/
theorem ordered_residues (p : Nat → Nat → Int) (N : Nat) (hN : 0 < N) (h00 : p 0 0 = 1) :
    iterRes p N eA = 0 ∧ iterRes p N eB = 1 :=
  ⟨res_w_after_z p N, by rw [res_z_after_w p N hN, h00]⟩

theorem neg_one_pow_ne_zero : ∀ b : Nat, (-1 : Int) ^ b ≠ 0
  | 0 => by decide
  | b + 1 => by
      rw [Int.pow_succ]; intro h
      have := neg_one_pow_ne_zero b
      rcases Int.mul_eq_zero.mp h with h | h
      · exact this h
      · omega

/-- The mixed divisor `z + w = 0`: `z + w` divides no monomial `z^a w^b`, so `G` has no
local Laurent expansion whose only polar divisors are the coordinate axes. -/
theorem mixed_divisor (a b : Nat) (P : Int → Int → Int) :
    ¬ (∀ z w : Int, z ^ a * w ^ b = (z + w) * P z w) := by
  intro h
  have := h 1 (-1)
  rw [Int.one_pow, Int.one_mul, show (1 : Int) + -1 = 0 by rfl, Int.zero_mul] at this
  exact neg_one_pow_ne_zero b this

/-- A single number cannot approximate both ordered targets `0` and `1` better than `1/2`. -/
theorem ordered_targets_gap (q : Rat) : 1 / 2 ≤ q ∨ 1 / 2 ≤ 1 - q := by grind

theorem ordered_targets_gap_int (a b : Int) : b ≤ 2 * a ∨ b ≤ 2 * (b - a) := by omega

/-! ## 4b. Simplex volume / gap count, discrete form (paper `ordered-scales.tex`)

Ordered scale tuples are parametrized by their nonnegative gaps.  The discrete
analogue of `vol S_m(T) = T^m/m!` is the number `C(T+m, m)` of gap tuples
`(d_0,…,d_{m-1}) ∈ ℕ^m` with `Σ d_i ≤ T`; fixing one gap `d_0 = u` leaves a
simplex of dimension `m-1` with `C(T-u+m-1, m-1)` points, the analogue of
`(T-u)^{m-1}/(m-1)!` in `(eq:orderedscale-gap)`. -/

/-- All gap tuples of length `m` with sum `≤ T`. -/
def tuples : Nat → Nat → List (List Nat)
  | 0, _ => [[]]
  | m + 1, T => (List.range (T + 1)).flatMap (fun u => (tuples m (T - u)).map (u :: ·))

theorem mem_tuples : ∀ (m T : Nat) (l : List Nat), l ∈ tuples m T ↔ l.length = m ∧ l.sum ≤ T
  | 0, T, l => by
      cases l with
      | nil => simp [tuples]
      | cons a l => simp [tuples]
  | m + 1, T, l => by
      simp only [tuples, List.mem_flatMap, List.mem_range, List.mem_map]
      constructor
      · rintro ⟨u, hu, l', hl', rfl⟩
        have := (mem_tuples m (T - u) l').1 hl'
        simp; omega
      · rintro ⟨hlen, hsum⟩
        cases l with
        | nil => simp at hlen
        | cons u l' =>
          simp at hlen hsum
          exact ⟨u, by omega, l', (mem_tuples m (T - u) l').2 ⟨hlen, by omega⟩, rfl⟩

def sumN (f : Nat → Nat) : Nat → Nat
  | 0 => 0
  | n + 1 => sumN f n + f n

theorem sumN_succ' (f : Nat → Nat) : ∀ n, sumN f (n + 1) = f 0 + sumN (fun i => f (i + 1)) n
  | 0 => by simp [sumN]
  | n + 1 => by
      rw [sumN, sumN_succ' f n, sumN]; omega

theorem sumN_congr {f g : Nat → Nat} : ∀ {n : Nat}, (∀ i, i < n → f i = g i) → sumN f n = sumN g n
  | 0, _ => rfl
  | n + 1, h => by
      simp only [sumN]
      rw [sumN_congr (fun i hi => h i (by omega)), h n (by omega)]

theorem sum_map_range (f : Nat → Nat) : ∀ n, ((List.range n).map f).sum = sumN f n
  | 0 => by simp [sumN]
  | n + 1 => by
      rw [List.range_succ, List.map_append, List.sum_append, sum_map_range f n]; simp [sumN]

/-- Hockey stick: `Σ_{u=0}^T C(T-u+m, m) = C(T+m+1, m+1)`. -/
theorem hockey_stick (m : Nat) : ∀ T, sumN (fun u => choose (T - u + m) m) (T + 1) = choose (T + (m + 1)) (m + 1)
  | 0 => by simp [sumN]
  | T + 1 => by
      rw [sumN_succ']
      rw [sumN_congr (g := fun u => choose (T - u + m) m) (fun i _ => by congr 1; omega),
        hockey_stick m T]
      rw [show T + 1 + (m + 1) = (T + m + 1) + 1 by omega, choose_succ_succ,
        show T + 1 - 0 + m = T + m + 1 by omega, show T + (m + 1) = T + m + 1 by omega]

/-- **Lattice-point count of the gap simplex**: `#{d ∈ ℕ^m : Σ d ≤ T} = C(T+m, m)`. -/
theorem length_tuples : ∀ m T, (tuples m T).length = choose (T + m) m
  | 0, T => by simp [tuples]
  | m + 1, T => by
      simp only [tuples, List.length_flatMap, List.length_map]
      rw [sum_map_range, sumN_congr (g := fun u => choose (T - u + m) m)
        (fun u _ => length_tuples m (T - u)), hockey_stick]

theorem sum_flatMap {α : Type} (L : List α) (g : α → List Int) :
    (L.flatMap g).sum = (L.map (fun a => (g a).sum)).sum := by
  induction L with
  | nil => simp
  | cons a L ih => simp [List.flatMap_cons, List.sum_append, ih]

theorem sum_map_const {α : Type} (L : List α) (c : Int) : (L.map (fun _ => c)).sum = L.length * c := by
  induction L with
  | nil => simp
  | cons a L ih => simp [ih]; grind

theorem sum_map_range_int (f : Nat → Int) : ∀ n, ((List.range n).map f).sum = sumR f n
  | 0 => by simp
  | n + 1 => by
      rw [List.range_succ, List.map_append, List.sum_append, sum_map_range_int f n]; simp

/-- **Gap identity** (discrete `(eq:orderedscale-gap)`): summing a weight of the first gap
over the simplex equals `Σ_u w(u) · #(simplex of dimension m with sum ≤ T-u)`. -/
theorem gap_identity (m T : Nat) (w : Nat → Int) :
    ((tuples (m + 1) T).map (fun l => w (l.headD 0))).sum
      = sumR (fun u => w u * (choose (T - u + m) m : Int)) (T + 1) := by
  simp only [tuples, List.map_flatMap, List.map_map]
  rw [sum_flatMap, ← sum_map_range_int]
  congr 1
  apply List.map_congr_left
  intro u _
  have : ((fun l => w (l.headD 0)) ∘ fun x => u :: x) = fun _ => w u := by funext l; rfl
  rw [this, sum_map_const, length_tuples]; grind

/-! ## 5. Chebyshev polynomials and the weight identity of Theorem `thm:cheb` -/

/-- Coefficients of `T_n` (index = power of `x`): `T_0 = 1`, `T_1 = x`,
`T_{n+2} = 2x T_{n+1} - T_n`. -/
def cheb : Nat → Nat → Int
  | 0 => fun k => if k = 0 then 1 else 0
  | 1 => fun k => if k = 1 then 1 else 0
  | n + 2 => fun k => 2 * (if k = 0 then 0 else cheb (n + 1) (k - 1)) - cheb n k

theorem cheb_deg_lead_aux : ∀ n : Nat,
    ((∀ k, n < k → cheb n k = 0) ∧ (1 ≤ n → cheb n n = 2 ^ (n - 1))) ∧
    ((∀ k, n + 1 < k → cheb (n + 1) k = 0) ∧ cheb (n + 1) (n + 1) = 2 ^ n)
  | 0 => by
      refine ⟨⟨fun k hk => ?_, fun h => absurd h (by decide)⟩, fun k hk => ?_, rfl⟩
      · simp [cheb]; omega
      · simp [cheb]; omega
  | n + 1 => by
      obtain ⟨_, ⟨h1, h2⟩⟩ := cheb_deg_lead_aux n
      obtain ⟨⟨h3, _⟩, _⟩ := cheb_deg_lead_aux n
      refine ⟨⟨h1, fun _ => h2⟩, fun k hk => ?_, ?_⟩
      · show 2 * (if k = 0 then 0 else cheb (n + 1) (k - 1)) - cheb n k = 0
        rw [if_neg (by omega), h1 (k - 1) (by omega), h3 k (by omega)]; rfl
      · show 2 * (if n + 2 = 0 then 0 else cheb (n + 1) (n + 2 - 1)) - cheb n (n + 2) = 2 ^ (n + 1)
        rw [if_neg (by omega), show n + 2 - 1 = n + 1 by omega, h2, h3 (n + 2) (by omega), Int.pow_succ]
        grind

/-- `deg T_n ≤ n`. -/
theorem cheb_degree (n k : Nat) (h : n < k) : cheb n k = 0 := (cheb_deg_lead_aux n).1.1 k h

/-- **Leading coefficient of `T_n` is `2^{n-1}`** for every `n ≥ 1`. -/
theorem cheb_leading (n : Nat) (h : 1 ≤ n) : cheb n n = 2 ^ (n - 1) := (cheb_deg_lead_aux n).1.2 h

/-- Coefficients of `Q_n(u) = T_n(2u - 1)` (Chebyshev polynomial on `[0,1]`):
`Q_0 = 1`, `Q_1 = 2u - 1`, `Q_{n+2} = 2(2u-1) Q_{n+1} - Q_n`. -/
def chebS : Nat → Nat → Int
  | 0 => fun k => if k = 0 then 1 else 0
  | 1 => fun k => if k = 0 then -1 else if k = 1 then 2 else 0
  | n + 2 => fun k => 4 * (if k = 0 then 0 else chebS (n + 1) (k - 1)) - 2 * chebS (n + 1) k - chebS n k

theorem chebS_deg_lead_aux : ∀ n : Nat,
    ((∀ k, n < k → chebS n k = 0) ∧ (1 ≤ n → chebS n n = 2 ^ (2 * n - 1))) ∧
    ((∀ k, n + 1 < k → chebS (n + 1) k = 0) ∧ chebS (n + 1) (n + 1) = 2 ^ (2 * n + 1))
  | 0 => by
      refine ⟨⟨fun k hk => ?_, fun h => absurd h (by decide)⟩, fun k hk => ?_, rfl⟩
      · simp [chebS]; omega
      · simp [chebS]; omega
  | n + 1 => by
      obtain ⟨_, ⟨h1, h2⟩⟩ := chebS_deg_lead_aux n
      obtain ⟨⟨h3, _⟩, _⟩ := chebS_deg_lead_aux n
      refine ⟨⟨h1, fun _ => by rw [h2]; congr 1⟩, fun k hk => ?_, ?_⟩
      · show 4 * (if k = 0 then 0 else chebS (n + 1) (k - 1)) - 2 * chebS (n + 1) k - chebS n k = 0
        rw [if_neg (by omega), h1 (k - 1) (by omega), h1 k (by omega), h3 k (by omega)]; rfl
      · show 4 * (if n + 2 = 0 then 0 else chebS (n + 1) (n + 2 - 1)) - 2 * chebS (n + 1) (n + 2)
          - chebS n (n + 2) = 2 ^ (2 * (n + 1) + 1)
        rw [if_neg (by omega), show n + 2 - 1 = n + 1 by omega, h2, h1 (n + 2) (by omega),
          h3 (n + 2) (by omega), show 2 * (n + 1) + 1 = (2 * n + 1) + 1 + 1 by omega,
          Int.pow_succ, Int.pow_succ]
        grind

/-- **Leading coefficient `2^{2m-1}` of `Q_m(u) = T_m(2u-1)`** for every `m ≥ 1`; with
`u = (t-L)/T` this is the constant `K_m = 2^{2m-1}/T^m` of `(eq:kappa)`. -/
theorem chebS_leading (m : Nat) (h : 1 ≤ m) : chebS m m = 2 ^ (2 * m - 1) := (chebS_deg_lead_aux m).1.2 h

theorem chebS_degree (m k : Nat) (h : m < k) : chebS m k = 0 := (chebS_deg_lead_aux m).1.1 k h

/-- Values: `T_n` by the recurrence (any integer argument). -/
def chebVal : Nat → Int → Int
  | 0, _ => 1
  | 1, x => x
  | n + 2, x => 2 * x * chebVal (n + 1) x - chebVal n x

def evalP (p : Nat → Int) (N : Nat) (x : Int) : Int := sumR (fun k => p k * x ^ k) N

theorem evalP_shift (p : Nat → Int) (N : Nat) (x : Int) :
    evalP (fun k => if k = 0 then 0 else p (k - 1)) (N + 1) x = x * evalP p N x := by
  unfold evalP
  rw [sumR_succ']
  have h0 : (fun k => if k = 0 then (0 : Int) else p (k - 1)) 0 * x ^ 0 = 0 := by simp
  rw [h0, Int.zero_add, ← sumR_mul_left]
  apply sumR_congr; intro i _
  simp only [Nat.add_one_ne_zero, if_false, Nat.add_sub_cancel, Int.pow_succ]; grind

theorem evalP_extend (p : Nat → Int) {N M : Nat} (h : N ≤ M) (hp : ∀ k, N ≤ k → p k = 0) (x : Int) :
    evalP p M x = evalP p N x := by
  unfold evalP; exact sumR_extend h (fun k h1 _ => by rw [hp k h1]; simp)

/-- The coefficient sequence `cheb n` really is `T_n`: its evaluation obeys the recurrence. -/
theorem cheb_eval : ∀ (n : Nat) (x : Int), evalP (cheb n) (n + 1) x = chebVal n x
  | 0, x => by simp [evalP, cheb, chebVal]
  | 1, x => by simp [evalP, cheb, chebVal, sumR, Int.pow_one]
  | n + 2, x => by
      have e : evalP (cheb (n + 2)) (n + 3) x
          = 2 * evalP (fun k => if k = 0 then 0 else cheb (n + 1) (k - 1)) (n + 2 + 1) x
            - evalP (cheb n) (n + 3) x := by
        unfold evalP; rw [← sumR_mul_left, ← sumR_sub]
        apply sumR_congr; intro k _; simp only [cheb]; grind
      rw [e, evalP_shift, cheb_eval (n + 1) x,
        evalP_extend (cheb n) (N := n + 1) (by omega) (fun k hk => cheb_degree n k (by omega)),
        cheb_eval n x]
      simp only [chebVal]; grind

/-- `chebS n` is `T_n(2u - 1)`. -/
theorem chebS_eval : ∀ (n : Nat) (u : Int), evalP (chebS n) (n + 1) u = chebVal n (2 * u - 1)
  | 0, u => by simp [evalP, chebS, chebVal]
  | 1, u => by simp [evalP, chebS, chebVal, sumR, Int.pow_one]; omega
  | n + 2, u => by
      have e : evalP (chebS (n + 2)) (n + 3) u
          = 4 * evalP (fun k => if k = 0 then 0 else chebS (n + 1) (k - 1)) (n + 2 + 1) u
            - 2 * evalP (chebS (n + 1)) (n + 3) u - evalP (chebS n) (n + 3) u := by
        unfold evalP; rw [← sumR_mul_left, ← sumR_mul_left, ← sumR_sub, ← sumR_sub]
        apply sumR_congr; intro k _; simp only [chebS]; grind
      rw [e, evalP_shift, chebS_eval (n + 1) u,
        evalP_extend (chebS (n + 1)) (N := n + 2) (by omega) (fun k hk => chebS_degree (n + 1) k (by omega)),
        chebS_eval (n + 1) u,
        evalP_extend (chebS n) (N := n + 1) (by omega) (fun k hk => chebS_degree n k (by omega)),
        chebS_eval n u]
      simp only [chebVal]; grind

/-- Endpoint nodes `x_0 = 1`, `x_m = -1`: `T_m(1) = 1`, `T_m(-1) = (-1)^m`, all `m`. -/
theorem chebVal_one_neg_one : ∀ n : Nat, chebVal n 1 = 1 ∧ chebVal n (-1) = (-1) ^ n
  | 0 => by simp [chebVal]
  | 1 => by simp [chebVal]
  | n + 2 => by
      obtain ⟨a1, b1⟩ := chebVal_one_neg_one (n + 1)
      obtain ⟨a0, b0⟩ := chebVal_one_neg_one n
      simp only [chebVal, a1, b1, a0, b0, Int.pow_succ]; constructor <;> grind

/-! ### Theorem `thm:cheb`, exact weights for `m = 1, 2, 3` (rational nodes) and `m = 4`

On `[-1,1]` (`L = -1`, `T = 2`): nodes `x_j = cos(jπ/m)`, Lagrange leading weights
`w_j = 1/Π_{k≠j}(x_j - x_k)`.  We check: (i) the `w_j` are those weights,
(ii) `Σ_j w_j P(x_j) = a_m` for every polynomial of degree `≤ m`,
(iii) `T_m(x_j) = (-1)^j`, (iv) `sign w_j = (-1)^j` and `Σ_j |w_j| = 2^{m-1}`;
then (v) on a general interval `[L, L+T]`, nodes `t_j = L + T(1+x_j)/2` and weights
`W_j = w_j (2/T)^m` reproduce `a_m` and `Σ_j (-1)^j W_j = 2^{2m-1}/T^m = K_m`. -/

theorem cheb_m1 :
    ((1 / 2 : Rat) * (1 - (-1)) = 1 ∧ (-1 / 2 : Rat) * ((-1) - 1) = 1) ∧
    (∀ c0 c1 : Rat, (1 / 2) * (c0 + c1 * 1) + (-1 / 2) * (c0 + c1 * (-1)) = c1) ∧
    ((1 : Rat) = 1 ∧ (-1 : Rat) = -1) ∧
    ((0 : Rat) < 1 / 2 ∧ (0 : Rat) < -(-1 / 2) ∧ (1 / 2 : Rat) - (-1 / 2) = 2 ^ 0) ∧
    (∀ L T c0 c1 : Rat, T ≠ 0 →
      (1 / 2) * (2 / T) * (c0 + c1 * (L + T)) + (-1 / 2) * (2 / T) * (c0 + c1 * L) = c1 ∧
      (1 / 2) * (2 / T) - (-1 / 2) * (2 / T) = 2 ^ 1 / T ^ 1) := by
  refine ⟨⟨by grind, by grind⟩, fun c0 c1 => by grind, ⟨rfl, rfl⟩, ⟨by grind, by grind, by grind⟩,
    fun L T c0 c1 hT => ⟨by grind, by grind⟩⟩

theorem cheb_m2 :
    ((1 / 2 : Rat) * ((1 - 0) * (1 - (-1))) = 1 ∧ (-1 : Rat) * ((0 - 1) * (0 - (-1))) = 1 ∧
      (1 / 2 : Rat) * (((-1) - 1) * ((-1) - 0)) = 1) ∧
    (∀ c0 c1 c2 : Rat, (1 / 2) * (c0 + c1 * 1 + c2 * 1 ^ 2) + (-1) * (c0 + c1 * 0 + c2 * 0 ^ 2)
        + (1 / 2) * (c0 + c1 * (-1) + c2 * (-1) ^ 2) = c2) ∧
    ((2 * 1 ^ 2 - 1 : Rat) = 1 ∧ (2 * 0 ^ 2 - 1 : Rat) = -1 ∧ (2 * (-1) ^ 2 - 1 : Rat) = 1) ∧
    ((0 : Rat) < 1 / 2 ∧ (0 : Rat) < -(-1) ∧ (0 : Rat) < 1 / 2 ∧
      (1 / 2 : Rat) - (-1) + 1 / 2 = 2 ^ 1) ∧
    (∀ L T c0 c1 c2 : Rat, T ≠ 0 →
      (1 / 2) * (2 / T) ^ 2 * (c0 + c1 * (L + T) + c2 * (L + T) ^ 2)
        + (-1) * (2 / T) ^ 2 * (c0 + c1 * (L + T / 2) + c2 * (L + T / 2) ^ 2)
        + (1 / 2) * (2 / T) ^ 2 * (c0 + c1 * L + c2 * L ^ 2) = c2 ∧
      (1 / 2) * (2 / T) ^ 2 - (-1) * (2 / T) ^ 2 + (1 / 2) * (2 / T) ^ 2 = 2 ^ 3 / T ^ 2) := by
  refine ⟨⟨by grind, by grind, by grind⟩, fun c0 c1 c2 => by grind, ⟨by grind, by grind, by grind⟩,
    ⟨by grind, by grind, by grind, by grind⟩, fun L T c0 c1 c2 hT => ⟨by grind, by grind⟩⟩

/-- `m = 3`: nodes `1, 1/2, -1/2, -1`, weights `2/3, -4/3, 4/3, -2/3`, `T_3 = 4x^3 - 3x`. -/
theorem cheb_m3 :
    ((2 / 3 : Rat) * ((1 - 1 / 2) * (1 + 1 / 2) * (1 + 1)) = 1 ∧
      (-4 / 3 : Rat) * ((1 / 2 - 1) * (1 / 2 + 1 / 2) * (1 / 2 + 1)) = 1 ∧
      (4 / 3 : Rat) * ((-1 / 2 - 1) * (-1 / 2 - 1 / 2) * (-1 / 2 + 1)) = 1 ∧
      (-2 / 3 : Rat) * ((-1 - 1) * (-1 - 1 / 2) * (-1 + 1 / 2)) = 1) ∧
    (∀ c0 c1 c2 c3 : Rat,
      (2 / 3) * (c0 + c1 + c2 + c3) + (-4 / 3) * (c0 + c1 * (1 / 2) + c2 * (1 / 2) ^ 2 + c3 * (1 / 2) ^ 3)
        + (4 / 3) * (c0 + c1 * (-1 / 2) + c2 * (-1 / 2) ^ 2 + c3 * (-1 / 2) ^ 3)
        + (-2 / 3) * (c0 - c1 + c2 - c3) = c3) ∧
    ((4 * 1 ^ 3 - 3 * 1 : Rat) = 1 ∧ (4 * (1 / 2) ^ 3 - 3 * (1 / 2) : Rat) = -1 ∧
      (4 * (-1 / 2) ^ 3 - 3 * (-1 / 2) : Rat) = 1 ∧ (4 * (-1) ^ 3 - 3 * (-1) : Rat) = -1) ∧
    ((0 : Rat) < 2 / 3 ∧ (0 : Rat) < 4 / 3 ∧
      (2 / 3 : Rat) - (-4 / 3) + 4 / 3 - (-2 / 3) = 2 ^ 2) ∧
    (∀ L T c0 c1 c2 c3 : Rat, T ≠ 0 →
      (2 / 3) * (2 / T) ^ 3 * (c0 + c1 * (L + T) + c2 * (L + T) ^ 2 + c3 * (L + T) ^ 3)
        + (-4 / 3) * (2 / T) ^ 3 * (c0 + c1 * (L + 3 * T / 4) + c2 * (L + 3 * T / 4) ^ 2 + c3 * (L + 3 * T / 4) ^ 3)
        + (4 / 3) * (2 / T) ^ 3 * (c0 + c1 * (L + T / 4) + c2 * (L + T / 4) ^ 2 + c3 * (L + T / 4) ^ 3)
        + (-2 / 3) * (2 / T) ^ 3 * (c0 + c1 * L + c2 * L ^ 2 + c3 * L ^ 3) = c3 ∧
      (2 / 3) * (2 / T) ^ 3 - (-4 / 3) * (2 / T) ^ 3 + (4 / 3) * (2 / T) ^ 3 - (-2 / 3) * (2 / T) ^ 3
        = 2 ^ 5 / T ^ 3) := by
  refine ⟨⟨by grind, by grind, by grind, by grind⟩, fun c0 c1 c2 c3 => by grind,
    ⟨by grind, by grind, by grind, by grind⟩, ⟨by grind, by grind, by grind⟩,
    fun L T c0 c1 c2 c3 hT => ⟨by grind, by grind⟩⟩

/-- `m = 4`: nodes `1, s, 0, -s, -1` with `s = cos(π/4)`, characterized only by
`2 s^2 = 1` (so the statement holds in `ℝ`, or any commutative ring with such `s`);
weights `1, -2, 2, -2, 1`, `T_4 = 8x^4 - 8x^2 + 1`, `Σ|w_j| = 8 = 2^3`. -/
theorem cheb_m4 {α : Type} [Lean.Grind.CommRing α] (s : α) (hs : 2 * s * s = 1) :
    ((1 : α) * ((1 - s) * (1 - 0) * (1 + s) * (1 + 1)) = 1 ∧
      (-2 : α) * ((s - 1) * (s - 0) * (s + s) * (s + 1)) = 1 ∧
      (2 : α) * ((0 - 1) * (0 - s) * (0 + s) * (0 + 1)) = 1 ∧
      (-2 : α) * ((-s - 1) * (-s - s) * (-s - 0) * (-s + 1)) = 1 ∧
      (1 : α) * ((-1 - 1) * (-1 - s) * (-1 - 0) * (-1 + s)) = 1) ∧
    (∀ c0 c1 c2 c3 c4 : α,
      1 * (c0 + c1 + c2 + c3 + c4) + (-2) * (c0 + c1 * s + c2 * s ^ 2 + c3 * s ^ 3 + c4 * s ^ 4)
        + 2 * c0 + (-2) * (c0 - c1 * s + c2 * s ^ 2 - c3 * s ^ 3 + c4 * s ^ 4)
        + 1 * (c0 - c1 + c2 - c3 + c4) = c4) ∧
    ((8 * 1 ^ 4 - 8 * 1 ^ 2 + 1 : α) = 1 ∧ 8 * s ^ 4 - 8 * s ^ 2 + 1 = -1 ∧
      (8 * 0 ^ 4 - 8 * 0 ^ 2 + 1 : α) = 1 ∧ 8 * (-s) ^ 4 - 8 * (-s) ^ 2 + 1 = -1 ∧
      (8 * (-1) ^ 4 - 8 * (-1) ^ 2 + 1 : α) = 1) ∧
    ((1 : Int) - (-2) + 2 - (-2) + 1 = 2 ^ 3) := by
  refine ⟨⟨by grind, by grind, by grind, by grind, by grind⟩, fun c0 c1 c2 c3 c4 => by grind,
    ⟨by grind, by grind, by grind, by grind, by grind⟩, by decide⟩

/-- `m = 4` on a general interval, over any field with `2 s^2 = 1`. -/
theorem cheb_m4_interval {α : Type} [Lean.Grind.Field α] (s : α) (hs : 2 * s * s = 1)
    (L T c0 c1 c2 c3 c4 : α) (hT : T ≠ 0) :
    let P : α → α := fun t => c0 + c1 * t + c2 * t ^ 2 + c3 * t ^ 3 + c4 * t ^ 4
    let t : α → α := fun x => L + T * (1 + x) / 2
    1 * (2 / T) ^ 4 * P (t 1) + (-2) * (2 / T) ^ 4 * P (t s) + 2 * (2 / T) ^ 4 * P (t 0)
      + (-2) * (2 / T) ^ 4 * P (t (-s)) + 1 * (2 / T) ^ 4 * P (t (-1)) = c4 ∧
    1 * (2 / T) ^ 4 - (-2) * (2 / T) ^ 4 + 2 * (2 / T) ^ 4 - (-2) * (2 / T) ^ 4 + 1 * (2 / T) ^ 4
      = 2 ^ 7 / T ^ 4 := by
  intro P t
  constructor
  · simp only [P, t]; grind
  · grind

end CelestialIR

#print axioms CelestialIR.binomial
#print axioms CelestialIR.choose_mul_choose
#print axioms CelestialIR.shear_kernel
#print axioms CelestialIR.shear_add
#print axioms CelestialIR.shear_zero
#print axioms CelestialIR.shear_inverse
#print axioms CelestialIR.shear_inverse'
#print axioms CelestialIR.shear_compensator
#print axioms CelestialIR.shear_top
#print axioms CelestialIR.shear_preserves_order
#print axioms CelestialIR.shear_smul
#print axioms CelestialIR.shear_integer_powers
#print axioms CelestialIR.shear_eval
#print axioms CelestialIR.jet_formula
#print axioms CelestialIR.jet_weight
#print axioms CelestialIR.jet_formula_rat
#print axioms CelestialIR.shearFin_add
#print axioms CelestialIR.shearFin_zero
#print axioms CelestialIR.shearFin_inverse
#print axioms CelestialIR.shearFin_compensator
#print axioms CelestialIR.shearFin_top
#print axioms CelestialIR.shearFin_two
#print axioms CelestialIR.shear_residue_moves
#print axioms CelestialIR.fdPow_eq
#print axioms CelestialIR.forward_difference_lt
#print axioms CelestialIR.forward_difference_eq
#print axioms CelestialIR.fdiff_affine_pow
#print axioms CelestialIR.fdiff_polynomial
#print axioms CelestialIR.bondi_norm
#print axioms CelestialIR.h2_domination
#print axioms CelestialIR.h2_gap
#print axioms CelestialIR.bondi_norm_nonneg
#print axioms CelestialIR.bondi_norm_eq_zero_iff
#print axioms CelestialIR.h2_spectral_bound
#print axioms CelestialIR.h2_harmonic_bound
#print axioms CelestialIR.eA_inverse
#print axioms CelestialIR.eB_inverse
#print axioms CelestialIR.eA_support
#print axioms CelestialIR.eB_support
#print axioms CelestialIR.res_w_after_z
#print axioms CelestialIR.res_z_after_w
#print axioms CelestialIR.ordered_residues
#print axioms CelestialIR.mixed_divisor
#print axioms CelestialIR.ordered_targets_gap
#print axioms CelestialIR.ordered_targets_gap_int
#print axioms CelestialIR.mem_tuples
#print axioms CelestialIR.hockey_stick
#print axioms CelestialIR.length_tuples
#print axioms CelestialIR.gap_identity
#print axioms CelestialIR.cheb_degree
#print axioms CelestialIR.cheb_leading
#print axioms CelestialIR.chebS_leading
#print axioms CelestialIR.chebS_degree
#print axioms CelestialIR.cheb_eval
#print axioms CelestialIR.chebS_eval
#print axioms CelestialIR.chebVal_one_neg_one
#print axioms CelestialIR.cheb_m1
#print axioms CelestialIR.cheb_m2
#print axioms CelestialIR.cheb_m3
#print axioms CelestialIR.cheb_m4
#print axioms CelestialIR.cheb_m4_interval
