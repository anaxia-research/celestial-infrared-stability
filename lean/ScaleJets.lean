import Std

/-!
Symbolic algebra of the neutral double-pole jet over the integers.

This uses Lean's standard library only. It proves the shear action, its
inverse, and cancellation by a reference-scale compensator for arbitrary
integer parameters. It does not formalize Mellin integration, meromorphic
continuation, connections, or the general real/complex analytic statement.
-/

namespace ScaleJets

/-- Coefficients are ordered (residue, double-pole coefficient). -/
def shift (s : Int) (v : Int × Int) : Int × Int :=
  (v.1 - s * v.2, v.2)

theorem shift_zero (v : Int × Int) : shift 0 v = v := by
  rcases v with ⟨a, b⟩
  simp [shift]

/-- The additive scale group acts by triangular shears. -/
theorem shift_group (s t : Int) (v : Int × Int) :
    shift t (shift s v) = shift (s + t) v := by
  rcases v with ⟨a, b⟩
  apply Prod.ext
  · change a - s * b - t * b = a - (s + t) * b
    rw [Int.add_mul]
    omega
  · rfl

theorem shift_inverse (s : Int) (v : Int × Int) :
    shift (-s) (shift s v) = v := by
  rw [shift_group]
  have cancel : s + -s = 0 := by omega
  rw [cancel]
  exact shift_zero v

/-- log(mu)' = log(mu) - s cancels the entire nilpotent jet action. -/
theorem reference_compensator (ell s : Int) (v : Int × Int) :
    shift (ell - s) (shift s v) = shift ell v := by
  rw [shift_group]
  congr 1
  omega

/-- Counterexample to invariance of the residue of a neutral double pole. -/
theorem residue_negative_control : (shift 5 (3, 2)).1 ≠ 3 := by
  decide

theorem leading_coefficient_invariant (s : Int) (v : Int × Int) :
    (shift s v).2 = v.2 := by
  rfl

end ScaleJets

#print axioms ScaleJets.shift_zero
#print axioms ScaleJets.shift_group
#print axioms ScaleJets.shift_inverse
#print axioms ScaleJets.reference_compensator
#print axioms ScaleJets.residue_negative_control
#print axioms ScaleJets.leading_coefficient_invariant
