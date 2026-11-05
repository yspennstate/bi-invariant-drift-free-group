import Mathlib.Data.Nat.Basic
import Lean.Elab.Tactic.Omega

-- The exponent comparison in the squared Catalan union-bound estimate.
-- n = k + 2*m, with at most a quarter of positions left unmatched.
theorem cancellation_union_bound_exponent
    (n k m : Nat) (hcount : k + 2 * m = n) (hsmall : 4 * k ≤ n) :
    3 * n ≤ 8 * m := by
  omega

theorem cancellation_pair_fraction
    (n k m : Nat) (hcount : k + 2 * m = n) (hsmall : 4 * k ≤ n) :
    3 * n ≤ 4 * (2 * m) := by
  omega
