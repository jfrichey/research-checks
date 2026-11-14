module

public import Mathlib.Data.Finset.Card
public import Mathlib.Data.Finset.SymmDiff
public import Mathlib.Tactic.Tauto

@[expose] public section

/-!
# Finite balloon processes under a common event order

The event list is arbitrary. For the geometric model it is the list of pairs
ordered by distance, using the same tie order in both inputs. No probability
or geometric premise is needed for discrepancy contraction.

Informal source: `balloon-process.tex`, Lemma 2.2 and Theorem 7.31.
-/

namespace Balloons

open scoped symmDiff

variable {α : Type*} [DecidableEq α]

/-- A proposed collision between two distinct vertices. -/
structure Event (α : Type*) where
  left : α
  right : α
  distinct : left ≠ right

/-- Accept an event precisely when both endpoints remain active. -/
def step (e : Event α) (A : Finset α) : Finset α :=
  if e.left ∈ A ∧ e.right ∈ A then (A.erase e.left).erase e.right else A

/-- Scan an ordered finite event list. -/
def run : List (Event α) → Finset α → Finset α
  | [], A => A
  | e :: es, A => run es (step e A)

/-- Number of active locations present in exactly one of two copies. -/
def discrepancy (A B : Finset α) : ℕ := (A ∆ B).card

@[simp] theorem run_nil (A : Finset α) : run [] A = A := rfl

@[simp] theorem run_cons (e : Event α) (es : List (Event α)) (A : Finset α) :
    run (e :: es) A = run es (step e A) := rfl

theorem run_append (es fs : List (Event α)) (A : Finset α) :
    run (es ++ fs) A = run fs (run es A) := by
  induction es generalizing A with
  | nil => rfl
  | cons e es ih => simpa only [List.cons_append, run_cons] using ih (step e A)

theorem step_subset (e : Event α) (A : Finset α) : step e A ⊆ A := by
  unfold step
  split
  · exact Finset.Subset.trans (Finset.erase_subset _ _) (Finset.erase_subset _ _)
  · exact Finset.Subset.refl _

theorem run_subset (es : List (Event α)) (A : Finset α) : run es A ⊆ A := by
  induction es generalizing A with
  | nil => exact Finset.Subset.refl _
  | cons e es ih => exact Finset.Subset.trans (ih _) (step_subset e A)

theorem step_mem_iff_of_ne (e : Event α) (A : Finset α) {x : α}
    (hx : x ≠ e.left) (hy : x ≠ e.right) : x ∈ step e A ↔ x ∈ A := by
  unfold step
  split <;> simp_all [Finset.mem_erase]

private def bit (A : Finset α) (x : α) : ℕ := if x ∈ A then 1 else 0

private theorem card_erase_add_bit (A : Finset α) (x : α) :
    (A.erase x).card + bit A x = A.card := by
  by_cases h : x ∈ A
  · simpa [bit, h] using Finset.card_erase_add_one h
  · simp [bit, h]

private theorem card_split_pair (A : Finset α) (e : Event α) :
    A.card = ((A.erase e.left).erase e.right).card + bit A e.left + bit A e.right := by
  have h₁ := card_erase_add_bit A e.left
  have h₂ := card_erase_add_bit (A.erase e.left) e.right
  have hbit : bit (A.erase e.left) e.right = bit A e.right := by
    simp [bit, Finset.mem_erase, e.distinct.symm]
  rw [hbit] at h₂
  omega

private theorem outside_step_symmDiff (e : Event α) (A B : Finset α) :
    (((step e A ∆ step e B).erase e.left).erase e.right) =
      (((A ∆ B).erase e.left).erase e.right) := by
  ext x
  by_cases hx : x = e.left
  · subst x
    simp
  by_cases hy : x = e.right
  · subst x
    simp
  simp [Finset.mem_erase, hx, hy, Finset.mem_symmDiff,
    step_mem_iff_of_ne e A hx hy, step_mem_iff_of_ne e B hx hy]

/-- The exact local identity: discrepancies either persist or annihilate in pairs. -/
theorem step_discrepancy_cases (e : Event α) (A B : Finset α) :
    discrepancy (step e A) (step e B) = discrepancy A B ∨
      discrepancy (step e A) (step e B) + 2 = discrepancy A B := by
  unfold discrepancy
  have hbefore := card_split_pair (A ∆ B) e
  have hafter := card_split_pair (step e A ∆ step e B) e
  rw [outside_step_symmDiff] at hafter
  by_cases ha : e.left ∈ A <;> by_cases hb : e.right ∈ A <;>
    by_cases hc : e.left ∈ B <;> by_cases hd : e.right ∈ B
  all_goals
    simp [bit, step, Finset.mem_symmDiff, Finset.mem_erase,
      ha, hb, hc, hd, e.distinct, e.distinct.symm] at hbefore hafter ⊢
  all_goals omega

/-- One collision cannot increase the number of active discrepancies. -/
theorem step_discrepancy_le (e : Event α) (A B : Finset α) :
    discrepancy (step e A) (step e B) ≤ discrepancy A B := by
  have h := step_discrepancy_cases e A B
  omega

/-- One collision preserves the parity of the discrepancy count. -/
theorem step_discrepancy_mod_two (e : Event α) (A B : Finset α) :
    discrepancy (step e A) (step e B) % 2 = discrepancy A B % 2 := by
  have h := step_discrepancy_cases e A B
  omega

/-- Contraction for every finite common event schedule. -/
theorem run_discrepancy_le (es : List (Event α)) (A B : Finset α) :
    discrepancy (run es A) (run es B) ≤ discrepancy A B := by
  induction es generalizing A B with
  | nil => exact le_rfl
  | cons e es ih => exact (ih _ _).trans (step_discrepancy_le e A B)

/-- Parity preservation for every finite common event schedule. -/
theorem run_discrepancy_mod_two (es : List (Event α)) (A B : Finset α) :
    discrepancy (run es A) (run es B) % 2 = discrepancy A B % 2 := by
  induction es generalizing A B with
  | nil => rfl
  | cons e es ih => exact (ih _ _).trans (step_discrepancy_mod_two e A B)

/-- Contraction between two stages, expressed as an initial and a further schedule. -/
theorem run_append_discrepancy_le (es fs : List (Event α)) (A B : Finset α) :
    discrepancy (run (es ++ fs) A) (run (es ++ fs) B) ≤
      discrepancy (run es A) (run es B) := by
  rw [run_append, run_append]
  exact run_discrepancy_le fs _ _

@[simp] theorem discrepancy_eq_zero_iff (A B : Finset α) :
    discrepancy A B = 0 ↔ A = B := by
  simp [discrepancy]

/-- Agreement at one stage is absorbing under every common continuation. -/
theorem agreement_absorbing (es fs : List (Event α)) (A B : Finset α)
    (h : run es A = run es B) : run (es ++ fs) A = run (es ++ fs) B := by
  rw [run_append, run_append, h]

/-- Exactly one discrepancy remains exactly one, including at terminal states. -/
theorem run_discrepancy_eq_one (es : List (Event α)) (A B : Finset α)
    (h : discrepancy A B = 1) : discrepancy (run es A) (run es B) = 1 := by
  have hle := run_discrepancy_le es A B
  have hmod := run_discrepancy_mod_two es A B
  omega

/-- Adding a previously absent vertex initially creates exactly one discrepancy. -/
theorem discrepancy_insert (A : Finset α) (z : α) (hz : z ∉ A) :
    discrepancy A (insert z A) = 1 := by
  have hs : A ∆ insert z A = {z} := by
    ext x
    simp only [Finset.mem_symmDiff, Finset.mem_insert, Finset.mem_singleton]
    aesop
  simp [discrepancy, hs]

/-- The single moving discrepancy lemma, for an arbitrary event schedule. -/
theorem run_insert_discrepancy (es : List (Event α)) (A : Finset α) (z : α)
    (hz : z ∉ A) : discrepancy (run es A) (run es (insert z A)) = 1 :=
  run_discrepancy_eq_one es A (insert z A) (discrepancy_insert A z hz)

/-- A common background cancels from the initial edit budget. -/
theorem discrepancy_union_le (C P Q : Finset α) :
    discrepancy (C ∪ P) (C ∪ Q) ≤ P.card + Q.card := by
  have hs : (C ∪ P) ∆ (C ∪ Q) ⊆ P ∪ Q := by
    intro x hx
    simp only [Finset.mem_symmDiff, Finset.mem_union] at hx ⊢
    tauto
  exact (Finset.card_le_card hs).trans (Finset.card_union_le P Q)

/-- At most `2*k` discrepancies from two sets of at most `k` replacement locations. -/
theorem resampling_discrepancy_le (es : List (Event α)) (C P Q : Finset α)
    (k : ℕ) (hp : P.card ≤ k) (hq : Q.card ≤ k) :
    discrepancy (run es (C ∪ P)) (run es (C ∪ Q)) ≤ 2 * k := by
  have h₁ := run_discrepancy_le es (C ∪ P) (C ∪ Q)
  have h₂ := discrepancy_union_le C P Q
  omega

/-- The total decrease is twice a nonnegative integer, with no asymptotic premise. -/
theorem run_discrepancy_pair_losses (es : List (Event α)) (A B : Finset α) :
    ∃ j : ℕ, discrepancy (run es A) (run es B) + 2 * j = discrepancy A B := by
  have hle := run_discrepancy_le es A B
  have hmod := run_discrepancy_mod_two es A B
  refine ⟨(discrepancy A B - discrepancy (run es A) (run es B)) / 2, ?_⟩
  omega

end Balloons
