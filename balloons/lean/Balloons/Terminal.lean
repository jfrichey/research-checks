module

public import Balloons.Finite

@[expose] public section

/-!
# Completion and parity of the finite process

A schedule covering all pairs in the initial configuration leaves one active
vertex for odd input size and none for even input size. Distances and ties do
not enter this conclusion. The coverage hypothesis is explicit.
-/

namespace Balloons

variable {α : Type*} [DecidableEq α]

/-- An accepted collision removes exactly two active vertices. -/
theorem step_card_of_accept (e : Event α) (A : Finset α)
    (h : e.left ∈ A ∧ e.right ∈ A) : (step e A).card + 2 = A.card := by
  have h₁ := Finset.card_erase_add_one h.1
  have hr : e.right ∈ A.erase e.left := by
    simp [Finset.mem_erase, e.distinct.symm, h.2]
  have h₂ := Finset.card_erase_add_one hr
  simp only [step, ite_eq_left h]
  omega

/-- A single event either removes two vertices or has no effect. -/
theorem step_card_cases (e : Event α) (A : Finset α) :
    (step e A).card = A.card ∨ (step e A).card + 2 = A.card := by
  by_cases h : e.left ∈ A ∧ e.right ∈ A
  · exact Or.inr (step_card_of_accept e A h)
  · exact Or.inl (by simp [step, h])

theorem step_card_mod_two (e : Event α) (A : Finset α) :
    (step e A).card % 2 = A.card % 2 := by
  have h := step_card_cases e A
  omega

theorem run_card_mod_two (es : List (Event α)) (A : Finset α) :
    (run es A).card % 2 = A.card % 2 := by
  induction es generalizing A with
  | nil => rfl
  | cons e es ih => exact (ih _).trans (step_card_mod_two e A)

/-- Immediately after an event, its two endpoints cannot both remain. -/
theorem step_not_both (e : Event α) (A : Finset α) :
    ¬(e.left ∈ step e A ∧ e.right ∈ step e A) := by
  by_cases h : e.left ∈ A ∧ e.right ∈ A
  · simp [step, h]
  · simp [step, h]

/-- An event that has been scanned excludes its two endpoints from co-survival. -/
theorem run_not_both_of_mem (es : List (Event α)) (A : Finset α)
    (e : Event α) (he : e ∈ es) :
    ¬(e.left ∈ run es A ∧ e.right ∈ run es A) := by
  induction es generalizing A with
  | nil => simp at he
  | cons f fs ih =>
    rcases List.mem_cons.mp he with h | h
    · subst f
      intro ha
      exact step_not_both e A ⟨run_subset fs _ ha.1, run_subset fs _ ha.2⟩
    · exact ih (step f A) h

/-- Every distinct pair of initially active vertices occurs, in either orientation. -/
def Covers (es : List (Event α)) (A : Finset α) : Prop :=
  ∀ u ∈ A, ∀ v ∈ A, u ≠ v →
    ∃ e ∈ es, (e.left = u ∧ e.right = v) ∨ (e.left = v ∧ e.right = u)

/-- A complete schedule leaves at most one active vertex. -/
theorem run_card_le_one (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) : (run es A).card ≤ 1 := by
  apply Finset.card_le_one.mpr
  intro u hu v hv
  by_contra hne
  obtain ⟨e, he, hpair⟩ := hc u (run_subset es A hu) v (run_subset es A hv) hne
  have hn := run_not_both_of_mem es A e he
  rcases hpair with ⟨hl, hr⟩ | ⟨hl, hr⟩
  · exact hn ⟨hl ▸ hu, hr ▸ hv⟩
  · exact hn ⟨hl ▸ hv, hr ▸ hu⟩

/-- Odd input size gives exactly one survivor. -/
theorem run_card_eq_one (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) (hodd : A.card % 2 = 1) : (run es A).card = 1 := by
  have hle := run_card_le_one es A hc
  have hmod := run_card_mod_two es A
  omega

/-- Even input size gives no survivors. -/
theorem run_eq_empty (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) (heven : A.card % 2 = 0) : run es A = ∅ := by
  have hle := run_card_le_one es A hc
  have hmod := run_card_mod_two es A
  apply Finset.card_eq_zero.mp
  omega

/-- The finite odd survivor is unique, with no choice of a favored vertex. -/
theorem exists_unique_survivor (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) (hodd : A.card % 2 = 1) : ∃! x, x ∈ run es A :=
  Finset.card_eq_one_iff_existsUnique.mp (run_card_eq_one es A hc hodd)

end Balloons
