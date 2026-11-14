module

public import Balloons.History

@[expose] public section

/-!
# Events involving initially absent vertices are irrelevant

Both coupled inputs may be run on a larger common ambient schedule. Filtering
that schedule to events supported by one input changes neither its active
configuration nor its accepted collision history.
-/

namespace Balloons

variable {α : Type*} [DecidableEq α]

/-- Retain only the events with both endpoints in an ambient input set. -/
def supportedEvents (S : Finset α) (es : List (Event α)) : List (Event α) :=
  es.filter fun e => decide (e.left ∈ S ∧ e.right ∈ S)

/-- Filtering away events outside an initial ambient set preserves every sub-input. -/
theorem run_supportedEvents (S A : Finset α) (hAS : A ⊆ S) (es : List (Event α)) :
    run (supportedEvents S es) A = run es A := by
  induction es generalizing A with
  | nil => rfl
  | cons e es ih =>
    by_cases h : e.left ∈ S ∧ e.right ∈ S
    · have ht : step e A ⊆ S := (step_subset e A).trans hAS
      simpa [supportedEvents, h] using ih (step e A) ht
    · have hA : ¬(e.left ∈ A ∧ e.right ∈ A) :=
        fun ha => h ⟨hAS ha.1, hAS ha.2⟩
      simpa [supportedEvents, h, step, hA] using ih A hAS

/-- The actual accepted history also ignores initially unsupported events. -/
theorem accepted_supportedEvents (S A : Finset α) (hAS : A ⊆ S)
    (es : List (Event α)) : accepted (supportedEvents S es) A = accepted es A := by
  induction es generalizing A with
  | nil => rfl
  | cons e es ih =>
    by_cases h : e.left ∈ S ∧ e.right ∈ S
    · have ht : step e A ⊆ S := (step_subset e A).trans hAS
      simp only [supportedEvents, List.filter_cons, decide_eq_true_eq, ite_eq_left h, accepted]
      have hi := ih (step e A) ht
      simp only [supportedEvents] at hi
      simp only [hi]
    · have hA : ¬(e.left ∈ A ∧ e.right ∈ A) :=
        fun ha => h ⟨hAS ha.1, hAS ha.2⟩
      simpa [supportedEvents, h, accepted, step, hA] using ih A hAS

/-- A common ambient schedule reduces to the input's own supported events. -/
theorem run_restrict_initial (A : Finset α) (es : List (Event α)) :
    run (supportedEvents A es) A = run es A :=
  run_supportedEvents A A (Finset.Subset.refl _) es

end Balloons
