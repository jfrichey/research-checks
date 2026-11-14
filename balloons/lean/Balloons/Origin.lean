module

public import Balloons.Terminal
public import Balloons.Relabel

@[expose] public section

/-!
# Following the actual discrepancy

`follow` transports a single insertion or deletion through the common event
schedule. Its output is a location, not merely a bound on the discrepancy
count. This supplies a finite origin map for later marked constructions.
-/

namespace Balloons

open scoped symmDiff

variable {α : Type*} [DecidableEq α]

/-- At an event involving the discrepancy, it transfers to the other common endpoint. -/
def transfer (e : Event α) (A : Finset α) (z : α) : α :=
  if z = e.left ∧ e.right ∈ A then e.right
  else if z = e.right ∧ e.left ∈ A then e.left
  else z

/-- The pointwise coupling identity for one event. -/
theorem step_toggle (e : Event α) (A : Finset α) (z : α) :
    step e (A ∆ {z}) = step e A ∆ {transfer e A z} := by
  ext x
  have hne := e.distinct
  by_cases hzL : z = e.left <;> by_cases hzR : z = e.right <;>
    by_cases hxL : x = e.left <;> by_cases hxR : x = e.right <;>
    by_cases hxZ : x = z <;>
    by_cases hL : e.left ∈ A <;> by_cases hR : e.right ∈ A
  all_goals simp_all [step, transfer, Finset.mem_symmDiff, Finset.mem_erase, eq_comm]

/-- Transport a discrepancy while the background itself follows its usual scan. -/
def follow : List (Event α) → Finset α → α → α
  | [], _, z => z
  | e :: es, A, z => follow es (step e A) (transfer e A z)

/-- A single changed origin remains a single changed location, with its identity retained. -/
theorem run_toggle (es : List (Event α)) (A : Finset α) (z : α) :
    run es (A ∆ {z}) = run es A ∆ {follow es A z} := by
  induction es generalizing A z with
  | nil => rfl
  | cons e es ih =>
    simp only [run_cons, follow, step_toggle]
    exact ih _ _

/-- The origin map composes exactly across two portions of the schedule. -/
theorem follow_append (es fs : List (Event α)) (A : Finset α) (z : α) :
    follow (es ++ fs) A z = follow fs (run es A) (follow es A z) := by
  induction es generalizing A z with
  | nil => rfl
  | cons e es ih => simpa only [List.cons_append, follow, run_cons] using ih _ _

private theorem toggle_eq_insert (A : Finset α) (z : α) (hz : z ∉ A) :
    A ∆ {z} = insert z A := by
  ext x
  simp only [Finset.mem_symmDiff, Finset.mem_singleton, Finset.mem_insert]
  aesop

/-- Inserting an absent point is the same exact pointwise coupling. -/
theorem run_insert_eq_toggle_follow (es : List (Event α)) (A : Finset α) (z : α)
    (hz : z ∉ A) : run es (insert z A) = run es A ∆ {follow es A z} := by
  rw [← toggle_eq_insert A z hz]
  exact run_toggle es A z

/-- The survivor after insertion into a complete even background is its transported origin. -/
theorem inserted_survivor_eq_follow (es : List (Event α)) (A : Finset α) (z : α)
    (hc : Covers es A) (heven : A.card % 2 = 0) (hz : z ∉ A) :
    run es (insert z A) = {follow es A z} := by
  rw [run_insert_eq_toggle_follow es A z hz, run_eq_empty es A hc heven]
  simp

/-- Agreement after two separate insertions is agreement of their transported origins. -/
theorem inserted_survivors_agree_iff (es : List (Event α)) (A : Finset α) (z w : α)
    (hc : Covers es A) (heven : A.card % 2 = 0) (hz : z ∉ A) (hw : w ∉ A) :
    run es (insert z A) = run es (insert w A) ↔ follow es A z = follow es A w := by
  rw [inserted_survivor_eq_follow es A z hc heven hz,
    inserted_survivor_eq_follow es A w hc heven hw]
  simp

/-- A transfer can only move to a background point or retain the original discrepancy. -/
theorem transfer_mem_insert (e : Event α) (A : Finset α) (z : α) :
    transfer e A z ∈ insert z A := by
  unfold transfer
  split_ifs with h₁ h₂
  · exact Finset.mem_insert_of_mem h₁.2
  · exact Finset.mem_insert_of_mem h₂.2
  · exact Finset.mem_insert_self _ _

/-- The transported origin is an original background vertex or the origin itself. -/
theorem follow_mem_insert (es : List (Event α)) (A : Finset α) (z : α) :
    follow es A z ∈ insert z A := by
  induction es generalizing A z with
  | nil => simp [follow]
  | cons e es ih =>
    rcases Finset.mem_insert.mp (ih (step e A) (transfer e A z)) with h | h
    · simp only [follow, h]
      exact transfer_mem_insert e A z
    · exact Finset.mem_insert_of_mem (step_subset e A h)

/-- The finite origin map has only old-vertex outputs and its possible unchanged origin. -/
theorem follow_eq_self_or_mem (es : List (Event α)) (A : Finset α) (z : α) :
    follow es A z = z ∨ follow es A z ∈ A :=
  Finset.mem_insert.mp (follow_mem_insert es A z)

variable {β : Type*} [DecidableEq β]

/-- Transporting one discrepancy commutes with injective relabeling. -/
theorem transfer_map (f : α ↪ β) (e : Event α) (A : Finset α) (z : α) :
    transfer (e.map f) (A.map f) (f z) = f (transfer e A z) := by
  by_cases h₁ : z = e.left ∧ e.right ∈ A <;>
    by_cases h₂ : z = e.right ∧ e.left ∈ A <;>
      simp [transfer, Event.map, f.injective.eq_iff, h₁, h₂, e.distinct.symm]

/-- The whole origin map commutes with injective relabeling. -/
theorem follow_map (f : α ↪ β) (es : List (Event α)) (A : Finset α) (z : α) :
    follow (es.map (Event.map f)) (A.map f) (f z) = f (follow es A z) := by
  induction es generalizing A z with
  | nil => rfl
  | cons e es ih =>
    simp only [List.map_cons, follow, step_map, transfer_map]
    exact ih _ _

end Balloons
