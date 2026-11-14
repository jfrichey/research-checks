module

public import Balloons.Finite

@[expose] public section

/-!
# Relabeling the finite coupling

Injective relabelings preserve the event scan and discrepancy counts. This is
the deterministic input for later exchangeability statements; no distributional
symmetry is assumed or concluded here.
-/

namespace Balloons

open scoped symmDiff

variable {α β : Type*} [DecidableEq α] [DecidableEq β]

/-- Transport the endpoints of an event along an embedding. -/
def Event.map (f : α ↪ β) (e : Event α) : Event β :=
  ⟨f e.left, f e.right, fun h => e.distinct (f.injective h)⟩

/-- One collision commutes with injective relabeling. -/
theorem step_map (f : α ↪ β) (e : Event α) (A : Finset α) :
    step (e.map f) (A.map f) = (step e A).map f := by
  unfold step Event.map
  simp only [Finset.mem_map']
  split <;> simp_all

/-- The entire event scan commutes with injective relabeling. -/
theorem run_map (f : α ↪ β) (es : List (Event α)) (A : Finset α) :
    run (es.map (Event.map f)) (A.map f) = (run es A).map f := by
  induction es generalizing A with
  | nil => rfl
  | cons e es ih =>
    simp only [List.map_cons, run_cons, step_map]
    exact ih (step e A)

/-- Injective relabeling preserves the active discrepancy count. -/
theorem discrepancy_map (f : α ↪ β) (A B : Finset α) :
    discrepancy (A.map f) (B.map f) = discrepancy A B := by
  have hm : (A ∆ B).map f = A.map f ∆ B.map f := by
    simp only [Finset.symmDiff_def, Finset.map_union, Finset.map_sdiff]
  unfold discrepancy
  rw [← hm, Finset.card_map]

/-- Relabeling both coupled copies leaves their discrepancy profile unchanged. -/
theorem run_map_discrepancy (f : α ↪ β) (es : List (Event α)) (A B : Finset α) :
    discrepancy (run (es.map (Event.map f)) (A.map f))
      (run (es.map (Event.map f)) (B.map f)) =
        discrepancy (run es A) (run es B) := by
  rw [run_map, run_map, discrepancy_map]

end Balloons
