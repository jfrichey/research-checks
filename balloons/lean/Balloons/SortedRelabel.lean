module
public import Balloons.Enumeration

@[expose] public section

/-! # Relabeling the complete distance-sorted process -/

namespace Balloons

variable {α β τ : Type*} [DecidableEq α] [DecidableEq β] [LinearOrder τ]

/-- Sorting preserves the exact ambient event membership. -/
theorem mem_timedEvents (S : Finset α) (p : Event α → τ) (e : Event α) :
    e ∈ timedEvents S p ↔ e.left ∈ S ∧ e.right ∈ S := by
  change e ∈ (allEvents S).mergeSort _ ↔ _
  rw [(List.mergeSort_perm _ _).mem_iff, mem_allEvents]

/-- Relabeling a sorted complete schedule preserves its event multiset. -/
theorem timedEvents_map_perm (f : α ↪ β) (S : Finset α) (p : Event β → τ) :
    ((timedEvents S (fun e => p (e.map f))).map (Event.map f)).Perm
      (timedEvents (S.map f) p) := by
  unfold timedEvents orderedEvents
  exact ((List.mergeSort_perm _ _).map _).trans
    ((allEvents_map_perm f S).trans (List.mergeSort_perm _ _).symm)

/-- The complete sorted metric scan is equivariant, including its spatial history,
when equal priorities in the ambient set can only represent one unordered pair. -/
theorem timed_run_history_map (f : α ↪ β) (S A : Finset α) (p : Event β → τ)
    (htie : ∀ e ∈ allEvents (S.map f), ∀ g ∈ allEvents (S.map f),
      p e = p g → e.vertices = g.vertices) :
    run (timedEvents (S.map f) p) (A.map f) =
      (run (timedEvents S (fun e => p (e.map f))) A).map f ∧
    history (timedEvents (S.map f) p) (A.map f) =
      (history (timedEvents S (fun e => p (e.map f))) A).image (fun edge => edge.map f) := by
  let es := timedEvents S (fun e => p (e.map f))
  let fs := timedEvents (S.map f) p
  have hp : (es.map (Event.map f)).Perm fs := timedEvents_map_perm f S p
  have he : (es.map (Event.map f)).Pairwise fun e g => p e ≤ p g := by
    rw [List.pairwise_map]
    exact timedEvents_sorted S _
  have hf : fs.Pairwise fun e g => p e ≤ p g := timedEvents_sorted _ _
  have hc : ∀ e ∈ es.map (Event.map f), ∀ g ∈ fs,
      p e = p g → e.vertices = g.vertices := by
    intro e he' g hg
    exact htie e ((mem_allEvents (S.map f) e).mpr
      ((mem_timedEvents (S.map f) p e).mp (hp.mem_iff.mp he')))
      g ((mem_allEvents (S.map f) g).mpr ((mem_timedEvents (S.map f) p g).mp hg))
  obtain ⟨hrun, hhist⟩ := run_history_eq_of_sorted_perm p _ _ hp he hf hc (A.map f)
  rw [run_map] at hrun
  rw [history_map] at hhist
  exact ⟨hrun.symm, hhist.symm⟩

end Balloons
