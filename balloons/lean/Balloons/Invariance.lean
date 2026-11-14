module

public import Balloons.HistoryDistance
public import Balloons.Time

@[expose] public section

/-!
# Independence of orientation and sorted tie order

Two sorted permutations have identical priority sequences. If a tied priority
can only represent the same unordered pair, corresponding events have the
same action and record the same spatial edge. This proves equality of the
whole spatial history as well as the final active configuration.
-/

namespace Balloons

variable {α τ : Type*} [DecidableEq α]

/-- Whether an event is accepted depends only on its unordered endpoints. -/
theorem event_active_iff (e : Event α) (A : Finset α) :
    e.vertices ⊆ A ↔ e.left ∈ A ∧ e.right ∈ A := by
  simp [Event.vertices, Finset.insert_subset_iff, Finset.singleton_subset_iff]

/-- The event action can be written entirely in terms of its unordered pair. -/
theorem step_vertices (e : Event α) (A : Finset α) :
    step e A = if e.vertices ⊆ A then A \ e.vertices else A := by
  simp only [event_active_iff]
  unfold step
  split_ifs with h
  · ext x
    simp only [Finset.mem_erase, Finset.mem_sdiff, Event.vertices,
      Finset.mem_insert, Finset.mem_singleton]
    tauto
  · rfl

/-- Reversing or otherwise representing the same spatial pair has the same action. -/
theorem step_eq_of_vertices_eq (e f : Event α) (A : Finset α)
    (h : e.vertices = f.vertices) : step e A = step f A := by
  rw [step_vertices, step_vertices, h]

/-- The spatial record of one step also depends only on the unordered pair. -/
theorem accepted_vertices_cons (e : Event α) (es : List (Event α)) (A : Finset α) :
    (accepted (e :: es) A).map Event.vertices =
      if e.vertices ⊆ A then e.vertices :: (accepted es (step e A)).map Event.vertices
      else (accepted es (step e A)).map Event.vertices := by
  simp only [event_active_iff]
  by_cases h : e.left ∈ A ∧ e.right ∈ A <;> simp [accepted, h]

/-- Equal priority sequences give identical configurations and spatial histories
when any corresponding tie represents the same unordered pair. -/
theorem run_history_eq_of_priority_sequences_eq (p : Event α → τ)
    (es fs : List (Event α)) (hseq : es.map p = fs.map p)
    (hcompat : ∀ e ∈ es, ∀ f ∈ fs, p e = p f → e.vertices = f.vertices) (A : Finset α) :
    run es A = run fs A ∧
      (accepted es A).map Event.vertices = (accepted fs A).map Event.vertices := by
  induction es generalizing fs A with
  | nil =>
    cases fs <;> simp_all [accepted]
  | cons e es ih =>
    cases fs with
    | nil => simp at hseq
    | cons f fs =>
      obtain ⟨hhead, htail⟩ := List.cons.inj (by simpa only [List.map_cons] using hseq)
      have hv := hcompat e (List.mem_cons_self) f (List.mem_cons_self) hhead
      have hs := step_eq_of_vertices_eq e f A hv
      have hc : ∀ e' ∈ es, ∀ f' ∈ fs, p e' = p f' → e'.vertices = f'.vertices := by
        intro e' he f' hf
        exact hcompat e' (List.mem_cons_of_mem e he) f' (List.mem_cons_of_mem f hf)
      have ht := ih fs htail hc (step f A)
      constructor
      · simpa only [run_cons, hs] using ht.1
      · simp only [accepted_vertices_cons, hv, hs, ht.2]

variable [LinearOrder τ]

omit [DecidableEq α] in
/-- Sorted permutations have identical priority sequences, even when priorities repeat. -/
theorem sorted_perm_priority_sequences_eq (p : Event α → τ) (es fs : List (Event α))
    (hp : es.Perm fs) (he : es.Pairwise fun e f => p e ≤ p f)
    (hf : fs.Pairwise fun e f => p e ≤ p f) : es.map p = fs.map p := by
  exact (hp.map p).eq_of_pairwise'
    (by simpa only [List.pairwise_map] using he)
    (by simpa only [List.pairwise_map] using hf)

/-- Sorting a fixed event collection gives the same trajectory and spatial history
whenever all priority ties are merely orientations of the same pair. -/
theorem run_history_eq_of_sorted_perm (p : Event α → τ) (es fs : List (Event α))
    (hp : es.Perm fs) (he : es.Pairwise fun e f => p e ≤ p f)
    (hf : fs.Pairwise fun e f => p e ≤ p f)
    (hcompat : ∀ e ∈ es, ∀ f ∈ fs, p e = p f → e.vertices = f.vertices) (A : Finset α) :
    run es A = run fs A ∧ history es A = history fs A := by
  obtain ⟨hrun, hhist⟩ := run_history_eq_of_priority_sequences_eq p es fs
    (sorted_perm_priority_sequences_eq p es fs hp he hf) hcompat A
  exact ⟨hrun, congrArg List.toFinset hhist⟩

/-- The same invariance holds at every separation cutoff, so it preserves the trajectory. -/
theorem cutoff_run_history_eq_of_sorted_perm (p : Event α → τ) (r : τ)
    (es fs : List (Event α)) (hp : es.Perm fs)
    (he : es.Pairwise fun e f => p e ≤ p f) (hf : fs.Pairwise fun e f => p e ≤ p f)
    (hcompat : ∀ e ∈ es, ∀ f ∈ fs, p e = p f → e.vertices = f.vertices) (A : Finset α) :
    run (eventsThrough p r es) A = run (eventsThrough p r fs) A ∧
      history (eventsThrough p r es) A = history (eventsThrough p r fs) A := by
  rw [eventsThrough_eq_filter p r es he, eventsThrough_eq_filter p r fs hf]
  exact run_history_eq_of_sorted_perm p _ _ (hp.filter _) (he.filter _) (hf.filter _)
    (fun e he' f hf' => hcompat e (List.mem_filter.mp he').1 f (List.mem_filter.mp hf').1) A

end Balloons
