module

public import Balloons.Time
public import Balloons.HistoryDistance
public import Mathlib.Topology.MetricSpace.Defs

@[expose] public section

/-!
# The finite metric balloon process

Pairs are scanned in increasing distance. The cutoff is separation, twice
physical time. A common ambient set gives all coupled inputs the same tie
order; no assumption that arbitrary ties give the same history is made.
-/

noncomputable section

namespace Balloons

variable {α : Type*} [MetricSpace α] [DecidableEq α]

/-- Separation at which the two growing balloons would meet. -/
def edgeLength (e : Event α) : ℝ := dist e.left e.right

/-- All ambient pairs, sorted in increasing separation. -/
def metricSchedule (S : Finset α) : List (Event α) := timedEvents S edgeLength

/-- All pair events scanned by a separation cutoff. -/
def eventsAt (S : Finset α) (r : ℝ) : List (Event α) :=
  eventsThrough edgeLength r (metricSchedule S)

/-- The active vertices in the metric process at separation `r`. -/
def activeAt (S A : Finset α) (r : ℝ) : Finset α := run (eventsAt S r) A

/-- Unordered spatial matching edges accepted by separation `r`. -/
def historyAt (S A : Finset α) (r : ℝ) : Finset (Finset α) := history (eventsAt S r) A

/-- The metric schedule is actually sorted by distance. -/
theorem metricSchedule_sorted (S : Finset α) :
    (metricSchedule S).Pairwise fun e f => edgeLength e ≤ edgeLength f :=
  timedEvents_sorted S edgeLength

/-- Every ambient pair is considered. -/
theorem metricSchedule_covers (S : Finset α) : Covers (metricSchedule S) S :=
  timedEvents_covers S edgeLength

/-- A pair is scanned by a cutoff exactly when its distance is at most that cutoff. -/
theorem mem_eventsAt (S : Finset α) (r : ℝ) (e : Event α) :
    e ∈ eventsAt S r ↔ e ∈ metricSchedule S ∧ dist e.left e.right ≤ r :=
  mem_eventsThrough edgeLength r (metricSchedule S) (metricSchedule_sorted S) e

/-- The metric process only removes points. -/
theorem activeAt_subset (S A : Finset α) (r : ℝ) : activeAt S A r ⊆ A :=
  run_subset _ A

/-- Any two remaining balloons are farther apart than the separation already scanned. -/
theorem activeAt_separated (S A : Finset α) (hAS : A ⊆ S) (r : ℝ)
    {u v : α} (hu : u ∈ activeAt S A r) (hv : v ∈ activeAt S A r) (hne : u ≠ v) :
    r < dist u v := by
  apply lt_of_not_ge
  intro hd
  obtain ⟨e, he, hp⟩ := metricSchedule_covers S u
    (hAS (activeAt_subset S A r hu)) v (hAS (activeAt_subset S A r hv)) hne
  have hlen : dist e.left e.right ≤ r := by
    rcases hp with ⟨hl, hr⟩ | ⟨hl, hr⟩
    · simpa only [hl, hr] using hd
    · simpa only [hl, hr, dist_comm v u] using hd
  have hn := run_not_both_of_mem (eventsAt S r) A e
    ((mem_eventsAt S r e).mpr ⟨he, hlen⟩)
  rcases hp with ⟨hl, hr⟩ | ⟨hl, hr⟩
  · exact hn ⟨hl ▸ hu, hr ▸ hv⟩
  · exact hn ⟨hl ▸ hv, hr ▸ hu⟩

/-- Distinct starting locations have no collision at separation zero. -/
theorem activeAt_zero (S A : Finset α) : activeAt S A 0 = A := by
  have he : eventsAt S 0 = [] := by
    apply List.eq_nil_iff_forall_not_mem.mpr
    intro e he
    have hd := ((mem_eventsAt S 0 e).mp he).2
    exact e.distinct (dist_eq_zero.mp (le_antisymm hd dist_nonneg))
  simp [activeAt, he]

/-- A cutoff bounding all input distances has considered every input pair. -/
theorem eventsAt_covers (S A : Finset α) (hAS : A ⊆ S) (r : ℝ)
    (hbound : ∀ u ∈ A, ∀ v ∈ A, dist u v ≤ r) : Covers (eventsAt S r) A := by
  intro u hu v hv huv
  obtain ⟨e, he, hp⟩ := metricSchedule_covers S u (hAS hu) v (hAS hv) huv
  refine ⟨e, (mem_eventsAt S r e).mpr ⟨he, ?_⟩, hp⟩
  rcases hp with ⟨hl, hr⟩ | ⟨hl, hr⟩
  · simpa only [hl, hr] using hbound u hu v hv
  · simpa only [hl, hr] using hbound v hv u hu

/-- A sufficiently large cutoff leaves exactly one balloon for an odd input. -/
theorem activeAt_odd_card (S A : Finset α) (hAS : A ⊆ S) (r : ℝ)
    (hbound : ∀ u ∈ A, ∀ v ∈ A, dist u v ≤ r) (hodd : A.card % 2 = 1) :
    (activeAt S A r).card = 1 :=
  run_card_eq_one _ A (eventsAt_covers S A hAS r hbound) hodd

/-- A sufficiently large cutoff removes every balloon for an even input. -/
theorem activeAt_even_empty (S A : Finset α) (hAS : A ⊆ S) (r : ℝ)
    (hbound : ∀ u ∈ A, ∀ v ∈ A, dist u v ≤ r) (heven : A.card % 2 = 0) :
    activeAt S A r = ∅ :=
  run_eq_empty _ A (eventsAt_covers S A hAS r hbound) heven

/-- The deterministic contraction theorem now holds at physical separation cutoffs. -/
theorem activeAt_discrepancy_le (S A B : Finset α) {r s : ℝ} (hrs : r ≤ s) :
    discrepancy (activeAt S A s) (activeAt S B s) ≤
      discrepancy (activeAt S A r) (activeAt S B r) :=
  activeThrough_discrepancy_le edgeLength hrs (metricSchedule S) A B

/-- Replacing at most `k` input points creates at most `2*k` active discrepancies. -/
theorem activeAt_resampling_bound (S C P Q : Finset α) (r : ℝ) (k : ℕ)
    (hp : P.card ≤ k) (hq : Q.card ≤ k) :
    discrepancy (activeAt S (C ∪ P) r) (activeAt S (C ∪ Q) r) ≤ 2 * k :=
  resampling_discrepancy_le _ C P Q k hp hq

/-- Every accepted pair has reached its collision separation. -/
theorem acceptedAt_length_le (S A : Finset α) (r : ℝ) (e : Event α)
    (he : e ∈ accepted (eventsAt S r) A) : dist e.left e.right ≤ r :=
  ((mem_eventsAt S r e).mp (accepted_mem_schedule _ A e he)).2

/-- Collision counting and the active count agree at every metric cutoff. -/
theorem historyAt_balance (S A : Finset α) (r : ℝ) :
    (activeAt S A r).card + 2 * (historyAt S A r).card = A.card :=
  history_card_balance _ A

end Balloons
