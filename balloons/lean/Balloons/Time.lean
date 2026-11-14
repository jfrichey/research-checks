module

public import Balloons.Schedule
public import Balloons.History

@[expose] public section

/-!
# Ordered collision times and cutoff configurations

Priorities take values in any linear order. The geometric process will use
pair distance as its priority. The definitions and cutoff contraction below
already apply to that choice; geometric distributional statements are separate.
-/

namespace Balloons

variable {α τ : Type*} [DecidableEq α] [LinearOrder τ]

/-- The initial part of a schedule whose priorities are at most the cutoff. -/
def eventsThrough (p : Event α → τ) (r : τ) (es : List (Event α)) : List (Event α) :=
  es.takeWhile fun e => decide (p e ≤ r)

omit [DecidableEq α] in
/-- Increasing a cutoff extends the scanned prefix. -/
theorem eventsThrough_prefix (p : Event α → τ) {r s : τ} (hrs : r ≤ s)
    (es : List (Event α)) :
    ∃ fs, eventsThrough p r es ++ fs = eventsThrough p s es := by
  induction es with
  | nil => exact ⟨[], rfl⟩
  | cons e es ih =>
    by_cases h : p e ≤ r
    · have hs : p e ≤ s := h.trans hrs
      obtain ⟨fs, hf⟩ := ih
      exact ⟨fs, by simpa [eventsThrough, h, hs] using congrArg (List.cons e) hf⟩
    · exact ⟨eventsThrough p s (e :: es), by simp [eventsThrough, h]⟩

omit [DecidableEq α] in
/-- A sorted schedule's prefix contains exactly the events below the cutoff. -/
theorem eventsThrough_eq_filter (p : Event α → τ) (r : τ) (es : List (Event α))
    (hsorted : es.Pairwise fun e f => p e ≤ p f) :
    eventsThrough p r es = es.filter (fun e => decide (p e ≤ r)) := by
  induction es with
  | nil => rfl
  | cons e es ih =>
    obtain ⟨hhead, htail⟩ := List.pairwise_cons.mp hsorted
    by_cases h : p e ≤ r
    · simpa [eventsThrough, h] using congrArg (List.cons e) (ih htail)
    · have hf : es.filter (fun f => decide (p f ≤ r)) = [] := by
        apply List.filter_eq_nil_iff.mpr
        intro f hmem
        simp only [decide_eq_true_eq]
        exact fun hfr => h ((hhead f hmem).trans hfr)
      simp [eventsThrough, h, hf]

omit [DecidableEq α] in
/-- Exact event membership at a sorted cutoff. -/
theorem mem_eventsThrough (p : Event α → τ) (r : τ) (es : List (Event α))
    (hsorted : es.Pairwise fun e f => p e ≤ p f) (e : Event α) :
    e ∈ eventsThrough p r es ↔ e ∈ es ∧ p e ≤ r := by
  rw [eventsThrough_eq_filter p r es hsorted]
  simp

noncomputable section

/-- Complete event schedule sorted by the supplied collision priority. -/
def timedEvents (S : Finset α) (p : Event α → τ) : List (Event α) :=
  orderedEvents S fun e f => decide (p e ≤ p f)

/-- The explicit schedule is nondecreasing in collision priority. -/
theorem timedEvents_sorted (S : Finset α) (p : Event α → τ) :
    (timedEvents S p).Pairwise fun e f => p e ≤ p f := by
  have ht : ∀ e f g : Event α,
      decide (p e ≤ p f) → decide (p f ≤ p g) → decide (p e ≤ p g) := by
    intro e f g h₁ h₂
    simp only [decide_eq_true_eq] at h₁ h₂ ⊢
    exact h₁.trans h₂
  have htotal : ∀ e f : Event α, decide (p e ≤ p f) || decide (p f ≤ p e) := by
    intro e f
    simpa using le_total (p e) (p f)
  simpa [timedEvents, orderedEvents] using
    List.pairwise_mergeSort ht htotal (allEvents S)

theorem timedEvents_covers (S : Finset α) (p : Event α → τ) :
    Covers (timedEvents S p) S := orderedEvents_covers S _

end

/-- Active configuration after scanning events up to a separation cutoff. -/
def activeThrough (p : Event α → τ) (r : τ) (es : List (Event α)) (A : Finset α) :
    Finset α := run (eventsThrough p r es) A

/-- Discrepancies contract between any two ordered cutoffs. -/
theorem activeThrough_discrepancy_le (p : Event α → τ) {r s : τ} (hrs : r ≤ s)
    (es : List (Event α)) (A B : Finset α) :
    discrepancy (activeThrough p s es A) (activeThrough p s es B) ≤
      discrepancy (activeThrough p r es A) (activeThrough p r es B) := by
  obtain ⟨fs, hfs⟩ := eventsThrough_prefix p hrs es
  simpa only [activeThrough, ← hfs] using
    run_append_discrepancy_le (eventsThrough p r es) fs A B

/-- Every cutoff obeys the initial discrepancy budget. -/
theorem activeThrough_initial_bound (p : Event α → τ) (r : τ)
    (es : List (Event α)) (A B : Finset α) :
    discrepancy (activeThrough p r es A) (activeThrough p r es B) ≤ discrepancy A B :=
  run_discrepancy_le _ A B

/-- One insertion has one active discrepancy at every cutoff. -/
theorem activeThrough_insert_discrepancy (p : Event α → τ) (r : τ)
    (es : List (Event α)) (A : Finset α) (z : α) (hz : z ∉ A) :
    discrepancy (activeThrough p r es A) (activeThrough p r es (insert z A)) = 1 :=
  run_insert_discrepancy _ A z hz

/-- The replacement budget holds simultaneously at all cutoffs. -/
theorem activeThrough_resampling_bound (p : Event α → τ) (r : τ)
    (es : List (Event α)) (C P Q : Finset α) (k : ℕ)
    (hp : P.card ≤ k) (hq : Q.card ≤ k) :
    discrepancy (activeThrough p r es (C ∪ P)) (activeThrough p r es (C ∪ Q)) ≤ 2 * k :=
  resampling_discrepancy_le _ C P Q k hp hq

end Balloons
