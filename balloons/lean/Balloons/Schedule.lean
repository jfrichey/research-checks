module

public import Balloons.Terminal
public import Mathlib.Data.List.Sort

@[expose] public section

/-!
# Explicit pair schedules

The finite model is not axiomatized: `allEvents` enumerates all distinct ordered
pairs, and sorting that list by a priority gives a complete greedy schedule.
Each unordered pair occurs twice; its second occurrence is ineffectual since
active vertices can only disappear.
-/

noncomputable section

namespace Balloons

variable {α : Type*} [DecidableEq α]

/-- Enumerate all ordered pairs of distinct vertices of the ambient finite set. -/
def allEvents (S : Finset α) : List (Event α) :=
  S.toList.flatMap fun u => S.toList.filterMap fun v =>
    if h : u ≠ v then some ⟨u, v, h⟩ else none

theorem allEvents_covers (S : Finset α) : Covers (allEvents S) S := by
  intro u hu v hv huv
  refine ⟨⟨u, v, huv⟩, ?_, Or.inl ⟨rfl, rfl⟩⟩
  apply List.mem_flatMap.mpr
  refine ⟨u, Finset.mem_toList.mpr hu, ?_⟩
  apply List.mem_filterMap.mpr
  exact ⟨v, Finset.mem_toList.mpr hv, by simp [huv]⟩

omit [DecidableEq α] in
/-- Coverage passes to any initial subset of the ambient configuration. -/
theorem Covers.mono {es : List (Event α)} {A S : Finset α}
    (h : Covers es S) (hAS : A ⊆ S) : Covers es A := by
  intro u hu v hv huv
  exact h u (hAS hu) v (hAS hv) huv

omit [DecidableEq α] in
/-- Reordering a schedule preserves its coverage, though it may change its matching. -/
theorem Covers.perm {es fs : List (Event α)} {A : Finset α}
    (h : Covers es A) (hp : es.Perm fs) : Covers fs A := by
  intro u hu v hv huv
  obtain ⟨e, he, hpair⟩ := h u hu v hv huv
  exact ⟨e, hp.mem_iff.mp he, hpair⟩

/-- Sort all events by a shared comparison; every event remains in the schedule. -/
def orderedEvents (S : Finset α) (le : Event α → Event α → Bool) : List (Event α) :=
  (allEvents S).mergeSort le

theorem orderedEvents_covers (S : Finset α) (le : Event α → Event α → Bool) :
    Covers (orderedEvents S le) S :=
  (allEvents_covers S).perm (List.mergeSort_perm _ _).symm

/-- The finite greedy process for any initial subset and any shared comparison. -/
def process (S A : Finset α) (le : Event α → Event α → Bool) : Finset α :=
  run (orderedEvents S le) A

theorem process_odd_card (S A : Finset α) (le : Event α → Event α → Bool)
    (hAS : A ⊆ S) (hodd : A.card % 2 = 1) : (process S A le).card = 1 :=
  run_card_eq_one _ A ((orderedEvents_covers S le).mono hAS) hodd

theorem process_even_empty (S A : Finset α) (le : Event α → Event α → Bool)
    (hAS : A ⊆ S) (heven : A.card % 2 = 0) : process S A le = ∅ :=
  run_eq_empty _ A ((orderedEvents_covers S le).mono hAS) heven

/-- Every finite prefix gives the same edit bound; no assumption on tie order is needed. -/
theorem prefix_discrepancy_le (es : List (Event α)) (A B : Finset α) (n : ℕ) :
    discrepancy (run (es.take n) A) (run (es.take n) B) ≤ discrepancy A B :=
  run_discrepancy_le _ A B

/-- The discrepancy count is nonincreasing as more events are scanned. -/
theorem prefix_discrepancy_antitone (es : List (Event α)) (A B : Finset α)
    {m n : ℕ} (hmn : m ≤ n) :
    discrepancy (run (es.take n) A) (run (es.take n) B) ≤
      discrepancy (run (es.take m) A) (run (es.take m) B) := by
  have hn : n = m + (n - m) := by omega
  rw [hn, List.take_add]
  exact run_append_discrepancy_le _ _ A B

end Balloons
