module

public import Balloons.Terminal
public import Balloons.Relabel

@[expose] public section

/-!
# Accepted collision histories

The history records the actual accepted events of the finite scan. Its
cardinality, endpoint exclusion, and matching property are proved directly
from that scan, not postulated as properties of a separate object.
-/

namespace Balloons

variable {α β : Type*} [DecidableEq α] [DecidableEq β]

/-- Accepted events, in the order in which they are scanned. -/
def accepted : List (Event α) → Finset α → List (Event α)
  | [], _ => []
  | e :: es, A =>
    if e.left ∈ A ∧ e.right ∈ A then e :: accepted es (step e A)
    else accepted es (step e A)

@[simp] theorem accepted_nil (A : Finset α) : accepted [] A = [] := rfl

theorem mem_accepted_cons (e f : Event α) (es : List (Event α)) (A : Finset α) :
    e ∈ accepted (f :: es) A ↔
      (e = f ∧ f.left ∈ A ∧ f.right ∈ A) ∨ e ∈ accepted es (step f A) := by
  by_cases h : f.left ∈ A ∧ f.right ∈ A
  · simp [accepted, h]
  · simp only [accepted, ite_eq_right h]
    tauto

/-- Every accepted event was in the proposed schedule. -/
theorem accepted_mem_schedule (es : List (Event α)) (A : Finset α)
    (e : Event α) (he : e ∈ accepted es A) : e ∈ es := by
  induction es generalizing A with
  | nil => simp at he
  | cons f fs ih =>
    rcases (mem_accepted_cons e f fs A).mp he with ⟨rfl, _⟩ | ht
    · exact List.mem_cons_self
    · exact List.mem_cons_of_mem f (ih _ ht)

/-- Both endpoints of an accepted edge were initially active. -/
theorem accepted_endpoints_mem (es : List (Event α)) (A : Finset α)
    (e : Event α) (he : e ∈ accepted es A) : e.left ∈ A ∧ e.right ∈ A := by
  induction es generalizing A with
  | nil => simp at he
  | cons f fs ih =>
    rcases (mem_accepted_cons e f fs A).mp he with ⟨rfl, h⟩ | ht
    · exact h
    · have h := ih _ ht
      exact ⟨step_subset f A h.1, step_subset f A h.2⟩

/-- Accepted endpoints are absent from the remaining active configuration. -/
theorem accepted_endpoints_not_active (es : List (Event α)) (A : Finset α)
    (e : Event α) (he : e ∈ accepted es A) :
    e.left ∉ run es A ∧ e.right ∉ run es A := by
  induction es generalizing A with
  | nil => simp at he
  | cons f fs ih =>
    rcases (mem_accepted_cons e f fs A).mp he with ⟨rfl, h⟩ | ht
    · have hl : e.left ∉ step e A := by simp [step, h]
      have hr : e.right ∉ step e A := by simp [step, h]
      exact ⟨fun hx => hl (run_subset fs _ hx), fun hx => hr (run_subset fs _ hx)⟩
    · exact ih _ ht

/-- Two vertices are removed per accepted edge, with no lost or duplicated vertex. -/
theorem accepted_length_balance (es : List (Event α)) (A : Finset α) :
    (run es A).card + 2 * (accepted es A).length = A.card := by
  induction es generalizing A with
  | nil => simp [accepted]
  | cons e es ih =>
    have ht := ih (step e A)
    by_cases h : e.left ∈ A ∧ e.right ∈ A
    · have hs := step_card_of_accept e A h
      simp only [run_cons, accepted, ite_eq_left h, List.length_cons]
      omega
    · simpa [run_cons, accepted, h, step] using ht

/-- Accepted edges of an odd complete input number exactly `(N-1)/2`. -/
theorem accepted_length_of_odd_complete (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) (hodd : A.card % 2 = 1) :
    (accepted es A).length = (A.card - 1) / 2 := by
  have h := accepted_length_balance es A
  rw [run_card_eq_one es A hc hodd] at h
  omega

/-- Accepted edges of an even complete input number exactly `N/2`. -/
theorem accepted_length_of_even_complete (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) (heven : A.card % 2 = 0) :
    (accepted es A).length = A.card / 2 := by
  have h := accepted_length_balance es A
  rw [run_eq_empty es A hc heven] at h
  simp only [Finset.card_empty] at h
  omega

/-- The endpoint sets of two proposed edges are disjoint. -/
def EndpointDisjoint (e f : Event α) : Prop :=
  e.left ≠ f.left ∧ e.left ≠ f.right ∧ e.right ≠ f.left ∧ e.right ≠ f.right

/-- No two accepted edges share an endpoint: the accepted history is a matching. -/
theorem accepted_pairwise_disjoint (es : List (Event α)) (A : Finset α) :
    (accepted es A).Pairwise EndpointDisjoint := by
  induction es generalizing A with
  | nil => simp [accepted]
  | cons e es ih =>
    by_cases h : e.left ∈ A ∧ e.right ∈ A
    · simp only [accepted, ite_eq_left h, List.pairwise_cons]
      refine ⟨?_, ih (step e A)⟩
      intro f hf
      have hp := accepted_endpoints_mem es (step e A) f hf
      simp only [step, ite_eq_left h, Finset.mem_erase] at hp
      unfold EndpointDisjoint
      aesop
    · simpa only [accepted, ite_eq_right h] using ih (step e A)

/-- Splitting the scan splits its history at the corresponding remaining state. -/
theorem accepted_append (es fs : List (Event α)) (A : Finset α) :
    accepted (es ++ fs) A = accepted es A ++ accepted fs (run es A) := by
  induction es generalizing A with
  | nil => simp [accepted]
  | cons e es ih =>
    by_cases h : e.left ∈ A ∧ e.right ∈ A <;>
      simp [accepted, run_cons, h, ih]

/-- Relabeling transports the actual accepted history. -/
theorem accepted_map (f : α ↪ β) (es : List (Event α)) (A : Finset α) :
    accepted (es.map (Event.map f)) (A.map f) = (accepted es A).map (Event.map f) := by
  induction es generalizing A with
  | nil => rfl
  | cons e es ih =>
    simp only [List.map_cons, accepted, Event.map, Finset.mem_map']
    change (if e.left ∈ A ∧ e.right ∈ A then
      e.map f :: accepted (es.map (Event.map f)) (step (e.map f) (A.map f))
      else accepted (es.map (Event.map f)) (step (e.map f) (A.map f))) = _
    rw [step_map, ih]
    split <;> simp_all

end Balloons
