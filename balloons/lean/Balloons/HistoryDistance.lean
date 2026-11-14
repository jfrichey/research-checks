module

public import Balloons.History

@[expose] public section

/-!
# Unordered spatial edges and history distance

Edges are two-element finite sets of physical vertices. Thus reversing an
event does not create a different spatial edge. The exact relation between
history distance and shared edges is proved from the finite process.
-/

namespace Balloons

open scoped symmDiff

variable {α : Type*} [DecidableEq α]

/-- The unordered spatial edge underlying an event. -/
def Event.vertices (e : Event α) : Finset α := {e.left, e.right}

@[simp] theorem Event.vertices_card (e : Event α) : e.vertices.card = 2 := by
  simp [Event.vertices, e.distinct]

/-- The accepted unordered edges form a set with no repeated edge. -/
theorem accepted_vertices_nodup (es : List (Event α)) (A : Finset α) :
    ((accepted es A).map Event.vertices).Nodup := by
  rw [List.nodup_iff_pairwise_ne, List.pairwise_map]
  apply (accepted_pairwise_disjoint es A).imp
  intro e f hdis heq
  have hm : e.left ∈ f.vertices := by
    rw [← heq]
    simp [Event.vertices]
  simp only [Event.vertices, Finset.mem_insert, Finset.mem_singleton] at hm
  exact hm.elim hdis.1 hdis.2.1

/-- The set of actual accepted unordered spatial edges. -/
def history (es : List (Event α)) (A : Finset α) : Finset (Finset α) :=
  ((accepted es A).map Event.vertices).toFinset

theorem history_card (es : List (Event α)) (A : Finset α) :
    (history es A).card = (accepted es A).length := by
  unfold history
  rw [List.toFinset_card_of_nodup (accepted_vertices_nodup es A), List.length_map]

/-- Every spatial edge of the history has exactly two endpoints. -/
theorem history_edge_card (es : List (Event α)) (A : Finset α) (edge : Finset α)
    (he : edge ∈ history es A) : edge.card = 2 := by
  obtain ⟨e, _, rfl⟩ := List.mem_map.mp (List.mem_toFinset.mp he)
  exact e.vertices_card

/-- Exact vertex accounting for the unordered spatial matching. -/
theorem history_card_balance (es : List (Event α)) (A : Finset α) :
    (run es A).card + 2 * (history es A).card = A.card := by
  rw [history_card]
  exact accepted_length_balance es A

/-- Symmetric-difference distance between two spatial histories. -/
def historyDistance (es : List (Event α)) (A B : Finset α) : ℕ :=
  (history es A ∆ history es B).card

/-- Number of spatial edges accepted in both copies. -/
def sharedEdges (es : List (Event α)) (A B : Finset α) : ℕ :=
  (history es A ∩ history es B).card

private theorem card_symmDiff_balance {β : Type*} [DecidableEq β] (S T : Finset β) :
    (S ∆ T).card + 2 * (S ∩ T).card = S.card + T.card := by
  have hd : Disjoint (S \ T) (T \ S) := by
    rw [Finset.disjoint_left]
    intro x hx hy
    simp only [Finset.mem_sdiff] at hx hy
    exact hx.2 hy.1
  rw [Finset.symmDiff_def, Finset.card_union_of_disjoint hd]
  have hs := Finset.card_sdiff_add_card_inter S T
  have ht := Finset.card_sdiff_add_card_inter T S
  rw [Finset.inter_comm T S] at ht
  omega

/-- Every common edge contributes twice to the sum of the two history sizes. -/
theorem historyDistance_balance (es : List (Event α)) (A B : Finset α) :
    historyDistance es A B + 2 * sharedEdges es A B =
      (history es A).card + (history es B).card :=
  card_symmDiff_balance _ _

/-- The exact overlap identity for equal-size complete odd inputs. -/
theorem historyDistance_odd_balance (es : List (Event α)) (A B : Finset α)
    (ha : Covers es A) (hb : Covers es B) (heq : A.card = B.card)
    (hodd : A.card % 2 = 1) :
    historyDistance es A B + 2 * sharedEdges es A B = A.card - 1 := by
  have hbodd : B.card % 2 = 1 := by omega
  have h₁ := history_card_balance es A
  have h₂ := history_card_balance es B
  rw [run_card_eq_one es A ha hodd] at h₁
  rw [run_card_eq_one es B hb hbodd] at h₂
  have hd := historyDistance_balance es A B
  omega

/-- In particular, `H = N - 1 - 2 J`, the finite identity at complete refresh. -/
theorem historyDistance_eq_sub_shared (es : List (Event α)) (A B : Finset α)
    (ha : Covers es A) (hb : Covers es B) (heq : A.card = B.card)
    (hodd : A.card % 2 = 1) :
    historyDistance es A B = A.card - 1 - 2 * sharedEdges es A B := by
  have h := historyDistance_odd_balance es A B ha hb heq hodd
  omega

end Balloons
