module

public import Balloons.Metric
public import Balloons.SortedRelabel
public import Mathlib.Data.List.Pairwise
public import Mathlib.Data.Set.Finite.Lemmas

@[expose] public section

/-!
# Greedy histories and the unique stable perfect matching

The matching representation is the existing finite set of unordered
two-vertex edges. Stability says that every omitted pair has an endpoint
already matched at no greater distance. For a perfect matching this is
equivalent to the absence of a pair preferred strictly by both endpoints.

The greedy scan satisfies stability even with ties. Distinct distances
between different unordered pairs make the blocking witness strictly
shorter and give uniqueness. Even cardinality supplies perfection.

Informal source: `balloon-process.tex`, Lemma `lem:stable`.
-/

noncomputable section

namespace Balloons

variable {α τ : Type*} [DecidableEq α]

/-- A perfect matching on `A`, in the same unordered-edge representation
as `history`: supported two-point edges, disjoint and covering every vertex. -/
structure PerfectMatchingOn (A : Finset α) (M : Finset (Finset α)) : Prop where
  edge_card : ∀ E ∈ M, E.card = 2
  support : ∀ E ∈ M, E ⊆ A
  disjoint : ∀ E ∈ M, ∀ F ∈ M, E ≠ F → Disjoint E F
  cover : ∀ x ∈ A, ∃ E ∈ M, x ∈ E

private theorem history_event (es : List (Event α)) (A E : Finset α)
    (hE : E ∈ history es A) : ∃ e ∈ accepted es A, e.vertices = E :=
  List.mem_map.mp (List.mem_toFinset.mp hE)

private theorem accepted_edge_mem_history (es : List (Event α)) (A : Finset α)
    (e : Event α) (he : e ∈ accepted es A) : e.vertices ∈ history es A :=
  List.mem_toFinset.mpr (List.mem_map.mpr ⟨e, he, rfl⟩)

private theorem matching_edge_event {A : Finset α} {M : Finset (Finset α)}
    (hM : PerfectMatchingOn A M) (E : Finset α) (hE : E ∈ M) :
    ∃ e ∈ allEvents A, e.vertices = E := by
  obtain ⟨a, b, hab, hE'⟩ := Finset.card_eq_two.mp (hM.edge_card E hE)
  let e : Event α := ⟨a, b, hab⟩
  have hv : e.vertices = E := hE'.symm
  have hs : e.vertices ⊆ A := by
    rw [hv]
    exact hM.support E hE
  exact ⟨e, (mem_allEvents A e).mpr ((event_active_iff e A).mp hs), hv⟩

/-- Every removed vertex belongs to an actual accepted edge. -/
theorem removed_vertex_accepted (es : List (Event α)) (A : Finset α)
    {x : α} (hx : x ∈ A) (hgone : x ∉ run es A) :
    ∃ e ∈ accepted es A, x ∈ e.vertices := by
  induction es generalizing A with
  | nil => exact (hgone hx).elim
  | cons e es ih =>
    by_cases hs : x ∈ step e A
    · obtain ⟨f, hf, hxf⟩ := ih (step e A) hs hgone
      exact ⟨f, (mem_accepted_cons f e es A).mpr (Or.inr hf), hxf⟩
    · have ha : e.left ∈ A ∧ e.right ∈ A := by
        by_contra h
        exact hs (by simpa [step, h] using hx)
      have hv : x ∈ e.vertices := by
        by_contra h
        apply hs
        rw [step_vertices, ite_eq_left ((event_active_iff e A).mpr ha)]
        exact Finset.mem_sdiff.mpr ⟨hx, h⟩
      exact ⟨e, (mem_accepted_cons e e es A).mpr (Or.inl ⟨rfl, ha⟩), hv⟩

/-- The actual history of any complete scan on an even input is a perfect matching. -/
theorem history_perfectMatchingOn (es : List (Event α)) (A : Finset α)
    (hc : Covers es A) (heven : A.card % 2 = 0) :
    PerfectMatchingOn A (history es A) := by
  refine ⟨history_edge_card es A, ?_, ?_, ?_⟩
  · intro E hE
    obtain ⟨e, he, rfl⟩ := history_event es A E hE
    exact (event_active_iff e A).mpr (accepted_endpoints_mem es A e he)
  · have hp : ((accepted es A).map Event.vertices).Pairwise Disjoint := by
      rw [List.pairwise_map]
      apply (accepted_pairwise_disjoint es A).imp
      intro e f h
      apply Finset.disjoint_left.mpr
      intro x hx hy
      simp only [Event.vertices, Finset.mem_insert, Finset.mem_singleton] at hx hy
      unfold EndpointDisjoint at h
      aesop
    intro E hE F hF hne
    exact hp.forall (List.mem_toFinset.mp hE) (List.mem_toFinset.mp hF) hne
  · intro x hx
    have hgone : x ∉ run es A := by
      rw [run_eq_empty es A hc heven]
      simp
    obtain ⟨e, he, hxe⟩ := removed_vertex_accepted es A hx hgone
    exact ⟨e.vertices, accepted_edge_mem_history es A e he, hxe⟩

variable [LinearOrder τ]

/-- Every proposed pair with initially active endpoints has an accepted
incident edge at no greater priority. This records the blocking edge itself. -/
theorem accepted_blocker_le (p : Event α → τ) (es : List (Event α)) (A : Finset α)
    (hsorted : es.Pairwise fun e f => p e ≤ p f)
    (e : Event α) (he : e ∈ es) (hA : e.vertices ⊆ A) :
    ∃ f ∈ accepted es A,
      (e.left ∈ f.vertices ∨ e.right ∈ f.vertices) ∧ p f ≤ p e := by
  induction es generalizing A with
  | nil => simp at he
  | cons g gs ih =>
    obtain ⟨hhead, htail⟩ := List.pairwise_cons.mp hsorted
    rcases List.mem_cons.mp he with rfl | he
    · have hg := (event_active_iff e A).mp hA
      exact ⟨e, (mem_accepted_cons e e gs A).mpr (Or.inl ⟨rfl, hg⟩),
        Or.inl (by simp [Event.vertices]), le_rfl⟩
    · by_cases hg : g.left ∈ A ∧ g.right ∈ A
      · by_cases ht : e.left ∈ g.vertices ∨ e.right ∈ g.vertices
        · exact ⟨g, (mem_accepted_cons g g gs A).mpr (Or.inl ⟨rfl, hg⟩),
            ht, hhead e he⟩
        · have hs : e.vertices ⊆ step g A := by
            intro x hx
            have hn : x ∉ g.vertices := by
              intro hxg
              simp only [Event.vertices, Finset.mem_insert, Finset.mem_singleton] at hx
              rcases hx with rfl | rfl
              · exact ht (Or.inl hxg)
              · exact ht (Or.inr hxg)
            rw [step_vertices, ite_eq_left ((event_active_iff g A).mpr hg)]
            exact Finset.mem_sdiff.mpr ⟨hA hx, hn⟩
          obtain ⟨f, hf, hb⟩ := ih (step g A) htail he hs
          exact ⟨f, (mem_accepted_cons f g gs A).mpr (Or.inr hf), hb⟩
      · have hs : step g A = A := by simp [step, hg]
        obtain ⟨f, hf, hb⟩ := ih A htail he hA
        exact ⟨f, (mem_accepted_cons f g gs A).mpr (Or.inr (hs.symm ▸ hf)), hb⟩

variable [MetricSpace α]

/-- Distance ties among pairs of `A` can only reverse the same unordered edge. -/
def DistinctPairDistances (A : Finset α) : Prop :=
  ∀ e ∈ allEvents A, ∀ f ∈ allEvents A,
    edgeLength e = edgeLength f → e.vertices = f.vertices

/-- On a perfect matching, this is the no-blocking-pair condition:
every omitted pair has an endpoint matched at no greater distance. -/
def StableOn (A : Finset α) (M : Finset (Finset α)) : Prop :=
  ∀ e ∈ allEvents A, e.vertices ∉ M →
    ∃ f : Event α, f.vertices ∈ M ∧
      (e.left ∈ f.vertices ∨ e.right ∈ f.vertices) ∧ edgeLength f ≤ edgeLength e

/-- Greedy histories have no blocking pair, with no parity or no-tie assumption. -/
theorem metric_history_stable (A : Finset α) :
    StableOn A (history (metricSchedule A) A) := by
  intro e he _
  have hA := (event_active_iff e A).mpr ((mem_allEvents A e).mp he)
  have he' : e ∈ metricSchedule A :=
    (mem_timedEvents A edgeLength e).mpr ((mem_allEvents A e).mp he)
  obtain ⟨f, hf, hb⟩ := accepted_blocker_le edgeLength (metricSchedule A) A
    (metricSchedule_sorted A) e he' hA
  exact ⟨f, accepted_edge_mem_history _ A f hf, hb⟩

/-- Under distinct pair distances, every unchosen unordered pair is blocked
by an actual accepted edge of strictly smaller length. -/
theorem rejected_metric_edge_has_shorter_blocker (A : Finset α)
    (htie : DistinctPairDistances A) (e : Event α) (he : e ∈ allEvents A)
    (hrejected : e.vertices ∉ history (metricSchedule A) A) :
    ∃ f ∈ accepted (metricSchedule A) A,
      (e.left ∈ f.vertices ∨ e.right ∈ f.vertices) ∧ edgeLength f < edgeLength e := by
  have he' : e ∈ metricSchedule A :=
    (mem_timedEvents A edgeLength e).mpr ((mem_allEvents A e).mp he)
  obtain ⟨f, hf, hshare, hle⟩ := accepted_blocker_le edgeLength (metricSchedule A) A
    (metricSchedule_sorted A) e he' ((event_active_iff e A).mpr ((mem_allEvents A e).mp he))
  have hfA : f ∈ allEvents A :=
    (mem_allEvents A f).mpr (accepted_endpoints_mem _ A f hf)
  refine ⟨f, hf, hshare, lt_of_le_of_ne hle ?_⟩
  intro hlen
  have hv := htie f hfA e he hlen
  exact hrejected (hv ▸ accepted_edge_mem_history _ A f hf)

/-- Two stable perfect matchings with distinct unordered pair distances agree.
The shortest edge in their symmetric difference would have a still shorter
incident edge in the opposite matching, a contradiction. -/
theorem stablePerfectMatching_unique (A : Finset α) (htie : DistinctPairDistances A)
    (M N : Finset (Finset α)) (hM : PerfectMatchingOn A M) (hN : PerfectMatchingOn A N)
    (hsM : StableOn A M) (hsN : StableOn A N) : M = N := by
  classical
  have descend (P Q : Finset (Finset α)) (hP : PerfectMatchingOn A P)
      (hQ : PerfectMatchingOn A Q) (hsQ : StableOn A Q)
      (e : Event α) (heA : e ∈ allEvents A) (heP : e.vertices ∈ P)
      (heQ : e.vertices ∉ Q) :
      ∃ f ∈ allEvents A, f.vertices ∈ Q ∧ f.vertices ∉ P ∧ edgeLength f < edgeLength e := by
    obtain ⟨f, hfQ, hshare, hle⟩ := hsQ e heA heQ
    have hfA : f ∈ allEvents A :=
      (mem_allEvents A f).mpr ((event_active_iff f A).mp (hQ.support _ hfQ))
    have hne : f.vertices ≠ e.vertices := fun h => heQ (h ▸ hfQ)
    have hlt : edgeLength f < edgeLength e :=
      lt_of_le_of_ne hle (fun h => hne (htie f hfA e heA h))
    have hfP : f.vertices ∉ P := by
      intro hf
      have hd := Finset.disjoint_left.mp (hP.disjoint _ heP _ hf hne.symm)
      rcases hshare with hl | hr
      · exact hd (by simp [Event.vertices]) hl
      · exact hd (by simp [Event.vertices]) hr
    exact ⟨f, hfA, hfQ, hfP, hlt⟩
  let D : Finset (Event α) := (allEvents A).toFinset.filter fun e =>
    (e.vertices ∈ M ∧ e.vertices ∉ N) ∨ (e.vertices ∈ N ∧ e.vertices ∉ M)
  have hD : D = ∅ := by
    by_contra hne
    obtain ⟨e, he, hmin⟩ := D.exists_min_image edgeLength
      (Finset.nonempty_iff_ne_empty.mpr hne)
    have heA : e ∈ allEvents A := List.mem_toFinset.mp (Finset.mem_filter.mp he).1
    rcases (Finset.mem_filter.mp he).2 with ⟨heM, heN⟩ | ⟨heN, heM⟩
    · obtain ⟨f, hfA, hfN, hfM, hlt⟩ := descend M N hM hN hsN e heA heM heN
      have hf : f ∈ D := Finset.mem_filter.mpr
        ⟨List.mem_toFinset.mpr hfA, Or.inr ⟨hfN, hfM⟩⟩
      exact (not_lt_of_ge (hmin f hf)) hlt
    · obtain ⟨f, hfA, hfM, hfN, hlt⟩ := descend N M hN hM hsM e heA heN heM
      have hf : f ∈ D := Finset.mem_filter.mpr
        ⟨List.mem_toFinset.mpr hfA, Or.inl ⟨hfM, hfN⟩⟩
      exact (not_lt_of_ge (hmin f hf)) hlt
  apply Finset.ext
  intro E
  constructor
  · intro hE
    by_contra hn
    obtain ⟨e, he, hv⟩ := matching_edge_event hM E hE
    have heD : e ∈ D := Finset.mem_filter.mpr
      ⟨List.mem_toFinset.mpr he, Or.inl ⟨hv.symm ▸ hE, hv.symm ▸ hn⟩⟩
    simp [hD] at heD
  · intro hE
    by_contra hn
    obtain ⟨e, he, hv⟩ := matching_edge_event hN E hE
    have heD : e ∈ D := Finset.mem_filter.mpr
      ⟨List.mem_toFinset.mpr he, Or.inr ⟨hv.symm ▸ hE, hv.symm ▸ hn⟩⟩
    simp [hD] at heD

/-- Stability characterizes the greedy perfect matching on an even finite
metric input with distinct distances between different unordered pairs. -/
theorem stablePerfectMatching_iff_eq_history (A : Finset α)
    (heven : A.card % 2 = 0) (htie : DistinctPairDistances A)
    (M : Finset (Finset α)) :
    (PerfectMatchingOn A M ∧ StableOn A M) ↔ M = history (metricSchedule A) A := by
  have hp := history_perfectMatchingOn (metricSchedule A) A (metricSchedule_covers A) heven
  have hs := metric_history_stable A
  constructor
  · rintro ⟨hM, hsM⟩
    exact stablePerfectMatching_unique A htie M _ hM hp hsM hs
  · rintro rfl
    exact ⟨hp, hs⟩

end Balloons
