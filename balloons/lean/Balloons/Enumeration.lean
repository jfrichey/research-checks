module
public import Balloons.Invariance
public import Balloons.Relabel
public import Mathlib.Data.List.Nodup

@[expose] public section

/-! # Complete pair enumeration and injective relabeling -/

namespace Balloons

variable {α β : Type*} [DecidableEq α] [DecidableEq β]

omit [DecidableEq α] in
@[ext] theorem Event.ext {e f : Event α} (hl : e.left = f.left) (hr : e.right = f.right) :
    e = f := by
  cases e
  cases f
  cases hl
  cases hr
  rfl

omit [DecidableEq α] [DecidableEq β] in
theorem Event.map_injective (f : α ↪ β) : Function.Injective (Event.map f) := by
  intro e g h
  apply Event.ext
  · exact f.injective (congrArg Event.left h)
  · exact f.injective (congrArg Event.right h)

private theorem mem_events_from (S : Finset α) (u : α) (e : Event α) :
    e ∈ S.toList.filterMap (fun v => if h : u ≠ v then some (Event.mk u v h) else none) ↔
      e.left = u ∧ e.right ∈ S := by
  constructor
  · intro h
    obtain ⟨v, hv, he⟩ := List.mem_filterMap.mp h
    by_cases huv : u ≠ v
    · simp only [dite_eq_left huv, Option.some.injEq] at he
      rw [← he]
      exact ⟨rfl, Finset.mem_toList.mp hv⟩
    · simp [huv] at he
  · rintro ⟨hu, hv⟩
    refine List.mem_filterMap.mpr ⟨e.right, Finset.mem_toList.mpr hv, ?_⟩
    have hn : u ≠ e.right := hu ▸ e.distinct
    simp only [dite_eq_left hn, Option.some.injEq]
    exact Event.ext hu.symm rfl

/-- The enumeration contains exactly the ordered distinct pairs of ambient vertices. -/
theorem mem_allEvents (S : Finset α) (e : Event α) :
    e ∈ allEvents S ↔ e.left ∈ S ∧ e.right ∈ S := by
  simp only [allEvents, List.mem_flatMap, mem_events_from, Finset.mem_toList]
  constructor
  · rintro ⟨u, hu, he, hr⟩
    exact ⟨he ▸ hu, hr⟩
  · intro he
    exact ⟨e.left, he.1, rfl, he.2⟩

/-- Each ordered pair is enumerated exactly once. -/
theorem allEvents_nodup (S : Finset α) : (allEvents S).Nodup := by
  apply List.nodup_flatMap.mpr
  constructor
  · intro u _
    apply List.Nodup.filterMap _ (Finset.nodup_toList S)
    intro v w e hv hw
    by_cases huv : u ≠ v
    · by_cases huw : u ≠ w
      · simp only [dite_eq_left huv, Option.mem_some_iff] at hv
        simp only [dite_eq_left huw, Option.mem_some_iff] at hw
        exact congrArg Event.right (hv.trans hw.symm)
      · simp [huw] at hw
    · simp [huv] at hv
  · apply (Finset.nodup_toList S).imp
    intro u v huv e he hf
    have h₁ := (mem_events_from S u e).mp he
    have h₂ := (mem_events_from S v e).mp hf
    exact huv (h₁.1.symm.trans h₂.1)

/-- Injective relabeling preserves the complete event enumeration up to permutation. -/
theorem allEvents_map_perm (f : α ↪ β) (S : Finset α) :
    ((allEvents S).map (Event.map f)).Perm (allEvents (S.map f)) := by
  apply (List.perm_ext_iff_of_nodup
    ((allEvents_nodup S).map (Event.map_injective f)) (allEvents_nodup (S.map f))).mpr
  intro e
  constructor
  · intro h
    obtain ⟨g, hg, rfl⟩ := List.mem_map.mp h
    have hg' := (mem_allEvents S g).mp hg
    exact (mem_allEvents (S.map f) _).mpr
      ⟨Finset.mem_map.mpr ⟨g.left, hg'.1, rfl⟩,
       Finset.mem_map.mpr ⟨g.right, hg'.2, rfl⟩⟩
  · intro h
    obtain ⟨hl, hr⟩ := (mem_allEvents (S.map f) e).mp h
    obtain ⟨u, hu, hfu⟩ := Finset.mem_map.mp hl
    obtain ⟨v, hv, hfv⟩ := Finset.mem_map.mp hr
    have huv : u ≠ v := by
      intro he
      exact e.distinct (hfu.symm.trans ((congrArg f he).trans hfv))
    refine List.mem_map.mpr ⟨⟨u, v, huv⟩, (mem_allEvents S _).mpr ⟨hu, hv⟩, ?_⟩
    exact Event.ext hfu hfv

@[simp] theorem Event.vertices_map (f : α ↪ β) (e : Event α) :
    (e.map f).vertices = e.vertices.map f := by
  simp only [Event.vertices, Event.map, Finset.map_insert, Finset.map_singleton]

/-- The unordered spatial history commutes with injective relabeling. -/
theorem history_map (f : α ↪ β) (es : List (Event α)) (A : Finset α) :
    history (es.map (Event.map f)) (A.map f) =
      (history es A).image (fun edge => edge.map f) := by
  ext edge
  simp only [history, accepted_map, List.mem_toFinset, List.mem_map, Finset.mem_image]
  constructor
  · rintro ⟨_, ⟨e, he, rfl⟩, rfl⟩
    exact ⟨e.vertices, ⟨e, he, rfl⟩, (Event.vertices_map f e).symm⟩
  · rintro ⟨_, ⟨e, he, rfl⟩, rfl⟩
    exact ⟨e.map f, ⟨e, he, rfl⟩, Event.vertices_map f e⟩

end Balloons
