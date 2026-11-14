module

public import Balloons.Measurable
public import Balloons.Basins

@[expose] public section

/-!
# A measurable insertion map for actual metric input

Old points have labels `some i` and the inserted origin has label `none`.
Every pair is sorted by its actual spatial distance. On an even background,
the transported label is the unique survivor. Its spatial location is a
measurable function of the inserted origin, so the basin formula applies
without a separate measurability hypothesis.
-/

noncomputable section

namespace Balloons

open MeasureTheory
open scoped ENNReal

variable {ι X : Type*} [Fintype ι] [DecidableEq ι]

/-- The old labels in an input with one distinguished insertion slot. -/
def backgroundLabels : Finset (Option ι) := Finset.univ.erase none

@[simp] theorem backgroundLabels_card : (backgroundLabels (ι := ι)).card = Fintype.card ι := by
  simp [backgroundLabels, Fintype.card_option]

@[simp] theorem insert_none_backgroundLabels :
    insert none (backgroundLabels (ι := ι)) = Finset.univ := by
  simp [backgroundLabels]

/-- Attach the fixed old locations and the new origin to their labels. -/
def augmentedLocation (x : ι → X) (z : X) : Option ι → X
  | none => z
  | some i => x i

variable [MetricSpace X]

/-- A complete schedule using actual pair distances of the augmented input. -/
def insertionEvents (x : ι → X) (z : X) : List (Event (Option ι)) :=
  timedEvents Finset.univ (fun e =>
    dist (augmentedLocation x z e.left) (augmentedLocation x z e.right))

/-- The transported insertion label. -/
def insertionLabel (x : ι → X) (z : X) : Option ι :=
  follow (insertionEvents x z) backgroundLabels none

/-- The spatial output of the finite insertion map. -/
def insertionLocation (x : ι → X) (z : X) : X :=
  augmentedLocation x z (insertionLabel x z)

/-- For an even background this map gives the unique final survivor label. -/
theorem run_insertion_eq_singleton (x : ι → X) (z : X)
    (heven : Fintype.card ι % 2 = 0) :
    run (insertionEvents x z) Finset.univ = {insertionLabel x z} := by
  have hc : Covers (insertionEvents x z) (backgroundLabels (ι := ι)) :=
    (timedEvents_covers Finset.univ _).mono (Finset.subset_univ _)
  simpa only [insert_none_backgroundLabels, insertionLabel] using
    inserted_survivor_eq_follow (insertionEvents x z) backgroundLabels none hc
      (by simpa using heven) (by simp [backgroundLabels])

/-- The finite set of distinct old spatial locations. -/
def backgroundLocations (x : ι → X) : Finset X := by
  classical
  exact Finset.univ.image x

/-- The insertion map outputs either its unchanged origin or an old location. -/
theorem insertionLocation_eq_self_or_mem (x : ι → X) (z : X) :
    insertionLocation x z = z ∨ insertionLocation x z ∈ backgroundLocations x := by
  classical
  cases h : insertionLabel x z with
  | none => exact Or.inl (by simp [insertionLocation, h, augmentedLocation])
  | some i =>
    right
    simp only [insertionLocation, h, augmentedLocation, backgroundLocations]
    exact Finset.mem_image.mpr ⟨i, Finset.mem_univ _, rfl⟩

variable [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X]

local instance : MeasurableSpace (Option ι) := ⊤
local instance : MeasurableSingletonClass (Option ι) := ⟨fun _ => trivial⟩

omit [Fintype ι] [DecidableEq ι] [MetricSpace X] [BorelSpace X]
    [SecondCountableTopology X] in
/-- Every fixed label has a measurable location as the origin varies. -/
theorem measurable_augmentedLocation (x : ι → X) (i : Option ι) :
    Measurable fun z => augmentedLocation x z i := by
  cases i
  · exact measurable_id
  · exact measurable_const

/-- Sorting the actual distance comparisons produces a measurable insertion label. -/
theorem measurable_insertionLabel (x : ι → X) : Measurable (insertionLabel x) :=
  measurable_priority_observable
    (fun z e => dist (augmentedLocation x z e.left) (augmentedLocation x z e.right))
    (fun _ => 0)
    (fun e => (measurable_augmentedLocation x e.left).dist
      (measurable_augmentedLocation x e.right)) measurable_const
    (fun d => follow (orderedEvents Finset.univ d.1) backgroundLabels none)

/-- The finite metric insertion map is measurable as a spatial map. -/
theorem measurable_insertionLocation (x : ι → X) : Measurable (insertionLocation x) := by
  have hm : Measurable fun p : X × Option ι => augmentedLocation x p.1 p.2 :=
    measurable_from_prod_countable_left (measurable_augmentedLocation x)
  exact hm.comp (measurable_id.prodMk (measurable_insertionLabel x))

omit [Fintype ι] [DecidableEq ι] [MetricSpace X] [BorelSpace X]
    [SecondCountableTopology X] in
/-- Each location is jointly measurable in the old sample and the new origin. -/
theorem measurable_augmentedLocation_joint (i : Option ι) :
    Measurable fun p : (ι → X) × X => augmentedLocation p.1 p.2 i := by
  cases i with
  | none => exact measurable_snd
  | some i => exact (measurable_pi_apply i).comp measurable_fst

/-- The transported label is jointly measurable in all input locations. -/
theorem measurable_insertionLabel_joint :
    Measurable fun p : (ι → X) × X => insertionLabel p.1 p.2 :=
  measurable_priority_observable
    (fun (p : (ι → X) × X) e => dist (augmentedLocation p.1 p.2 e.left)
      (augmentedLocation p.1 p.2 e.right)) (fun _ => 0)
    (fun e => (measurable_augmentedLocation_joint e.left).dist
      (measurable_augmentedLocation_joint e.right)) measurable_const
    (fun d => follow (orderedEvents Finset.univ d.1) backgroundLabels none)

/-- The spatial insertion map is jointly measurable in all input locations. -/
theorem measurable_insertionLocation_joint :
    Measurable fun p : (ι → X) × X => insertionLocation p.1 p.2 := by
  have hm : Measurable fun q : ((ι → X) × X) × Option ι =>
      augmentedLocation q.1.1 q.1.2 q.2 :=
    measurable_from_prod_countable_left measurable_augmentedLocation_joint
  exact hm.comp (measurable_id.prodMk measurable_insertionLabel_joint)

/-- The finite metric insertion map satisfies the basin identity with no
unproved measurability premise. For even backgrounds it is survivor agreement. -/
theorem insertionLocation_agreement_eq_basin_sq (μ : Measure X) [SFinite μ]
    [NullSingletonClass μ] (x : ι → X) :
    (μ.prod μ) {p : X × X | insertionLocation x p.1 = insertionLocation x p.2} =
      ∑ a ∈ backgroundLocations x, basinMass μ (insertionLocation x) a ^ 2 :=
  agreement_eq_sum_basin_sq μ (backgroundLocations x) (insertionLocation x)
    (measurable_insertionLocation x) (insertionLocation_eq_self_or_mem x)

/-- Averaging over any independent random background gives the exact finite
agreement identity as an expectation of squared basin masses. -/
theorem random_background_agreement_eq_basin_integral
    (ρ : Measure (ι → X)) (μ : Measure X) [SFinite μ] [NullSingletonClass μ] :
    (ρ.prod (μ.prod μ)) {p : (ι → X) × (X × X) |
      insertionLocation p.1 p.2.1 = insertionLocation p.1 p.2.2} =
      ∫⁻ x, ∑ a ∈ backgroundLocations x, basinMass μ (insertionLocation x) a ^ 2 ∂ρ := by
  have h₁ : Measurable fun p : (ι → X) × (X × X) =>
      insertionLocation p.1 p.2.1 :=
    measurable_insertionLocation_joint.comp
      (measurable_fst.prodMk (measurable_fst.comp measurable_snd))
  have h₂ : Measurable fun p : (ι → X) × (X × X) =>
      insertionLocation p.1 p.2.2 :=
    measurable_insertionLocation_joint.comp
      (measurable_fst.prodMk (measurable_snd.comp measurable_snd))
  rw [Measure.prod_apply (measurableSet_eq_fun h₁ h₂)]
  apply lintegral_congr
  intro x
  exact insertionLocation_agreement_eq_basin_sq μ x

/-- The disagreement probability is the complement of the basin-square expectation. -/
theorem random_background_disagreement_eq_basin_integral
    (ρ : Measure (ι → X)) (μ : Measure X) [IsProbabilityMeasure ρ]
    [IsProbabilityMeasure μ] [NullSingletonClass μ] :
    (ρ.prod (μ.prod μ)) {p : (ι → X) × (X × X) |
      insertionLocation p.1 p.2.1 ≠ insertionLocation p.1 p.2.2} = 1 -
      ∫⁻ x, ∑ a ∈ backgroundLocations x, basinMass μ (insertionLocation x) a ^ 2 ∂ρ := by
  have h₁ : Measurable fun p : (ι → X) × (X × X) =>
      insertionLocation p.1 p.2.1 :=
    measurable_insertionLocation_joint.comp
      (measurable_fst.prodMk (measurable_fst.comp measurable_snd))
  have h₂ : Measurable fun p : (ι → X) × (X × X) =>
      insertionLocation p.1 p.2.2 :=
    measurable_insertionLocation_joint.comp
      (measurable_fst.prodMk (measurable_snd.comp measurable_snd))
  have hc := measure_compl (μ := ρ.prod (μ.prod μ)) (measurableSet_eq_fun h₁ h₂)
    (measure_ne_top _ _)
  rw [measure_univ, random_background_agreement_eq_basin_integral ρ μ] at hc
  exact hc

end Balloons
