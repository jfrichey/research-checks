module

public import Balloons.InsertionLaw

@[expose] public section

/-! # Comparing label and position noise for one resampling -/

noncomputable section

namespace Balloons

open MeasureTheory
open scoped ENNReal

variable {ι X : Type*} [Fintype ι] [DecidableEq ι]
  [MetricSpace X] [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X]

omit [DecidableEq ι] [BorelSpace X] in
/-- A common iid background and two independent origins form a larger iid sample. -/
theorem two_origins_preserving (μ : Measure X) [IsProbabilityMeasure μ] :
    MeasurePreserving (fun p : (ι → X) × (X × X) =>
      augmentedLocation (augmentedLocation p.1 p.2.1) p.2.2)
      ((Measure.pi fun _ : ι => μ).prod (μ.prod μ))
      (Measure.pi fun _ : Option (Option ι) => μ) :=
  ((augmented_input_preserving (ι := Option ι) μ).comp
    ((augmented_input_preserving (ι := ι) μ).prod (MeasurePreserving.id μ))).comp
      (MeasurePreserving.symm MeasurableEquiv.prodAssoc
        (measurePreserving_prodAssoc (Measure.pi fun _ : ι => μ) μ μ))

omit [DecidableEq ι] in
/-- All old points and both possible new values are distinct almost surely. -/
theorem two_origins_injective_ae (μ : Measure X) [IsProbabilityMeasure μ]
    [NullSingletonClass μ] :
    ∀ᵐ p : (ι → X) × (X × X) ∂(Measure.pi (fun _ : ι => μ)).prod (μ.prod μ),
      Function.Injective (augmentedLocation (augmentedLocation p.1 p.2.1) p.2.2) :=
  (two_origins_preserving μ).quasiMeasurePreserving.ae (iid_sample_injective_ae μ)

omit [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X] in
/-- When all physical locations are distinct, position agreement forces label agreement. -/
theorem insertion_location_eq_implies_label_eq (x : ι → X) (z w : X)
    (hinj : Function.Injective (augmentedLocation (augmentedLocation x z) w))
    (heq : insertionLocation x z = insertionLocation x w) :
    insertionLabel x z = insertionLabel x w := by
  have hright : insertionLocation x w =
      augmentedLocation (augmentedLocation x z) w
        (Option.map Option.some (insertionLabel x w)) := by
    cases h : insertionLabel x w <;> simp [insertionLocation, h, augmentedLocation]
  have hindex : some (insertionLabel x z) = Option.map Option.some (insertionLabel x w) :=
    hinj (heq.trans hright)
  cases hw : insertionLabel x w with
  | none => simp [hw] at hindex
  | some i => simpa only [hw, Option.map_some, Option.some.injEq] using hindex

omit [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X] in
/-- If labels agree but positions differ, the common winning label must be the refreshed one. -/
theorem position_disagree_implies_label_or_inserted (x : ι → X) (z w : X)
    (hpos : insertionLocation x z ≠ insertionLocation x w) :
    insertionLabel x z ≠ insertionLabel x w ∨ insertionLabel x z = none := by
  by_cases hl : insertionLabel x z = insertionLabel x w
  · right
    cases hz : insertionLabel x z with
    | none => rfl
    | some i =>
      exfalso
      apply hpos
      simp only [insertionLocation, ← hl, hz, augmentedLocation]
  · exact Or.inl hl

/-- The original one-coordinate label-disagreement probability for `N=2*m+1`. -/
def survivorLabelNoise (m : ℕ) : ℝ≥0∞ :=
  oneResamplingLaw m {p | insertionLabel p.1 p.2.1 ≠ insertionLabel p.1 p.2.2}

/-- With continuous circle inputs, different surviving labels have different positions. -/
theorem survivorLabelNoise_le_positionNoise (m : ℕ) :
    survivorLabelNoise m ≤ survivorPositionNoise m := by
  have ha : ∀ᵐ p ∂oneResamplingLaw m,
      Function.Injective (augmentedLocation (augmentedLocation p.1 p.2.1) p.2.2) := by
    simpa only [oneResamplingLaw, circleSampleLaw] using
      (two_origins_injective_ae (ι := Fin (2 * m)) circleLaw)
  apply measure_mono_ae
  filter_upwards [ha] with p hp
  exact fun hl hpos => hl (insertion_location_eq_implies_label_eq p.1 p.2.1 p.2.2 hp hpos)

/-- The probability that the first copy is won by its refreshed coordinate is `1/N`. -/
theorem circle_refreshed_label_survives (m : ℕ) :
    oneResamplingLaw m {p | insertionLabel p.1 p.2.1 = none} =
      ((2 * m + 1 : ℕ) : ℝ≥0∞)⁻¹ := by
  simpa only [oneResamplingLaw, circleSampleLaw, Fintype.card_option, Fintype.card_fin] using
    resampling_first_inserted_survives (ι := Fin (2 * m)) circleLaw circle_sphere_null (by simp)

/-- The two finite noise observables differ by at most `1/N`, with a fixed sign. -/
theorem survivorPositionNoise_le_labelNoise_add (m : ℕ) :
    survivorPositionNoise m ≤ survivorLabelNoise m + ((2 * m + 1 : ℕ) : ℝ≥0∞)⁻¹ := by
  have hsub : {p | insertionLocation p.1 p.2.1 ≠ insertionLocation p.1 p.2.2} ⊆
      ({p | insertionLabel p.1 p.2.1 ≠ insertionLabel p.1 p.2.2} ∪
        {p | insertionLabel p.1 p.2.1 = none} :
        Set ((Fin (2 * m) → UnitAddCircle) × (UnitAddCircle × UnitAddCircle))) := by
    intro p hp
    exact position_disagree_implies_label_or_inserted p.1 p.2.1 p.2.2 hp
  calc
    survivorPositionNoise m ≤ oneResamplingLaw m
        ({p | insertionLabel p.1 p.2.1 ≠ insertionLabel p.1 p.2.2} ∪
          {p | insertionLabel p.1 p.2.1 = none}) := measure_mono hsub
    _ ≤ survivorLabelNoise m + oneResamplingLaw m
        {p | insertionLabel p.1 p.2.1 = none} := measure_union_le _ _
    _ = _ := by rw [circle_refreshed_label_survives]

end Balloons
