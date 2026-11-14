module

public import Balloons.UniformWinner

@[expose] public section

/-! # The iid law of an inserted point and its surviving label -/

noncomputable section

namespace Balloons

open MeasureTheory
open scoped ENNReal

variable {ι X : Type*} [Fintype ι] [DecidableEq ι]
  [MetricSpace X] [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X]

omit [DecidableEq ι] in
private theorem option_card_odd (heven : Fintype.card ι % 2 = 0) :
    Fintype.card (Option ι) % 2 = 1 := by
  rw [Fintype.card_option]
  omega

omit [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X] in
/-- The complete-process survivor agrees with the transported insertion label. -/
theorem winner_augmented_eq_insertionLabel (heven : Fintype.card ι % 2 = 0)
    (x : ι → X) (z : X) :
    winner (by rw [Fintype.card_option]; omega) (augmentedLocation x z) = insertionLabel x z := by
  apply Finset.singleton_inj.mp
  exact (labeledSurvivors_eq_singleton (option_card_odd heven) (augmentedLocation x z)).symm.trans
    (run_insertion_eq_singleton x z heven)

omit [DecidableEq ι] [BorelSpace X] in
/-- Adding one independent point to an iid background gives the full augmented iid law. -/
theorem augmented_input_preserving (μ : Measure X) [IsProbabilityMeasure μ] :
    MeasurePreserving (fun p : (ι → X) × X => augmentedLocation p.1 p.2)
      ((Measure.pi fun _ : ι => μ).prod μ) (Measure.pi fun _ : Option ι => μ) := by
  have he : (fun p : (ι → X) × X => augmentedLocation p.1 p.2) =
      (MeasurableEquiv.piOptionEquivProd (fun _ : Option ι => X)).symm := by
    funext p i
    cases i <;> rfl
  rw [he]
  exact ⟨(MeasurableEquiv.piOptionEquivProd (fun _ : Option ι => X)).symm.measurable,
    Measure.pi_map_piOptionEquivProd (fun _ : Option ι => μ)⟩

local instance : MeasurableSpace (Option ι) := ⊤
local instance : MeasurableSingletonClass (Option ι) := ⟨fun _ => trivial⟩

/-- The inserted label, like every old label, survives with probability exactly `1/(N+1)`. -/
theorem insertion_label_uniform (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0)
    (heven : Fintype.card ι % 2 = 0) (i : Option ι) :
    ((Measure.pi fun _ : ι => μ).prod μ) {p | insertionLabel p.1 p.2 = i} =
      (Fintype.card (Option ι) : ℝ≥0∞)⁻¹ := by
  have hm : MeasurableSet {x : Option ι → X | winner (option_card_odd heven) x = i} :=
    measurableSet_eq_fun (measurable_winner _) measurable_const
  calc
    _ = ((Measure.pi fun _ : ι => μ).prod μ)
        ((fun p : (ι → X) × X => augmentedLocation p.1 p.2) ⁻¹'
          {x | winner (option_card_odd heven) x = i}) := by
      congr 1
      ext p
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, winner_augmented_eq_insertionLabel heven]
    _ = _ := ((augmented_input_preserving (ι := ι) μ).measure_preimage hm.nullMeasurableSet).trans
      (iid_winner_uniform μ hsphere (option_card_odd heven) i)

omit [DecidableEq ι] [BorelSpace X] in
/-- Keeping either value of the refreshed coordinate gives its usual augmented iid law. -/
theorem resampling_first_input_preserving (μ : Measure X) [IsProbabilityMeasure μ] :
    MeasurePreserving (fun p : (ι → X) × (X × X) => (p.1, p.2.1))
      ((Measure.pi fun _ : ι => μ).prod (μ.prod μ)) ((Measure.pi fun _ : ι => μ).prod μ) :=
  (MeasurePreserving.id _).prod measurePreserving_fst

/-- The probability that the first copy is won by the refreshed label is its reciprocal size. -/
theorem resampling_first_inserted_survives (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0)
    (heven : Fintype.card ι % 2 = 0) :
    ((Measure.pi fun _ : ι => μ).prod (μ.prod μ))
      {p | insertionLabel p.1 p.2.1 = none} = (Fintype.card (Option ι) : ℝ≥0∞)⁻¹ := by
  have hm : MeasurableSet {p : (ι → X) × X | insertionLabel p.1 p.2 = none} :=
    measurableSet_eq_fun measurable_insertionLabel_joint measurable_const
  exact ((resampling_first_input_preserving (ι := ι) μ).measure_preimage hm.nullMeasurableSet).trans
    (insertion_label_uniform μ hsphere heven none)

end Balloons
