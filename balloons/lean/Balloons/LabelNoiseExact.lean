module

public import Balloons.LabelNoise

@[expose] public section

/-!
# The exact finite correction from position noise to label noise

For a background `x`, `insertionAvailableMass μ x` is the mass of origins
whose inserted label survives. Independent origins both retain that label
with conditional probability equal to its square. On the almost-sure event
of joint injectivity, this is exactly the event contributing to position
disagreement but not label disagreement.

The circle identities hold for every finite odd size `N = 2*m+1`.
-/

noncomputable section

namespace Balloons

open MeasureTheory
open scoped ENNReal

variable {ι X : Type*} [Fintype ι] [DecidableEq ι] [MetricSpace X]

/-- With distinct joint locations, the extra position-disagreement event is
precisely that both copies are won by the inserted label. -/
theorem position_disagree_iff_label_disagree_or_both_inserted (x : ι → X) (z w : X)
    (hinj : Function.Injective (augmentedLocation (augmentedLocation x z) w)) :
    insertionLocation x z ≠ insertionLocation x w ↔
      insertionLabel x z ≠ insertionLabel x w ∨
        (insertionLabel x z = none ∧ insertionLabel x w = none) := by
  constructor
  · intro hpos
    rcases position_disagree_implies_label_or_inserted x z w hpos with hl | hz
    · exact Or.inl hl
    · by_cases hl : insertionLabel x z = insertionLabel x w
      · exact Or.inr ⟨hz, hl.symm.trans hz⟩
      · exact Or.inl hl
  · rintro (hl | ⟨hz, hw⟩) heq
    · exact hl (insertion_location_eq_implies_label_eq x z w hinj heq)
    · have hzw : z = w := by
        simpa only [insertionLocation, hz, hw, augmentedLocation] using heq
      have hi : (some none : Option (Option ι)) = none := hinj hzw
      cases hi

variable [MeasurableSpace X]

/-- The origin mass for which the transported insertion label is `none`;
this is the available mass `Theta` on an even background. -/
def insertionAvailableMass (μ : Measure X) (x : ι → X) : ℝ≥0∞ :=
  μ {z | insertionLabel x z = none}

/-- For a fixed background, two independent origins both retain their
inserted labels with mass equal to the square of the available mass. -/
theorem both_inserted_eq_available_mass_sq (μ : Measure X) [SFinite μ] (x : ι → X) :
    (μ.prod μ) {p : X × X |
      insertionLabel x p.1 = none ∧ insertionLabel x p.2 = none} =
        insertionAvailableMass μ x ^ 2 := by
  change (μ.prod μ)
    ({z | insertionLabel x z = none} ×ˢ {z | insertionLabel x z = none}) =
      μ {z | insertionLabel x z = none} ^ 2
  rw [Measure.prod_prod, pow_two]

variable [BorelSpace X] [SecondCountableTopology X]

local instance : MeasurableSpace (Option ι) := ⊤
local instance : MeasurableSingletonClass (Option ι) := ⟨fun _ => trivial⟩

/-- The available mass is measurable in all background locations. -/
theorem measurable_insertionAvailableMass (μ : Measure X) [SFinite μ] :
    Measurable (insertionAvailableMass (ι := ι) μ) :=
  measurable_measure_prodMk_left
    (measurableSet_eq_fun measurable_insertionLabel_joint measurable_const)

/-- Product integration over independent origins gives the mean square
available mass for any independent random background. -/
theorem random_background_both_inserted_eq_available_integral
    (ρ : Measure (ι → X)) (μ : Measure X) [SFinite μ] :
    (ρ.prod (μ.prod μ)) {p : (ι → X) × (X × X) |
      insertionLabel p.1 p.2.1 = none ∧ insertionLabel p.1 p.2.2 = none} =
        ∫⁻ x, insertionAvailableMass μ x ^ 2 ∂ρ := by
  have h₁ : Measurable fun p : (ι → X) × (X × X) =>
      insertionLabel p.1 p.2.1 :=
    measurable_insertionLabel_joint.comp
      (measurable_fst.prodMk (measurable_fst.comp measurable_snd))
  have h₂ : Measurable fun p : (ι → X) × (X × X) =>
      insertionLabel p.1 p.2.2 :=
    measurable_insertionLabel_joint.comp
      (measurable_fst.prodMk (measurable_snd.comp measurable_snd))
  have hm : MeasurableSet {p : (ι → X) × (X × X) |
      insertionLabel p.1 p.2.1 = none ∧ insertionLabel p.1 p.2.2 = none} :=
    (measurableSet_eq_fun h₁ measurable_const).inter
      (measurableSet_eq_fun h₂ measurable_const)
  rw [Measure.prod_apply hm]
  apply lintegral_congr
  intro x
  exact both_inserted_eq_available_mass_sq μ x

/-- Exact finite correction: `p_pos = p_lab + E[Theta^2]` for `N = 2*m+1`. -/
theorem survivorPositionNoise_eq_labelNoise_add_available_integral (m : ℕ) :
    survivorPositionNoise m = survivorLabelNoise m +
      ∫⁻ x, insertionAvailableMass circleLaw x ^ 2 ∂circleSampleLaw (2 * m) := by
  let L : Set ((Fin (2 * m) → UnitAddCircle) × (UnitAddCircle × UnitAddCircle)) :=
    {p | insertionLabel p.1 p.2.1 ≠ insertionLabel p.1 p.2.2}
  let I : Set ((Fin (2 * m) → UnitAddCircle) × (UnitAddCircle × UnitAddCircle)) :=
    {p | insertionLabel p.1 p.2.1 = none ∧ insertionLabel p.1 p.2.2 = none}
  have h₁ : Measurable fun p :
      (Fin (2 * m) → UnitAddCircle) × (UnitAddCircle × UnitAddCircle) =>
      insertionLabel p.1 p.2.1 :=
    measurable_insertionLabel_joint.comp
      (measurable_fst.prodMk (measurable_fst.comp measurable_snd))
  have h₂ : Measurable fun p :
      (Fin (2 * m) → UnitAddCircle) × (UnitAddCircle × UnitAddCircle) =>
      insertionLabel p.1 p.2.2 :=
    measurable_insertionLabel_joint.comp
      (measurable_fst.prodMk (measurable_snd.comp measurable_snd))
  have hmI : MeasurableSet I :=
    (measurableSet_eq_fun h₁ measurable_const).inter
      (measurableSet_eq_fun h₂ measurable_const)
  have ha : ∀ᵐ p ∂oneResamplingLaw m,
      Function.Injective (augmentedLocation (augmentedLocation p.1 p.2.1) p.2.2) := by
    simpa only [oneResamplingLaw, circleSampleLaw] using
      (two_origins_injective_ae (ι := Fin (2 * m)) circleLaw)
  have he : {p | insertionLocation p.1 p.2.1 ≠ insertionLocation p.1 p.2.2} =ᵐ[
      oneResamplingLaw m] L ∪ I := by
    filter_upwards [ha] with p hp
    exact propext (position_disagree_iff_label_disagree_or_both_inserted
      p.1 p.2.1 p.2.2 hp)
  have hd : Disjoint L I := by
    apply Set.disjoint_left.mpr
    intro p hp hi
    exact hp (hi.1.trans hi.2.symm)
  calc
    survivorPositionNoise m = oneResamplingLaw m (L ∪ I) := measure_congr he
    _ = survivorLabelNoise m + oneResamplingLaw m I := measure_union hd hmI
    _ = _ := by
      congr 1
      exact random_background_both_inserted_eq_available_integral
        (circleSampleLaw (2 * m)) circleLaw

/-- The label-disagreement probability, old-basin agreement mass, and
available-mass correction sum to one. -/
theorem survivorLabelNoise_add_basin_integral_add_available_integral (m : ℕ) :
    survivorLabelNoise m +
      (∫⁻ x, ∑ a ∈ backgroundLocations x,
        basinMass circleLaw (insertionLocation x) a ^ 2 ∂circleSampleLaw (2 * m)) +
      (∫⁻ x, insertionAvailableMass circleLaw x ^ 2 ∂circleSampleLaw (2 * m)) = 1 := by
  have hb : (∫⁻ x, ∑ a ∈ backgroundLocations x,
      basinMass circleLaw (insertionLocation x) a ^ 2 ∂circleSampleLaw (2 * m)) ≤ 1 := by
    rw [← circle_agreement_eq_basin_integral m]
    exact prob_le_one
  rw [add_right_comm, ← survivorPositionNoise_eq_labelNoise_add_available_integral,
    survivorPositionNoise_eq_basin_integral, tsub_add_cancel_of_le hb]

/-- The exact finite label identity, with subtraction in `ENNReal`:
`p_lab = 1 - E[sum_old basin^2] - E[Theta^2]`. -/
theorem survivorLabelNoise_eq_basin_integral_sub_available_integral (m : ℕ) :
    survivorLabelNoise m = 1 -
      (∫⁻ x, ∑ a ∈ backgroundLocations x,
        basinMass circleLaw (insertionLocation x) a ^ 2 ∂circleSampleLaw (2 * m)) -
      (∫⁻ x, insertionAvailableMass circleLaw x ^ 2 ∂circleSampleLaw (2 * m)) := by
  have hfinite : (∫⁻ x, insertionAvailableMass circleLaw x ^ 2
      ∂circleSampleLaw (2 * m)) ≠ ⊤ := by
    rw [← random_background_both_inserted_eq_available_integral
      (circleSampleLaw (2 * m)) circleLaw]
    exact measure_ne_top _ _
  calc
    survivorLabelNoise m = survivorPositionNoise m -
        (∫⁻ x, insertionAvailableMass circleLaw x ^ 2 ∂circleSampleLaw (2 * m)) :=
      ENNReal.eq_sub_of_add_eq hfinite
        (survivorPositionNoise_eq_labelNoise_add_available_integral m).symm
    _ = _ := by rw [survivorPositionNoise_eq_basin_integral]

end Balloons
