module

public import Balloons.Insertion
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Periodic
public import Mathlib.MeasureTheory.Constructions.Pi

@[expose] public section

/-!
# Uniform circle inputs and the finite survivor-noise identity

The unit additive circle has circumference one and diameter one half.
The probability model uses normalized Haar measure and independent samples.
Ties are resolved by the fixed label enumeration; their almost-sure absence
and invariance under changing that enumeration are separate assertions.
-/

noncomputable section

namespace Balloons

open MeasureTheory
open scoped ENNReal

local instance : DecidableEq UnitAddCircle := Classical.decEq _

/-- Uniform probability measure on the circle of circumference one. -/
def circleLaw : Measure UnitAddCircle := volume

instance circleLaw_probability : IsProbabilityMeasure circleLaw :=
  ⟨UnitAddCircle.measure_univ⟩

instance circleLaw_diffuse : NullSingletonClass circleLaw := by
  constructor
  intro x
  change (volume : Measure UnitAddCircle) {x} = 0
  simpa using (AddCircle.volume_closedBall (T := 1) (x := x) 0)

/-- The circle's metric has diameter one half. -/
theorem circle_dist_le_half (x y : UnitAddCircle) : dist x y ≤ (1 : ℝ) / 2 := by
  rw [dist_eq_norm]
  simpa using (AddCircle.norm_le_half_period (p := 1) (x := x - y) (by norm_num))

/-- Separation one half completes every odd circle input. -/
theorem circle_activeAt_odd_card (S A : Finset UnitAddCircle) (hAS : A ⊆ S)
    (hodd : A.card % 2 = 1) : (activeAt S A ((1 : ℝ) / 2)).card = 1 := by
  classical
  exact activeAt_odd_card S A hAS _ (fun u _ v _ => circle_dist_le_half u v) hodd

/-- Separation one half completes every even circle input. -/
theorem circle_activeAt_even_empty (S A : Finset UnitAddCircle) (hAS : A ⊆ S)
    (heven : A.card % 2 = 0) : activeAt S A ((1 : ℝ) / 2) = ∅ := by
  classical
  exact activeAt_even_empty S A hAS _ (fun u _ v _ => circle_dist_le_half u v) heven

/-- An iid sample of uniform circle locations. -/
def circleSampleLaw (n : ℕ) : Measure (Fin n → UnitAddCircle) :=
  Measure.pi fun _ => circleLaw

instance circleSampleLaw_probability (n : ℕ) : IsProbabilityMeasure (circleSampleLaw n) :=
  inferInstanceAs (IsProbabilityMeasure (Measure.pi fun _ : Fin n => circleLaw))

/-- An even common background and two independent possible values of one coordinate. -/
def oneResamplingLaw (m : ℕ) :
    Measure ((Fin (2 * m) → UnitAddCircle) × (UnitAddCircle × UnitAddCircle)) :=
  (circleSampleLaw (2 * m)).prod (circleLaw.prod circleLaw)

instance oneResamplingLaw_probability (m : ℕ) : IsProbabilityMeasure (oneResamplingLaw m) :=
  inferInstanceAs (IsProbabilityMeasure
    ((circleSampleLaw (2 * m)).prod (circleLaw.prod circleLaw)))

/-- The spatial disagreement probability for one resampling in an odd sample of size `2*m+1`. -/
def survivorPositionNoise (m : ℕ) : ℝ≥0∞ :=
  oneResamplingLaw m {p |
    insertionLocation p.1 p.2.1 ≠ insertionLocation p.1 p.2.2}

/-- Both outputs used in the probability model are actual odd-input survivors. -/
theorem circle_insertion_is_survivor (m : ℕ) (x : Fin (2 * m) → UnitAddCircle)
    (z : UnitAddCircle) :
    run (insertionEvents x z) Finset.univ = {insertionLabel x z} :=
  run_insertion_eq_singleton x z (by simp)

/-- Exact finite agreement probability for the iid circle model. -/
theorem circle_agreement_eq_basin_integral (m : ℕ) :
    oneResamplingLaw m {p |
      insertionLocation p.1 p.2.1 = insertionLocation p.1 p.2.2} =
      ∫⁻ x, ∑ a ∈ backgroundLocations x,
        basinMass circleLaw (insertionLocation x) a ^ 2 ∂circleSampleLaw (2 * m) :=
  random_background_agreement_eq_basin_integral (circleSampleLaw (2 * m)) circleLaw

/-- The finite spatial noise probability is one minus the expected basin-square sum. -/
theorem survivorPositionNoise_eq_basin_integral (m : ℕ) :
    survivorPositionNoise m = 1 -
      ∫⁻ x, ∑ a ∈ backgroundLocations x,
        basinMass circleLaw (insertionLocation x) a ^ 2 ∂circleSampleLaw (2 * m) :=
  random_background_disagreement_eq_basin_integral (circleSampleLaw (2 * m)) circleLaw

end Balloons
