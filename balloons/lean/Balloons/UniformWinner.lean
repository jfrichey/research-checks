module

public import Balloons.Winner

@[expose] public section

/-! # Uniformity of the surviving label under iid sampling -/

noncomputable section

namespace Balloons

open MeasureTheory
open scoped ENNReal

variable {ι X : Type*} [Fintype ι] [DecidableEq ι] [MeasurableSpace X]

omit [DecidableEq ι] in
/-- Relabeling iid coordinates preserves their joint law. -/
theorem iid_precompose_preserving (μ : Measure X) [IsProbabilityMeasure μ] (σ : ι ≃ ι) :
    MeasurePreserving (fun x : ι → X => x ∘ σ)
      (Measure.pi fun _ : ι => μ) (Measure.pi fun _ : ι => μ) := by
  have he : (fun x : ι → X => x ∘ σ) =
      MeasurableEquiv.piCongrLeft (fun _ : ι => X) σ.symm := by
    funext x i
    simp [MeasurableEquiv.coe_piCongrLeft, Equiv.piCongrLeft_apply]
  rw [he]
  exact measurePreserving_piCongrLeft (fun _ : ι => μ) σ.symm

variable [MetricSpace X] [BorelSpace X] [SecondCountableTopology X]
  [MeasurableSpace ι] [MeasurableSingletonClass ι]

/-- All labels related by a permutation have equal survival probability. -/
theorem winner_mass_permutation (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0)
    (hodd : Fintype.card ι % 2 = 1) (σ : ι ≃ ι) (i : ι) :
    (Measure.pi fun _ : ι => μ) {x | winner hodd x = σ i} =
      (Measure.pi fun _ : ι => μ) {x | winner hodd x = i} := by
  let P := Measure.pi fun _ : ι => μ
  have hm : MeasurableSet {x : ι → X | winner hodd x = i} :=
    measurableSet_eq_fun (measurable_winner hodd) measurable_const
  calc
    _ = P ((fun x : ι → X => x ∘ σ) ⁻¹' {x | winner hodd x = i}) := by
      apply measure_congr
      filter_upwards [iid_no_distance_ties_ae μ hsphere] with x hx
      have hw := winner_relabel hodd σ x hx
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, ← hw, σ.injective.eq_iff]
    _ = _ := (iid_precompose_preserving μ σ).measure_preimage hm.nullMeasurableSet

/-- Each label in an odd iid input survives with probability exactly the reciprocal input size. -/
theorem iid_winner_uniform (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0)
    (hodd : Fintype.card ι % 2 = 1) (i : ι) :
    (Measure.pi fun _ : ι => μ) {x | winner hodd x = i} = (Fintype.card ι : ℝ≥0∞)⁻¹ := by
  let P := Measure.pi fun _ : ι => μ
  have heq : ∀ j : ι, P {x | winner hodd x = j} = P {x | winner hodd x = i} := by
    intro j
    simpa only [Equiv.swap_apply_right] using
      (winner_mass_permutation μ hsphere hodd (Equiv.swap i j) j).symm
  have hsum : (∑ j : ι, P {x | winner hodd x = j}) = 1 := by
    have hpre : ∀ j : ι, winner (X := X) hodd ⁻¹' {j} = {x | winner hodd x = j} := by
      intro j
      ext x
      simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_ofPred_eq]
    have h := sum_measure_preimage_singleton (μ := P) (f := winner (X := X) hodd)
      Finset.univ (fun j _ => (measurable_winner hodd) (measurableSet_singleton j))
    simpa only [hpre, Finset.coe_univ, Set.preimage_univ, measure_univ] using h
  have hmul : P {x | winner hodd x = i} * (Fintype.card ι : ℝ≥0∞) = 1 := by
    simpa only [heq, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_comm] using hsum
  exact ENNReal.eq_inv_of_mul_eq_one_left hmul

/-- The original odd uniform-circle input has a uniform surviving coordinate label. -/
theorem circle_winner_uniform (n : ℕ) (hodd : n % 2 = 1) (i : Fin n) :
    circleSampleLaw n {x | winner (ι := Fin n) (by simpa using hodd) x = i} = (n : ℝ≥0∞)⁻¹ := by
  simpa [circleSampleLaw] using iid_winner_uniform circleLaw circle_sphere_null (by simpa using hodd) i

end Balloons
