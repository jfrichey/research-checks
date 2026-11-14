module

public import Balloons.Circle
public import Mathlib.Probability.Independence.Basic

@[expose] public section

/-! # Distinct points in the finite iid input -/

noncomputable section

namespace Balloons

open MeasureTheory ProbabilityTheory

variable {ι X : Type*} [Fintype ι] [MeasurableSpace X]

/-- Two distinct iid coordinates have exactly the product law. -/
theorem iid_pair_law (μ : Measure X) [IsProbabilityMeasure μ] {i j : ι} (hij : i ≠ j) :
    (Measure.pi fun _ : ι => μ).map (fun x => (x i, x j)) = μ.prod μ := by
  have hi : iIndepFun (fun i (x : ι → X) => x i) (Measure.pi fun _ : ι => μ) :=
    iIndepFun_pi (fun _ => aemeasurable_id)
  rw [(hi.indepFun hij).map_prod_eq_prod_map_map
    (measurable_pi_apply i).aemeasurable (measurable_pi_apply j).aemeasurable,
    (measurePreserving_eval (fun _ : ι => μ) i).map_eq,
    (measurePreserving_eval (fun _ : ι => μ) j).map_eq]

variable [MeasurableEq X]

/-- Distinct coordinates of a diffuse iid sample coincide on a null event. -/
theorem iid_pair_coincidence_null (μ : Measure X) [IsProbabilityMeasure μ]
    [NullSingletonClass μ] {i j : ι} (hij : i ≠ j) :
    (Measure.pi fun _ : ι => μ) {x | x i = x j} = 0 := by
  have hd := prod_diagonal_null μ
  rw [← iid_pair_law μ hij,
    Measure.map_apply ((measurable_pi_apply i).prodMk (measurable_pi_apply j))
      (measurableSet_eq_fun measurable_fst measurable_snd)] at hd
  exact hd

/-- A finite diffuse iid sample consists of distinct spatial locations almost surely. -/
theorem iid_sample_injective_ae (μ : Measure X) [IsProbabilityMeasure μ]
    [NullSingletonClass μ] :
    ∀ᵐ x : ι → X ∂Measure.pi (fun _ : ι => μ), Function.Injective x := by
  classical
  have h : ∀ i j : ι, ∀ᵐ x : ι → X ∂Measure.pi (fun _ : ι => μ), x i = x j → i = j := by
    intro i j
    by_cases hij : i = j
    · exact Filter.Eventually.of_forall (fun _ _ => hij)
    · have hn : ∀ᵐ x : ι → X ∂Measure.pi (fun _ : ι => μ), x i ≠ x j := by
        rw [ae_iff]
        simpa using iid_pair_coincidence_null μ hij
      filter_upwards [hn] with x hx
      exact fun hxy => (hx hxy).elim
  filter_upwards [Filter.eventually_all.mpr (fun i => Filter.eventually_all.mpr (h i))] with x hx
  exact fun i j hxy => hx i j hxy

/-- In particular the finite circle model has no repeated starting location. -/
theorem circle_sample_injective_ae (n : ℕ) :
    ∀ᵐ x ∂circleSampleLaw n, Function.Injective x :=
  iid_sample_injective_ae circleLaw

end Balloons
