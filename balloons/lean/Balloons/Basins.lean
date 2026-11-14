module

public import Balloons.Origin
public import Mathlib.MeasureTheory.Measure.Prod

@[expose] public section

/-!
# Agreement as a sum of squared insertion-basin masses

The finite origin map can output an old vertex or the unchanged origin.
For independent diffuse origins, agreement away from the old vertices is
confined to a null diagonal. Thus only the old-vertex basins contribute.
-/

noncomputable section

namespace Balloons

open MeasureTheory Set
open scoped ENNReal

variable {α : Type*} [MeasurableSpace α] [MeasurableEq α]

/-- Two independent diffuse points agree with probability zero. -/
theorem prod_diagonal_null (μ : Measure α) [SFinite μ] [NullSingletonClass μ] :
    (μ.prod μ) {p : α × α | p.1 = p.2} = 0 := by
  rw [Measure.prod_apply (measurableSet_eq_fun measurable_fst measurable_snd)]
  have hsec : ∀ x : α, Prod.mk x ⁻¹' {p : α × α | p.1 = p.2} = {x} := by
    intro x
    ext y
    simp [eq_comm]
  simp only [hsec, measure_singleton, lintegral_zero]

/-- The mass of origins transported to one specified old vertex. -/
def basinMass (μ : Measure α) (f : α → α) (a : α) : ℝ≥0∞ := μ {z | f z = a}

/-- Independent origins agree exactly through squared old-vertex basin masses. -/
theorem agreement_eq_sum_basin_sq (μ : Measure α) [SFinite μ] [NullSingletonClass μ]
    (A : Finset α) (f : α → α) (hf : Measurable f)
    (hrange : ∀ z, f z = z ∨ f z ∈ A) :
    (μ.prod μ) {p : α × α | f p.1 = f p.2} = ∑ a ∈ A, basinMass μ f a ^ 2 := by
  classical
  let B : α → Set α := fun a => {z | f z = a}
  let U : Set (α × α) := ⋃ a ∈ A, B a ×ˢ B a
  let C : Set (α × α) := {p | f p.1 = f p.2}
  let D : Set (α × α) := {p | p.1 = p.2}
  have hUC : U ⊆ C := by
    rintro p hp
    obtain ⟨a, _, ha⟩ := mem_iUnion₂.mp hp
    exact ha.1.trans ha.2.symm
  have hCU : C ⊆ U ∪ D := by
    intro p hp
    change f p.1 = f p.2 at hp
    by_cases ha : f p.1 ∈ A
    · exact Or.inl (mem_iUnion₂.mpr ⟨f p.1, ha, rfl, hp.symm⟩)
    · have hx : f p.1 = p.1 := (hrange p.1).resolve_right ha
      have hy : f p.2 = p.2 := (hrange p.2).resolve_right (by simpa only [← hp] using ha)
      exact Or.inr (hx.symm.trans (hp.trans hy))
  have hmass : (μ.prod μ) C = (μ.prod μ) U := by
    apply le_antisymm _ (measure_mono hUC)
    calc
      (μ.prod μ) C ≤ (μ.prod μ) (U ∪ D) := measure_mono hCU
      _ ≤ (μ.prod μ) U + (μ.prod μ) D := measure_union_le _ _
      _ = (μ.prod μ) U := by rw [prod_diagonal_null μ, add_zero]
  have hdis : PairwiseDisjoint (↑A : Set α) (fun a => B a ×ˢ B a) := by
    intro a _ b _ hab
    apply Set.disjoint_left.mpr
    intro p hpa hpb
    exact hab (hpa.1.symm.trans hpb.1)
  have hm : ∀ a ∈ A, MeasurableSet (B a ×ˢ B a) := by
    intro a _
    exact (measurableSet_eq_fun hf measurable_const).prod
      (measurableSet_eq_fun hf measurable_const)
  change (μ.prod μ) C = _
  rw [hmass, show (μ.prod μ) U = ∑ a ∈ A, (μ.prod μ) (B a ×ˢ B a) from
    measure_biUnion_finset hdis hm]
  apply Finset.sum_congr rfl
  intro a _
  simp only [Measure.prod_prod, basinMass, B, pow_two]

variable [DecidableEq α]

/-- The basin formula applies to the actual transported origin whenever its map is measurable.
The event schedule may depend on the inserted origin. -/
theorem follow_agreement_eq_sum_basin_sq (μ : Measure α) [SFinite μ]
    [NullSingletonClass μ] (A : Finset α) (es : α → List (Event α))
    (hf : Measurable fun z => follow (es z) A z) :
    (μ.prod μ) {p : α × α | follow (es p.1) A p.1 = follow (es p.2) A p.2} =
      ∑ a ∈ A, basinMass μ (fun z => follow (es z) A z) a ^ 2 :=
  agreement_eq_sum_basin_sq μ A _ hf (fun z => follow_eq_self_or_mem (es z) A z)

/-- For an even background the basin identity gives actual survivor agreement.
The null event that an inserted origin was already present contributes nothing. -/
theorem inserted_survivor_agreement_eq_sum_basin_sq (μ : Measure α) [SFinite μ]
    [NullSingletonClass μ] (A : Finset α) (es : α → List (Event α))
    (hc : ∀ z, Covers (es z) A) (heven : A.card % 2 = 0)
    (hf : Measurable fun z => follow (es z) A z) :
    (μ.prod μ) {p : α × α |
      run (es p.1) (insert p.1 A) = run (es p.2) (insert p.2 A)} =
      ∑ a ∈ A, basinMass μ (fun z => follow (es z) A z) a ^ 2 := by
  have ha : ∀ᵐ z ∂μ, z ∉ A := A.countable_toSet.ae_notMem μ
  have h₁ : ∀ᵐ p ∂μ.prod μ, p.1 ∉ A := Measure.quasiMeasurePreserving_fst.ae ha
  have h₂ : ∀ᵐ p ∂μ.prod μ, p.2 ∉ A := Measure.quasiMeasurePreserving_snd.ae ha
  calc
    _ = (μ.prod μ) {p : α × α |
        follow (es p.1) A p.1 = follow (es p.2) A p.2} := by
      apply measure_congr
      filter_upwards [h₁, h₂] with p hp₁ hp₂
      simp only [
        inserted_survivor_eq_follow (es p.1) A p.1 (hc p.1) heven hp₁,
        inserted_survivor_eq_follow (es p.2) A p.2 (hc p.2) heven hp₂,
        Finset.singleton_inj]
    _ = _ := follow_agreement_eq_sum_basin_sq μ A es hf

end Balloons
