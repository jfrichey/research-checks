module

public import Balloons.Sampling
public import Mathlib.MeasureTheory.Group.AddCircle

@[expose] public section

/-! # Excluding distance ties by conditioning on one fresh coordinate -/

noncomputable section

namespace Balloons

open MeasureTheory ProbabilityTheory

variable {Ω X ι : Type*} [MeasurableSpace Ω] [MeasurableSpace X]
  [MetricSpace X] [BorelSpace X] [SecondCountableTopology X]

/-- An independent point avoids a sphere whose center and radius depend on the other inputs. -/
theorem random_sphere_null (ν : Measure Ω) (μ : Measure X) [SFinite μ]
    (c : Ω → X) (r : Ω → ℝ) (hc : Measurable c) (hr : Measurable r)
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0) :
    (ν.prod μ) {p : Ω × X | dist p.2 (c p.1) = r p.1} = 0 := by
  have hm : MeasurableSet {p : Ω × X | dist p.2 (c p.1) = r p.1} :=
    measurableSet_eq_fun
      (measurable_snd.dist (hc.comp measurable_fst)) (hr.comp measurable_fst)
  rw [Measure.prod_apply hm]
  have hsec : ∀ ω : Ω, Prod.mk ω ⁻¹'
      {p : Ω × X | dist p.2 (c p.1) = r p.1} = Metric.sphere (c ω) (r ω) :=
    fun _ => rfl
  simp only [hsec, hsphere, lintegral_zero]

variable [Fintype ι]

omit [MetricSpace X] [BorelSpace X] [SecondCountableTopology X] in
/-- Each iid coordinate is independent of the complete vector of remaining coordinates. -/
theorem iid_coordinate_independent_rest (μ : Measure X) [IsProbabilityMeasure μ] (i : ι) :
    IndepFun (fun x : ι → X => x i)
      (fun x : ι → X => fun j : {j : ι // j ≠ i} => x j)
      (Measure.pi fun _ : ι => μ) := by
  have hi : iIndepFun (fun i (x : ι → X) => x i) (Measure.pi fun _ : ι => μ) :=
    iIndepFun_pi (fun _ => aemeasurable_id)
  have hd : Disjoint ({i} : Set ι) {j | j ≠ i} := by
    rw [Set.disjoint_left]
    intro j hj hnj
    exact hnj (Set.mem_singleton_iff.mp hj)
  have h := hi.indepFun_set ({i} : Set ι) {j | j ≠ i} hd (fun j => measurable_pi_apply j)
  exact h.comp (measurable_pi_apply ⟨i, Set.mem_singleton i⟩) measurable_id

/-- A distance involving a fresh coordinate cannot equal a distance between other coordinates. -/
theorem iid_distance_tie_null (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0) (i j k l : ι)
    (hij : i ≠ j) (hik : i ≠ k) (hil : i ≠ l) :
    (Measure.pi fun _ : ι => μ) {x | dist (x i) (x j) = dist (x k) (x l)} = 0 := by
  let P := Measure.pi fun _ : ι => μ
  let target : (ι → X) → X × ℝ := fun x => (x j, dist (x k) (x l))
  have ht : Measurable target :=
    (measurable_pi_apply j).prodMk ((measurable_pi_apply k).dist (measurable_pi_apply l))
  have hs : Measurable fun y : {a : ι // a ≠ i} → X =>
      (y ⟨j, hij.symm⟩, dist (y ⟨k, hik.symm⟩) (y ⟨l, hil.symm⟩)) :=
    (measurable_pi_apply _).prodMk ((measurable_pi_apply _).dist (measurable_pi_apply _))
  have hi : IndepFun target (fun x : ι → X => x i) P :=
    ((iid_coordinate_independent_rest μ i).comp measurable_id hs).symm
  have hmap : P.map (fun x => (target x, x i)) = (P.map target).prod μ := by
    rw [hi.map_prod_eq_prod_map_map ht.aemeasurable (measurable_pi_apply i).aemeasurable,
      (measurePreserving_eval (fun _ : ι => μ) i).map_eq]
  have hn := random_sphere_null (P.map target) μ Prod.fst Prod.snd
    measurable_fst measurable_snd hsphere
  have hm : MeasurableSet {p : (X × ℝ) × X | dist p.2 p.1.1 = p.1.2} :=
    measurableSet_eq_fun
      (measurable_snd.dist (measurable_fst.comp measurable_fst))
      (measurable_snd.comp measurable_fst)
  rw [← hmap, Measure.map_apply (ht.prodMk (measurable_pi_apply i)) hm] at hn
  exact hn

/-- Every sphere in the circle has uniform measure zero, including the antipodal sphere. -/
theorem circle_sphere_null (x : UnitAddCircle) (r : ℝ) :
    circleLaw (Metric.sphere x r) = 0 := by
  have h := (ae_eq_set.mp (AddCircle.closedBall_ae_eq_ball (T := 1) (x := x) (ε := r))).1
  simpa only [circleLaw, Metric.closedBall_sdiff_ball] using h

/-- The circle iid input has no tie between a distance with a fresh endpoint and another pair. -/
theorem circle_distance_tie_null (n : ℕ) (i j k l : Fin n)
    (hij : i ≠ j) (hik : i ≠ k) (hil : i ≠ l) :
    circleSampleLaw n {x | dist (x i) (x j) = dist (x k) (x l)} = 0 :=
  iid_distance_tie_null circleLaw circle_sphere_null i j k l hij hik hil

private theorem fresh_endpoint {κ : Type*} [DecidableEq κ] (e f : Event κ)
    (hne : e.vertices ≠ f.vertices) :
    (e.left ≠ f.left ∧ e.left ≠ f.right) ∨ (e.right ≠ f.left ∧ e.right ≠ f.right) := by
  by_cases hL : e.left = f.left ∨ e.left = f.right
  · by_cases hR : e.right = f.left ∨ e.right = f.right
    · have he := e.distinct
      exfalso
      apply hne
      rcases hL with hL | hL <;> rcases hR with hR | hR
      all_goals simp_all [Event.vertices, Finset.pair_comm]
    · exact Or.inr ⟨fun h => hR (Or.inl h), fun h => hR (Or.inr h)⟩
  · exact Or.inl ⟨fun h => hL (Or.inl h), fun h => hL (Or.inr h)⟩

variable [DecidableEq ι]

/-- Distinct unordered pairs have unequal distances almost surely. -/
theorem iid_event_tie_null (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0) (e f : Event ι)
    (hne : e.vertices ≠ f.vertices) :
    (Measure.pi fun _ : ι => μ)
      {x | dist (x e.left) (x e.right) = dist (x f.left) (x f.right)} = 0 := by
  rcases fresh_endpoint e f hne with h | h
  · exact iid_distance_tie_null μ hsphere e.left e.right f.left f.right e.distinct h.1 h.2
  · have heq : {x : ι → X | dist (x e.left) (x e.right) = dist (x f.left) (x f.right)} =
        {x : ι → X | dist (x e.right) (x e.left) = dist (x f.left) (x f.right)} := by
      ext x
      simp only [Set.mem_ofPred_eq]
      rw [dist_comm (x e.left) (x e.right)]
    rw [heq]
    exact iid_distance_tie_null μ hsphere e.right e.left f.left f.right
      e.distinct.symm h.1 h.2

/-- Simultaneously, every distance tie in a finite iid input is only a reversal of one pair. -/
theorem iid_no_distance_ties_ae (μ : Measure X) [IsProbabilityMeasure μ]
    (hsphere : ∀ a t, μ (Metric.sphere a t) = 0) :
    ∀ᵐ x : ι → X ∂Measure.pi (fun _ : ι => μ), ∀ e f : Event ι,
      dist (x e.left) (x e.right) = dist (x f.left) (x f.right) → e.vertices = f.vertices := by
  classical
  have h : ∀ e f : Event ι, ∀ᵐ x : ι → X ∂Measure.pi (fun _ : ι => μ),
      dist (x e.left) (x e.right) = dist (x f.left) (x f.right) → e.vertices = f.vertices := by
    intro e f
    by_cases he : e.vertices = f.vertices
    · exact Filter.Eventually.of_forall (fun _ _ => he)
    · have hn : ∀ᵐ x : ι → X ∂Measure.pi (fun _ : ι => μ),
          dist (x e.left) (x e.right) ≠ dist (x f.left) (x f.right) := by
        rw [ae_iff]
        simpa using iid_event_tie_null μ hsphere e f he
      filter_upwards [hn] with x hx
      exact fun hxy => (hx hxy).elim
  exact Filter.eventually_all.mpr (fun e => Filter.eventually_all.mpr (h e))

/-- The unit-circle iid model has no ties between distinct unordered pairs. -/
theorem circle_no_distance_ties_ae (n : ℕ) :
    ∀ᵐ x ∂circleSampleLaw n, ∀ e f : Event (Fin n),
      dist (x e.left) (x e.right) = dist (x f.left) (x f.right) → e.vertices = f.vertices :=
  iid_no_distance_ties_ae circleLaw circle_sphere_null

end Balloons
