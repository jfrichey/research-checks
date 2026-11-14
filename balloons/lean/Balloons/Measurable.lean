module

public import Balloons.Metric
public import Balloons.Origin
public import Mathlib.MeasureTheory.Constructions.BorelSpace.Metric
public import Mathlib.MeasureTheory.Constructions.BorelSpace.Order

@[expose] public section

/-!
# Measurable finite balloon observables

For finitely many labels the scan depends only on finitely many comparisons
of pair priorities and cutoff comparisons. This makes its configurations,
histories, and transported origins measurable without imposing an unproved
measurability assumption on the sorting algorithm.
-/

noncomputable section

namespace Balloons

instance Event.finite {ι : Type*} [Finite ι] : Finite (Event ι) :=
  Finite.of_injective (fun e : Event ι => (e.left, e.right)) (by
    intro e f h
    cases e with
    | mk a b hab =>
      cases f with
      | mk c d hcd =>
        obtain ⟨rfl, rfl⟩ := Prod.mk.inj h
        rfl)

variable {Ω ι β : Type*} [MeasurableSpace Ω] [Finite ι]

private theorem measurable_decide_le {f g : Ω → ℝ}
    (hf : Measurable f) (hg : Measurable g) :
    Measurable fun ω => decide (f ω ≤ g ω) := by
  classical
  exact (measurable_of_finite (fun p : Prop => decide p)).comp
    (measurableSet_setOfPred.mp (measurableSet_le hf hg))

/-- Any observable of the finite comparison table is measurable. -/
theorem measurable_priority_observable [MeasurableSpace β]
    (p : Ω → Event ι → ℝ) (r : Ω → ℝ)
    (hp : ∀ e, Measurable fun ω => p ω e) (hr : Measurable r)
    (F : ((Event ι → Event ι → Bool) × (Event ι → Bool)) → β) :
    Measurable fun ω => F
      (fun e f => decide (p ω e ≤ p ω f), fun e => decide (p ω e ≤ r ω)) := by
  apply (measurable_of_finite F).comp
  exact (measurable_pi_iff.mpr fun e => measurable_pi_iff.mpr fun f =>
    measurable_decide_le (hp e) (hp f)).prodMk
      (measurable_pi_iff.mpr fun e => measurable_decide_le (hp e) hr)

variable [DecidableEq ι]

/-- The active labeled configuration is measurable, including a measurable random cutoff. -/
theorem measurable_activeThrough_timed (S A : Finset ι)
    (p : Ω → Event ι → ℝ) (r : Ω → ℝ)
    (hp : ∀ e, Measurable fun ω => p ω e) (hr : Measurable r) :
    Measurable fun ω => activeThrough (p ω) (r ω) (timedEvents S (p ω)) A :=
  measurable_priority_observable p r hp hr
    (fun d => run ((orderedEvents S d.1).takeWhile d.2) A)

/-- The full matching history up to a measurable cutoff is measurable. -/
theorem measurable_historyThrough_timed (S A : Finset ι)
    (p : Ω → Event ι → ℝ) (r : Ω → ℝ)
    (hp : ∀ e, Measurable fun ω => p ω e) (hr : Measurable r) :
    Measurable fun ω => history (eventsThrough (p ω) (r ω) (timedEvents S (p ω))) A :=
  measurable_priority_observable p r hp hr
    (fun d => history ((orderedEvents S d.1).takeWhile d.2) A)

/-- A fixed insertion origin has a measurable transported label. -/
theorem measurable_followThrough_timed [MeasurableSpace ι] (S A : Finset ι) (z : ι)
    (p : Ω → Event ι → ℝ) (r : Ω → ℝ)
    (hp : ∀ e, Measurable fun ω => p ω e) (hr : Measurable r) :
    Measurable fun ω => follow (eventsThrough (p ω) (r ω) (timedEvents S (p ω))) A z :=
  measurable_priority_observable p r hp hr
    (fun d => follow ((orderedEvents S d.1).takeWhile d.2) A z)

variable {X : Type*} [MetricSpace X] [MeasurableSpace X]
  [BorelSpace X] [SecondCountableTopology X]

omit [Finite ι] [DecidableEq ι] in
/-- Distances between measurable input locations supply measurable priorities. -/
theorem measurable_location_priority (x : Ω → ι → X)
    (hx : ∀ i, Measurable fun ω => x ω i) (e : Event ι) :
    Measurable fun ω => dist (x ω e.left) (x ω e.right) :=
  (hx e.left).dist (hx e.right)

/-- The labeled metric balloon configuration is a random variable. -/
theorem measurable_metric_active (S A : Finset ι) (x : Ω → ι → X)
    (hx : ∀ i, Measurable fun ω => x ω i) (r : Ω → ℝ) (hr : Measurable r) :
    Measurable fun ω =>
      activeThrough (fun e => dist (x ω e.left) (x ω e.right)) (r ω)
        (timedEvents S (fun e => dist (x ω e.left) (x ω e.right))) A :=
  measurable_activeThrough_timed S A _ r (measurable_location_priority x hx) hr

/-- The labeled collision history of measurable metric input is a random variable. -/
theorem measurable_metric_history (S A : Finset ι) (x : Ω → ι → X)
    (hx : ∀ i, Measurable fun ω => x ω i) (r : Ω → ℝ) (hr : Measurable r) :
    Measurable fun ω => history
      (eventsThrough (fun e => dist (x ω e.left) (x ω e.right)) (r ω)
        (timedEvents S (fun e => dist (x ω e.left) (x ω e.right)))) A :=
  measurable_historyThrough_timed S A _ r (measurable_location_priority x hx) hr

end Balloons
