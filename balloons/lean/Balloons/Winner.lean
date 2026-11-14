module

public import Balloons.SortedRelabel
public import Balloons.Ties

@[expose] public section

/-! # The measurable surviving label and its equivariance -/

noncomputable section

namespace Balloons

open MeasureTheory

variable {ι X : Type*} [Fintype ι] [DecidableEq ι] [MetricSpace X]

/-- Surviving labels in the complete distance-sorted input. -/
def labeledSurvivors (x : ι → X) : Finset ι :=
  run (timedEvents Finset.univ (fun e => dist (x e.left) (x e.right))) Finset.univ

/-- Every odd labeled metric input has a single surviving label. -/
theorem labeledSurvivors_card (x : ι → X) (hodd : Fintype.card ι % 2 = 1) :
    (labeledSurvivors x).card = 1 :=
  run_card_eq_one _ _ (timedEvents_covers _ _) (by simpa using hodd)

/-- Select the unique survivor; uniqueness makes the choice immaterial. -/
def winner (hodd : Fintype.card ι % 2 = 1) (x : ι → X) : ι :=
  Classical.choose (Finset.card_eq_one.mp (labeledSurvivors_card x hodd))

theorem labeledSurvivors_eq_singleton (hodd : Fintype.card ι % 2 = 1) (x : ι → X) :
    labeledSurvivors x = {winner hodd x} :=
  Classical.choose_spec (Finset.card_eq_one.mp (labeledSurvivors_card x hodd))

/-- The spatially ordered process respects every relabeling of a no-tie input. -/
theorem labeledSurvivors_relabel (σ : ι ≃ ι) (x : ι → X)
    (htie : ∀ e f : Event ι, dist (x e.left) (x e.right) = dist (x f.left) (x f.right) →
      e.vertices = f.vertices) :
    labeledSurvivors x = (labeledSurvivors (x ∘ σ)).map σ.toEmbedding := by
  have h := (timed_run_history_map σ.toEmbedding Finset.univ Finset.univ
    (fun e => dist (x e.left) (x e.right)) (fun e _ f _ => htie e f)).1
  simpa [labeledSurvivors, Event.map, Function.comp_def] using h

/-- The winning label transforms by the same permutation as the input labels. -/
theorem winner_relabel (hodd : Fintype.card ι % 2 = 1) (σ : ι ≃ ι) (x : ι → X)
    (htie : ∀ e f : Event ι, dist (x e.left) (x e.right) = dist (x f.left) (x f.right) →
      e.vertices = f.vertices) : σ (winner hodd (x ∘ σ)) = winner hodd x := by
  have h := labeledSurvivors_relabel σ x htie
  rw [labeledSurvivors_eq_singleton, labeledSurvivors_eq_singleton, Finset.map_singleton] at h
  exact (Finset.singleton_inj.mp h).symm

variable [MeasurableSpace X] [BorelSpace X] [SecondCountableTopology X]

/-- The complete finite survivor configuration is measurable. -/
theorem measurable_labeledSurvivors : Measurable (labeledSurvivors (ι := ι) (X := X)) :=
  measurable_priority_observable
    (fun (x : ι → X) e => dist (x e.left) (x e.right)) (fun _ => 0)
    (fun e => (measurable_pi_apply e.left).dist (measurable_pi_apply e.right)) measurable_const
    (fun d => run (orderedEvents Finset.univ d.1) Finset.univ)

/-- Selecting the unique survivor is measurable; no measurable-choice premise is assumed. -/
theorem measurable_winner [MeasurableSpace ι] (hodd : Fintype.card ι % 2 = 1) :
    Measurable (winner (X := X) hodd) := by
  intro s _
  have hs : winner (X := X) hodd ⁻¹' s =
      ⋃ i ∈ s, {x : ι → X | i ∈ labeledSurvivors x} := by
    ext x
    simp only [Set.mem_preimage, Set.mem_iUnion, Set.mem_ofPred_eq,
      labeledSurvivors_eq_singleton hodd, Finset.mem_singleton]
    aesop
  rw [hs]
  exact MeasurableSet.iUnion fun i => MeasurableSet.iUnion fun _ =>
    ((measurable_finset_mem i).comp measurable_labeledSurvivors).setOf

end Balloons
