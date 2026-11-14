# Response renewal and backward creation laws

Selected mathematical record supporting the reader. Working proof status and all open sign conditions are retained. Reviews attributed to research agents are AI reviews; independent human review remains outstanding.

## 1. Responses on sleeping configurations

Identify a stable configuration with its sleeping support B. For a terminal
payoff f, let `F_f(C)` be its expectation after stabilizing configuration C,
and define

`g_f(B,x) = F_f(sleep B + active x) - f(B)`, with `g_f(B,s)=0`.

When x is vacant in B, put `t_x(B)=f(B union {x})-f(B)`. Conditioning the
new particle's first instruction gives

`g_f(B,x) = p_x t_x(B) + q_x sum_z K(x,z) g_f(B,z)`.             (1)

For x in B define a subprobability kernel `W_x(B,C)` as follows. Remove
the sleeper at x, force one particle to take a K-step from x, and stabilize
among the remaining sleepers. Discard the whole history if any particle
visits x, including a forced self-loop. Original-sink exits are permitted.
Otherwise record its final support C. Thus `x notin C` and `|C|<=|B|`.

**First-return identity.** For x in B,

`g_f(B,x) = sum_z K(x,z) g_f(B,z)`
`            + sum_C W_x(B,C) [g_f(C,x)-t_x(C)]`.               (2)

**Proof.** An active addition at the occupied x produces two active
particles; force one to jump and hold the other at x. Compare this with
one forced jump added to the original sleeping B. Run the jumped particle
and any descendants outside x in both systems. At a first return to x,
both configurations agree, with two active particles there. Those branches
cancel. On a no-return branch ending at C, the difference is
`F_f(sleep C+active x)-f(C union {x})=g_f(C,x)-t_x(C)`.
Fresh stacks at this stopping line justify the remaining expectations.
This argument allows arbitrary other sleepers and uses no monotonicity.

## 2. The nonnegative operator

Fix a common-particle bound N. For nonempty B with `|B|<=N`, let
`Z_B(z,t)`, t in B, be the probability that one particle started at z first
hits B at t before sleeping outside B or exiting. Set `Z_B(t,u)=1{t=u}`
for t,u in B. Equivalently, on `U=V\B`,

`Z_B|_U = (I-diag(q_U)K_UU)^(-1) diag(q_U)K_UB`.

For the homogeneous versions of (1)–(2), (1) eliminates all vacant marked
positions: `g(B,z)=sum_(t in B) Z_B(z,t) g(B,t)`. Define T on occupied
marked states `(B,x)`, x in B, by

`T((B,x),(C,t)) = 1{C=B} sum_z K(x,z) Z_B(z,t)`
`                   + W_x(B,C) Z_C(x,t)`.                    (3)

Empty C contributes zero. Every entry is nonnegative and T is triangular
by particle count: it never moves from B to a larger C. Its row sums need
not be at most one; on CB's four-site example with N=2, one exceeds 1.10.
The following weighted contraction is therefore substantive.

**Theorem.** Each diagonal count-k block `T_k` has spectral radius below
one. Consequently `I-T` is invertible with nonnegative inverse.

**Proof.** Use the payoff `f_k(S)=1{|S|<k}`. At an occupied marked state
with `|B|=k`, put

`h_k(B,x)=P(|Stab(sleep B+active x)|<k)`.

At this count, the vacant-site creation term in (1) is zero. Applying
(2), and then eliminating vacant marks with Z, gives the exact identity

`h_k = T_k h_k + b_k`,                                      (4)

where

`b_k(B,x)=sum_(|C|=k-1) W_x(B,C)`
`                      P(|Stab(sleep C+active x)|<k)`.

Indeed, when `|C|=k`, the bracket in (2) is the corresponding h-response.
When `|C|=k-1`, subtracting `t_x(C)=-1` leaves the displayed probability.
When `|C|<=k-2`, particle conservation makes that bracket zero.

Every entry of b is strictly positive. To see this, choose a simple
positive-probability K-path from x to s, with internal vertices
`z_1,...,z_m`, and take `C=B\{x}`. A forced particle can follow this path,
leaving C unchanged and never returning to x: at an initially occupied
internal vertex leave one sleeper and send the other particle onwards;
at a vacant vertex jump onwards. The event has probability

`alpha(B,x)=K(x,z_1) product_i [r_i K(z_i,z_(i+1))] > 0`,

where `z_(m+1)=s`, `r_i=p_(z_i)` on C and `r_i=q_(z_i)` elsewhere.
For a direct sink edge, alpha is `K(x,s)`. Thus
`W_x(B,C)>=alpha(B,x)`. An ordinary active addition at the vacant x can
first jump and then realize the same event, giving the second probability
at least `q_x alpha(B,x)`. Therefore

`b_k(B,x) >= q_x alpha(B,x)^2 > 0`.

Equation (4) and `0<=h_k<=1` imply `h_k>=b_k>0`. With
`epsilon_k=min_(|B|=k,x in B) q_x alpha(B,x)^2>0`,

`T_k h_k <= (1-epsilon_k) h_k`.

The weighted maximum norm proves `rho(T_k)<=1-epsilon_k<1`.
There are finitely many blocks, so `rho(T)<1` and
`(I-T)^(-1)=sum_(j>=0) T^j` is finite and nonnegative. This bound is for
the specified finite system; it is not uniform in volume or parameters.

## 3. The killed complete-block identity

Let P map binary all-active supports A to their final sleeping supports B,
including the empty support, and let `E=P^(-1)` be CB's complete-block
inverse. Restricting both matrices to supports of size at most N preserves
the inverse identity, because neither matrix increases particle count.

For a vertex x, let `N_x(A,B)` be the subprobability of starting from
all-active A, never visiting x, and stabilizing at B. Its row is zero if
x is in A. The identity needed below is

`(E N_x)(B,C) = 1{B=C}`                         if x is outside B,
`(E N_x)(B,C) = -(q_x/p_x) W_x(B,C)`             if x is in B.  (5)

**Proof.** Use active-addition operators killed at the first visit to x.
Write them as `a_hat_z`, with `a_hat_x=0`, `a_hat_s=I`, and set
`j_hat_z=sum_w K(z,w)a_hat_w`,
`b_hat_z=(a_hat_z-q_z j_hat_z)/p_z`.
The killed operators commute: stabilization and the event of visiting x
are Abelian. They satisfy `a_hat_z^2=a_hat_z j_hat_z`.

If x is outside B, adding each `b_hat_z`, z in B, starting empty simply
creates B, by the empty-site first-instruction identity. If x is in B,
apply `b_hat_x` last. Since `b_hat_x=-(q_x/p_x)j_hat_x`, this produces
`-(q_x/p_x)W_x(B,.)`. Expanding the same product into simultaneous initial
jumps and then IDLA gives `E N_x`, exactly as in CB. Any IDLA visit to x
leaves a particle there, so subsequent no-visit stabilization discards that
term. The selected particles' prescribed initial jumps precede IDLA;
successive ordinary forced-jump updates would not give E. This proves (5).

## 4. Invertibility on all formal signed-renewal coordinates

Let D have arbitrary entries `(A,x)`, with `1<=|A|<=N` and x outside A,
and suppose `D=M D`, for the signed first-contact matrix formed by unreading
sleepers, reducing multiplicities by IDLA, and eliminating occupied sources.
Extend `D(empty,x)=0`, `D(A,s)=0`, and, at occupied sources, extend by
`D(A,x)=sum_z K(x,z)D(A,z)` until the marked walk reaches a vacancy or s.
Sink accessibility makes this extension unique, including A=V.

For each marked site x, view `D_x` as a column on common supports and put
`g_x=E D_x`, so `D_x=P g_x`. The mixed-common-state extension generated by
unreading and IDLA is exactly its stabilization law applied to g_x:
both obey the same linear first-instruction identities, and agree on every
all-active input. This assertion concerns formal responses, not only actual
gradients of a single terminal payoff.

Let `H_A(x,y)` be R's first-sleep law for the extra walker, allowing sleep
only outside A. First contact followed by stabilizing the common state
then gives the factorization

`(M D)(A,x)=sum_(y outside A) H_A(x,y) [(P-N_y)g_y](A)`.        (6)

Inverting this one-particle first-sleep equation, and using the occupied
source extension, gives, for every A,x,

`q_x P(g_x-sum_z K(x,z)g_z) = -p_x N_x g_x`.

Multiplying by E and applying (5) yields precisely

`g(B,x)=q_x sum_z K(x,z)g(B,z)`                  if x is outside B,
`g(B,x)=sum_z K(x,z)g(B,z)+sum_C W_x(B,C)g(C,x)` if x is in B.

Also `g(empty,x)=0`. Section 2 reduces these equations to `g=Tg` at
occupied marks. Its contraction theorem forces g=0, hence D=0. A finite
square matrix with trivial kernel is invertible, proving the stated claim
for `I-M`. This does **not** assert `rho(M)<1`, entrywise positivity of
`(I-M)^(-1)`, or convergence of the old signed Neumann series.

## 5. What remains for occupation positivity

For target v, take `f(B)=1{v in B}`. Then `t_x(C)=1{x=v}` whenever x is
vacant in C. Equation (1) has a positive target-creation term. Equation
(2), however, contains the negative term
`-1{x=v} sum_C W_x(B,C)` when the mark starts at an occupied target.
The positive inverse of `I-T` alone does not remove that sign.

The desired all-active increment is
`D(A,x)=sum_B P(A,B)g_f(B,x)`. A proof must use this creation law to compare
the positive and negative terms. The new theorem permits a convergent
nonnegative history expansion of the transformed response, while leaving
its target gain/loss balance as the substantive open question. It also
removes the formerly separate general invertibility gap in the AF program.

## 6. A last-arrival deletion identity for complete creation laws

Let `mu_eta(B)` be the final-support law of any finite all-active input eta,
including multiplicities. Let `N_eta^v(C)` be its subprobability law of
ending at C without ever visiting v. For every C not containing v,

`mu_eta(C) = N_eta^v(C)`
`             + (q_v/p_v) sum_(B contains v) mu_eta(B) W_v(B,C)`.       (7)

For binary inputs this follows at once by multiplying (5) by P. There is
also a direct full-book proof which explains the role of all-active input.
On a completed history ending with v occupied, change the **last** target
record's sleep bit from one to zero, retaining its potential arrow d.
Follow the original history while holding the newly emitted packet d.
This leaves the old final support with v removed, and a forced K-step
packet from v. Restrict to its complete continuation avoiding v. Its law
is precisely W: choose a canonical scheduler which never previews an
unused arrow, and condition only on the actually used bits and arrows of
the completed creation history. The unspent terminal resident arrows and
unused book suffixes then retain their independent laws.

The target's last index stays unchanged on this event. Conversely, for
any history which visits v but ends with v empty, changing its last target
jump to a sleep suppresses one packet and leaves v occupied; the suppressed
packet's later continuation avoids v, as in AX §2. The two operations are
inverse on full books. The changed bit gives the weight ratio `q_v/p_v`.
Partitioning into finite histories proves (7). An initially sleeping,
never-visited v has no such creator bit, so this proof does not extend to
arbitrary sleeping inputs.

Put `a_eta=P_eta(no visit to v)` and
`w_v(B)=sum_C W_v(B,C)`. Writing `Q_eta=mu_eta(v occupied)`, summing (7) gives

`1-a_eta = Q_eta [1+(q_v/p_v) E_mu_eta(w_v(B) | v occupied)]`.          (8)

Here `q_v=1-p_v` is the jump probability, while Q is the occupation
probability. The function w_v is decreasing under inclusion of the
sleeping background `B\{v}`: extra residents can be delayed until their
first contact, as in AX §1, so they can only increase the visited set.
Thus an order on the conditioned creation backgrounds would be sufficient
for (8), together with the known visit monotonicity. **That stronger
conditional order is false.**

An exact example has target 1, lower input `{1,3}`, and upper input
`{1,2,3}`, with jump rows (last column is the sink)

```
0: (0, 1/2, 0,   0,   1/2)
1: (0, 0,   1/3, 1/3, 1/3)
2: (1/3, 0, 0,   1/3, 1/3)
3: (1/2, 0, 0,   0,   1/2)
p: (2/3, 1/2, 9/10, 9/10).
```

Adding site 2 increases both occupation at 1 and joint occupation of 1,3,
but the increment of `P(3 occupied | 1 occupied)` is

`-567216996528315050522131684863826450 /`
` 4778182373573871099409781457610716184503 < 0`.

Therefore a coupling of all target-occupied creation backgrounds by
inclusion cannot supply the missing argument. This does not refute the
specific scalar comparison involving w_v or the occupation conjecture.
When eta already contains v, `a_eta=0`, and decreasing the conditional
expectation of w_v in (8) is exactly the original occupation inequality;
restating that scalar comparison is not a proof.


## 7. Deletion residuals intertwine with active additions

Updated 2026-10-15 UTC. **The identity and cone invariance below are
proved. The proposed occupation inequality on the cone is OPEN.**

Let A_y be the stochastic matrix which adds one active particle at y to
a sleeping support and stabilizes. For a marked site z, let
Ahat_y^z be the corresponding subprobability matrix on supports avoiding
z, killed on any visit to z; Ahat_z^z=0. Define the signed matrix L_z,
with rows all stable supports B and columns C avoiding z, by

L_z(B,C) = 1{B=C}                       when z is outside B,
L_z(B,C) = -[(1-p_z)/p_z] W_z(B,C)       when z is in B.

For a measure mu, its deletion residual is nu_z=mu L_z. Equation (7)
says that when mu is an all-active stabilization law, nu_z is exactly
the unnormalized no-visit creation law. In particular nu_z>=0.

**Intertwining theorem.** For every y,z,

A_y L_z = L_z Ahat_y^z.                                      (9)

**Proof.** Let P be the invertible binary all-active stabilization matrix.
For each binary input A, equation (7), also valid with multiplicities,
identifies the row of P A_y L_z with the no-z-visit final law from
active A+delta_y. Stabilizing A first and then adding y, killing on a visit
to z, gives that same law: stabilization and its visited set are Abelian,
and the unused suffixes are fresh. This is the corresponding row of
P L_z Ahat_y^z. Cancel P. If y=z both sides vanish. No positivity of
P inverse is used.

Consequently the cone

C = {mu>=0 : mu L_z>=0 for every z}

is preserved by every active-addition matrix A_y. Indeed mu A_y>=0
and (mu A_y)L_z=(mu L_z)Ahat_y^z>=0. This is invariance of the full cone;
the subset C_N supported on at most N sleepers maps to C_(N+1), not
necessarily to itself. Every finite all-active creation law belongs to C.
No conditioning or division by a survival probability occurs.

**Candidate at two common particles.** For mu in C_2, every source x
and target v, is mu g_v(.,x)>=0? Here
g_v(B,x)=P(v occupied after adding active x to sleeping B)-1{v in B}.
This would imply the desired two-common-particle theorem. It asks for
more than the original conjecture, because membership in C is only a
necessary creation-law condition; equivalence is not asserted.

There is an exact finite dual formulation. A sufficient, and by finite
cone duality necessary, certificate consists of nonnegative functions
h_z(C) on supports C avoiding z and a nonnegative remainder r(B) such that

g_v(B,x) = r(B) + sum_z (L_z h_z)(B),  |B|<=2.                (10)

It yields mu g=mu r+sum_z (mu L_z)h_z>=0. Thus the unresolved task is
an arbitrary-graph construction of h and r from the backward histories,
rather than choosing a favorable history without an inverse. Loop
cancellation may run through residuals whose terminal support contains v;
they cannot simply be discarded.

**Finite evidence, not a proof.** The exact rational simplex calculation
in verify_creation_cone.py checks both primal feasibility and dual
inequalities. It tests all source/target pairs for the archived kernels
and the AX forced-return, RN conditional-order, CR uniform-routing and
CX creator-capacity examples, at common bound two. The earlier in-memory
run also tested ten seeded three/four-site kernels (seed 20261015);
none had a negative cone response. The maintained script reproduces that
screen with --screen. No general conclusion follows from these cases.

The target's deletion constraint alone is insufficient. On the archived
two-site path it leaves negative responses even at count two. Keeping
only residual coordinates whose final supports avoid the target also
fails on the archived four-site graph. A normalized feasible measure for
that relaxation has masses

empty:8/175, {2}:4/105, {1}:4/175, {1,2}:1/105,
{0}:16/175, {0,3}:128/175, {0,2}:4/105, {0,1}:4/175,

and its response for source 0, target 3 is
-45855894649/22349040532000. This is a failure of the relaxed constraints,
not of C, an actual all-active law, or occupation monotonicity.


### One sleeper and a strict enlargement of creation laws

**Proposition.** The cone sign in (10) holds for every mu in C_1.
Only the target's deletion inequality is needed.

Let rho be the stabilization law of one active particle at v. Its mass
rho({v}) is positive. Since it certainly visits v, rho L_v=0. The
constraint mu L_v>=0 therefore gives

mu = [mu({v})/rho({v})] rho + sigma,

where sigma>=0 is supported on target-free configurations. On that support
g_v(B,x)>=0. Also rho g_v is the proved one-common-particle increment
[R §4], so mu g_v>=0. The same proof gives an explicit certificate:
h_v(C)=g_v(C,x) for C avoiding v, all other h_z=0; its only possibly
nonzero remainder is at {v}, and equals the one-common increment divided
by rho({v}).

The cone is strictly larger than nonnegative mixtures of all-active
creation laws. On the symmetric two-site path with half sink leakage
at each site and p=(1/2,1/2), the law

mu(empty)=5/21, mu({0})=mu({1})=8/21, mu({0,1})=0

has every deletion residual nonnegative, but its unique binary all-active
input expansion mu E has empty-input coefficient -1/7. Surplus IDLA
therefore cannot supply a nonnegative arbitrary-multiplicity mixture
either. This confirms that C_2 positivity is an additional conjecture.

## 8. A normalized backward chain which retains the original roots

**Theorems below are proved.** They do not establish occupation order.

Let pi=P(V,.) be the stabilization law from one active particle at every
vertex. It is stationary for every active addition A_x: perform surplus
IDLA on V+delta_x, sending the surplus to the accessible sink, then
stabilize the remaining all-active V with fresh suffixes. Also pi(B)>0
for every stable support B. Hold one active at each desired B-site, route
all other particles to the sink using jumps, then sleep the retained
particles. There is a finite positive-probability legal realization.
In particular pi L_z=0 for every z.

Write ell_z=(1-p_z)/p_z. For C not containing z define

R_z(C,B)=ell_z pi(B) W_z(B,C)/pi(C),  B containing z.          (11)

**Normalization and termination.** Equation pi L_z=0 says that every
row of R_z sums to one. Thus (11) is a genuine backward conditional law
of a complete no-return deletion avalanche; its selected histories are
not assigned their unconditioned probabilities. It can restore several
sleepers and can move existing ones. It is not inclusion-monotone.

A positive transition satisfies |B|>=|C|, by particle conservation in W.
Moreover R_z(C,C union {z})>0. To see this, force the removed z-particle
along a simple positive-probability sink path avoiding z after departure,
leaving each encountered C-resident asleep again. This is the same
positive no-return path used in the contraction proof in Section 2.

There are finitely many choices (C,z). The minimum of these strict
one-particle-increase probabilities is positive. Consequently any adapted
rule which selects a missing z at every step reaches V almost surely,
with finite expected time: each step has a uniformly positive conditional
probability of increasing the particle count. This is a bound for the
specified finite system only. Stopping on exceeding a prescribed count
also gives finite termination. Cycles at a fixed particle count do not
require an assumed spectral sign.

For any nonnegative measure mu, put h_mu(C)=mu(C)/pi(C). Then

(mu L_z)(C)=pi(C)[h_mu(C)-sum_B R_z(C,B)h_mu(B)].              (12)

Thus C is exactly the cone of nonnegative functions excessive for every
permitted backward insertion, expressed as densities against pi.

**The entire root-exit law is independent of the rule.** Fix any root
set R. Run backward insertions at missing roots until the state contains
R, without a particle-count cap. For each B containing R let
H_R(C,B) be the exit probability at B from C. These probabilities are
the same for every adapted choice of missing roots.

Here is a finite-dimensional proof. For a signed row mu, write mu=alpha P.
Since P L_z=N_z and N_z has zero rows at inputs containing z, while its
square block on inputs avoiding z has full rank, mu L_z=0 exactly when
alpha is supported on inputs containing z. The block's rank follows
directly from invertibility of P and the identity rows of L_z on supports
avoiding z. Consequently the measures mu_D with D containing R form a
basis of the common kernel of all L_z, z in R. Densities against pi
give 2^(|V|-|R|) independent functions harmonic for every missing-root
choice.

For any fixed policy the absorbing states are precisely the supports
containing R, and absorption has finite mean. Its harmonic space has
that same dimension. The common harmonic space therefore equals the
policy's whole harmonic space: every boundary function has a unique
extension harmonic for every choice. Bounded optional stopping proves
the asserted rule independence, including adapted rules. This is
independence of exit laws for posterior resampling; it is distinct from
the fixed-book restoration theorem in Section 11.

**Original-root representation.** Fix a binary all-active input A of
size k. Starting from stable C, repeatedly select any missing original
root z in A\C and apply R_z. Stop on reaching A or exceeding k particles.
Let H_A(C) be the probability of reaching A; set H_A(C)=0 when |C|>k.

For the actual creation law mu_A, equation (7) gives mu_A L_z=0 for
every z in A. Its density h_A is therefore harmonic for each permitted
choice of a missing root, and is zero above count k. The preceding
termination argument and bounded optional stopping give

mu_A(C)=pi(C) [mu_A(A)/pi(A)] H_A(C).                         (13)

The positive factor mu_A(A)/pi(A) is fixed by normalization. The empty
input is the immediate-stop case. Formula (13) holds for every adapted
selection rule; in particular the hitting probability H_A is independent
of that rule. It retains the original common roots, unlike replacing their
creation law by an arbitrary sleeping background. General multiplicities
remain covered by the saved surplus-IDLA reduction.

There is a stronger comparison, still insufficient for occupation order.
For any mu in C and any binary input A, apply the A-root stopping rule
to h_mu=mu/pi, stopping at A or above |A|. This function is excessive
for every insertion, so optional stopping gives
h_mu(C)>=H_A(C)h_mu(A). Combining with (13) yields

mu(C)>=[mu(A)/mu_A(A)] mu_A(C), for every stable C.          (14)

For mu=mu_D this holds for every all-active D, without a nesting
assumption. Summing gives mu_D(A)<=mu_A(A), and the likelihood ratio's
minimum is attained at A. The remainder after subtracting this component
is nonnegative, but need not belong to C, so repeated peeling is not
justified. The comparison does not determine a target marginal: the
remaining mass may have either target bit.

The exact verifier checks (11)--(13) for all binary inputs of the archived
two/three/four-site kernels, using both smallest-missing-root and
largest-missing-root rules. The optional-stopping argument proves the
arbitrary-graph result; the finite checks only audit its implementation.

## 9. Completing target-return cycles before accounting for credit

There is a smaller sufficient form of the cone certificate (10).
Fix source x and target v, write g(B)=g_v(B,x), and let O be the supports
containing v with |B|<=2. Eliminate target-free coordinates first by setting

r_v(B)=g(B)+ell_v sum_(D avoids v) W_v(B,D)g(D), B in O.       (15)

The singleton value r_v({v}) is nonnegative by Section 7. It can be
negative at a pair {v,w}. For z!=v and C in O avoiding z, define the column

M(B;z,C)=1{B=C}-ell_z 1{z in B} W_z(B,C)
  -ell_v ell_z sum_(D avoids v, z in D) W_v(B,D)W_z(D,C).     (16)

A sufficient certificate is a>=0 with r_v-Ma>=0. Indeed, put h_z(C)=a(z,C)
on these target-occupied coordinates, with h_z zero elsewhere for z!=v,
and set

h_v(D)=g(D)-sum_(z!=v)(L_z h_z)(D), D avoiding v.

Each non-target summand here is nonpositive, so h_v>=0. The remaining
target-occupied coefficients are exactly r_v-Ma; target-free remainders
vanish. This proves the sufficiency without selecting or renormalizing
favorable histories. Necessity of this restricted certificate is not
asserted.

The normalized interpretation is particularly useful. From C containing
v, apply R_z for a missing z; if the resulting state lacks v, immediately
apply R_v. Call this complete transition T_z^v. Then

pi(B) M(B;z,C)=pi(C)[1{B=C}-T_z^v(C,B)].

It is a stochastic backward step ending with the target occupied.
Particle count never decreases, and there is a positive chance of
increasing it. At common bound two, killing on reaching three leaves
a chain on {v} and the pairs {v,w}: the remaining state is just one
off-target sleeper. Arbitrary adapted routing still terminates. Thus
normalization and fixed-count loop convergence are proved; the remaining
question is the capacity to transport the signed forcing (15).

**Finite evidence and failed rules.** Exact controlled certificates (16)
exist on the archived and four named obstruction examples. A uniform
choice of the missing z fails: on AX's forced-return graph, source and
target 0 give uniform potential -8/525 at {0,1}. It also fails on CR's
irreducible four-site graph, so irreducibility does not repair this rule.
Ordinary-walk Green weights, sleep-odds Green weights, and first-sleep
routing also fail on named examples. These are failures of specified
policies, not of the controlled-chain formulation or monotonicity.

A more structured source of credit is the proved one-common-particle
coupling. Partition its positive target rewards by the last vertex of
the common walk, including both a terminal sleep and a final jump to the
sink. After division by P_v({v}), these nonnegative credits sum exactly
to r_v({v}). Assigning each such credit to its marked column at {v},
accounting for the column's self-return factor, leaves feasible pair
transport on the four named obstruction examples. Keeping only sleeping
endpoints loses necessary exit credit. Uniform routing of the subsequent
pair transport still fails on CR. The verifier retains this particular
last-vertex proposal; its general capacity and a canonical subsequent
routing rule remain OPEN.

## 10. The last-arrival weighting does not order backgrounds

For an input which initially includes v, define
m_v(B)=1+ell_v sum_C W_v(B,C) on supports containing v. Equation (8)
makes rho_eta(B)=mu_eta(B)m_v(B) a probability measure. It is the
pushforward obtained by changing the last target jump, when present,
into a sleep; complete deletion histories supply the additional mass.

This is a different law from conditioning on final target occupation.
Nevertheless its occupied-background domination also fails. In Section
6's four-site example, with v=1, lower roots {1,3}, and extra root 2,

rho_upper(3 occupied)-rho_lower(3 occupied)
 = -5116827833076768435777/211074938042774429635582400 < 0.

The actual target occupation increment is positive. Consequently the
backward normalization cannot be completed by requiring inclusion of
these weighted backgrounds. The kernel and root-hitting representation
in Section 8 retain their validity.

## 11. Canonical root restoration on complete books

**The statements in this section are proved.** They concern complete
record books, not a pathwise realization of the posterior kernels R_z.
Repeated posterior resampling and repeated changes to one fixed book
must not be identified.

Use the buffer records from [SX](stack-exploration.md):
at arrival i to z, read a sleep bit b_z(i) and a potential jump arrow
d_z(i). An old resident emits its reserved arrow; the new arrival stays
when b_z(i)=1, and otherwise emits its own arrow. Emitted arrows at each
vertex are therefore a prefix of its arrow sequence. For an all-active
input eta, let m_z be the number emitted at z and n_z the arrivals. If
I_z(m) counts arrows entering z among these emitted prefixes, then

n_z=eta_z+I_z(m),  m_z=n_z-b_z(n_z),  b_z(0)=0.              (17)

The local output function n-b_z(n) is nondecreasing. Hence the usual
least-action argument for these monotone local processors gives:
adding input increases emission counts, and changing any sleep bits
from zero to one decreases them. In particular all such sleep changes
preserve finite stabilization of a specified stabilizable book.

Fix a set of roots R contained in the support of eta. If a root z is
absent after stabilization, its last arrival exists and has bit zero.
Change that bit to one and restabilize the original input. Repeat with
any missing root. Every step decreases the emission vector
coordinatewise, strictly at the selected root: the old vector no longer
satisfies (17) there. There are at most sum_z m_z initial steps.

At termination every root is occupied. Every changed index is at or
beyond the corresponding final terminal index, since arrivals decrease.
Reset changes beyond that index to their original values. This does
not affect the final history. The resulting book differs from the
original only by forcing its final terminal root bits to one.

**Canonical restoration theorem.** The final used prefixes and this
reset book are independent of all choices of missing roots. The final
emission vector is the greatest vector obtainable by forcing only
terminal root sleeps and ending with all of R occupied. For fixed R,
it is nondecreasing in eta and nonincreasing in the original sleep
bits, with the arrows held fixed.

Here is the comparison proof, including the maximality assertion.
Take any candidate stabilization from a smaller input eta' and a book
with more original sleep bits. Suppose it ends with all roots occupied,
and only terminal root sleeps were additionally forced. Write m' for
its emissions and t_z=m'_z+1 for its terminal root indices. Initially
the current original process emits at least m', by local monotonicity.
Inductively, assume all sleep changes already made by the restoration
algorithm are at indices at least t_z. These changes cannot affect
the candidate: they are either at its already-sleeping terminal root
or beyond its used prefix. Local monotonicity again gives m>=m'.
If the algorithm now selects a missing root z, its changed index is

n_z=eta_z+I_z(m) >= eta'_z+I_z(m')=m'_z+1=t_z.

This maintains the induction. For the same input and original book,
each terminal outcome thus dominates every other candidate, proving
uniqueness and maximality. Comparing two inputs or two sleep-bit fields
in the same argument gives the stated monotonicities. Unused changes
before resetting can depend on the order; no confluence of those
discarded changes is asserted.

**Exact finite fibres.** Denote this canonical map by F_(eta,R), and
let E_(eta,R) be the event that the original stabilization occupies R.
It fixes every book in E. For a book omega in E, its inverse candidates
are exactly the at most 2^|R| books omega^S obtained by changing the
terminal root bits in a subset S of R from one to zero. Retain the
candidate precisely when F_(eta,R)(omega^S)=omega; call the resulting
family I(omega). No other preimages are possible because all discarded
changes have been reset.

The family I is downward closed. Indeed, starting from omega and
decreasing more of these bits increases the canonical emission vector,
which is already at least that of omega. If the vector for a larger
subset equals the original one, the vector for every smaller subset
does too. Its terminal indices then agree, so its reset book is omega.

For independent ARW books put ell_z=(1-p_z)/p_z. The pushforward of
product book measure under F_(eta,R) has density, relative to original
book measure on E,

Z_(eta,R)(omega)=sum_(S in I(omega)) product_(z in S) ell_z. (18)

To justify the variable indices, partition by the finite tuple of
terminal indices and the flipped subset. On each piece the inverse is
a fixed finite coordinate flip, with probability ratio product ell_z.
The pieces enumerate exactly the preimages just described. Almost-sure
finite stabilization also holds after these finitely many modifications;
null exceptions can be discarded. Consequently

1 <= Z_(eta,R) <= product_(z in R) 1/p_z,
P_eta(R occupied) >= product_(z in R) p_z.                 (19)

The latter bound permits arbitrary additional active input, but is
weaker than comparison with its actual lower-input occupation law.

For two roots a,b the polynomial (18) is

1 + alpha ell_a + beta ell_b + gamma ell_a ell_b,
alpha,beta,gamma in {0,1}, gamma<=alpha beta.

Thus the two indicators of inverse flips are negatively correlated
conditional on the complete canonical image. They are **not** the final
occupation indicators. The admissible family need not factor: take
roots {0,2}, initial terminal records at 0 and 2 equal to (sleep,1),
and outside records at 1 equal to (sleep,0), then (jump,2).
Give all later root records sleep bits and sink arrows. Releasing
either root alone can be undone; releasing both sends particles back
to both roots and creates a larger root-occupied prefix. The family is
{empty,{0},{2}} and Z=1+ell_0+ell_2.

The exact audit in verify_root_restoration.py checks order independence,
input and bit monotonicity, and the density identity on all 1,296
two-site books with two variable records per site and fixed sleeping
sink tails; the density uses zero jump odds at these deterministic
tails. It also checks the nonfactorizing three-site witness and 600
seeded finite books. These audit the general proof; they do not turn
emission monotonicity into occupation monotonicity.

**Pointwise fibre positivity fails.** Apply restoration first to the
lower input and then to the upper input, using the common roots in
both. Upper restoration preserves the restored lower prefix, since a
missing upper root must have a strictly later terminal index. The joint
map has at most 2^(2|R|) inverse candidates, but its weighted target
increment can be negative in a single fibre.

A small exact witness has common roots {1,2}, extra source 0 and target
0. Site 2 is isolated from the others and jumps to the sink. The canonical
book has records (bit,arrow)

0: (1,sink), (0,0), (0,1);
1: (1,0), (0,sink);
2: (1,sink),

followed by sleeping records with sink arrows. Both root terminal
indices are one in both canonical stabilizations. All four subsets
of these two root bits are admissible inverse flips. Flipping root 1
changes the target increment from +1 to -1; flipping root 2 does not
change it. The fibre's signed weight is
(1-ell_1)(1+ell_2), negative for p_1<1/2. Only the finitely used records
are prescribed; all continuations are immaterial. Positive K-probability
for the displayed arrows gives a genuine ARW cylinder. This refutes
conditioning on a fixed joint canonical book as a sufficient positivity
argument, even when the extra source is outside the common support.
A proof must also move records between fibres.


The restriction to small sleep probabilities is unnecessary. With one
common root 1, extra source 0 and target 1, take the canonical records

0: (0,1), (0,sink);
1: (0,1), (1,0), (1,sink).

Both canonical terminal supports are {1}. Among flips of the two
terminal-root records (1,2),(1,3), precisely the empty flip and the
flip of (1,3) return to this joint canonical book. Their target
increments are 0 and -1, respectively. The fibre therefore has
weight -ell_1<0 for every 0<p_1<1, including p_1=1/2.
All other continuations are immaterial. An independent second common
root with a sleeping sink record multiplies this by 1+ell_2, giving
the same obstruction at N2. The root-restoration verifier checks the
two admissible inverses exactly.

## 12. Positive entrance measures and the source-weighted remaining sign

For the root set R and a boundary support B containing R, define

Gamma_R(B,C)=pi(C) H_R(C,B)/pi(B).                          (20)

Each row is nonnegative, supported on |C|<=|B|, and equals the identity
on boundary columns C containing R. Its density against pi is harmonic
for all missing roots in R. Write m_R(B)=sum_C Gamma_R(B,C) and let
nu_(R,B)=Gamma_R(B,.)/m_R(B). These are genuine probability measures.

Every all-active input eta containing R has the exact positive expansion

mu_eta(C)=sum_(B contains R) mu_eta(B) Gamma_R(B,C).         (21)

Indeed, mu_eta/pi is common harmonic, so its boundary values determine
it by Section 8. Thus the normalized mixture weights are
mu_eta(B)m_R(B), and sum to one. For the binary input R itself, the only
possible boundary support is R; hence nu_(R,R)=mu_R.

For one extra particle at x outside R, only B=R and B=R union {y}
can contribute. For a target v put
d_y(v)=nu_(R,R union {y})(v occupied)-mu_R(v occupied).
The desired increment is exactly

D_R(x;v)=sum_(y outside R) mu_(R+x)(R union {y})
                         m_R(R union {y}) d_y(v).          (22)

This isolates the actual source-dependent boundary weights. Removing
those weights and demanding each d_y(v)>=0 is **false**, even for two
common roots and a target among them.

An exact obstruction is the directed four-cycle: from z jump to z+1
modulo four with probability 31/32 and to the sink with probability 1/32.
Take p=(1/2,1/2,1/8,1/8), R={0,2}, and target v=0. Then

d_3(0)=
-139677858771237304857855520692679096496067445789542055936 /
 20105157511516783257569470544841023362289343683299598529833 < 0.

The actual add-one increments at both possible new sources 1 and 3 are
positive. The other entrance component compensates the negative one
with the weights in (22). This is not a counterexample to monotonicity.

The entrance measures can be computed without knowing pi. Let Q_R have
rows and columns indexed by supports containing R, with
Q_R(D,B)=P(D,B). The common-harmonic-space argument proves Q_R invertible,
and Gamma_R=Q_R^(-1) P_R, where P_R consists of those rows of P. Count
triangularity permits restriction to supports of size at most |R|+1
when computing (22). The exact verifier checks nonnegativity, every
missing-root harmonic equation, and all boundary columns for the
archived two/three/four-site kernels. It also verifies the four-cycle
obstruction and the positive weighted sums (22).

A further tempting reduction also fails. Even the one-extra block of
Q_R need not be a Green matrix for a killed Markov chain on extra-particle
locations, after any positive diagonal rescaling. With p_z=1/2 on five
vertices, let each jump row be uniform on the following destinations
(s is the original sink):

0:{1,3,s}, 1:{2,s}, 2:{0,3,s}, 3:{4,s}, 4:{0,s}.

For R={3}, restrict Q_R to R and its one-site enlargements. Then

Q_R^(-1)({1,3},{3,4})=1/80636>0.

A killed-chain Green inverse has nonpositive off-diagonal entries;
positive row/column scalings do not change this sign. Adding an
independent sixth vertex, initially active, with p=1/2 and a certain
sink jump gives two common roots {3,5} and the corresponding entry
1/40318. The verifier retains the exact calculation. This rules out
collapsing the source-weight computation to a Markov chain whose only
state is the possible extra-particle location.

## 13. A positive collision lift and the trace/killing identity

**The reduction and identity below are proved; their occupation sign
consequence remains OPEN.** They retain a creator location in addition
to the moving pair, rather than the source-only state ruled out above.

For one ordinary ARW walker, write S(t,u) for its first-sleep law,
e_t=1-sum_u S(t,u) for its sink probability, and
h_t(u)=S(t,u)/S(u,u) for its probability of visiting u. Superscript D
means killing at the first visit to D: S^D and e^D count sleeps and
original-sink exits before that visit. Let a_t^D(y), y in D, be the
probability of first reaching D at y before sleeping or exiting.
Visits to a marked starting vertex count immediately.

The two-creator collision kernel is

C_2(t,u)=S(t,u)h_t(u).

Its row sums are strictly below one. Indeed, run the first active
particle to sleep/loss, then the second until it reaches the first
sleeping endpoint or stops. A contact at u restarts two active
particles at u with fresh suffixes. This gives a positive renewal for
the complete two-particle law, including its killed variants. For
target occupation the one-common increment d_t=Q_v(2 delta_t)-S(t,v)
satisfies

d=C_2 d+g,  g_t=S(t,v)[1-h_t(v)]>=0.                       (23)

For w fixed, let H_w(t,u), u!=w, be the first-sleep law when sleep is
disabled at w, and E_w(t) its sink probability. Explicitly,

H_w(t,u)=S(t,u)+S(t,w)S(w,u)/(1-S(w,w)),
E_w(t)=e_t+S(t,w)e_w/(1-S(w,w)).

Let B(t,w) denote the expected nonnegative terminal payoff f from two
active particles at t and one sleeper at w, with t!=w. Also let U_w be
the payoff from two active particles at w and V_(t,w) the payoff from
one active particle at t plus a sleeper at w. Both are two-particle
quantities supplied by the preceding collision renewal.

There is a strictly substochastic kernel on the n(n-1) states (t,w).
For each u!=w add the following two transitions, allowing coincidences
between their destinations:

T_contact((t,w),(u,w)) += C_2^w(t,u),
T_contact((t,w),(w,u)) += X_w(t,u),

where C_2^w is the two-creator collision kernel killed on visits to w,
and

X_w(t,u)=h_t(w)H_w(t,u)+S^w(t,u)a_t^{u,w}(w).               (24)

The nonnegative reward is

r_f(t,w)=h_t(w)E_w(t)U_w + e_t^w V_(t,w)
 + sum_(u!=w) S^w(t,u)[e_t^{u,w} f({u,w})
             + sum_(z outside {u,w}) S^{u,w}(t,z) f({u,w,z})].

Then B=r_f+T_contact B. All coefficients are probabilities of complete
stopped-walk histories; none are conditioned and assigned their old
weights.

**Proof and termination.** Explore one of the active t-particles until
it hits w, sleeps, or exits, holding the other initial packet. A first
hit at w gives two active at w and one active at t. Move the extra
t-particle until sleep outside w or exit; this gives the H_w,E_w branch.
An exit before hitting w leaves the two-particle state V_(t,w).
A sleep at u leaves one active at t and sleepers u,w. Move that particle
until it hits {u,w}, sleeps, or exits. Hits produce B(u,w) or B(w,u);
no-hit endings give the displayed stable payoffs. The pending-delivery
buffer schedule is legal, with unused arrows and suffixes fresh.

Each row has deficit at least h_t(w)E_w(t)+e_t^w>0. If w is reachable
from t, the first term is positive; otherwise the second is positive
by original-sink accessibility. Thus the finite inverse
(I-T_contact)^(-1) is nonnegative. Any three-active input is recovered
by first exploring one root to sleep/loss and applying the same
two-explorer decomposition to the remaining roots; all nonterminal
contacts reduce to these B states. Only ordinary-walk probabilities
and the two-particle renewal are needed.

There is a useful exact connection to the backwards creation problem.
Erase visits to collision state w from C_2, defining its trace kernel on
t,u!=w by

C_hat_2^w(t,u)=C_2(t,u)
             + C_2(t,w)C_2(w,u)/(1-C_2(w,w)).

Then

X_w(t,u) h_w(u)=C_hat_2^w(t,u)-C_2^w(t,u)>=0.              (25)

The identity remains valid when h_w(u)=0; division then is not allowed.
For a direct algebraic proof, use
S^w(t,u)=S(t,u)-h_t(w)S(w,u) and
a_t^{u,w}(w)=h_t(w)-[S^w(t,u)/S^w(u,u)]h_u(w)
in (24). Multiplication by h_w(u) gives (25).
Probabilistically, tracing erases a renewal collision at w, whereas
physical killing excludes visits to w by either creator walk. Their
difference is exactly the term which changes the protected sleeper
and the active-pair location.

This lift does not by itself prove the sign. If one subtracts the
baseline U_t from B(t,w), the resulting forcing can be negative: on
the symmetric two-site half-sink kernel with p=(1/2,1/2), its entry at
(t,w,v)=(0,1,0) is -21/520. Thus a positive resolvent alone is again
insufficient. The new identity gives an explicit way to track that
forcing across creator collisions rather than discarding its history.

The exact verifier verify_contact_renewal.py reconstructs the complete
two- and three-particle laws and checks (25) against raw ARW absorption
on all eight maintained kernels. These are audits, not a proof of
nonnegative add-one or sleeper-insertion gains.

### Localizing the remaining debt at the protected target

There is a further **proved identity**, but its required flow inequality
is OPEN. Regard the one-common gains for all targets as a matrix

D_1(t,v)=Q_v(2 delta_t)-S(t,v)
       =(I-C_2)^(-1)(S-C_2)(t,v).

This is strictly substochastic: S-C_2 is nonnegative, and
D_1 1=1-(I-C_2)^(-1)e<1. It is the absorption kernel of the
collision chain, with terminal rewards S-C_2 and sink losses e.

For a contact state (t,w), set b=h_t(w). Define the n-column matrices

A((t,w),u)=b(1-b) 1{u=w},
J=(I-T_contact)^(-1) A,
L((t,w),u)=b 1{u=w}+J((t,w),u).

All their entries are nonnegative. Let beta((t,w),v) be the occupation
gain from adding a sleeper at w to two active particles at t. Then

beta((t,w),v)=1{v=w}-[L(I-D_1)]((t,w),v).                 (26)

Equivalently, if
k((t,w),v)=(1-b)1{v=w}+b D_1(w,v) is the proved gain from
adding that sleeper to one active particle at t, then

beta=k+J(D_1-I).                                         (27)

Here is the cancellation behind the identity. First subtract
U_t(v)+b D_1(w,v) from the contact payoff B(t,w;v). Substitution in
(24), using (25) and (I-C_2)D_1=S-C_2, leaves forcing

b(1-b)D_1(w,v)-X_w(t,v)[1-h_w(v)],                        v!=w,

and, when v=w,

b(1-b)D_1(w,w)+(1-b)^2
  -sum_(u!=w) C_2^w(t,u)[1-h_u(w)].

Subtracting the additional boundary term (1-b)1{v=w}
cancels the negative cross-transition term in the first line; in the
second it leaves b(1-b)[D_1(w,w)-1]. Thus, for gamma=beta-k,

(I-T_contact)gamma=A(D_1-I),

which proves (27) and (26). In particular the remaining negative
forcing is supported only where the protected sleeper is the target.
No termwise positivity of the previous forcing was used.

Formula (26) asks whether L is a possible occupation flow, before a
randomized stopping time, for the D_1-chain started at w. Precisely the
needed inequality is

L(I-D_1)<=delta_w, entrywise.                            (28)

If (28) holds, beta is its nonnegative terminal stopping measure,
including the possibility of sink loss. Conversely, since D_1 is
transient and L>=0, that inequality permits such a randomized stopping
rule: at each site stop or continue in the proportions given by the
terminal mass and the flow L. This is a proposed interpretation, not
a construction of the required rule from creator histories.

The verifier checks (26)--(27), strict substochasticity of D_1, and
the raw occupation laws on the eight maintained kernels. An in-memory
screen of 100 four/five-site kernels (seed 2026101507) found no negative
coincident-root sleeper gain. This is finite evidence only. The following
extension retains distinct roots and the actual source; its stopping-flow
inequality remains OPEN.

### Retaining two distinct roots and the actual extra source

The preceding flow identity extends to arbitrary active roots a,b.
They may coincide. For w outside R={a,b}, write b_a=h_a(w), b_b=h_b(w),
and, for u!=w, put

M_ab^w(u)=S^w(a,u)S^w(b,u)/S^w(u,u),
X_ab^w(u)=b_a H_w(b,u)+S^w(a,u)a_b^{u,w}(w).

The second expression is symmetric in a,b. Indeed, with
A_z(u)=S^w(z,u), it equals

b_a A_b(u)+b_b A_a(u)-M_ab^w(u)h_u(w)
 +b_a b_b S(w,u)/(1-S(w,w)).

Let P_ab^w be the nonnegative contact row which assigns M_ab^w(u)
to (u,w), and X_ab^w(u) to (w,u). These are first-contact
probabilities for the same two-explorer schedule as above. Define

L_ab^w=[b_a+b_b-b_a b_b] delta_w+P_ab^w J.                 (29)

This is nonnegative and satisfies the exact identity

Q_v(active a+active b+sleep w)-Q_v(active a+active b)
 =1{v=w}-[L_ab^w(I-D_1)](v).                             (30)

To check the algebra, the bilinear version of (25) is

X_ab^w(u)h_w(u)
 =S(a,u)S(b,u)/S(u,u)
  +[S(a,w)S(b,w)/S(w,w)]C_2(w,u)/(1-C_2(w,w))
  -M_ab^w(u).

Insert this in the two-explorer decomposition and use
Q_v(a+b)=S(a,v)+S(b,v)-S(a,v)h_b(v)
 +sum_u S(a,u)h_b(u)D_1(u,v). The terms with no contact
leave [b_a+b_b-b_a b_b]delta_w in the flow; the contact terms
leave P_ab^w J. This proves (29)--(30). For a=b=t, (29) reduces
to L=h_t(w)delta_w+J by J=A+T_contact J.

For an added active source x, let rho_x^R(w) be its first-sleep law
with sleep disabled on R, and original-sink exits discarded.
Exploring this particle first leaves the common roots unchanged.
Consequently the full N2 increment is

D_(a+b)(x;.)=rho_x^R-Lambda_(a,b,x)(I-D_1),
Lambda_(a,b,x)=sum_(w outside R) rho_x^R(w)L_ab^w.          (31)

Thus the actual source weights are retained explicitly. The N2 sign
is precisely the assertion that this specified nonnegative Lambda is
a valid partial occupation flow of the D_1-chain started with rho.
Nonnegativity of Lambda alone does not establish that assertion.
The exact source-flow verifier checks (29)--(31) against raw ARW
absorption, including distinct common roots and multiplicities.

### A causal charging rule already fails on a directed path

The stronger sleeper bounds (28) and (30) remain OPEN; negative
entrance-component gains from Section 12 do not refute them.
There is, however, an exact obstruction to charging the terms of
(29) chronologically at their corresponding contact-chain occurrences.

Augment T_contact by an initial state i_w for each w, with row
Q(i_w,c)=P_ab^w(c), and initialize with eta(i_w)=rho_x^R(w).
Give i_w reward f_w=h_a(w)+h_b(w)-h_a(w)h_b(w), and contact (t,w)
reward h_t(w)[1-h_t(w)]. If G_Q=(I-Q)^(-1), the canonical charge
measure on these lifted states is d=eta G_Q diag(reward).
Its projection to protected vertices is Lambda.

For any transient substochastic Q, a measure d is realizable as a
randomized stopping distribution from eta exactly when

d G_Q <= eta G_Q.                                       (31a)

Necessity counts expected future visits after stopping. Conversely,
u=(eta-d)G_Q>=0 gives arrival measure u+d=eta+uQ.
Stop at z with probability d(z)/(d(z)+u(z)); transience gives
precisely d. Equivalently d H_A<=eta H_A for every set A,
where H_A is the probability of visiting A. Singleton A already suffice.
By contrast, selecting a visited state after seeing the complete path
has Hall cuts d(A)<=eta H_A. That selection need not be a stopping time.

Take 0->1->2->s, common roots 0,1, extra source w=2, and arbitrary
interior p. Write q_i=1-p_i. Then rho=p_2 delta_2 and D_1(2,.)=0,
so a chronological D_1 realization permits at most one charge.
Only i=i_2 and c=(1,2) have positive canonical charges, with

f=q_1(1+q_0 p_1),  Q(i,c)=q_0 p_1,
Q(c,c)=p_1,  reward(c)=q_1 p_1,
d(i)=p_2 f,  d(c)=p_2 q_0 p_1^2.

Since H_{c}(i)=q_0 p_1 and H_{c}(c)=1, the causal cut fails by

d H_{c}-eta H_{c}=p_2 q_0^2 q_1 p_1^2>0.                 (31b)

At p_i=1/2 the demand is 9/64 and capacity 8/64.
Nevertheless Lambda=p_2(1-p_0 p_1)delta_2 and the actual
increment is p_2 p_0 p_1 delta_2>=0. Complete-path Hall cuts also
pass here. Thus the obstruction is charging the direct f_w term
before resolving the later contact history; delaying or redistributing
that charge remains possible.

A useful special case of the stronger sleeper inequality is proved:
if w jumps only to itself or the sink, then
beta_ab^w=P_(a+b)(no particle visits w) delta_w>=0.
Stabilize outside w with w treated as a sink. Both systems have identical
outside dynamics. If any particle enters w, both have positive active
populations there; every such population has final occupation probability
p_w/[1-(1-p_w)K(w,w)], so their contributions cancel. On no entrance,
only the added sleeper remains. This allows arbitrary cycles elsewhere.


### A terminal pair gives an explicit stopping flow

There is a **proved special case**, allowing any finite all-active
common input eta with eta(w)=0. Suppose B={a,w} has no jump exits
except to the original sink. Put f=P_eta(some particle visits w).
Stabilize outside B while holding every common arrival in B, then
perform surplus IDLA on those arrivals together with the common inputs
initially in B. Let zeta be the probability that its support fills B.
The sleeper-insertion gain is

beta_eta^w=(1-f)delta_w+(f-zeta)D_1(w,.),
0<=zeta<=f<=1.                                          (31c)

Indeed, outside histories agree in the two copies. Replay the common
surplus topplings in the copy with the added sleeper at w; they remain
legal, and a common arrival at w wakes that sleeper. Conditional on
this common exploration, unused suffixes are fresh. If its IDLA support
is empty, {a}, {w}, or B, the remaining gains are respectively
delta_w, (1-b)delta_w+bD_1(w,.), D_1(w,.), and zero, where b=h_a(w).

For the last case, retain one active everywhere in B and drain the
upper surplus to the sink; both marginal continuations have the same
fresh full-input law. This averages over the surplus route and exit
forest. If the first three support probabilities are q_empty,q_a,q_w,
then f=zeta+q_w+bq_a. This proves (31c), including
beta_eta^w(w)>=P_eta(no visit to w). Self-loops are allowed.

For two common particles, possibly at distinct roots, (26) or (30)
identifies the specified flow, with G_D=(I-D_1)^(-1), as

L=zeta G_D(w,.)+(f-zeta)delta_w.                          (31d)

Its explicit D_1-chain rule stops immediately with probability 1-f,
takes one step then stops with probability f-zeta, and runs until
natural killing with probability zeta. The corresponding departure
measures are zero, delta_w, and G_D(w,.).

A positive smaller-resolvent calculation also applies. For any absorbing
boundary B, on its complement E put U=S^B,
C=U^(odot2)diag(U_ii^(-1)), and G=(I-C)^(-1).
Let H_z be the lone ARW particle's probability of entering B at z
before sleeping or reaching the original sink, and c=1-sum_(z in B)H_z.
Then gamma_z=U^(-1)H_z>=0 and nu=U^(-1)c>=1. Explicitly,
gamma_z=diag(p)^(-1)diag(1-p)K_Ez and
nu=1+diag(p)^(-1)diag(1-p)K_Es.
For the one-entry reward g_z=2H_z odot c-C H_z, (32) gives

C^(-1)g_z>=gamma_z odot c+(nu-1)odot H_z>=0.             (31e)

To check this, sum the endpoint bound (32) over both endpoint
indices in H_z=U gamma_z and c=U nu. The terms with repeated endpoint
i contribute equality, those with exactly one endpoint i give the
displayed two products, and all others are nonnegative.

For B={a,w} and two coincident roots outside it, put
kappa_z=K(z,B\\{z})/[1-K(z,z)] and s_z=1-kappa_z.
The outside collision renewal gives one-arrival mass Gg_z,
same-site two-arrival mass G(H_z^2), and different-site mass
2G(H_a H_w). Therefore

zeta=G[kappa_a H_a^2+2H_a H_w+kappa_w H_w^2],
f-zeta=b[Gg_a+s_a G(H_a^2)]+Gg_w+s_w G(H_w^2).

The general flow inequality remains open when a third reachable
vertex prevents the two-site full-input reset. A separate AI audit
verified the preprocessing, reset, and all displayed identities.

### Localization to the target's strongly connected component

The main all-active target question reduces to the target's strongly
connected component B in the positive-jump graph. If target occupation
for K killed on exit from B is increasing for every finite all-active
input on B, then the same target inequality holds on the full graph.

Let A consist of vertices outside B which can reach B. There are no
edges from B to A, or from V\\(A union B) into A union B; exits from B
never return. Stabilize A first, holding inputs and arrivals in B.
For ordered all-active inputs eta<=eta', common books on A give
ordered emission counts m_z^eta<=m_z^eta' by Section 11.
The emitted destination lists are prefixes of the same arrows, so

xi_eta(b)=eta(b)+sum_(z in A)sum_(j<=m_z^eta)1{d_z(j)=b}

is coordinatewise at most xi_eta'(b). Conditional on the entire
upstream books, these vectors are fixed and the books in B remain
independent and fresh. Applying the assumed local inequality and
averaging proves the reduction; no upstream occupation comparison
is needed.

In particular the conjecture is proved for any target whose component
has at most two vertices, with unrestricted upstream components.
For the local two-site assertion, run common surplus IDLA and hold the
extra particle. Empty residual support is immediate, a singleton uses
N1, and full support uses the surplus-to-sink reset above. All steps
allow multiplicities and self-loops. This is also a consequence of the
two-site calculation in [DI Section 7](drain-inventory.md#7-one-off-target-site-with-nonzero-sleep-counts).
Sink accessibility ensures the local surplus drain terminates.
This localization was separately audited; it does not resolve an
irreducible component of three or more vertices.

## 14. A positive decomposition from Hadamard potential inequalities

This section uses a known theorem: the entrywise square of an inverse
M-matrix is again an inverse M-matrix; for a potential matrix the
square is again a potential. See Dellacherie, Martinez and San Martin,
[Hadamard functions of inverse M-Matrices, Theorems 2.2--2.3](https://arxiv.org/abs/math/0601688).
The following algebraic consequences and their application are proved
here. They do not imply the occupation sign.

Let U be a nonnegative inverse M-matrix and V=U^{odot 2}. For column
indices i,j define
c_ij=V^(-1)(U(:,i) odot U(:,j)). Then c_ij>=0. Also, for i!=j,

c_ij(i)>=U(i,j)/(2U(i,i)).                                (32)

For the latter bound, square the vector
U(:,j)-[U(i,j)/U(i,i)]U(:,i), whose i-th entry is zero.
Applying the i-th row of V^(-1), with nonpositive off-diagonal
entries, gives c^2-2c c_ij(i)<=0. For a coordinate k different
from i,j, choose the analogous linear combination vanishing at k;
its transformed square gives c_ij(k)>=0. If a required entry
vanishes, the identity V^(-1)U(:,i)^2=e_i shows that each
corresponding off-diagonal contribution vanishes. The diagonal
case c_ii=e_i is immediate. Transposition gives the same statement
for products of rows.

There is also the four-point inequality

U(a,i)U(b,j)+U(a,j)U(b,i)
 >=2[U(a,:) odot U(b,:)] V^(-1)
                         [U(:,i) odot U(:,j)].            (33)

Here is a direct consequence of the same square theorem. For
nonnegative vectors r,s, the triangular block matrix

[ alpha, r^T U, r^T U s;
  0,     U,     U s;
  0,     0,     beta ]

is inverse M for alpha,beta>0: its inverse has off-diagonal blocks
-r^T/alpha, -s/beta, and zero in the corner. Apply the square
theorem and inspect the inverse's corner. It gives

(r^T U s)^2 >= (r^T U)^{odot 2} V^(-1)(U s)^{odot 2}.

Take r=e_a+e_b and s=e_i+e_j. Terms with a repeated row or column
cancel by V^(-1)V=I, leaving (33). Thus this inequality does not
require reversibility or symmetry of the jump kernel.

For U=S, define lambda_ab=[S(a,:) odot S(b,:)](S^{odot 2})^(-1).
It is a subprobability row. Nonnegativity follows from (32) applied
to the transpose. To see its mass is at most one, use the preceding
block inequality with r supported on a,b and s=S^(-1)1>=0.
The first-sleep matrix is a potential, so this s is nonnegative.

There is an exact positive decomposition of complete two-particle laws:

mu_(a+b)=sum_t lambda_ab(t) mu_(2t)+zeta_ab,  zeta_ab>=0.  (34)

Indeed, lambda_ab C_2=S(a,:) odot h_b. Subtract the corresponding
coincident renewal equations. The residual empty-support weight is
e_a e_b-sum_t lambda_ab(t)e_t^2. Its singleton-u weight is

e_a S(b,u)+e_b S(a,u)
 -2 sum_t lambda_ab(t)e_t S(t,u),

and its {u,z}-weight is the defect in (33), with U=S.
These are nonnegative by (33), since e=S[S^(-1)e] and
S^(-1)e=S^(-1)1-1>=0. Their total mass is
1-sum_t lambda_ab(t).

The killed version U=S^w also improves the contact representation.
Let C^w=U^{odot 2}diag(U(u,u)^(-1)), b=h_.(w), and e=e_.^w.
Both b and e are U-excessive: U^(-1)b,U^(-1)e>=0.
Premultiplying each fixed-w contact block by (C^w)^(-1)
leaves nonnegative cross-transition rates, nonnegative terminal
rewards, and nonnegative row deficits. For example, the cross block is

2 diag(U(u,u))(U^{odot 2})^(-1)diag(b)U-diag(b)
 +[(C^w)^(-1)b^2] S(w,:)/(1-S(w,w)).

Each term, including the displayed difference, is nonnegative by
(32) and expansion of b into nonnegative U-columns. Terminal
two-sleeper rewards have the same form with e in place of b;
three-sleeper rewards use (32) for two distinct columns, and the
remaining rewards are nonnegative products of U-excessive functions.
The rate (C^w)^(-1)[b(1-b)] is nonnegative for the same reason.
Thus the transformed contact system is a finite killed continuous-time
Markov generator. Its inverse is the old positive resolvent times
the nonnegative block diagonal matrix C^w.

### Compatible merging weights for complete creator groups

The preceding inequalities extend to any integer m>=1. Put
V_m=S^{odot m}, and for an ordered list of roots a_1,...,a_m define

lambda_(a_1,...,a_m)
 =[S(a_1,:) odot ... odot S(a_m,:)]V_m^(-1).              (35)

These rows are nonnegative and have mass at most one. Here is a proof
which also supplies the needed polynomial inequalities. For
M=V_m^(-1), the off-diagonal entries are nonpositive. Holder's inequality
with weights -M(i,j)/M(i,i), j!=i, gives

M [S(:,u_1) odot ... odot S(:,u_m)] >=0,

because each column's m-th power has M-image a coordinate vector.
Apply the same argument to the transpose for nonnegativity in (35).
Writing nu=V_m^(-1)1>=0, another application of Holder gives
sum_j nu_j product_r S(a_r,j)<=1, proving the mass bound.

For any endpoint list u_1,...,u_m, the symmetrized independent endpoint
weight dominates the mixture from (35):

(1/m!) sum_permutations sigma product_r S(a_r,u_sigma(r))
 >=sum_t lambda_(a_1,...,a_m)(t) product_r S(t,u_r).        (36)

Indeed, g=V_m^(-1)product_r S(:,u_r)>=0. Holder bounds the
right-hand side by product_(r,j) S(a_r,u_j)^(1/m);
the arithmetic-geometric mean inequality bounds the left-hand side
from below by the same quantity. The argument also allows any
endpoint column which is a nonnegative S-potential, in particular
the original-sink column e=S(S^(-1)1-1).

These merging kernels are compatible with group size. If groups of
sizes r,s,k occupy a,b,c, respectively, then

sum_t lambda_(a^r,b^s)(t) lambda_(t^(r+s),c^k)(u)
 =lambda_(a^r,b^s,c^k)(u),                               (37)

where a^r denotes r copies of a. This follows by inserting V_(r+s)
between the two factors. Thus the probability and final location of
a successful full merge do not depend on the binary merge order.
No independence of the failed branches is asserted.

This has an exact operator consequence relevant to the history
exchange. Let A_z be active addition and B_z the signed sleep-insertion
operator from CB. They commute and satisfy

A_z=e_z I+sum_u S(z,u)B_u.

For commuting variables z_u put ell_a(z)=e_a+sum_u S(a,u)z_u.
By (36), every coefficient of

R_(a_1,...,a_m)(z)
 =product_r ell_(a_r)(z)
  -sum_t lambda_(a_1,...,a_m)(t) ell_t(z)^m

is nonnegative; its pure degree-m powers z_u^m vanish. For m=2,
substitution z_u=B_u gives exactly the positive residual law in (34).
For m=3 the polynomial can still contain B_u^2 B_w and B_u^2.
Repeated signed sleep insertion is not a positive transition, so
coefficient positivity alone does not prove the ARW sign.

Crucially, the two orders are related coefficientwise by the exact
nonnegative identity

R_(a,b,c)=R_(a,b)ell_c
          +sum_t lambda_(a,b)(t)R_(t,t,c).               (38)

This follows from (37). The verifier checks degree two and three,
including original-sink terms, all multiplicities, and (37)--(38).
The current use to investigate is exchanging these complete
degree-three blocks while resolving their repeated-sleeper terms
through the contact renewal; the necessary target inequality remains
OPEN.

### A positive reference model for the signed sleep operators

There is a useful parameter regime in which the signed operators
themselves are stochastic. If p_z>=1/2, then B_z is a Markov kernel.
For a sleeping support C containing z, put
J_z=sum_y K(z,y)A_y, with A_sink=I. The complete-block identities give

delta_C B_z
 =[(2p_z-1)/p_z] delta_C J_z
  +[(1-p_z)/p_z] delta_(C without z) J_z^2.               (39)

Both coefficients are nonnegative and sum to one. On supports avoiding
z, B_z simply inserts its sleeper. For (39), use A_z B_z=A_z^2 and
A_z^2=A_z J_z, then expand A_z=p_z B_z+(1-p_z)J_z.
Thus if every p_z>=1/2, every R_m(B) is a nonnegative kernel with
row mass 1-sum_t lambda_m(t). Equations (35)--(38) are then actual
positive kernel decompositions, on arbitrary initial sleeping
backgrounds. This is not an assertion that those kernels preserve
occupation order.

This also supplies positive creator decompositions for general p.
Choose p_z^*=max(p_z,1/2), and construct the weights and polynomials
at p^*. Let T=E_(p^*)P_p be the saved stochastic parameter-transfer
kernel from CB; the jump odds at p^* are no larger than those at p.
It fixes the empty state and intertwines active additions:

A_z^(p^*) T=T A_z^p.

Right-multiplying the positive reference identities by T gives, for
any list of active roots of length m,

mu_(a_1+...+a_m)^p
 =sum_t lambda_(a_1,...,a_m)^(p^*)(t) mu_(m t)^p
  +delta_empty R_(a_1,...,a_m)^(p^*)(B^(p^*)) T.           (40)

The residual is nonnegative, and the compatible refinement (38)
survives with the actual additions A^p. This uses reference weights,
not an unproved claim that every native low-p polynomial is positive
after substitution. Nor does the transfer preserve an asserted
occupation order: that sign still has to be proved with the final
activation history retained. The verifier audits (40) at degrees two
and three for the maintained kernels.


### The reference operators are a capacity-one sandpile

**Proved representation; it does not assert sandpile monotonicity.**
Suppose all p_z>=1/2. Set
Q(z,u)=sum_y K(z,y)S(y,u) for physical u, and
Q(z,sink)=K(z,sink)+sum_y K(z,y)e_y. Then Q is stochastic,
its sink is accessible, and J_z=Q(z,sink)I+sum_u Q(z,u)B_u.
The operator identity A_z^2=A_z J_z gives

B_z^2=alpha_z B_z J_z+(1-alpha_z)J_z^2,
alpha_z=(2p_z-1)/p_z in [0,1].                           (41)

Thus B_z is addition and stabilization in the following sandpile.
A stable site holds at most one particle. When at least two are
present at z, emit one particle with probability alpha_z, and two
with probability 1-alpha_z; each emitted particle takes an independent
Q-step. Repeat until stable. Sink accessibility gives finite expected
stabilization: the sum of the particles' expected Q-hitting times
drops in conditional expectation by the number emitted at a toppling.
The bounded absorption solution is unique. Products of the stochastic
B_z satisfy these same equations by (41), and give the correct
squarefree stable boundary laws, which proves the representation.

Moreover A_a=e_a I+sum_u S(a,u)B_u. Consequently an all-active ARW
input has the law obtained by independently depositing each initial
particle with endpoint law S(a,.) (or deleting it with probability e_a),
then stabilizing this sandpile. At p=1/2 every toppling emits two,
Q_V=2S-I, and A_z=(B_z+J_z)/2. This is the usual stochastic-sandpile
toppling rule; its active/sleeping half-toppling description is discussed
by Campailla and Forien,
[Sections 2.3–2.5](https://arxiv.org/abs/2507.01572).
Here the effective routing kernel Q need not be a simple random walk.

The potential structure of S is essential to this proposed use.
For a general two-emission sandpile, depositing each source by
S=(I+Q_V)/2 need not give binary add-one monotonicity. On three sites,
with the last column the sink, take

Q = [[0,   1/3, 1/2, 1/6],
     [0,   1/3, 1/3, 1/3],
     [1/2, 0,   0,   1/2]].

For lower roots {0,2}, extra source 1 and target 2, the occupation
increment is -2273/434166. Exact finite-state absorption verifies it.
This is **not** an ARW counterexample: S^(-1)(1,0)=6/19>0, so S
is not an inverse M-matrix. It rules out dropping the ARW potential
condition and applying a general sandpile monotonicity assertion.

### Merge order does not give consistent creator partitions

The successful merging weights (35) cannot in general be realized
as one random partition whose restrictions to all subsets have the
corresponding merging probabilities. On the symmetric two-site
half-sink kernel with p=1/2, S has diagonal 8/15 and off-diagonal 2/15.
For three labelled creators at locations 0,0,1,

sum lambda_(0,0)=1,  sum lambda_(0,1)=8/17,
sum lambda_(0,0,1)=4/13.

If such a partition existed, the two creators at 0 would be together
almost surely. Either of their pairwise merges with the creator at 1
would therefore be exactly the full-merge event. But 8/17!=4/13.
This does not contradict (37)--(38), which concern group-size-aware
merge order. Those identities must not be used as a consistent
coalescing-walk coupling when changing which creator is extra.

**Limits of this decomposition.** The residual zeta_ab in (34)
need not have nonnegative add-one response. On the symmetric two-site
half-sink kernel with p=(1/2,1/2), a=0,b=1,x=v=0, its response is
-58/1105. Also the killed mixing weights need not be bounded by the
unrestricted weights: on the three-site half-sink path, killing site
2 increases lambda_(0,1)(1) by 1/1088. Finally, preconditioning the
contact forcing does not make it nonnegative. These exact checks
rule out componentwise positivity and a simple nesting argument;
the positive decomposition must still be used with global accounting.


The exact verifier verify_creator_potentials.py checks the source-flow
identities, positive decompositions, transformed rates, and stated
obstructions on all eight maintained kernels.

## 15. Cancelling routes inside a complete root block

**The identities and latent-law comparison below are proved.** They do
not prove an occupation comparison. Put r_z=(1-p_z)/p_z,
F_z=A_z-J_z, a_C=product_(z in C) A_z, M_C=I-K_CC, and
d_C=det M_C. Empty products and determinants equal one. Sink
accessibility gives d_C>0 and M_C^(-1)>=0. The commuting additions
satisfy A_z F_z=0, and B_z=A_z+r_z F_z.

For T subset V let J_z^out=sum_(y outside T, including s) K(z,y)A_y.
Then, with determinants in the commutative addition algebra,

product_(z in T) F_z
 = det(diag(A_T) M_T-diag(J_T^out))
 = sum_(U subset T) (-1)^|T\U| d_U a_U
                    product_(z in T\U) J_z^out.             (42)

Here is a direct proof, requiring no diagonalization. The two matrices
diag(A_T)M_T-diag(J^out) and M_T diag(A_T)-diag(J^out)
have the same determinant by principal-minor expansion. Row z of the
second is
F_z e_z + sum_(y in T, y!=z) K(z,y)A_y(e_z-e_y).
Expand by row multilinearity. A selected directed cycle makes the
corresponding row vectors dependent. An acyclic nonempty arrow selection
has a path ending at a row which selected F_u e_u; its last arrow
contributes A_u, so that term contains A_u F_u=0. Only the selection
with no arrows survives. Principal-minor expansion gives the last
expression in (42).

Fix R and C subset R. For z outside C define the positive exit operator

H_z^(C,R)=sum_(y outside R, including s)
 [K(z,y)+K(z,C)M_C^(-1)K(C,y)] A_y.

The bracket sums routes through C before exiting it. On multiplication
by a_C, the equations A_c=J_c hold for c in C. To use them without
cancelling a possibly singular operator, work in the quotient by
Ann(a_C). Eliminating C in that quotient uses only the scalar inverse
M_C^(-1). Apply (42) to R\C for the resulting effective kernel.
The Schur determinant identity gives its principal minor on D\C as
d_D/d_C. Multiplying back by a_C and expanding product B yields

product_(z in R) B_z
 = sum_(D subset R) (-1)^|R\D| d_D a_D
   sum_(C subset D) [r_(R\C)/d_C]
                    product_(z in R\D) H_z^(C,R),           (43)

where r_U=product_(z in U)r_z. Every weight and exit kernel inside
the alternating sum is nonnegative. This cancels all traffic inside R
before releasing particles outside it.

Let E=P^(-1) be the binary stabilization inverse. Applying (43) to
the empty support and using surplus IDLA proves

E(R,R)=Z_R=d_R sum_(C subset R) r_(R\C)/d_C.               (44)

Indeed, a term with D properly contained in R either loses a particle
to the sink or starts at a site outside R. IDLA cannot remove such an
initially occupied outside site, so the term cannot yield binary R.
This argument only identifies the diagonal; it does not assert signs
for every entry in a partial row of E.

For R=V all exit operators are scalars
kappa_z^C=K(z,s)+K(z,C)M_C^(-1)K(C,s). Hence

E(V,D)=(-1)^|V\D| d_D sum_(C subset D)
       [r_(V\C)/d_C] product_(z in V\D) kappa_z^C,
pi(V)=1/Z_V.                                              (45)

The full-row signs are weak signs: coefficients can vanish, for example
E(V,empty)=0 on 0->1->s. The reciprocal follows from EP=I and
P(D,V)=0 for D!=V. For proper R, generally P(R,R)!=1/Z_R.

There is an attractive auxiliary subset law

nu_R(C)=[r_(R\C)/d_C] / sum_(D subset R) r_(R\D)/d_D.       (46)

Its conditional insertion odds at z, given the other occupied sites C,
are 1/[r_z sigma_z(C)], where
sigma_z(C)=1-K(z,z)-K(z,C)M_C^(-1)K(C,z)>0.
The subtracted return-path sum increases with C. Thus the odds increase,
and common-uniform single-site heat-bath updates couple nu_R below
the restriction of nu_S whenever R subset S. The same coupling gives
monotonicity when the r_z decrease.

This is **not the ARW occupation law**. On the two-site half-sink
kernel at p=1/2, in support order empty,{0},{1},V,
nu_V=(9,9,9,12)/39 whereas pi=(7,10,10,12)/39.
Only their full-support masses agree. No representation transferring
(46)'s order to the target response has been proved.

The creator verifier checks (43) for all root sets of size at most
three on the eight maintained kernels, and checks all 38 inverse
diagonals and 38 full-row entries in the five archived matrices.
The general arguments above received a separate AI proof audit;
this is not a novelty or user-review claim.

## 16. Fresh continuations and the full-root cycle-popping fibre

**The results below are proved.** They identify complete-book measures
and compare joint retention of the original roots. General single-target
add-one occupation order remains open.

Fix roots R contained in the support of an all-active input eta.
Let omega stabilize at B containing R, with emission vector m^0.
Remove the first m_z^0 records at each site. At an occupied site the
first remaining record has bit one and contains the resident's unread
arrow; at a vacant site it is entirely unread. Call this shifted book U_B.
Starting one active at each site of B in U_B immediately sleeps them.

For every terminal-root flip S subset R, its inverse stabilization in
omega equals the original prefix followed by the corresponding shifted
inverse process. To replay the prefix legally, hold the newly emitted
terminal-root packets until all of the original prefix has been processed.
Abelianness then gives the same raw endpoint and the emission increment.

The admissible family also shifts exactly. The original unflipped book
is a root-occupied candidate, so canonical restoration of its inverse
cannot fall below m^0. Candidates above m^0 correspond in both
directions to shifted candidates: a terminal root index m_z^0+a_z+1
becomes a_z+1. Maximality in Section 11 therefore shows that admissibility
depends only on B and U_B, not eta or its previous creation history.

Reveal a deterministic forward history using bits when needed and arrows
only when emitted. Conditional on any finite completed history ending at
B, all unread resident arrows and subsequent bits/arrows have their
original independent laws. The first bits at occupied sites are fixed
to one in U_B; their probabilities are already included in the creation
history. Write I_R(B,U_B) for its admissible flip family and C_S for
the raw terminal support after flip S. Define

G_R(B,C)=E [sum_(S in I_R(B,U_B)) r_S 1{C_S=C}].           (47)

Disintegrating the finite-fibre change of variables gives
mu_eta(C)=sum_(B contains R) mu_eta(B)G_R(B,C).
For binary eta=D containing R, invertibility of Q_R from Section 12
therefore identifies G_R=Gamma_R. This gives the existing positive
entrance measures a direct fresh-record construction. In particular,

m_R(B)=E sum_(S in I_R(B,U_B)) r_S,
m_R(R)=1/P(R,R),  G_R(R,C)=mu_R(C)/P(R,R).

These expectations use tails conditional on an original forward history.
A canonical restoration image biases the tails by its fibre density;
they are not independent under that image law.

### Entrance masses decrease under boundary enlargement

**Theorem.** For arbitrary roots R and B containing R,

m_R(B union {y})<=m_R(B), y outside B.                    (47a)

Here and below r_S=product_(z in S)(1-p_z)/p_z. Couple the fresh books
U^-=U_B and U^+=U_(B union {y}) as follows. Keep every arrow at its
original emission index, and keep every bit except those at y. At y
insert a leading sleep bit, shifting only the bit sequence:

b_y^+(1)=1, b_y^+(k+1)=b_y^-(k), k>=1.

Both books have exactly the fresh laws used in (47). The insertion
does not shift arrows. More generally, a common base book for U_R,
with this insertion at each site of B\R, couples all B containing R
simultaneously.

Fix S subset R and flip its first root bits in both books. Let u^-
be the lower canonical restored emission vector and use its final
reset book as the lower comparison book. Its terminal root indices are
t_r=1+I_r(u^-)=u_r^-+1. We claim that the upper restoration ends with
u^+>=u^- coordinatewise.

For k incoming arrows at y, the lower and upper local emission functions
(including their respective initial inputs zero and one) satisfy

Phi_y^-(0)=Phi_y^+(0)=0,
Phi_y^+(k)=k+1-b_y^-(k)=Phi_y^-(k)+1, k>=1.

At roots, all lower restoration edits were at indices at least t_r.
Inductively suppose the same holds for all edits made so far in the
upper restoration; initially there are none. The root books agree below
t_r. At t_r the lower reset bit is one, so the current upper processor
emits at least as much. Other nonroot processors agree, apart from the
comparison at y just proved.

For precision this is a comparison on the used lower prefixes, not on
the complete edited root books. Write T^\pm(v)_z=Phi_z^\pm(I_z(v)), with
the initial input included in Phi. The monotone iteration w^(0)=0,
w^(n+1)=T^-(w^(n)) increases to u^-. Its root arrival indices never
exceed t_r, so T^+(w^(n))>=T^-(w^(n)). If v is the current upper
stabilization vector, then T^+(v)=v. Induction gives w^(n)<=v and hence
u^-<=v. If the upper algorithm next restores a missing root r, its
edit index is 1+I_r(v)>=1+I_r(u^-)=t_r. This maintains the invariant
through termination and proves u^+>=u^-.

The unflipped fresh books stabilize with zero emissions. A flip S is
admissible precisely when its canonical restored emission vector is
zero. Indeed, zero emissions put every terminal root index at one,
restore all flipped first bits, and reset every later edit, reproducing
the entire unflipped book. The converse is immediate. Therefore

I_R(B union {y},U^+) subset I_R(B,U^-)                    (47b)

pointwise. Summing the nonnegative weights r_S and taking expectations
in (47) proves (47a). All arguments take place on the common event of
finite stabilization of the finitely modified books. Accessible sinks
give this event probability one. Self-loops, heterogeneous p, and zero
jump probabilities are allowed. Later root bits have been retained;
no arrow-only criterion or fresh-tail claim after restoration is used.

**Joint-root consequence.** For every finite all-active eta>=1_R,

P_eta(R occupied)>=P_(1_R)(R occupied)=1/m_R(R).          (47c)

Indeed, the boundary columns of Gamma_R form the identity and (21)
gives

1=sum_(B contains R)mu_eta(B)m_R(B)
 <=m_R(R)sum_(B contains R)mu_eta(B).

Equivalently the exact nonnegative slack is

m_R(R)P_eta(R occupied)-1
 =sum_(B contains R)mu_eta(B)[m_R(R)-m_R(B)].

This compares arbitrary larger inputs with the binary-root baseline;
it does not compare two arbitrary larger inputs or each root separately.
In particular it is consistent with Section 12's negative individual
entrance gain. The proof and its restricted-prefix comparison received
a separate AI audit. The maintained root-restoration verifier checks the
simultaneous coupling, edit-index invariant, whole-book admissibility,
and the saved nonfactorizing and later-root-bit witnesses, as well as
exact entrance masses and normalized expansions on the archived kernels.
Those finite checks supplement the general argument above.

**Component factorization and strictness.** Let D range over the strongly
connected components of the positive-edge jump graph. In each local D
process, kill every jump leaving D. Then the fresh fibre factors pathwise:

I_R(B,U) = product_D I^D_(R intersect D)(B intersect D,U_D),
m_R(B) = product_D m^D_(R intersect D)(B intersect D).     (47d)

The first equality identifies a flip set with its component restrictions;
root-free factors contain just the empty flip and have mass one. To prove
it, use the zero-emission criterion above. If the local canonical vector
in D is nonzero, first force the flipped first root bits in all strict
ancestor components to one. They now emit nothing. Restore D before its
descendants: its evolution is exactly local, and descendant edits cannot
affect it. Sleep-bit monotonicity shows that the original global canonical
vector dominates this nonzero local vector. Conversely, if all local
vectors vanish but the global vector does not, choose its first nonzero
component in topological order. It has no incoming emissions. Restricting
the final reset book to this component gives a nonzero local terminal-root
candidate, contradicting local maximality. Independence of the component
books then gives the expectation factorization in (47d).

For y outside B, the inequality in (47a) is strict exactly when the
component of y meets R. Necessity follows from (47d). For sufficiency,
choose r in R in that component, a simple path r to y, and a simple
path y to r. Flip only the first bit of r. Prescribe a single packet
following these paths, assigning arrows in emission order. At other
initially occupied sites require each subsequent arrival bit to be one;
at initially vacant intermediate sites require it to be zero. Require
the first lower y-bit and the second r-bit to be one. The paths may
overlap; their successive visits prescribe different records consistently.
Below, the packet stops at vacant y, leaving only r missing; restoring
r's first bit returns the whole book and gives zero emissions. Above,
the inserted y-sleeper passes its reserved arrow onward, and the packet
returns to r. Every root is already occupied, so restoration stops with
nonzero emissions. This finite cylinder has positive probability; no
unused tails are constrained. Thus its flip {r}, of weight r_r>0,
belongs only to the lower family.

Put U=(union of components meeting R) without R, and let Anc(U) be all
vertices from which U is reachable, including paths of length zero.
Then the joint-root comparison has the sharp form

P_eta(R occupied)>P_(1_R)(R occupied)
 iff supp(eta-1_R) intersects Anc(U).                    (47e)

Indeed, m_R(B)<m_R(R) exactly when B meets U. If no extra origin belongs
to Anc(U), that ancestor set initially has exactly |R intersect Anc(U)|
particles and has no incoming edges from its complement. A terminal
boundary support containing R cannot also occupy U, so the slack in
(47c) is zero. Conversely, route one extra particle to some y in U,
hold one active particle at every root, send all other extras to the
sink, and finally sleep the retained R union {y}. A finite legal
positive-probability realization exists, even with root multiplicities.
Its strictly smaller entrance mass makes the slack positive. This
factorization and strictness proof received a separate AI audit. The
maintained verifier checks 120 finite-book component factorizations,
four exact mass products, and four raw source comparisons on a fixed
condensation graph, including self-loops and repeated roots.

### Full roots reduce to ordinary cycle popping

Now take eta=1_V and R=V. A forward book ending at V has lost no
particles. From its unread arrow tails let I_hat(m) count arrivals
from finite prefixes m>=0. Call m balanced when I_hat(m)=m;
balance excludes sink arrows. Let L be the union of the supports of
all such balanced prefixes. Then

I(omega)={S subset V: S intersect L is empty},
Z_book(omega)=product_(z outside L)(1+r_z).                (48)

**Proof of the fibre identity.** If restoration of omega^S has emissions
u different from m^0, then u>=m^0 and both terminal supports equal V.
Conservation makes a=u-m^0 a nonzero balanced prefix. Its support must
meet S: otherwise every flipped old terminal bit was restored at its
original index and all other changes are beyond the old history, which
would reproduce m^0.

Conversely let a balanced m meet S, and put M=m^0+m. Starting from
omega^S, temporarily force each bit at M_z+1 to one. Conservation makes
M a stabilizing upper bound, so least action gives actual emissions
u<=M. Restoring all old terminal bits gives the original prefix;
sleep-bit monotonicity also gives u>=m^0. No prefix through M contains
a sink arrow, so no particle is lost and the stable outcome is V.
Also u!=m^0, since at a site of S intersect support(m) the old terminal
bit remains zero.

Every temporary edit at M_z+1 is either beyond the actual used bits
or at the actual terminal index u_z+1. Reset the former. The latter
are precisely the allowed terminal-root edits of Section 11. Thus
this is a strictly larger root-occupied candidate, and canonical
maximality excludes S from the original fibre. This proves (48);
no future sleep bit entered the criterion.

For completeness, balanced prefixes are exactly finite cycle-popping
prefixes. Given a nonzero balanced remaining prefix, every selected arrow
stays inside its positive-count support. Its current first arrows contain
a directed cycle. Pop that cycle and repeat; the residual remains balanced.

Cycle popping terminates almost surely. Every positive-probability arrow
occurs infinitely often in its iid stack. If some vertices were popped
infinitely often, that nonempty set would be closed under positive
successors and hence would reach the sink, a contradiction. This
regular-arrow event survives finite shifts and also ensures stabilization
of the finitely modified books used above.

Let M* be a completed cycle-popping vector, whose remaining first-arrow
graph is acyclic. Any legal popping execution stays below M*: if its
next cycle reached a vertex already exhausted relative to M*, balance
of the residual prefix forces every predecessor on the cycle to be
exhausted too, leaving a cycle in the terminal graph. Any other completed
vector below M* must equal it, since a nonzero balanced residual would
contain a current cycle. Thus counts are order-independent and
L=support(M*).

Conditional on the original forward history ending at V, these arrow
tails are iid K, independently of p. Fibre normalization and (45) give

1=pi(V) E product_(z outside L)(1+r_z).

Comparing coefficients for all positive r yields

P(L avoids S)=d_V/d_(V\S),  S subset V.                   (49)

This also identifies L as the vertex range of intensity-one Markov
loops. Indeed, put independent Poisson counts on each rooted closed
word of length k with intensity (1/k) times its K-step product.
The total intensity is sum_(k>=1) tr(K_VV^k)/k=-log d_V.
The probability that all loops avoid S is d_V/d_(V\S).
These avoidance probabilities determine the random subset's law.
Only equality of vertex-range laws is asserted, not that the raw
popped cycles themselves are a Poisson family. Wilson's algorithm and
its connection with Markov loops are established background; see
[Chang--Le Jan, Section 5.2](https://arxiv.org/html/1402.1064#S5.SS2).

For any R, loops confined to R similarly give
Z_R=E product_(z in R outside the loop range)(1+r_z).
The attractive law (46) is consequently an auxiliary loop-weighted
subset law, not an ARW occupation distribution.

### Proper-root limitation and audits

For a proper root set, even later bits at the roots can affect
admissibility. Take R={0,1}, input R, and records

0: (1,1), (0,2);
1: (1,1), (0,2), (b,0);
2: (1,2), (0,0), (0,1),

with all unspecified records (1,sink). If b=1 the admissible flips
are empty and {0}; if b=0 only empty is admissible. Thus (48)'s
arrow-only characterization cannot be extended to proper roots.

The root-restoration verifier checks fresh-tail replay on 600 seeded
histories, including multiplicities, with 1,866 inverse candidates;
it checks (48) on 900 further restored full-root books, and checks the
proper-root witness exactly. These supplement the existing exhaustive
1,296-book audit. A separate AI audit checked the general arguments,
including the distinction between forward-history and image conditioning.

## 17. A positive creator law from cycles inside the input set

**The construction and identities below are proved.** The canonical
restoration used here kills particles on leaving the initial root set R.
It is not the proper-root restoration in the full graph from Section 16.

Treat every outside destination as a labelled absorbing exit, and apply
full-root restoration to input 1_R in this killed system. Under original
product-book measure, the event E_kill that killed stabilization ends at
R has probability 1/Z_R, by (45) for the killed graph. Conditional on
an original forward history in E_kill, there was no exit; its shifted
root arrow tails are iid K, including outside labels, and the outside
books remain untouched and iid.
Let L_R(U) be the cycle-popping range of these root tails with all exits
absorbing. For S subset R, set the first root bits to zero on S and
one on R\S; subsequent root bits and all outside books are fresh.
Let Y_R(U,S) be the full-graph terminal support, with all outside returns
permitted. The exact identity is

mu_R(D)=Z_R^(-1) E sum_(S subset R, S avoids L_R(U))
                    r_S 1{Y_R(U,S)=D}.                   (50)

Apply the killed-system fibre change of variables to the observable
given by full-graph stabilization. A canonical killed prefix has no
outside exit, so it is also a legal full-graph prefix. After terminal
flips, replay that prefix while holding the released root packets;
the full endpoint then depends only on U and S. Integrating the
original canonical-prefix history yields its probability 1/Z_R.
The change of variables applies to any bounded terminal-support payoff;
it need not concern the killed process. The removed prefix would still
matter for total odometers or recorded creator histories.

Since P(L_R avoids S)=d_R/d_(R\S), (50) gives a positive sampler:
choose C subset R with law nu_R from (46), condition the root arrow
tails on L_R avoiding R\C, fix their first bits to 1_C, and stabilize
in the full graph with fresh future bits and outside books. Thus nu_R
is an auxiliary input to an exact creator construction, although
it remains different from the final occupation law. The conditional
arrow environments change with R, and no order of the final outputs
has been established.

For two roots a,b with K(a,a)=K(b,b)=0, every nonempty S has
|R\S|<=1, hence no nonzero popped prefix can avoid S.
Its conditioning is just exclusion of the mutual first-arrow pair
a->b,b->a. Writing c=K(a,b)K(b,a),

Z_R=1+(1-c)(r_a+r_b+r_a r_b),
mu_R=[delta_R+sum_(empty!=S subset R) r_S
 E(1{first arrows are not mutual} delta_(Y_R(U,S)))]/Z_R.  (51)

Only the first arrows are conditioned; a retained root's reserved arrow
must be kept until emitted. All later walks and outside returns remain
part of Y_R. Augmented raw-ARW absorption retaining these frozen arrows
checks (51) in 271 law entries across the maintained kernels.
A separate finite-book audit checks killed-prefix replay in the full
graph on 600 books and 2,143 admissible inverse candidates.

### Independent forest and balanced loop prefixes

There is an exact way to sample the conditioned root arrows without
rejection. For a finite collection alpha of local arrow words, write
w(alpha) for its K-product. Call alpha balanced if its incoming and
outgoing counts agree at every site; all its arrows then stay inside
its support. Let F assign one arrow at each root and have no directed
cycle inside R. Outside destination labels are retained.

The event that ordinary cycle popping removes exactly alpha and
leaves F is precisely the cylinder prescribing alpha followed by F
at each site. Every balanced prefix is poppable by Section 16, and
the acyclic terminal F makes that execution complete. Conversely,
confluence gives this unique decomposition. Thus

P(popped words=alpha, terminal forest=F)=w(alpha)w(F).

Cycle inclusion-exclusion gives sum_F w(F)=d_R. Almost-sure termination
and Tonelli then give sum_alpha w(alpha)=1/d_R. The popped words and
forest are independent, with laws d_R w(alpha) and w(F)/d_R.
Arrows after F are fresh iid tails. Conditional on the popped support
being contained in C, only the prefix law changes: it becomes
d_C w(alpha) on balanced words supported in C. The forest remains
the forest on R, independently of C.

Consequently (50) has this positive sampler. Draw C with law nu_R.
Independently draw balanced words on C with law d_C w(alpha), an
exit forest F on R with law w(F)/d_R, and fresh arrow tails. Concatenate
alpha_z,F_z,tail_z, taking alpha_z empty outside C. Prescribe the bit
at record index 1 to be 1_C(z), and sample all later bits freshly.
Then stabilize binary R in the full graph. The popped words must
remain in the reconstructed stacks; cycle popping did not perform
physical ARW moves, and the ARW need not consume the whole prefix.
In particular the prescribed bit stays at index 1, not at the forest
arrow if a prefix was prepended.

The balanced words can themselves be generated from the intensity-one
loop soup on C. Let N_zy be their edge counts and t_zy in [0,1].
Applying the prefix normalization to the substochastic matrix
K_CC entrywise multiplied by t gives

E product_(z,y) t_zy^N_zy
 =d_C/det(I-K_CC entrywise-multiplied-by t).                (52)

The Poisson loop-soup edge counts have the same generating function
by the absolutely convergent trace-log series. Given a balanced count
table n, every local outgoing ordering is allowed and has the same
weight. Hence its probability is

d_C product_z [n_z! / product_y n_zy!] product_(z,y) K_zy^n_zy.

Conditional on the counts, independently shuffle each site's outgoing
multiset uniformly. This produces precisely the required balanced words.
Raw popped cycles are not asserted to form a Poisson family.

### The latent Gibbs coupling and its limitation

Put H_R=Z_R/d_R. For a finite loop multiset ell with range L(ell),
let a(ell) be the product of its Poisson intensity powers divided by
multiplicity factorials, without the exponential normalizer. The joint
law of C and its loop soup is

Q_R(C,ell)=H_R^(-1) r_(R\C) a(ell) 1{L(ell) subset C}.     (53)

Therefore loops conditional on C form the ordinary soup restricted to C;
conditional on the loops, C contains their range and includes each
remaining site independently with probability p_z.

This gives an explicit monotone coupling for nested root domains.
Use a common soup on the larger domain, retain only loops entirely in
each current subset, and then use common Bernoulli bits for the subset
updates. Inclusion is preserved. Every new disagreement inside the
common domain belongs to an excess loop meeting a previous disagreement;
the added sites are additional seeds. A common empty soup, of probability
d_larger>0 per update, removes dependence on the previous subsets.
Thus a stationary coupling can be constructed from a finite backward
sequence of updates. This is a normalized latent dependency mechanism,
not yet a coupling of the ARW continuations.

There is also a useful conditioning formulation. In original iid books
let C_0={z in R:b_z(1)=1} and
A_R={L_R(U) subset C_0}. Its probability is
P(A_R)=(product_(z in R)p_z) Z_R. Formula (50) says that conditioning
on A_R leaves the entire terminal sleeping-support law mu_R unchanged.
The loop/forest construction samples this conditioned complete-book law;
all later bits and outside books remain fresh.

Attractiveness of C alone is insufficient. For R={0,1},
K(0,1)=K(0,s)=1/2, K(1,s)=1, and p_0=p_1=1/2,
there are no internal cycles. The target-1 continuation probability
given C=empty is 1/4, whereas given C={0} it is zero.
This refutes monotonicity of the continuation kernel in its latent
subset, not add-one monotonicity. The comparison must also use the
new root's forest and subsequent activity.


Fixing the sampled exit forest is insufficient even after averaging the
latent subset and all loop prefixes. Take common roots R={0,1}, with
K(0,1)=K(0,s)=K(1,0)=K(1,s)=1/2, extra source 2 jumping certainly
to 0, and all p=1/2. No root jumps to 2. Condition the sampler on the
exit forest F(0)=1,F(1)=s.

The latent weights of empty,{0},{1},R are (3,3,3,4)/13.
For a proper latent subset there is no balanced prefix. At target 0
the respective lower continuation probabilities are
S(1,0),1,S(1,0), with S(1,0)=2/15; for C=R the probability is one.
Thus the lower conditional marginal is 3/5.

If the extra from 2 jumps to 0, surplus IDLA on R follows every
prepended mutual cycle and then 0->1->s. It consumes all prescribed
root arrows and leaves fresh suffixes and one active at each root.
The upper root law is therefore the unconditioned two-site stationary
law, whose target-0 marginal is 22/39. If the extra sleeps at 2, the
root outcome is unchanged. The conditional increment is consequently
(1/2)(22/39-3/5)=-7/390.

The unconditioned increment is zero, since the common two-site law is
stationary and there is no return to 2. This obstruction is for the
sampler's latent exit forest, not an ARW last-exit forest or a fixed
complete-book fibre. A successful comparison must exchange exit forests
as well as subsets or loop prefixes.

Separate AI audits checked the replay, disintegration, forest factorization,
loop edge counts, and Gibbs coupling. The main occupation and N2
conjectures remain OPEN.


## 18. Erasing active cycles preserves a common family of inputs

**The following normalization theorem is proved.** At each root z in R,
inspect its current first record (bit, arrow). Call a directed cycle
inside R active if at least one of its current bits is zero. Delete
one entire record at each vertex of any active cycle, and repeat.
Let D_R omega be the resulting book and N its vector of deleted counts.

Every execution is an ordinary arrow-cycle-popping execution, hence
terminates on the regular-stack event of Section 16. Available distinct
cycles are vertex-disjoint, so their deletions commute. Termination and
this local diamond property give a unique deleted prefix and final book.

For every finite all-active eta>=1_R, simultaneously,

T_eta(omega)=T_eta(D_R omega),
u_eta(omega)=N+u_eta(D_R omega),                           (54)

where T is terminal support and u counts emitted arrows. To verify one
cycle deletion, deliver one initial input at each cycle vertex, holding
all other inputs. Every zero bit emits its first arrow. Deliver these
first-record arrows around the cycle, continuing whenever delivery
releases another first-record arrow, and holding arrows emitted from
second records. Since the cycle contains a zero, every old arrow is
emitted exactly once and each cycle vertex receives exactly one replacement.
Multiple zero bits merely start distinct segments of this round.
The resulting residents and held packets are exactly those obtained
from one input per vertex in the shifted book. All operations are legal
atomic buffer operations; no sleep test ignores a physically present
particle. Deferred extra inputs, including repeats at roots, can now
be processed. Iteration and abelianness prove (54).

Let B_R be the event that every directed cycle in the initial root
arrow graph consists entirely of bit-one vertices. Then

(D_R)_* P=P(. | B_R).                                    (55)

Here is the full preimage accounting. Let W_R be the finite tuples of
local record words that admit an active-cycle-popping order exhausting
exactly those words. Whether a word tuple belongs to W_R depends only
on its own records, not the terminal book. For every w in W_R and
y in B_R, concatenation wy has a complete legal execution deleting w;
confluence identifies this with the canonical execution. Conversely,
every actual deleted prefix belongs to W_R. Count each word tuple once,
not each available popping order. Its product record weight q(w) factors
from the terminal book, and the corresponding inverse classes are disjoint.
Thus, for nonnegative measurable F,

E F(D_R omega)=(sum_(w in W_R) q(w)) E[1_(B_R) F(omega)].

Since b_R=P(B_R)>=product_R p_z>0, putting F=1 identifies the sum
as 1/b_R and proves (55). Equations (54)--(55) imply that B_R is
independent of the terminal-support vector of any finite family
eta_j>=1_R. This does not assert independence of total odometers.

The normalized first-arrow law is explicit. If a assigns one destination
to each root, write cyc_R(a) for its cyclic vertices, and let kappa_z
include all full-graph destinations and the sink. Its weight is

b_R^(-1) product_(z in R) kappa_z(a_z)
            product_(z in cyc_R(a)) p_z.                 (56)

Conditional on these arrows, cyclic vertices have first bit one;
the remaining first bits are independent Bernoulli(p_z), and later
records and outside books are fresh. An adapted block may therefore
be normalized when its records are still fresh and every root has a
reserved input. Selecting a block after previewing its records needs
additional accounting.

The earlier conditioning event A_R in Section 17 has the same
simultaneous-family extension. Apply its killed-system fibre identity
to a payoff of (T_eta_j)_j. Replay the common killed prefix while
holding each eta_j-1_R and the terminal packets released by its inverse
flips. Every full continuation then uses the same shifted tails and
flip set, with no residual dependence on the removed prefix. Hence

Law((T_eta_j)_j | A_R)=Law((T_eta_j)_j), eta_j>=1_R.        (57)

No equality of A_R and B_R is needed or asserted.

The B_R first-sleep field is not attractive in general. On three sites
with every off-diagonal K entry 12/25 and sink probability 1/25,
let h(C) be the probability that every current arrow cycle is contained
in C. Then d_R=1369/15625, h({0})=d_R,
h({0,1})=h({0,2})=4969/15625, and h(R)=1. Consequently

h({0,1})h({0,2})-h({0})h(R)=3300336/244140625>0.

The independent Bernoulli weights cancel from the lattice-condition
test, so this is a strict failure of log-supermodularity for every
interior sleep vector. In contrast, nu_R in Section 15 is attractive.
Neither conditional normalization proves an occupation comparison.

The root verifier checks 1,600 support and emission-count replays on
800 seeded books, with two opposite cycle-selection orders, for both
binary R and R plus one extra input. It also checks the lattice witness
exactly. A separate AI audit verified the legal replay, confluence,
preimage multiplicities, and both simultaneous-family statements.


### A cycle containing one missing initial input

The following local identity is **proved** and retains the joint lower
and upper endpoints. Fix a directed first-arrow cycle
x=z_0->z_1->...->z_(k-1)->x, with k>=2 and old bits b_i not all one.
Let U delete its first record at every cycle vertex. The upper input
1_C+xi has the same endpoint on U by (54). Here C is the cycle's vertex
set and all additional inputs xi may be deferred. The lower input
1_(C\\{x})+xi has these exact boundary descriptions on U:

- If b_0=0 and another bit is zero, keep active inputs at C\\{x}.
- If b_0=1, let j>=1 be the first zero after x. Keep reserved sleepers
  at z_0,...,z_(j-1), each with its old arrow to the next cycle vertex;
  leave a hole at z_j and active inputs at z_(j+1),...,z_(k-1).
- If b_0=0 and every other bit is one, keep the reserved sleepers at
  C\\{x} and prepend the still-unread old record (0,z_1) at x.

Add xi in every case. To prove the classification, process old
first-record arrows while holding packets emitted from replacement
records. In the first case, activity reaches x, passes through its zero,
and completes the round. In the second, activity from the first zero
propagates forward to x and stops at its one bit, leaving the displayed
sleeping arc. In the third no initial input encounters a zero, so no
old arrow is emitted. Comparing with initialization on U justifies
these descriptions using atomic buffer operations.

The third case can instead be kept as a dormant component: it retains
C\\{x} asleep until the first external/deferred arrival at y in C.
If no arrival occurs, those sleepers remain. At the first arrival,
replace the component by fresh active input 1_(C\\{x})+delta_y on U,
holding other pending packets. That arrival triggers the entire old
cycle. This is a one-time activation rule; for y!=x it puts two inputs
at y and none at x. The orientation disappears from this continuation
rule but remains in its probability weight.

Put P_*=product_(i>=1)p_(z_i) and
a_gamma=product_i K(z_i,z_(i+1)), with z_k=z_0.
After summing the old bits, the three types have unnormalized weights

w_rot=(1-p_x)(1-P_*),
w_j=[product_(i<j)p_(z_i)](1-p_(z_j)), 1<=j<k,
w_dorm=(1-p_x)P_*.

Their sum is 1-product_C p. Multiply by a_gamma for the unconditioned
first-cycle cylinder; divide by 1-product_C p conditional on that
specified cycle being active. The identity holds pointwise in U,
including for any finite family of deferred inputs. Individual patterns
retain their full factors a_gamma product_i p_i^b_i(1-p_i)^(1-b_i);
unused bits and reserved arrows cannot be dropped.

A moved hole with at most one reserved sleeper does not close this
identity. If k>=3 and only b_x=0, with no extra input, the lower endpoint
is deterministically C\\{x}. On a U-cylinder of zero/sink records,
any state with at most one reserved sleeper yields at most one sleeper.
Even for k=2 the unread zero matters: take old records
x:(0,y), y:(1,x), fresh U_x(1)=(1,s), U_y(1)=(0,s), and remaining
records zero/sink. With no extra input the lower endpoint is {y};
with delta_x it is {x}. Forgetting the old zero at x while keeping
only the reserved sleeper at y gives {x,y} in the latter experiment.

For a self-loop the old active record is (0,x). It is erasable for
every input, even one missing x: it is unused until an arrival, and
then makes exactly one self-jump before exposing U. The emission
correction is one precisely when x is visited; the cylinder weight
is (1-p_x)K(x,x).

The root verifier checks the joint boundary identity on 1,120 seeded
books for cycle lengths two through five, and checks both closure
obstructions. This is a local exact transformation, not a global
comparison: iterating it requires the reserved arcs and activation rules.


## 19. Higher set inequalities for the inverse diagonals

**Proved:** for fixed accessible K and positive odds r, there are
nonnegative weights w_S, depending only on K,r, such that

log Z_R=sum_(i in R) log(1+r_i)-sum_(empty!=S subset R) w_S. (58)

In particular all mixed set differences of log Z of order at least two
are nonpositive. For distinct x,y outside A this gives

Z_(A+x) Z_(A+y)>=Z_A Z_(A+x+y).                          (59)

Here is a formal power-series proof with convergence. Write
p_i=(1+r_i)^(-1), q_i=1-p_i, and introduce vertex variables t_i.
Set f(t)=det(I-diag(t)K)^(-1). Let D_i set t_i=0 and
T_i=p_i I+q_i D_i. For any series g with constant term one put
h_i(g)=1-D_i g/g.

Every h_i(f) has nonnegative coefficients. If J=V\\{i}, the Schur
complement gives

h_i(f)=t_i K_ii+t_i K_iJ
       (I-diag(t_J)K_JJ)^(-1) diag(t_J)K_Ji.

This property survives T_j, since direct algebra gives

h_i(T_j g)=D_j h_i(g)
 +p_j [h_i(g)-D_j h_i(g)]/[1-q_j h_j(g)].                 (60)

The numerator series are nonnegative, and the denominator has a
nonnegative geometric expansion. This includes i=j. Since
T_j g=g[1-q_j h_j(g)], it follows that

H(t)=log f(t)-log[(product_i T_i)f(t)]
    =sum_(alpha!=0) a_alpha t^alpha,  a_alpha>=0.

All expansions converge at t=1. Initially this follows from rho(K)<1
and the Schur complement, with 0<=h_i(f)(1)<1. Inductively the geometric
denominator in (60) has q_j h_j(g)(1)<1; its convergent value equals
1-D_i(T_j g)(1)/(T_j g)(1)<1. The ratio uses finite positive
combinations of reciprocal principal determinants. This also justifies
the logarithmic series and coefficient summation.

At t=1_R, f=1/d_R and (product_i T_i)f
=sum_(C subset R) p_C q_(R\\C)/d_C. Hence
H(1_R)=-log(p_R Z_R). Taking
w_S=sum_(support(alpha)=S) a_alpha proves (58).
More explicitly, for nonempty B disjoint from A, with |B|>=2,

Delta_B log Z_A=-sum_(B subset S subset A union B) w_S<=0.

No reversibility or positive direct-sink probability is required.
The exact creator verifier checks all 128 mixed differences of order
at least two on the eight maintained kernels, comparing rational
products of Z rather than rounded logarithms. These checks supplement
the general series argument; no occupation sign follows from (59) alone.


### Independent thinning of ordinary Markov loops

The series argument gives a stronger **proved construction**. Introduce
independent directed-edge variables t_ij, put B_ij=K_ij t_ij,
f=det(I-B)^(-1), and let D_i delete all outgoing variables at i.
The Schur proof of (60) applies unchanged. Deletions also preserve
coefficient positivity, since h_i(D_j g)=D_j h_i(g). Telescoping
log g-log D_i g=-log(1-h_i(g)) therefore gives

0<=_coeff log g<=_coeff log f,
g=(product_i T_i)f.                                     (61)

Convergence holds on a polydisc larger than the unit polydisc: choose
a>1 with a rho(K)<1, and repeat (60)'s convergence induction at
t_ij=a. Write ell_n=[t^n]log f and c_n=[t^n]log g for nonzero
edge-count tables n. Then 0<=c_n<=ell_n and

g(t)/g(1)=exp(sum_(n!=0) c_n(t^n-1)).                     (62)

The trace-log formula identifies ell_n with the ordinary Markov-loop
intensity aggregated by edge counts. Sample that finite-intensity
Poisson soup and independently retain every loop of count vector n
with probability c_n/ell_n. Its aggregate counts have (62)'s law.
Zero-intensity classes are absent. Every retained atom admits a
single Eulerian circuit; no disconnected count vector can occur
as an atom, because its ell_n is zero.

This reproduces the creator's full joint law of subset and balanced
words, not merely a marginal subset law. Indeed, if alpha is balanced
with support L and K-weight w(alpha), then

g_R(t)=sum_alpha p_L w(alpha)t^N(alpha).

Given the retained counts, independently shuffle each site's outgoing
multiset uniformly. The resulting word tuple has probability
p_L w(alpha)/g_R(1). Include L in C and include every site of R\\L
independently with probability p_i. Since g_R(1)=p_R Z_R/d_R,

P(C,alpha)=(d_R/Z_R) r_(R\\C) w(alpha) 1{L subset C}
         =nu_R(C) d_C w(alpha).                          (63)

Append an independent exit forest and fresh tails, prescribe the
first root bits to be 1_C, and perform the full continuation exactly
as in Section 17. The retained Poisson atoms are an auxiliary
decomposition of the aggregate counts; they need not equal the
original latent loop decomposition or the raw popped cycles.

The coefficients are local: c_n and ell_n depend only on the
vertices and edges appearing in n, and their ratio depends only on
n and the corresponding sleep probabilities, since the common
K-monomial cancels. Moreover, if u=product_(i!=j) T_i f, then

partial_(p_j) log g=h_j(u)/[1-(1-p_j)h_j(u)]>=_coeff0.

One ordinary soup with common retention uniforms therefore couples
all nested root domains and ordered sleep vectors. Retain whole loops
contained in each domain, rather than edge fragments of crossing loops.
Every additional loop for R subset T meets T\\R. Common independent
Bernoulli uniforms then give C_R subset C_T. For T=R+x, if x is absent
from C_T there is no additional loop, and C_T=C_R.

Independent priorities on outgoing edge occurrences give uniform local
shuffles in which the smaller words are subwords of the larger words.
They need not be initial prefixes. The exit forests change with the
root domain, and the bit prescribed at record one must stay there.
No order of the full continuations follows from this coupling.

The rejected loops give a positive meaning to (58):
w_S is their total intensity with vertex range exactly S, and
p_R Z_R is the probability that no rejected loop is wholly contained
in R. For positive odds, w_S>0 precisely when S is the range of a
positive-weight closed walk. Thus (59) is strict exactly when x and y
are mutually reachable inside A union {x,y}.
Separate AI audits checked both the series theorem and the joint sampler.


Even the common retained loops cannot be held fixed in a conditional
occupation comparison. Use Section 17's two half-sink roots 0,1 and
extra source 2->0, with all p=1/2. Condition on a nonempty retained
prefix on R={0,1}. It consists of m>=1 mutual cycles, forces C_R=R,
and makes the lower target-0 occupation one. There is no loop through 2,
so the nested soup adds no loop.

If the extra sleeps at 2, the root outcome is unchanged. Otherwise
surplus IDLA from 0 consumes all m mutual cycles and follows the
sampled forest to the sink. The three root forests are equiprobable.
For F_0=1,F_1=s both forest arrows are consumed; for F_0=s the
remaining F_1 is an independent ordinary K-arrow. In either case,
after averaging the forests, the root arrows are iid again.
Both roots were visited, so their remaining bits are fresh too.
The upper root continuation is therefore the ordinary full-input law,
with target-0 occupation 22/39. The conditional increment is

(1/2)(22/39-1)=-17/78,

whereas the unconditional increment is zero. The creator verifier
checks the finite-prefix absorption at m=1,2 exactly; the argument
covers every m>=1. Thus an occupation coupling must allow common
loop sectors to be exchanged or resampled, not only add new loops.

### Released inputs and negative dependence of their loss markers

In the Section 17 sampler on Q, let
E_x={x notin C, F_x=s}, where s is the original sink, and set
a_x=r_x K(x,s), a_S=product_(x in S)a_x.
For every bounded terminal-support payoff f,

E_Q[1_(E_S) f(T_Q)]
 =a_S [Z_(Q\\S)/Z_Q] mu_(Q\\S)(f), E_S=intersection_(x in S)E_x. (64)

To prove this, the joint (C,alpha,F) weight is
r_(Q\\C) w(alpha)w(F)/Z_Q. On E_S, C avoids S, so the balanced
prefix at each x in S is empty and its first record is (0,s).
Execute those records, losing precisely the corresponding inputs.
Deleting their forest arrows gives bijectively an exit forest on Q\\S,
including labelled arrows back into S. Its weight factors as
a_S Z_(Q\\S)/Z_Q times the normalized smaller sampler weight.
Every deleted root now has a fresh full-graph book and may later be
revisited. A jump to a labelled outside site instead of the original
sink would not give this identity.

In particular P_Q(E_S)=a_S Z_(Q\\S)/Z_Q, and conditioning on E_S
gives exactly mu_(Q\\S), while preserving the remaining markers.
The simultaneous-family identity (57) extends this to eta>=1_Q,
replacing the endpoint laws by mu_eta and mu_(eta-1_S).
Thus, for X_v=1{v in T_eta},

Cov(1_(E_S),X_v)=P(E_S)[q_v(eta-1_S)-q_v(eta)].           (65)

For positive K(x,s), the desired add-one sign is exactly a
nonpositive marker--occupation covariance. Markers with K(x,s)=0
are identically zero, so that particular reformulation is vacuous.

The markers themselves have pairwise negative correlation under
every positive-probability partial assignment, including prescribed
absences. Here is a consequence of the full hierarchy (58), stronger
than bare log-submodularity. For a set B of absent markers and U disjoint
from B define

Phi_B(U)=sum_(D subset B)(-1)^|D| a_D Z_(U union (B\\D))>0.

Positivity follows from its interpretation as Z_(U union B) times the
probability that no B-marker occurs; choosing C=U union B already
has positive weight. Its recursion is
Phi_(B+b)(U)=Phi_B(U+b)-a_b Phi_B(U).

This recursion preserves nonpositive log-Mobius coefficients of
orders at least two. For a positive function F with that property,
put t(U)=a_b F(U)/F(U+b)<1. Every nonconstant Mobius coefficient
of log t=log a_b-Delta_b log F is nonnegative. In the Boolean
multilinear algebra, products of nonnegative-coefficient functions
remain nonnegative; hence every t^m has nonnegative coefficients.
The convergent expansion

log[F(U+b)-a_b F(U)]
 =log F(U+b)-sum_(m>=1)t(U)^m/m

has nonpositive coefficients of orders at least two. The case a_b=0
is immediate. Induction applies this to every Phi_B.

Given G=E_A intersect intersection_(b in B) E_b^c, let
R=Q\\(A union B). Inclusion-exclusion gives
P(E_S | G)=a_S Phi_B(R\\S)/Phi_B(R).
Log-submodularity of Phi_B gives conditional pairwise negative
correlation; sequential conditioning also gives both cylinder bounds

P(E_S | G)<=product_(i in S)P(E_i | G),
P(intersection_(i in S)E_i^c | G)<=product_(i in S)P(E_i^c | G). (66)

These are marker statements, not the missing sign in (65).
For example, on the symmetric two-site half-sink kernel at p=1/2,
P(E_0)=4/13, q_1({0,1})=22/39, and
P(X_1=1 | E_0)=8/15, giving covariance -8/845.
But after also fixing C={1},
P(X_1=1 | E_0,C)=1 and P(X_1=1 | E_0^c,C)=8/15,
so the conditional covariance is positive, 14/135.
A released root may also be reoccupied:
P(X_0=1 | E_0)=2/15. These prevent an independence shortcut.
A separate AI audit verified the deletion identity and conditional
negative-dependence argument; the occupation covariance remains OPEN.


## 20. Periodic loops and the sleep bits in the balanced prefixes

**The following identities and sampler are proved.** The coefficients c_N in (62) admit
a physical interpretation, including heterogeneous sleep probabilities
and nonprimitive walks. Let W_N be the set of rooted closed vertex
words with edge-count table N, each vertex word counted once, and put
m=sum_ij N_ij. From a word w, list each site's outgoing destinations
in their order of occurrence, attach independent Bernoulli(p_i) bits,
and repeat each site's complete finite record word periodically.
Start one active particle at every visited vertex in this closed book.

If K^N=product_ij K_ij^N_ij>0, the identity is

c_N=(K^N/m)sum_(w in W_N)P_p(periodic full input stabilizes),
c_N/ell_N=average_(w in W_N)P_p(periodic full input stabilizes). (67)

The auxiliary periodic system has no sink. A successful stabilization
therefore ends with a sleeper at every visited vertex. This is not
the assertion that all periodic books stabilize.

Here are the algebraic and combinatorial ingredients of the argument.
Let H be the monoid of finite balanced local record words, under
sitewise concatenation, with record weights
K_ij t_ij p_i^b(1-p_i)^(1-b).
It is a trace monoid on bit-decorated directed simple cycles:
every balanced book is exhausted by deleting current cycles, and
two available cycles are disjoint. All deletion orders differ only
by exchanges of disjoint cycles.

Let M consist of the empty book and the balanced books whose first
record at every used vertex has bit one. Then
sum_(h in H)wt(h)=f and sum_(h in M)wt(h)=g.
The submonoid M is itself a trace monoid on its nonempty indecomposables,
with disjoint-support commutations. To check uniqueness, suppose
h=aA=bB are two factorizations with a,b indecomposable in M.
The trace-prefix diagram gives a=de, b=df, A=fz, B=ez,
where e and f have disjoint supports. Here d is a balanced left divisor
of a; e and f are balanced left divisors of B and A respectively.
Thus d,e,f lie in M. If d is nonempty, indecomposability forces
a=b=d. Otherwise a=e and b=f are disjoint. Also z belongs to M:
at each vertex at least one of e,f is empty, so its first bit is
exposed in ez or fz. Induction gives precisely the trace relations.
Every M-indecomposable has connected support.

The external trace facts used below are that a connected conjugacy
class has logarithmic weight 1/r, where r is its maximal repetition
exponent, and that Lyndon traces give a complete ordered factorization.
Both statements, for arbitrary commutation graphs, are explicit in
[Duchamp--Flouret--Luque, Sections 2.1--2.2](https://www.mat.univie.ac.at/~slc/wpapers/s54duflolu.pdf).
The graded alphabets here are finite in each bounded record degree,
so finite-degree truncation extends those facts to the M alphabet.

An H-conjugacy class meeting M meets it in exactly one M-conjugacy
class. Indeed, transpositions from h to k in M can be realized by
successively deleting balanced prefixes of h repeated periodically.
Their combined prefix u satisfies hu=uk and divides h^q for a
sufficiently large q. Thus u lies in M. To see directly that this
one-sided equality gives transpositions inside M, apply the prefix
diagram to hu=uk: h=de, u=df=fz, k=ez, with all factors in M.
If d is empty, cancellation gives z empty and h=k. Otherwise,
transposing de to ed leaves (ed)f=fk with a shorter conjugator.
Induction completes the claim. The same argument works in H.
Repetition exponents agree: if h=t^r in H and h is in M, its left
divisor t belongs to M. Connectedness agrees because every
M-indecomposable has connected support.

Moreover, the periodic book h^infinity stabilizes full input precisely
when its H-conjugacy class meets M. A completed stabilization has a
balanced emission prefix and all next bits one; deleting its cycles
moves the period to an M representative by transpositions.
Conversely, a balanced prefix exposing all-one next bits is a
stabilizing upper bound. Least action bounds every legal emission
count by it and forces termination. These cycle deletions identify
possible cuts; they are not claimed to be legal ARW moves themselves.

It follows that log g is the sum of wt(C)/r(C) over connected
H-conjugacy classes C whose periodic books stabilize. To convert
this to the rooted-walk sum (67), call a heap a pyramid if it has
one maximal cycle letter gamma. For a connected class C of record
degree m,

sum_(h in C, h a pyramid)|gamma(h)|=m/r(C).               (68)

One proof marks a cycle occurrence, with weight its record length.
Its principal ideal is a pyramid, giving E F=P F in the completed
heap algebra, where E multiplies degree and
P=sum_pyramids |gamma(h)|h. Summing over one conjugacy class is
a cyclic trace, so E log F has the same class sum as F^(-1)E F.
The class logarithm 1/r gives (68).

Finally, rooted closed words with their record bits are in bijection
with pairs (h,v), where h is a pyramid and v lies on its unique
last-arrow cycle. A rooted Euler word has exactly such a last-arrow
graph. Conversely, following the local words from v uses every
record: an unused balanced remainder would contain another
last-arrow cycle avoiding the exhausted v. Thus no edge labels
or additional choice of a root occurrence are introduced.
Stabilization is constant on H-conjugacy classes. Equation (68)
therefore makes the total stabilizing rooted-word weight m c_N/K^N,
which proves (67).

The repetition exponent is that of the decorated heap: random bits
can break a repetition of the uncoloured walk. Ordinary cyclic
re-rooting of one decorated vertex word need not preserve stabilization;
the averaging and class accounting above are essential.

The creator verifier checks (67) by exact repeated-state detection
on 5,792 decorated rooted words in seven classes. The two heterogeneous
examples at p=(1/3,2/5,3/7) give 1504/11025 on the complete directed
triangle and 6152/33075 on the doubled bowtie. Two separate AI audits checked the structural proof, rooting factors,
and its interface with the creator law.

### A canonical sampler retaining the paired bits and arrows

Fix one ordered alphabet of M-indecomposables on a finite ambient
domain, and its standard Lyndon complete factorization. Let L range
over its representatives of primitive connected conjugacy classes,
and put a_L=wt(L). The factor order is the theorem's Lyndon order;
an arbitrary ordering of arbitrary class representatives is insufficient.

Start with independent rooted Markov-loop atoms, each word w having
intensity K(w)/|w|. Attach independent sleep bits to every occurrence
and retain exactly the periodically stabilizing books. Mapping an
accepted decorated pyramid to its H-conjugacy class gives class
intensity wt(C)/r(C), by (68). The unique M-class identification
and preserved repetition exponent therefore give independent counts

A_(L,k)~Poisson(a_L^k/k), k>=1.
N_L=sum_(k>=1)k A_(L,k) has P(N_L=n)=(1-a_L)a_L^n.        (69)

Form the canonical nonincreasing product h=product_L L^N_L.
Complete factorization and (69) give

P(h)=wt(h)/g_R(1), h in M_R,
g_R(1)=product_(L supported in R)(1-a_L)^(-1).             (70)

Only finitely many accepted atoms occur almost surely, since their
total intensity is log g_R(1)<=-log d_R<infinity. Thus the product
is finite; the countable alphabet causes no infinite sorting issue.

This operation rearranges paired records by conjugacy and canonical
factorization; it does not simply copy the accepted rooted word.
For example, a one-vertex accepted periodic word 01 must be moved
to an M representative such as 10 before it is a valid used prefix.
Given the arrow projection alpha of h, (70) forces each first used
bit to one and leaves every later prefix bit independently Bernoulli(p).
This follows by dividing its product record weight by p_support(alpha).

Independently sample the exit forest F and Bernoulli(p_i) bits B_i,
then append the record (B_i,F_i) and fresh tails at every root.
Put C=support(h) union {i:B_i=1}. At used sites B_i is a fresh later
bit; at unused sites it is precisely the first bit 1_C(i).
This gives exactly Section 17's full creator-book law, including all
prefix bits, and requires no resampling of the accepted records' bits.

For R subset T, the M-indecomposables supported in R are precisely
the restrictions of the ambient alphabet: a factorization cannot
introduce records at a new vertex. Restricting the single global
Lyndon factorization therefore remains complete. Under one decorated
loop soup, acceptance and primitive class labels are unchanged for
loops contained in R. Deleting from h_T every whole Lyndon block
whose support is not contained in R leaves exactly h_R.

Consequently the smaller local record words are subsequences of the
larger ones, preserving bits and arrows together. Common B_i also
give C_R subset C_T. This is a coupling of balanced record prefixes,
not of initial prefixes of the completed books or their endpoints.
The forests still change with the domain, and the negative common-loop
example in Section 19 still applies. A second independent AI audit
checked this sampler, including arbitrary commutation graphs, powers,
infinite graded alphabets, and domain consistency.



### A coherent coupling by initial record prefixes

There is a stronger **proved coupling**, different from restriction of the
Lyndon soup. For T=R union {x}, call an M_T-indecomposable old if its
support is contained in R, and new otherwise. In the finite heap of
h in M_T, take every occurrence all of whose predecessors, including
itself, are old. This is the greatest all-old order ideal. Its trace
alpha is the greatest M_R left divisor of h, and h=alpha beta uniquely.

Every minimal occurrence of beta is new. Conversely, that condition
characterizes the decomposition: every occurrence of beta has a new
ancestor and cannot belong to an all-old initial ideal. Since every new
letter contains x, they are mutually dependent. Thus a nonempty beta
has exactly one minimal M-letter, which contains x.

Let Q_(R,x) be the set of these beta, including the empty word. The
weight-preserving bijection M_R times Q_(R,x) -> M_T gives

sum_(beta in Q_(R,x)) wt(beta)=g_T(1)/g_R(1),
P(alpha=a,beta=b)=[wt(a)/g_R(1)] [g_R(1)wt(b)/g_T(1)].     (70a)

In particular the two factors are independent. This argument is on
finite heaps and needs neither an alphabet ordering nor an infinite
sorting convention. Nonnegative summation and finiteness of g_T
justify the probabilistic normalization.

Write pi_R(h) for the greatest old left divisor. If R subset S, then
pi_R(pi_S(h))=pi_R(h): every R-left divisor is an S-left divisor.
Consequently one ambient M-word gives simultaneous correct marginals,
with h_R an actual initial record prefix of h_S at every site,
preserving the paired bits and arrows. This projection can regroup the
common Lyndon factors; it is not deletion of the new factors from the
ordered soup. The separately appended forest records still sit at
different boundaries, so the completed books are not nested prefixes.

The distinction from ordinary forest promotion is substantive.
Take R={a} and append to the record (1,a) at a the cycle records
(0,x) at a and (1,a) at x. The complete word is an M-indecomposable
and belongs to Q_(R,x), although its ordinary first-arrow cycle is
the self-loop at a. Its ordinary H-decomposition is a nonempty old
prefix followed by an x-cycle whose first a-bit is zero.
In general beta=delta omega under ordinary demotion, with
delta in M_R and omega's arrows in W_x from Section 21;
omega need not belong to M. Thus h=(alpha delta)omega describes
precisely how ordinary promotion absorbs old stable records.

Marking each paired record type by a commuting variable gives an
additional identity. The inventory generating function of beta is
[g_R(1)/g_T(1)]g_T(z)/g_R(z). By the decorated-class logarithm, this
is exactly the inventory law of the accepted periodic-loop soup
restricted to loops visiting x. It does not identify its ordered word
with raw popped cycles. If no positive loop visits x, beta is empty
almost surely, and Section 19's common-loop obstruction remains.

An independent AI derivation checked the heap factorization and its
coherence. The root verifier enumerates 334 paired balanced words through
degree four, checks 9,018 nested projections, and checks the graded
partition-function quotient at heterogeneous sleep probabilities.

## 21. Forest promotion and resampling common creator blocks

**The following kernels are proved.** Put T=R union {x}. Let W_x
contain the empty local word tuple and every nonempty balanced tuple
beta whose support S contains x and whose first arrows on S\\{x}
form a directed tree toward x. Including the first arrow at x,
the first-arrow graph on S has exactly one cycle, containing x.
For an exit forest F on T define H_R(beta,F) by taking the first
beta-arrow at each used old root, and F otherwise.

Pop iid arrow cycles inside R first, obtaining alpha, then continue
popping inside T, obtaining beta and F_T. The exact cylinder law is

P(alpha_R=alpha,beta=b,F_T=F)=w(alpha)w(b)w(F),             (71)

for alpha balanced on R, b in W_x, and F acyclic on T, and zero
otherwise. Indeed, alpha beta F is a complete T-popping cylinder.
Its first-stage stopping condition is that H_R(beta,F) is acyclic
on R. Balance keeps all first beta-arrows inside S; therefore this
condition is exactly beta in W_x. A cycle avoiding S would already
be a cycle of F. There is no popping-order multiplicity.

Consequently alpha, beta, and F_T are independent, with laws

d_R w(alpha), (d_T/d_R)w(beta), and w(F_T)/d_T, respectively,
F_R=H_R(beta,F_T).                                      (72)

Independence applies to the entire word beta, not merely its counts.

There is a concrete sampler for beta. Reverse each local outgoing
word, without reversing any arrow. Its last arrows away from x now
form a tree toward x, so following these reversed words from x
uses every arrow and returns to x. Balance prevents getting stuck
elsewhere; an unused balanced remainder would be closed under the
last-arrow tree while excluding x, which is impossible.
Conversely, a closed walk based at x has last-departure arrows
forming a tree toward x. Reversing its local departure lists gives
a unique element of W_x. Thus

sum_(beta in W_x) w(beta)t^N(beta)
 =(I-K_T(t))^(-1)(x,x)
 =det(I-K_R(t))/det(I-K_T(t)).                            (73)

Run an ordinary K-walk from x until its first exit from T, keep the
segment through its last visit to x, and reverse its local departure
lists. A fixed closed prefix omega is selected with probability
w(omega) sigma, where
sigma=P_x(exit T before a positive-time return to x)=d_T/d_R.
This proves the sampler directly; the last-return time is not
misused as a stopping time. Empty words, self-loops, proper visited
supports, and distinct labelled outside exits are all allowed.

Equation (73), after normalization, is also the count generating
function of ordinary Markov loops contained in T and visiting x.
It does not identify the raw second-stage popped cycles with a
Poisson family. Given F_T, draw beta by (72) independently to recover
F_R. Conversely, for a specified positive-weight lower forest f,

Pi_f(beta,F)=[w(beta)w(F)/w(f)]
  1{beta in W_x, F acyclic on T, H_R(beta,F)=f}            (74)

is normalized and is exactly continued cycle popping from f with
fresh suffixes and a fresh x-stack.

### Promotion in the positive creator law

Write the creator's joint subset/arrow-prefix/forest weight as
S_A(C,gamma,F)=r_(A\\C)w(gamma)w(F)/Z_A, with
support(gamma) subset C and gamma balanced.
Draw an R-creator state (C,alpha,f), promote f by (74), and draw a
fresh Bernoulli(p_x) bit. Let C^+ add x precisely when this bit is one.
Accept exactly when support(beta) subset C^+, and output
(C^+,gamma=alpha beta,F). The accepted weight is

[p_x/Z_R] r_(T\\C^+)w(alpha)w(beta)w(F).

Staged cycle popping gives a unique alpha beta decomposition, hence

P(accept)=p_x Z_T/Z_R in [d_T/d_R,1],
Law(output | accept)=S_T.                               (75)

An accepted proposal biases its common R-state; it cannot be used
as an unconditioned lower sample. Keeping the first lower proposal
and redrawing complete proposals until acceptance gives correct
marginals, but no occupation comparison.

There are also exact updates which exchange common loops. Decompose
an upper creator prefix gamma=alpha beta by staged popping and set
S=support(beta) intersect R. Then

Law(C intersect R,alpha | beta,F)
 =Law_R(C_R,alpha_R | S subset C_R).                      (76)

This concerns the latent pair, not the demoted forest's law.
It can be sampled by setting latent probabilities p_i^S=1 on S
and p_i^S=p_i elsewhere, drawing the Section 19 thinned-loop arrow
prefix with these parameters, and including S and its support in C.
Other old vertices have independent original Bernoulli(p_i) bits.
Include x certainly if beta is nonempty, and by Bernoulli(p_x)
otherwise. The physical later sleep bits always retain the original p:
the modified parameters are only for this conditional latent sampler.

Conversely, conditional on C,alpha,F, beta is empty if x notin C;
otherwise its law is

P(beta | C,alpha,F)=[d_C/d_(C\\{x})]w(beta),
 beta in W_x, support(beta) subset C.                    (77)

This is ordinary promotion on C\\{x} subset C. The empty-word
probability is at least d_T/d_R>0. Equations (76)--(77) are normalized
alternating resampling steps; beta empty permits a complete refresh
of the common latent pair. They remain valid in Section 19's example
where no loop can visit the added source, so ordinary promotion alone
leaves the losing common-loop sector fixed. Normalization alone
does not provide the missing target sign.

These cylinder, last-return, acceptance, and conditional-law arguments
received an independent AI audit. The remaining task is a joint
record/forest exchange using these kernels with controlled occupation.


### A surplus walk and the exact remaining books

The following **proved reset kernel** retains the entire consumed walk
and its labelled exit. Begin with the R-creator book, equivalently an iid
book conditioned on A_R, and inject a surplus at z in R while holding
one active particle at every root. Run only that surplus until it exits R.
The source z is fixed independently of the creator book; independent
randomization of z is also allowed.

Let pi be the path from z to its outside destination in the sampled
forest F_R, P its vertices inside R, and U=R\\P. The surplus consumes
precisely the forest arrows on pi and a balanced prefix epsilon of alpha.
There is a unique decomposition alpha=epsilon delta, where delta is
balanced and supported in U, and the last epsilon-arrows on U are
acyclic. No fresh arrow after F_R is used. The walk's chronological
loop-erasure is pi, and its visited set is D=P union support(epsilon).

To prove this, on reaching a successive forest-path vertex, follow the
remaining balanced arrows until that vertex's queue is exhausted.
Balance prevents exhaustion elsewhere, so this is a closed excursion
leaving a balanced remainder. Then take its forest arrow. An exhausted
vertex has no balanced incoming arrows remaining, so the forest path
is never revisited. Last-departure times exclude a last-arrow cycle
avoiding P. Conversely, these conditions reconstruct all the excursions
and the path, proving uniqueness. A visited vertex outside P may retain
a nonempty delta-word.

For A subset U define q_i=1-p_i and

Phi(U,A)=d_U sum_(B subset A)
             p_B q_(A\\B)/d_((U\\A) union B).             (78)

This is the iid probability that the U-cycle-popping range meets A only
where the current first bit is one. Hence Phi(U,empty)=1 and
Phi(R,R)=p_R Z_R. For a fixed finite surplus walk omega, let
A=R\\D be its unvisited set. Its creator-law probability is

P(omega)=w(omega) p_support(epsilon)
                     Phi(U,A)/Phi(R,R).                 (79)

Conditional on omega, the remaining books at P are iid, and those at U
are iid conditioned on

L_U intersect A subset {i in A: current first bit at i is one}. (80)

All remaining bits at visited vertices are fresh. At unvisited vertices
the original first-bit flags are retained; every later bit is fresh.
Treating P as absorbing in (80) is only an auxiliary cycle-popping test:
physical continuation still uses every root in R and permits returns.

Here is a cylinder proof and explicit conditional sampler. Write
F_R=pi union f, where f is an exit forest on U. The original creator
cylinder has weight

[p_C q_(R\\C)/Phi(R,R)]
                 w(pi)w(epsilon)w(delta)w(f),

with support(epsilon) union support(delta) subset C. Summing its consumed
first bits forces those on support(epsilon) to one and leaves the bits on
P\\support(epsilon) independently Bernoulli(p). The surviving flags
C_A=C intersect A have conditional distribution

P(C_A=B | omega)
 =p_B q_(A\\B) d_U/[d_J Phi(U,A)],
 J=(U\\A) union B.                                      (81)

Given B, independently sample delta with law d_J w(delta), supported
in J, and f with law w(f)/d_U. Append delta,f and fresh tails, retaining
first bits 1_B on A and fresh bits elsewhere. The cylinder factorization
proves (79)--(81), including independence of the fresh books on P.

Thus an ordinary K-walk from z until exit R, accepted with probability
p_support(epsilon) Phi(U,A), samples omega; its total acceptance
probability is Phi(R,R). The residual kernel depends on the sets P,D
only, even when the entire route and exit label are retained.
The exit label itself has the ordinary harmonic-measure marginal
[(I-K_RR)^(-1)K_(R,e)]_z: it follows the independently sampled forest
path, whose marginal is the loop-erased ordinary walk.

There is also a complete refresh of the residual arrow forest.
Pop U first, removing delta, then pop inside R. The additional balanced
word beta has first-arrow graph acyclic on U and law
(d_R/d_U)w(beta). Conditional on the consumed walk, it is independent
of delta and the surviving flags. The final forest F' has law w(F')/d_R,
independently of all of these and of the consumed walk. This follows
from the promotion cylinder w(beta)w(F')/d_U. The new prefix delta beta
need not have first bit one at every used vertex.

Two cases give an entire iid reset, including all remaining bits.
If D=R, the set A is empty. If |R|=2 and there are no self-loops at
the roots, U has at most one vertex and no cycle, so (80) is vacuous
for every consumed walk. In both cases the whole remaining root book
is iid, independently of the walk and its exit label.
For general R a reset cannot be asserted: with a->sink surely,
K_bb=K_b,sink=1/2 and p_b=1/2, a surplus from a leaves b's first bit
equal to one with probability 2/3, not 1/2.

The no-self-loop hypothesis can be arranged at the level of endpoint
laws. Delete self-jumps and replace
K_ij by K_ij/(1-K_ii) for j!=i (including the sink), and replace
p_i by p_i/[1-(1-p_i)K_ii]. For a lone active, sum its geometric
self-jumps; with surplus, discard the self-jumps that do not change
the configuration. Both yield these same effective parameters.
This is a change of model preserving stabilization laws, not an
identification of the original and transformed record books.

The root verifier checks the path, balanced decomposition, and residual
U-prefix on 400 seeded books. The creator verifier checks (79)'s total
normalization, labelled exit marginals, and (81) on 150 terminal walk
classes and 198 flag weights, using exact rational arithmetic.
A separate AI audit checked the conditional laws and the forest refresh.
They condition on the consumed walk, not additionally on the original
forest F_R. The regenerated forest F' is a different random object;
its path from z need not equal pi.
These identities do not give a conditional occupation sign; in particular
they do not contradict the negative retained-loop sector in Section 19.


### Root reset does not erase the lower outside path

The reset still leaves an essential **unread path in the lower copy**.
Let eta be the extra particle's path outside R up to its first entrance
at z in R, before sleeping, and let mu_R^eta be the lower continuation
with those zero-bit records retained. For two roots without self-loops,
averaging every root loop sector and forest gives the exact contact gain

Delta_v(eta)=sum_e H_R(z,e)mu_(R+e)(v)-mu_R^eta(v),        (81a)

where H_R is ordinary labelled harmonic measure and R+sink means R.
The upper's consumed outside records and root records leave fresh
suffixes; the lower has consumed neither. All later outside returns
remain in the continuations. Replacing mu_R^eta by mu_R is unjustified.

Even this fully averaged contact gain can be negative. At p=1/2 take
a->b,x equally, b->a surely, and x->a,sink equally. For roots {a,b},
target a, and contact record eta=(0,a) at x,

mu_R^eta(a)=10382/12675,
sum_e H_R(a,e)mu_(R+e)(a)=70459/93795,
Delta_a(eta)=-10613/156325.

The common empty/nonempty retained-loop sectors have weights 4/5,1/5;
both are already included. The contact event has probability 1/4.
The probability-weighted no-contact credit is 120119/1875900, and
the actual unconditional gain is positive, 4414/93795. The creator
verifier independently checks the fixed-record law, all three initial
outside branches, and their weighted sum by exact raw absorption.

Thus the remaining reset-based inequality must retain credit from
no-contact histories as well as the lower unread outside paths. This
is not a counterexample to monotonicity. If every exit from R goes
straight to the original sink, contact gain is zero and the extra's
first-sleep law outside R gives the nonnegative total increment, as
already covered by the terminal-region argument.

## 22. Addition characters and a full-input activation field

**The algebraic formulas below are proved. Their coefficient sign is
a new OPEN candidate.** Continue to write F_i=A_i-J_i, a_R=product_R A_i,
and r_i=(1-p_i)/p_i. Initially assume K(i,s)>0 at every vertex.
For C subset V define ordinary-walk probabilities

h_C(z)=P_z(original sink before visiting C),
j_C(i)=K(i,s)+sum_(z outside C) K(i,z)h_C(z), i in C.

Thus h_C=0 on C and is strictly positive elsewhere, and j_C(i)>0.
These h_C differ from Section 13's walks killed by sleeping.
Put R=V\\C and

U_C=a_R product_(i in C)F_i,
lambda_C=(-1)^|C| product_(i in R)h_C(i) product_(i in C)j_C(i).
                                                               (82)

Then U_C A_z=h_C(z)U_C, U_C U_D=0 for C!=D, and
U_C^2=lambda_C U_C. Indeed, A_iF_i=0 kills additions at C.
At R the relation a_R A_z=a_R J_z, followed by inversion of
I-K_RR, gives the harmonic eigenvalues h_C. Orthogonality follows
from an A_iF_i factor, and squaring gives lambda_C.

The normalized U_C/lambda_C form a complete family of projectors.
For a proof without assuming diagonalizability, let D_z be the
surplus-IDLA addition matrix on binary supports. The saved reduction
gives P A_z=D_z P. Ordered by cardinality, D_z is upper triangular;
its diagonal at support V\\C is h_C(z). Thus U_C is nonzero, having
diagonal lambda_C in this similar representation. There are 2^|V|
nonzero orthogonal projectors on a space of that dimension, proving

sum_C U_C/lambda_C=I.                                    (83)

For a terminal target set T, let f_T(B)=1{T subset B}, and define

b_T(C)=delta_empty U_C f_T,
W_T(t)=delta_empty product_i(A_i+t_iF_i) f_T
      =sum_C b_T(C)t_C.                                  (84)

The target occupation coefficient is b_{ {v} }(C).
For 0<=t_i<=r_i, the factor A_i+t_iF_i is
(1-t_i/r_i)A_i+(t_i/r_i)B_i. Expanding and applying the commuting
sleep-insertion factors first proves that W_T(t) is the probability
T is occupied when every vertex starts with one particle, independently
left asleep with probability t_i/r_i and otherwise active.

In particular, with m_T(S) denoting that probability when precisely
S is initially asleep,

Theta_T(C):=r_C b_T(C)
 =sum_(S subset C)(-1)^|C\\S| m_T(S),
sum_C Theta_T(C)=1.                                      (85)

All identities (84)--(85) remain valid with only sink accessibility;
they do not divide by possibly zero lambda_C.

The same algebra gives a surplus identity without using (83).
For extra all-active input eta, initially asleep set S, and one initial
particle at every vertex, let m_T(S;eta) denote the terminal probability.
Since U_C A_z=h_C(z)U_C, including degenerate characters,

m_T(S;eta)=sum_(C subset S) Theta_T(C)
                             product_z h_C(z)^eta(z).    (85a)

There is also a direct surplus-IDLA proof: hold the initial particle
at each vertex and move all surplus particles to the sink. Their
independent ordinary walks erase precisely the initial sleep labels
at visited sites and leave fresh remaining books. Thus
m_T(S;eta)=E[m_T(S without the extra walks' range)].
This does not assert that m_T is increasing in S.

The **general candidate sign is false**, already on three sites, by the
counterexample below. On any class where Theta_T(C)>=0 for all C,T is
proved, it gives a probability distribution on an auxiliary set H_T such that

m_T(S)=P(H_T subset S).                                  (86)

Then (85a) would further identify m_T(S;eta) with the probability
that H_T is contained in S and avoided by independent ordinary walks
from eta. The random set would be independent of those walks.

Such a stochastic backward construction is established below for full
occupation and for graphs with strongly connected components of size
at most two. A proper-target coefficient can be negative on a three-site
strongly connected component, so that component-size theorem is sharp.
This concerns a full initial particle configuration with varying
activation labels, not arbitrary sleeping-background attractiveness.
It also does not automatically prove add-one monotonicity.
For example, (83) gives for an all-active eta

q_v(eta)=sum_C [b_{ {v} }(C)/lambda_C]
                         product_z h_C(z)^eta(z),        (87)

with 0^0=1. Even nonnegative b would leave alternating contributions
in this expression. Zero direct-sink probabilities can coalesce the
characters; (87) must not then be used by dividing by zero.

There is a precise first consequence if the singleton-C case is proved.
For C={x}, put sigma_x=j_{ {x} }(x)>0. Eliminating V\\{x} in the
addition algebra gives

b_{ {v} }({x})=sigma_x[
       q_v(1_V)-q_v(1_(V\\{x}))].                        (88)

Thus that case is exactly the full-minus-one comparison. On three
vertices it supplies the binary N2 comparison. The arbitrary-graph off-diagonal sign
remains OPEN; Section 25 proves the three-site case by an exact positive
certificate. The added target itself is treated below.


### Two positive cases and the remaining singleton cancellation

For T=V, the entire dual row is **proved positive** and already present
in Section 15's latent law. Let Z_V(t)=d_V sum_C t_C/d_(V\C), with
the physical odds r fixed. Since P_r(A,V)=0 unless A=V and
P_r(V,V)=1/Z_V(r), the parameter-transfer identity gives

W_V(t)=Z_V(t)/Z_V(r),
Theta_V(C)=r_C d_V/[d_(V\C) Z_V(r)]=nu_V(V\C).            (89)

Thus the dual set for full occupation is the complement of the latent
subset, on arbitrary directed graphs. Formula (85a) already supplies
its positive walk-avoidance construction with incoming particles.

For the singleton C={x}, let mu=mu_(V\x), pi=mu_V, and let nu_x(B)
be the unrenormalized creation law from active V\x on histories never
visiting x. Holding the particle at x gives the exact vector identity

r_x sigma_x(pi-mu)
 =sum_(B subset V\x) nu_x(B)[delta_(B+x)-delta_B A_x].     (90)

On a first contact with x, its active/sleeping alternatives cancel
because both now have two active particles there. No-contact histories
give the right side. The left side is
mu(B_x-A_x)=r_x mu F_x=r_x sigma_x(pi-mu).

For target x every summand in (90) is nonnegative: it is the probability
that adding active x to sleeping B leaves x vacant. Hence

b_{ {x} }({x})>=K(x,s) product_i p_i.                     (91)

In fact b_{ {x} }({x})>0 under mere sink accessibility. First let every
other initial particle sleep immediately. Choose a simple positive
x-to-sink path, let x jump, and route every particle awakened on that
path forward along it. A finite positive-probability cylinder empties x
without returning there. Together with sigma_x>0, this proves strict
full-minus-one monotonicity at the added site.

For v!=x the right side of (90) has both signs. On the two-site
half-sink example, the B=empty summand for target v!=x is -1/30.
The averaging over the no-visit creation law is therefore essential.

That law has an exact positive description when h_i=h_{ {x} }(i)>0
for every i!=x. On R=V\x define

K^h(i,j)=K(i,j)h_j/h_i,  K^h(i,s)=K(i,s)/h_i,

retaining the original p_i, and let pi_R^h be its full-input stationary
law. Then

nu_x(B)=[h_R/h_B] pi_R^h(B),  B subset R.                 (92)

The likelihood ratio of a no-x history telescopes over particle jumps:
the transformed history weight is the old weight times h_B/h_R.
The transformed walk is stochastic and sink-accessible.
Also nu_x Ahat_y^x=h_y nu_x for y in R, by the deletion intertwining
and mu A_y=h_y mu+(1-h_y)pi. Formula (92) does not assert occupation
order under this tilt. The eigenmeasure identity remains valid when
some h_i=0; the displayed Doob representation requires positive h.
The exact verifier checks (90) in 122 terminal-law
coordinates, (92) in 45 entries, and strict positivity in all 13 archived
source choices.

### Full occupation and labelled boundary outputs

The full-target dual also describes all exit labels. Kill on a labelled
absorbing boundary E, and choose weights t_e in [0,1]. Define

h_C^t(z)=E_z[t_(exit); exit before visiting C].

For initial one particle everywhere, initially asleep precisely on S,
and extra all-active input eta, the identity is

E[1{final support=V} product_e t_e^(exit count at e)]
 =sum_(C subset S) Theta_V(C) product_z h_C^t(z)^eta(z).   (93)

Here Theta_V is (89), independent of the boundary weights.
For a proof with h=h_empty^t>0, perform the Doob transform with
internal probabilities K(i,j)h_j/h_i and exit probabilities
K(i,e)t_e/h_i. Internal principal determinants, hence Z and Theta_V,
are unchanged by diagonal similarity. The weighted history ratio is
h^(1+eta)/h^(final); on full occupation this is h^eta.
Apply (85a), noting that h(z) times the transformed avoidance probability
equals h_C^t(z). Boundary weights which give zero h follow by continuity.

Equivalently, the subprobability law of exits together with full
occupation is a mixture of independent ordinary walks from eta,
required to avoid the sampled C and with their exit labels retained.
These are exit-count multisets, not the identities or order of physical
ARW exits. Conditional on full occupation, the C weights are proportional
to 1{C subset S} Theta_V(C) product_z h_C^1(z)^eta(z); the walks in
each positive component are conditioned to avoid C.
This is an identity for a fixed incoming configuration, not an assertion
that adaptive exterior feedback can be substituted without accounting.
Exact weighted raw absorption checks 168 mixed-input/incoming comparisons
on two- and three-site kernels with two exit labels, including zero
boundary weights. A separate AI audit verified (89), (92), and (93),
including zero-weight limits and the conditional mixture weights.

**Evidence and limits.** All 122 one-target coefficients, and all 404
joint-target coefficients including empty targets, are nonnegative on
the five archived kernels. The projector equations were audited in
19,592 matrix entries; three kernels have complete nondegenerate
projectors, while the other two have 14 zero-lambda sectors.
A separate exact screen found no negative one-target coefficient on
25 three-site and five four-site kernels, seed 2026101612.
Direct raw absorption from every mixed full initial state on one
five-site kernel, seed 2026101615, independently checks all 160
coefficients in (85) against the projector computation and finds no
negative coefficient. These are finite evidence only.

The ARW restriction is substantive. For the general two-emission
sandpile routing Q displayed in Section 14, define the active insertion
by the deposition matrix S=(I+Q)/2 and a sleeping insertion by direct
deposition. Its coefficient Theta_{ {2} }({1}) equals
-9092/4124577. Thus commutation and the analogous quadratic addition
identity alone cannot establish this candidate; the nonnegative
original walk kernel, equivalently the relevant potential structure,
must enter.


### A three-site obstruction to general field positivity

**Exact counterexample, not an occupation-monotonicity counterexample.**
Take the irreducible kernel, with the last column the sink,

K=[[0,0,1/2,1/2],
   [1/2,0,0,1/2],
   [15/16,3/64,0,1/64]],
p_0=p_1=p_2=1/1001.

For target 0 and marked set {1,2}, exact absorption from the four full
mixed inputs gives

Theta_{ {0} }({1,2})
 =m_0({1,2})-m_0({1})-m_0({2})+m_0(empty)
 =-4885292907192310219717609728515625
   /21076002824735985616104392981206548077528 < 0.

All direct-sink probabilities are positive. Thus neither zero boundary
rates nor a reducible graph explains the failure. In particular there
is no probability distribution H_0 with m_0(S)=P(H_0 subset S) for this
kernel. The full-target law and the size-two SCC theorem remain valid.
The singleton marked coefficients on three sites remain nonnegative by
Section 25; it is the proposed positivity of every higher coefficient
that fails.

The maintained `proper_target_field_obstruction` audit computes the
four laws directly with exact rational ARW absorption. A separate
compact calculation using Section 25's one-particle and pair resolvents
gives the same number. For this fixed kernel and common sleep probability
p, that calculation also gives
Theta_{ {0} }({1,2})=-(555/1157632)p+O(p^2) as p decreases to zero.
The exact finite counterexample alone suffices to refute the general
field-positivity proposal.


The failure has an explicit first-order mechanism. For loopless K put
G=(I-K)^(-1), g=G_00, sigma=1/g, and, for {i,j}={1,2},

h_i=(k_i+K_ij k_j)/(1-K_12 K_21),
u_i=(K_0i+K_0j K_ji)h_i(1-h_i),
u_12=sigma-k_0-sum_(i=1,2)K_0i[
         k_i^2+h_j(h_i^2-k_i^2)],
chi_i=k_i+k_0 K_i0/(1-K_0i K_i0).

Here h_i is the ordinary sink-before-0 probability. The quantity u_i
is the increase in target-hitting probability after the forced departure
from 0 when dormant resident i is added at zero sleep probabilities;
u_12 is the increase with both residents. Thus 0<=u_i<=u_12.
For the latter formula, after first jumping to i the two active walks
both avoid j and 0 with probability k_i^2. On the remaining avoidance
mass h_i^2-k_i^2, resident j is activated and contributes h_j.

For p_i=epsilon rho_i, with fixed rho_i>=0,

Theta_{ {0} }({1,2})=epsilon rho_0 L+O(epsilon^2),
L=g[u_12-chi_1u_2-chi_2u_1].                            (88a)

Here is a derivation from Section 25's pair renewal, independent of the
large scalar positivity certificate. First take common odds t=(1-p)/p,
write d=det(I-K), d_i=det(I-K_(V without i)), D=diag(d_i), and
kappa_i=k_i+K_i0 k_0. The exact renewal matrix is

N=t^2 D(I-K)+t[3D-dG]+I.

Let beta be the full-input pair vector and gamma^i the vector from the
input missing i. Set
y=d beta-d_1 kappa_1 gamma^1-d_2 kappa_2 gamma^2.
Since N gamma^i=e_i and N beta=b/Z_V, expansion gives

y=t^(-2)k_0 e_0
  +t^(-3)[-k+k_0(G diag(G_ii)^(-1)G-3G)e_0]
  +O(t^(-4)).

Indeed, N_2^(-1)=G D^(-1), b/Z_V=Dk/d+
t^(-1)D(K-I)k/d+O(t^(-2)), and
k-kappa_1e_1-kappa_2e_2=k_0(I-K)e_0.
For i=1,2 the killed-Green identity gives
y_i=-chi_i t^(-3)+O(t^(-4)). The single-particle terms in the two-mark
field cancel because
d-d_1 kappa_1-d_2 kappa_2+kappa_1 kappa_2=0.
The remaining occupation coefficients tend to g u_2, g u_1 and g u_12;
also d/Z_V=t^(-3)+O(t^(-4)). This proves (88a) on the common-p diagonal.
The finite absorption law is analytic near p=0, and every target-0
moment here vanishes identically when p_0=0. Hence all other first
derivatives vanish, proving the stated heterogeneous expansion.

Each chi_i lies in [0,1]: chi_i<=k_i+K_i0<=1. Their sum need not be at
most one. In fact

L/g=(1-chi_1-chi_2)u_12
     +chi_1(u_12-u_2)+chi_2(u_12-u_1).

Thus a negative slope requires chi_1+chi_2>1; the two escape weights
cannot be counted as disjoint portions of one probability budget.

There is a sharper necessary support condition:

L<0 implies
K_10 K_20 k_0(K_01 K_12 k_2+K_02 K_21 k_1)>0.             (88b)

The small signed polynomial audit clears the positive denominator
d_0^2 d_1 d_2 d. Its 231-term numerator has an exactly invertible
row-homogeneous expansion with 2,700 terms; all 58 negative terms require
the displayed edges and sink factors. This proves (88b), including
boundary kernels. A simple directed three-cycle therefore has
nonnegative rare-sleep slope. The minimum irreducible support permitting
a negative slope is a three-cycle with a shortcut back to the target.

For K_02=K_10=1/2, K_20=E/(E+4), K_21=3/(E+4), E>=0,

L(E)=(-2E^2+77E+360)/[8(E+4)(E+8)(2E+13)].

It is negative exactly when E>(77+sqrt(8809))/4; E=60 gives the saved
counterexample and L=-555/1157632. The maintained
verify_three_site_rare_sleep.py checks the small polynomial identities,
the family formula and the support certificate using SymPy, without
parameter sampling or rebuilding the large occupation certificate.


## 23. Backward sets on acyclic graphs and two-site components

**The following dual is proved for a directed acyclic internal graph.**
Self-loops may first be eliminated as in Section 21. There is one base
particle at every vertex, independently asleep with probability s_i and
otherwise active; extra particles eta are active. Include every killing
destination among the labelled absorbing exits.

For target set T and exit weights t_e in [0,1], let F_(U,eta)(s;t)
be the joint expectation of final occupation of every target in U and
the exit-count weight. Eliminate vertices in reverse topological order.
At a sinkward vertex u, put ell=sum_e K(u,e)t_e, where successors
already eliminated are now weighted boundary vertices, and q_u=1-p_u.
Define

(alpha_u,beta_u)=
  (p_u+q_u ell, q_u(1-ell)), u notin T;
  (p_u, q_u),               u in T.                     (94)

Then, writing F_- for the upstream problem after deleting u,

F_(U,eta)
 =alpha_u ell^eta(u) F_-(t_u=ell)
  +s_u beta_u 1{eta(u)=0} F_-(t_u=0).                    (95)

For the proof, process the physical system in forward topological order.
Its incoming number N at u is determined by upstream data, independently
of u's local book and initial flag. Except for an untouched initial
sleeper with N+eta(u)=0, the site has N+eta(u)+1 active particles.
It makes N+eta(u) forced jumps and uses a fresh final sleep coin.
That coin is independent of the emitted destinations. For u notin T
the local joint weight is

(p_u+q_u ell)ell^(N+eta(u))
 +s_u q_u(1-ell)1{N+eta(u)=0};

for u in T it is

p_u ell^(N+eta(u))+s_u q_u1{N+eta(u)=0}.

This proves (95), including the case ell=0 with 0^0=1.

In the no-extra recursion, alpha_u,beta_u are nonnegative and sum to
one. Choose these two branches backward, mark u on the second branch,
and assign its boundary weight zero; on the first assign ell.
Let H be the resulting marked set. The construction is independent of
the initial flags s, and gives

F_(U,0)(s;t)=E s_H=sum_C Theta_(T,t)(C)s_C,
Theta_(T,t)(C)=P(H=C)>=0.                                (96)

Given the final marked set C, backward harmonic substitution identifies
the assigned boundary weight at each vertex with h_C^t, the weighted
ordinary-walk exit probability before hitting C. Hence the same
normalized branch law, for every eta, gives

F_(U,eta)(s;t)=sum_C Theta_(T,t)(C)s_C
                              product_i h_C^t(i)^eta(i). (97)

The extra factors are killing weights, not a reason to renormalize
the branch law. At t=1, this proves the general-field candidate in
Section 22 for DAGs, for every joint target set and every extra input.

There is a single graphical realization for all target sets when t=1.
Independently choose A_i~Bernoulli(q_i), and an ordinary K-walk W_i
from each vertex to the absorbing boundary. On the vertices with A_i=1,
draw an auxiliary arrow i->j whenever W_i visits j after time zero.
Let H_T consist of the ancestors of T intersect {i:A_i=1} in this graph,
including those target vertices themselves.

Expose vertices backward. A target is marked with probability q_i.
A non-target is marked when its coin is one and its independent walk
hits the already marked downstream set, with probability q_i(1-ell).
These are exactly (94)'s second-branch probabilities. Thus this graph
samples the backward sets and

H_(T union T')=H_T union H_T'.                            (98)

Independently sample ordinary walks Z_(i,k), k>=1, from each i.
For deterministic initial sleeping set S and extra configuration eta,
declare v occupied exactly when H_{ {v} } is contained in S and every
Z_(i,k), k<=eta(i), avoids it. Equations (97)--(98) give every joint
occupation probability, hence the whole terminal occupied-set law.
One common graph and walk family therefore works for all S,eta.

This is a distributional representation, not an identification with the
physical ARW paths. The DAG assumption supplies the independence of
incoming counts from the local book and permits the backward substitution.
Simply reusing these independent coins on a cyclic graph is wrong:
it would give stationary occupation p_v, whereas the two-site half-sink
example has stationary occupation 22/39 at p_v=1/2.
A cyclic extension must change the backward construction; positivity of
individual proposed dual rows would not itself supply the coherent
random-graph realization.

An independent AI audit checked the recursion, extra-particle factors,
and common-graph representation. The creator verifier checks 5,568
weighted mixed-input/extra/target identities by raw absorption on
three- and four-site DAGs, and enumerates 1,088 independent coin/walk
draws to check the graphical law and its union identity exactly.

### Boundary weights as random additional targets

There is a general **proved closure identity**. For an arbitrary finite
block, let h=h_empty^t be its weighted ordinary exit potential, initially
assumed positive. Perform the Doob transform with internal probabilities
K(i,j)h_j/h_i and exit probabilities K(i,e)t_e/h_i.
Independently include each i outside the requested target set T in D
with probability 1-h_i. For full mixed initial input and extras eta,

F_(T,eta)(s;t)=h^eta E_D[m_(T union D)^(K^h)(s;eta)].      (99)

Indeed, weighted history telescoping gives h^(1+eta)/h_B at terminal
support B. On {T subset B}, the factor h_V/h_B equals P(D subset B).
This proves (99) directly, with no normalization of a survival event.

Apply (85a) in the transformed block and use
h_i h_C^(K^h)(i)=h_C^t(i). Then

Theta_(T,t)(C)=E_D[Theta_(T union D)^(K^h)(C)],
F_(T,eta)(s;t)=sum_C Theta_(T,t)(C)s_C (h_C^t)^eta.         (100)

Thus positivity of all joint-target rows, for the relevant Doob
transforms, implies positivity of every weighted-exit row. These are
convex mixtures of target rows, and sum to one. Zero h is handled by
continuity of the finite field coefficients. In an irreducible block
with h identically zero, the weighted no-extra event permits no exit
and hence forces full occupation; its row is exactly Theta_V.

### Two-vertex blocks and composition across components

**The full-wake dual coefficients are nonnegative on every graph whose
strongly connected components have at most two vertices.**
This extends the field result beyond DAGs. It still concerns full mixed
base inputs, and does not settle the main arbitrary-graph conjecture.

First consider a two-vertex block without self-loops. Write
K(0,1)=a, K(1,0)=b, with remaining probability sent to absorbing exits,
and assume ab<1. Put q_i=1-p_i and

c=ab, d=1-c, D_s=1-c q_0 q_1,
D_f=1-c(1-p_0 p_1), u=1-b q_1.

For target {0}, the four coefficients are

Theta(empty)=p_0/D_s+p_0 p_1 c q_0 u/(D_s D_f),
Theta({0})=q_0 d u(1-c q_1)/(D_s D_f),
Theta({1})=p_0 q_0 q_1 c d u/(D_s D_f),
Theta({0,1})=q_0 q_1[b(1-a)D_f+c q_0 d u]/(D_s D_f).     (101)

Every term is nonnegative, including a or b equal to zero or one when
ab<1. Interchanging the vertices gives target {1}; the full target uses
(89), and the empty target has only its empty-set coefficient.

Here is a derivation. Write delta=pi({0,1})=p_0p_1/D_f,
beta_i=pi({i}), and r_i=q_i/p_i. The stationary deletion identities
pi L_i=0 give

beta_0=r_1[b p_0 beta_1+
               ((1-b)+b(1-a)p_0)delta],
beta_1=r_0[a p_1 beta_0+
               ((1-a)+a(1-b)p_1)delta].

Their determinant is D_s>0. Solving gives pi f_0 equal to the first
line of (101). The singleton-C coefficients follow from (88), with
sigma_0=sigma_1=d and the one-particle values
S(0,0)=p_0/D_s, S(1,0)=b q_1p_0/D_s. Finally,

F_0 F_1=d A_0 A_1-(1-b)A_0-(1-a)A_1+(1-a)(1-b)I

gives the last line of (101). This proves the coefficients rather than
inferring their signs from numerical examples.

To compose blocks, process the condensation DAG in reverse order.
All successors of the current component are weighted absorbing labels.
Use (100) to sample its marked subset C from the normalized positive row
for the local target set. Record the factor s_C, and assign to each
vertex in that component the boundary weight h_C^t(i) for incoming
particles from upstream. Conditional on upstream data, the numbers
of incoming particles are independent of this component's records,
because no directed edge returns from a later component.

This is the block form of (95). Induction therefore proves a normalized
positive field for the whole graph. The factors from any specified
extra input again equal the global ordinary-walk avoidance potential
to the power eta. Positive Doob transforms preserve internal SCC
structure; self-loops may be removed using Section 21. Singleton
components use (94), and two-vertex components use (101), completing
the theorem.

This constructs a backward set for each fixed target set, simultaneously
for all initial flags and extra inputs. Union-coherence across different
target sets on cyclic graphs is not asserted; the common graphical
construction in (98) remains a DAG result.

An independent AI audit checked the polynomial identities, boundary
mixture, zero-weight cases, and component composition.
The exact verifier checks all 384 two-site field coefficients on 24
rational kernels, including zero transition probabilities and a site
without a direct sink edge. It also checks 4,608 weighted
mixed-input/extra/target identities on a four-site graph with two
nontrivial SCCs, including vanishing boundary weights.

### Two finite coherent laws and a conditioning obstruction

The stronger union-coherence requirement (98) is **exactly feasible** on
the symmetric two-site half-sink kernel at p=1/2 and the archived
heterogeneous three-site kernel rn.CASES[3]. The verifier stores rational
laws on tuples (H_v): 10 positive atoms on two sites and 50 on three.
It checks normalization and all 68 nonempty-target constraints
P(union_(v in T) H_v=C)=Theta_T(C), independently computing the right
side by raw mixed-input absorption. These two finite certificates give
no general cyclic construction.

One tempting construction already fails. Reuse Section 17's full-root
creator books unchanged when some initial particles are asleep. On the
two-site half-sink example, the latent C has masses 4/13 for the full
set and 3/13 for each other subset. For proper C the balanced prefix is
empty; the three acyclic exit forests have equal probability, the first
bits are 1_C, and subsequent records are fresh. For C=V all first bits
are one, so every mixed full input immediately ends full.

Starting with precisely site 0 asleep, this artificial sampler gives
terminal masses, in order empty,{0},{1},{0,1},

(2/39, 76/195, 4/195, 7/13).

The actual mixed-input law is

(5/78, 68/195, 19/390, 7/13).

Full occupation agrees, but proper targets do not. In particular the
artificial coefficient Theta_{ {1} }({0}) is -1/195, whereas the physical
coefficient is 3/130. The creator-book conditioning preserves the saved
all-active family, not arbitrary initial sleep labels. The exact
verifier checks both laws with first-record-augmented raw absorption.

### Why a homogeneous target-set Markov chain is insufficient

There is a scoped **proved obstruction** even in the two-site half-sink
example. Regard Theta as a matrix with rows the requested target T and
columns the resulting dependency set C, both ordered empty,{0},{1},V.
It is

Theta =
[1,     0,      0,      0;
 22/39, 21/65,  3/130,  7/78;
 22/39, 3/130,  21/65,  7/78;
 4/13,  3/13,   3/13,   3/13].

Its distinct positive eigenvalues are 1,1/2,3/10,1/13. Any real logarithm
commutes with Theta and preserves its one-dimensional real eigenspaces;
it is therefore the unique spectral logarithm. Exact projectors give

(log Theta)(V,empty)=(16 log 2-5 log 13)/11<0,

since 2^16<13^5. Consequently Theta is not exp(tG) for any Markov
generator G and t>0.

Nor can it be the absorption kernel of a homogeneous subset chain which
stops at state T with rate kappa_T and outputs that current state.
Such a generator would necessarily satisfy

G Theta=diag(kappa)(Theta-I),
G=diag(kappa)(I-Theta^(-1)),

while

I-Theta^(-1) =
[0,   0,     0,     0;
 8/3, -11/3, -4/3,  7/3;
 8/3, -4/3,  -11/3, 7/3;
 -4,  6,     6,     -8].

Positive probability of output V requires kappa_V>0, contradicting
the negative transition rate G(V,empty)=-4 kappa_V. The verifier checks
the inverse and logarithm's spectral coefficients exactly.

These statements exclude two specific homogeneous continuous-time
descriptions on target sets alone. They do not exclude the finite
coherent laws above, a discrete transition with kernel Theta, a
time-dependent construction, or a dual with additional state.

### An exact marked-site deletion with a signed target kernel

The first-contact cancellation extends to arbitrary mixed flags.
Fix x and put R=V\\{x}, h_i=P_i(sink before x)>0 on R. Let Khat be
the no-x Doob kernel from (92), retaining p_i. For stable B subset R put

L_(x,T)(B)=1{T subset B union {x}}-(A_x f_T)(B),
kappa_(x,T)(D)=sum_(B subset D) (-1)^|D\\B| L_(x,T)(B)/h_B.

Then the following identities are **proved**:

partial_(s_x) W_T^V(s)
 =h_R sum_(D subset R) kappa_(x,T)(D) W_D^Khat(s_R),       (102)

Theta_T^V(C union {x})
 =h_R sum_(D subset R) kappa_(x,T)(D) Theta_D^Khat(C),
 C subset R.                                             (103)

Hold x's resident untouched while stabilizing the mixed input on R.
A first contact with x makes its active and sleeping alternatives
identical, so those histories cancel. On no-contact histories ending
at B their difference is L_(x,T)(B). The no-contact law is
(h_R/h_B) times the mixed full-input law for Khat, by history
telescoping, including the initial sleeping particles. Expanding
L(B)/h_B=sum_(D subset B)kappa(D) proves (102), and comparing field
coefficients proves (103). With zero h the undivided no-contact identity
remains valid; these divided formulas must not be used directly.

This is not a positive target-mixture recursion. For T={v}, v!=x,

kappa_(x,{v})(empty)=-S(x,v),

where S=(I-diag(q)K)^(-1)diag(p) is the one-particle first-sleep kernel.
It is negative whenever x can reach v. With positive direct-sink
probabilities the smaller polynomials W_D^Khat are linearly independent:
the nonzero orthogonal rows delta_empty U_C are independent by (82),
their nonvanishing following from the positive full-target coefficient,
and the target monomials form a payoff basis. Thus no alternative
nonnegative mixture of these same reduced target polynomials represents
the marked branch.

Even allowing their weighted boundary versions does not always help.
On the two-site half-sink example, x=0 and T={1} give

partial_(s_0) W_{ {1} }(s)
 =3/130+(7/78)s_1
 =-(1/15)W_empty^Khat+(7/39)W_{ {1} }^Khat.                (104)

The remaining site's fields with exit weight ell in [0,1] are
F_empty=(1+ell)/2+(1-ell)s_1/2 and F_{ {1} }=(1+s_1)/2.
Every nonnegative mixture has linear coefficient at most its constant
coefficient, whereas 7/78>3/130 in (104).

There can also be a negative coefficient at a nonempty target.
On the three-cycle with each clockwise/sink jump of probability 1/2
and all p_i=1/2, x=0,T={2} give

kappa_(0,{2})({1})=-54646/725571,
Theta_{ {2} }({0,1})=12211/2042348>0.

The verifier checks (103) for every C on these two kernels. The identity
retains A_x from the full graph, so it is not a closed smaller-graph
algorithm. These obstructions leave open a recursion which marks several
vertices together or keeps additional path/particle state.

## 24. Excursion dynamics and a graphical dual before the first loss

**The following kernels and full-target dual are proved.** They use
physical iid instructions, not creator-conditioned books.

### A finite occupation/activity kernel

At a configuration with at most one particle per site, write B for the
occupied support and A subset B for its active sites. Select i in A.
With probability p_i it sleeps, replacing A by A\\{i}. Otherwise take
an ordinary K-walk from i, starting with its first jump, until it hits
i again, a vacancy in V\\B, or the sink. At each intermediate occupied
site hold one active resident and move the surplus; a sleeper encountered
there wakes. This is a legal ARW schedule. If the endpoint is j in V,
replace B by (B\\{i}) union {j}; at a sink endpoint replace it by B\\{i}.
Remove i from A and activate every visited internal vertex, including
the endpoint when it lies in V.

The unused relevant records are fresh after every such block. The
procedure only reads instructions needed for the selected walk; forced
departures skip unsuccessful sleeps. Untouched sleepers retain their
flags and unread records. Thus this defines an exact Markov kernel on
the 3^|V| vacant/sleeping/active configurations. It retains all later
returns to a vacancy. The excursion verifier independently checks 1,198
terminal-law entries and 61 walk blocks on the archived kernels.

It also has a direct termination bound. Let h_i be ordinary-walk
expected time to the original sink, and put Phi(B)=sum_(i in B)h_i.
A jump block of length tau changes its conditional expectation by
-E tau<=-1, since only the starting and ending occupied sites differ.
Hence every nonterminal block decreases expected Phi by at least
q_min=min_i(1-p_i)>0. Optional stopping gives expected block count at
most Phi(B_initial)/q_min.

### The full-support activity process and its inverse images

When B=V, a jump block either returns to its source or loses a particle
at the sink. A returned excursion from i, visiting D outside i, changes

A -> A union D.

Any sink loss precludes full terminal occupation. Use independent
recovery marks at rate p_i and excursion marks at rate q_i, ignoring
marks whose source is inactive. These rates define auxiliary block time,
not the physical continuous-time ARW clock.

The returned-excursion/recovery process is additive in A. For the
avoidance test 1{A intersect H=empty}, reversing the graphical maps gives

recovery at i:             H -> H\\{i};
returned excursion i,D:    add i if D intersects H;
sink excursion from i:    add i unconditionally.

The last map records the constraint that a loss mark may only occur at
an inactive source. Thus the backward subset chain has death rate p_i
at i in H and, for i outside H, birth rate q_i e_i(H), where

e_i(H)=P_i(hit H or the sink before the first return to i)
      =d_(V\\H)/d_(V\\(H union {i})).                     (105)

The ratio follows from the diagonal cofactor identity
e_i(H)=1/[(I-K_(V\\H,V\\H))^(-1)]_(i,i).
A self-loop is a return at time one with D empty, so no self-loop
reduction is needed.

All birth and death rates are positive under the standing hypotheses.
The chain is irreducible and reversible with stationary law

w(H)=d_V r_H/[Z_V(r)d_(V\\H)]=Theta_V(H),                 (106)

because w(H+i)/w(H)=r_i e_i(H), the birth/death rate ratio.
This gives a direct graphical explanation of the determinant law (89).

Start the backward chain empty. Reversing the graphical constraints
over an interval of auxiliary length t gives

P_(V\\S)(no loss by t)=P_empty(H_t subset S).              (107)

To pass to infinite time, ignore the failure marks in the forward
activity process. In each unit interval, no excursion marks and at least
one recovery at every vertex have probability
exp(-sum_i q_i) product_i(1-exp(-p_i))>0, uniformly over its state.
It therefore becomes empty almost surely. No loss forever is equivalent
to reaching complete sleep before any loss. Letting t tend to infinity
in (107) recovers the full-target dual. Independent extra particles
are incorporated by the surplus walks in (85a).

A separate AI audit checked the legal scheduling, freshness, inverse
maps, self-loops, determinant ratio, and the limiting argument.
The exact excursion verifier checks 404 generator-duality entries,
404 detailed-balance entries and all 38 full mixed-input probabilities
on the five archived kernels, against raw ARW absorption.

### Why the continuation after a loss needs additional information

For a proper target, histories after the first sink loss still contribute:
a site vacated by that loss can later be filled. The general binary
kernel above retains its location and the excursion's activated sites.

One cannot assume that the full-wake coefficient sign survives a
vacancy. At p_i=1/2 take 0->2 surely, 1->2 surely, and
2->0,sink equally. With initial occupied support {0,1}, target 0,
the all-active probability is 43/70; making precisely site 1 initially
asleep gives 4/7. Their difference is -3/70. This is an exact
activation-field obstruction for a partial base, not a counterexample
to all-active add-one monotonicity.

### Exact first-loss currents and a normalized potential

For a fixed target set T, abbreviate theta(H)=Theta_T(H), and let
theta^(i)(H) be the analogous coefficient with initial support
V\\{i}, on the original graph. The initially vacant i remains
available for later returns. For i outside H set

h_H(i)=P_i(sink before H),
gamma_i(H)=P_i(sink before H or the first return to i).

Excursion renewal at i gives gamma_i(H)=e_i(H)h_H(i).       (108)

Forcing the first instruction at i yields the **exact edge identity**

q_i e_i(H) theta(H)-p_i theta(H+i)
 =q_i gamma_i(H) theta^(i)(H)=:j_i(H).                    (109)

To prove it, extract the coefficient of the other initial flags s_H.
The sleep branch contributes p_i[theta(H)+theta(H+i)].
A returned excursion contributes only if it avoids H, with coefficient
q_i(1-e_i(H))theta(H). A sink excursion also must avoid H, giving
q_i gamma_i(H)theta^(i)(H). Any visited flag in H has been erased, so
its mixed difference vanishes. Rearrangement gives (109).

Let Q be (105)'s birth/death generator, acting on functions. The
first-loss source at sleeping set S is

ell_T(S)=sum_(i outside S) q_i E_i[
  1{sink before first return i} m_T^(i)(S\\range(excursion))].

Its multilinear extension is
sum_i(1-s_i)sum_(H subset V\\i)j_i(H)s_H. Its coefficient b_T(H)
therefore satisfies

b_T(H)=sum_(i outside H)j_i(H)
       -sum_(i in H)j_i(H\\{i}),
Q^transpose theta=-b_T,  sum_H theta(H)=1.                (110)

The flux is signed because the one-hole coefficients can have either
sign. In particular sum_H b_T(H)=0; irreducibility of Q leaves no
separate solvability issue.

There is more information than divergence alone. Let pi(H)=Theta_V(H)
be (106), and u(H)=theta(H)/pi(H). Detailed balance converts (109) to

u(H)-u(H+i)
 =h_H(i)theta^(i)(H)/pi(H)=:g_i(H).                       (111)

Hence g_i(H)+g_j(H+i)=g_j(H)+g_i(H+j). Define G(empty)=0 and
G(H+i)-G(H)=g_i(H); this is independent of the insertion order.
Normalization reconstructs the whole row and its exact remaining sign:

theta(H)=pi(H)[1+E_pi G-G(H)],
theta>=0  iff  max_H G(H)-E_pi G<=1.                      (112)

The potential-height bound is **false in general**, by Section 22's
negative three-site coefficient. The identity remains an exact criterion
on restricted classes. Reversibility and zero circulation do not
establish the bound. If h_H(i)=0, the edge has zero
current, giving theta(H+i)=r_i e_i(H)theta(H) without a sign assumption.

The preceding one-vacancy example gives an explicitly negative current:
for T={0}, H={1}, i=2,
theta^(2)({1})=-3/70, e_2({1})=gamma_2({1})=1/2,
pi({1})=1/5, and j_2({1})=-3/280.
Indeed theta({1})=0 and theta({1,2})=3/140, so u increases across
that edge. This does not contradict nonnegative full-input coefficients.

The exact verifier checks 722 edge identities, 484 square circulations,
and 404 potential reconstructions against raw full and one-hole laws
on the archived kernels.

### Circulation renewal on a directed cycle

A special case isolates a precise remaining comparison. Let
K(i,i+1 mod n)=a_i and K(i,sink)=1-a_i, with c=product_i a_i<1.
For S subset R subset V, let m_(R,T)(S) be the target-T probability
from one particle at each R-site, asleep precisely on S. The graph
remains the original cycle, retaining returns through V\\R. Write

Theta_(R,T)(D)=sum_(S subset D)(-1)^|D\\S|m_(R,T)(S),
mu_(R,T)=m_(R,T)(empty).

For D nonempty and i outside D, let j be the first clockwise site of D
after i and let alpha_(i,D) be the product of a_z on that arc,
including i and excluding j. The following identities are **proved**:

Theta_(V,T)(D union {i})
 =r_i[Theta_(V,T)(D)
       -(1-alpha_(i,D))Theta_(V\\{i},T)(D)],              (113)

Theta_(V,T)({i})
 =r_i(1-c)[mu_(V,T)-mu_(V\\{i},T)].                       (113a)

Apply the binary excursion kernel with i active and take the finite
difference over D in the other initial flags. The sleep branch has
coefficient p_i[Theta_(V,T)(D)+Theta_(V,T)(D+i)].
Once a jump excursion reaches any site of D, that initial flag has
been erased and its D-coefficient vanishes. The remaining branches lose
the particle before reaching j, with total probability 1-alpha_(i,D);
their coefficient is Theta_(V\\i,T)(D). This proves (113).
When D is empty, keep the complete-circulation branch c mu_(V,T),
which proves (113a).

In particular, a leak-free arc gives Theta(D+i)=r_i Theta(D).
For a leaking arc, positivity still requires the comparison inside
(113)'s brackets. The one-hole coefficient need not be nonnegative;
the identity is not an induction proving the general cycle sign.

There is also a normalized scalar first-loss renewal. Fix the testing
order 0,...,n-1, put E_i={0,...,i-1}, x_k=i+k mod n, and define
L_(i,k)={x_1,...,x_k},
beta_(i,k)=(product_(h<k)a_(x_h))(1-a_(x_k)).
Testing sleep coins until the first jump gives

[1-c(1-p_V)]mu_(V,T)
 =p_V+sum_i q_i p_(E_i) sum_(k=0)^(n-1)
       beta_(i,k) m_(V\\{i},T)(E_i\\L_(i,k)).              (114)

All displayed mixture weights are nonnegative and sum to
1-c(1-p_V)>0. A complete circulation returns the full input to the
all-active state with fresh suffixes. These physical weights do not
remain termwise positive after arbitrary flag differences.

The exact excursion verifier checks 192 marked-field identities and
16 first-loss renewals on two heterogeneous three-cycles, one with
a single leaking edge.

## 25. A three-site pair renewal and its occupation numerator

**The reduction below is proved. A positive polynomial certificate then
proves the three-site occupation comparison for all parameters.**
Let V={0,1,2}. Remove self-loops as in Section 21, using the resulting
K,p throughout. In particular r'_i=r_i(1-K_ii). Put

k_i=K(i,sink), R_i=V\\{i}, d_i=det(I-K_(R_i,R_i)),
d=det(I-K), sigma_i=d/d_i,
Z_i=1+d_i(r_j+r_k+r_j r_k), {j,k}=R_i,
Z_V=d sum_(C subset V) r_(V\\C)/det(I-K_CC).

Let pi be the full all-active law and mu_x the all-active law from R_x.
Write beta_i=pi(R_i), gamma_i=mu_x(R_i), using **exact support** masses.
Define a three-by-three nonnegative matrix Rcal with zero diagonal by

Rcal_(i,j)=(r_i/Z_i)
 [K_ij(1+d_i r_k)+K_ik K_kj], {i,j,k}=V,

z_i=(r_i/Z_i)[sigma_i+(Z_i-1)k_i
       +d_i(r_k K_ij k_j+r_j K_ik k_k)].

Then

(I-Rcal) beta=z/Z_V,
(I-Rcal) gamma=e_x/Z_x.                                  (115)

For a derivation, apply the last-arrival identity at i to terminal
support R_i. A pair predecessor R_j contributes r_i W_i(R_j,R_i).
If its forced first jump lands at the empty j, the remaining two-site
full-output probability is (1+d_i r_k)/Z_i. If it lands at the
sleeper k, force another jump; full output requires a jump to j,
then a full active pair, giving K_kj/Z_i. These yield Rcal_(i,j).

For the full predecessor, let h_j=P_j(sink before i) on R_i. The
weighted two-site full-output formula gives, from sleeping R_i plus
one active at j, probability (h_j+d_i r_k k_j)/Z_i of avoiding i
and ending full in R_i. Thus

W_i(V,R_i)=k_i+
 [K_ij h_j+K_ik h_k+d_i(r_k K_ij k_j+r_j K_ik k_k)]/Z_i.

Use sigma_i=k_i+K_ij h_j+K_ik h_k and pi(V)=1/Z_V to obtain z.
For the lower law, the only no-i-visit source is i=x; its R_x mass is
1/Z_x on the killed two-site domain. This is generally not gamma_x.

The inverse G=(I-Rcal)^(-1) exists and is nonnegative without an
occupation assumption. Indeed, z_i>=r_i sigma_i/Z_i>0, and the actual
full-input pair vector satisfies beta_i>=z_i/Z_V>0 and Rcal beta<beta.
Consequently rho(Rcal)<1. No bound on its ordinary row sums is asserted.

### Extracting occupation with nonnegative avoidance coefficients

Fix v!=x. Let w_v(B) be the no-return probability after removing v
from sleeping B, forcing one K-jump from v, and stabilizing while
killing on any visit to v. Put

S=(I-diag(q)K)^(-1)diag(p),
a_v=0,
a_i=r_v S_vv[w_v({v})-w_v(R_i)], i!=v,
a_*=r_v S_vv[w_v({v})-w_v(V)].

Adding sleeping residents can only increase the visited set, hence
0<=a_i<=a_*. For any all-active creation law mu whose input contains v,
the scalar deletion identity and
S_vv=1/[1+r_v w_v({v})] give

mu(f_v)=S_vv+sum_i a_i mu(R_i)+a_* mu(V).

Substitution of (115) therefore yields

Delta_(x,v):=pi(f_v)-mu_x(f_v)
 =[a_*+a^transpose Gz]/Z_V
   -a^transpose Ge_x/Z_x.                               (116)

Every avoidance coefficient is explicitly a two-site quantity. Let
R_v={j,k}, and let theta(C) be its weighted empty-target field from
Section 23, with boundary weights zero at v and one at the original
sink. Set

h=(I-K_(R_v,R_v))^(-1) k_(R_v),
u=(I-diag(q_(R_v))K_(R_v,R_v))^(-1)
       [p_(R_v)+diag(q_(R_v))k_(R_v)].

Then

w_v({v})=k_v+sum_(j in R_v)K_vj u_j,

w_v({v,j})=k_v+K_vk[theta(empty)+theta({j})]
                 +K_vj[k_j u_j+K_jk theta(empty)],

w_v(V)=k_v+sum_(j in R_v)K_vj[
          theta(empty)h_j+theta(R_v\\{j})k_j].             (117)

The second formula separates the first jump to an empty or occupied
site. The third is the weighted full-base scattering formula with
one extra particle. Section 23's continuous extension handles zero h.
There is no full-graph active-insertion operator hidden in (117).

### The Schur complement

Put J=V\\{x}, M=(I-Rcal_(J,J))^(-1), and define

tau=1-Rcal_(x,J) M Rcal_(J,x)>0,
A=a_x+a_J^transpose M Rcal_(J,x),
B=a_*+a_J^transpose M z_J,
ell=Z_V/Z_x-z_x-Rcal_(x,J)M z_J.

Here A,B are nonnegative, and

Delta_(x,v)=[tau B-A ell]/(Z_V tau),
Theta_{ {v} }({x})
 =[r_x sigma_x/(Z_V tau)][tau B-A ell].                  (118)

Thus the occupation comparison is exactly tau B>=A ell.
It is automatic for ell<=0; ell can be positive. The reductions and
positivity of their ingredients alone do not prove this inequality.

A separate AI audit checked the deletion entries, source normalization,
transience, and occupation extraction, including kernels with zero
direct-sink probabilities. The maintained three-site renewal verifier
checks the pair masses, explicit avoidance formulas and all 18
off-diagonal gains on three specified heterogeneous kernels. These are
identity audits; the general sign is supplied by the certificate below.

### A target-independent pair-deficit capacity

There is a stronger **computer-assisted theorem** which supplies a
target-independent occupation proof:

pi(V)>=sum_i [mu_x(R_i)-pi(R_i)]_+.                       (119)

Thus the upper law's full-support mass can cover every deficit in its
exact pair masses, with the two occupied sites preserved. This concerns
the pair/triple layer; it is not an occupied-set coupling of the whole law.

For the certificate relabel x=1. Put P_d=d_0d_1d_2,
Zbar=P_d Z_V, N=diag(Z_i)(I-Rcal), b=diag(Z_i)z,
bbar=P_d b, and xi=bbar-Zbar e_1. For a family J of missing-site indices,

pi(V)+sum_(i in J)(beta_i-gamma_i)
 =T_J/[Zbar det N],
T_J=P_d det N+1_J^transpose adj(N)xi=d_1 H_J.             (120)

The independent pair-deficit verifier reconstructs the system,
checks the adjugate orientation with column-replacement determinants,
and verifies the exact division by d_1. Each H_J then has a positive
expansion in jump and sink probabilities:

H_J=sum_(alpha,nu) c^J_(alpha,nu) r^alpha
       product_(i!=j)K_ij^nu_ij
       product_i k_i^(m^J_(alpha,i)-sum_(j!=i)nu_ij).       (121)

All stored coefficients are positive integers and all exponents are
nonnegative. Exact inversion of the row substitutions recovers every
original coefficient.

| Family J | Terms in H_J | Sleep-odds blocks | Positive coefficients |
|---|---:|---:|---:|
| {1} | 1,725 | 24 | 44,017 |
| {1,2} | 4,058 | 33 | 83,375 |
| V | 6,415 | 42 | 117,062 |

These polynomial identities hold on the whole kernel simplex, including
zero direct-sink probabilities. All denominator factors in (120) are
positive by the proved renewal transience and sink accessibility.

It remains to justify why these three families prove the positive-part
bound. Set u=(gamma-beta)/pi(V), and U=V\\{x}. From (115),

u_U=(I-Rcal_(U,U))^(-1)[Rcal_(U,x)u_x-z_U].

Thus u_x<=0 implies u_U<=0. Every nonempty positive-deficit family
contains x. On three vertices its size is one, two, or three; relabelling
reduces it to one of the three certified families. Applying (120) to
that family proves (119).

For v!=x, combine (116), 0<=a_i<=a_*, and (119):

Delta_(x,v)
 >=a_*[pi(V)-sum_i(gamma_i-beta_i)_+]>=0.

The added target is covered by the general diagonal argument. This
proves the three-site full-minus-one theorem without needing a
target-dependent bound on the individual avoidance coefficients.

The maintained verifier is `verify_three_site_pair_deficit.py`
(SymPy 1.14.0); it regenerates and checks the saved
[coefficient artifact](three-site-pair-deficit-certificate.json.gz).
Its SHA256 is
`c6e4487ba10102d45006515e1063233077e76cf3f58ec0fc153c28ae46a9bb67`.

A simpler localization of the debt is false. On
K_01=K_12=1, K_20=K_2s=1/2,
p=(1/3,2/5,3/7), x=2, the deficit vector gamma-beta is
(420,333,2886)/45103, positive at all three pairs, while
pi(V)=4876/45103. The total capacity bound, not a single deficient
pair, is the relevant statement.

### A scalar threshold for the three-site capacity

The pair certificate has a compact Schur reduction. Continue with x=1,
and put u=(gamma-beta)/pi(V), so Nu=Z_V e_1-b. Write A_ij=-N_ij>=0
for i!=j and define

D=Z_0Z_2-A_02A_20,
v_0=(Z_2A_01+A_02A_21)/D,
v_2=(A_20A_01+Z_0A_21)/D,
w_0=(Z_2b_0+A_02b_2)/D,
w_2=(A_20b_0+Z_0b_2)/D,
tau=Z_1-A_10v_0-A_12v_2,
ell=Z_V-b_1-A_10w_0-A_12w_2.

The M-matrix property gives D,tau>0 and v_i,w_i>=0. Then

u_1=ell/tau, u_i=v_i u_1-w_i for i=0,2.

Consequently sum_i(u_i)_+<=1 is exactly the single threshold

ell/tau<=min_(B subset {0,2})
              (1+sum_(i in B)w_i)/(1+sum_(i in B)v_i).

For u_1<=0 all deficits are nonpositive; for u_1>0 maximize over
the subset of positive off-source deficits. The four cuts have

T_({1} union B)=P_d D[tau(1+sum_B w_i)-ell(1+sum_B v_i)],
det N=D tau.

This reduction does not supply the sign from transience alone. For
abstract renewal data R_01=R_21=4/5, other entries zero,
z=(1,1,11)/10, and source injection (17/20)e_1, the deficit vector
(I-R)^(-1)[(17/20)e_1-z]=(1/2,3/4,-1/2). Singleton-source and
full-sum cuts pass, but the two-positive-deficit sum is 5/4. These
are abstract data, not an ARW counterexample.

Some certificate factors have a short walk interpretation. Put

G_01=K_01+K_02K_21, G_10=K_10+K_12K_20,
L_0=k_0+K_02k_2, L_2=k_2+K_20k_0,
Q_i=K_(i,1)(1-K_(1,i))+K_(i,j)K_(j,1)
    +K_(1,i)K_(i,1)K_(i,j)K_(j,i), {i,j}={0,2}.

At r_2=0 both H_{ {1} } and H_{ {1,2} } equal
r_0 d_2 G_01G_10[1+r_1(k_1+K_12k_2)]Z_2.
Their highest total-odds-degree terms, and that of H_V, are respectively

d_0d_2d[d_0K_10Q_0+d_2K_12Q_2]r_0^2r_1^2r_2^2,
d_0^2d_2^2d L_2 r_0^2r_1^2r_2^3,
d_0^2d_2^2d r_0^2r_1^2r_2^2(L_0r_0+L_2r_2).

Hence along r_i=t rho_i with fixed rho_i>0 and t tending to infinity,
u_i=-t rho_i h_i+O(1) for i=0,2, where h_i=L_i/d_1, while

u_1 -> q_*=k_1+K_10 h_0/d_2+K_12 h_2/d_0<=1,
1-q_*=K_10Q_0/(d_1d_2)+K_12Q_2/(d_1d_0).

The probability q_* counts immediate sink escape from 1, or a first
jump to i, any number of direct i->1->i round trips, and then sink
absorption before another visit to 1. These disjoint ordinary-walk
events explain its bound. This is an asymptotic explanation of one
cut, not a replacement for the finite-odds capacity certificate.
`structural_audit` in the pair-deficit verifier checks these factors,
the generic Schur identities, and the abstract normalization obstruction.

### An independent scalar certificate and strictness

The scalar certificate below independently verifies the occupation sign
and identifies when the off-diagonal comparison is strict.

Relabel the target 0 and the added site 1. Define

P_d=d_0 d_1 d_2, Zbar=P_d Z_V,
N=diag(Z_i)(I-Rcal), b=diag(Z_i)z, bbar=P_d b,
M_1=I+diag(r)(I-K), D_one=det M_1, u=adj(M_1)e_0.

For {i,j,k}=V set

L_ii=r_j r_k(K_ji k_k+K_ki k_j+K_ji K_ki k_i),
L_ij=-r_k[(1+r_j)k_k+K_kj k_j],
c_i=r_j r_k[r_i k_j k_k+
            (k_j+K_ji k_i)(k_k+K_ki k_i)],
xi=D_one(1-e_0)-Lu.

The complete two-sleep and three-sleep inverse formulas (43)--(45) give,
writing Q_i=mu_i(f_0),

N^transpose Q=1-e_0-L S_(.,0),
Z_V pi(f_0)=1+b^transpose Q-c^transpose S_(.,0).

Thus there is an exactly specified polynomial H(K,r) satisfying

r_0 d_1 H
 =P_d det(N)[D_one-c^transpose u]
   +xi^transpose adj(N)[bbar-Zbar e_1],                   (122)

and the desired increments are

pi(f_0)-mu_1(f_0)=r_0 d_1 H/[D_one Zbar det(N)],
Theta_{ {0} }({1})=r_0 r_1 d H/[D_one Zbar det(N)].        (123)

All formulas here substitute k_i=1-sum_(j!=i)K_ij and K_ii=0.
Exact polynomial division in (122) has zero remainder. A separate
derivation from the killed two-site fields (117) gives the same
polynomial identically, not merely at sampled parameter choices.

The positive grouping uses jump rates relative to the sink:
for k_i>0 write A_ij=K_ij/k_i and t_i=1+sum_(j!=i)A_ij,
so K_ij=A_ij/t_i and k_i=1/t_i. The certificate is

H(K(A),r)=sum_alpha r^alpha
 P_alpha(A)/(t_0^m_(alpha,0)t_1^m_(alpha,1)t_2^m_(alpha,2)),
alpha in {0,1,2} x {0,1,2,3} x {0,1,2,3,4}.             (124)

Each of the 60 P_alpha has nonnegative integer coefficients.
The clearing degrees are at most (5,7,6), respectively. Altogether
there are 330,334 nonzero coefficients, each strictly positive.
The reduced H has 23,312 monomials before the rate substitution.

The maintained symbolic verifier constructs H from the compact formulas,
checks its equality with the independent complete-block numerator,
computes every rate coefficient exactly, and inverts all 60 substitutions
back to the original coefficients. This is a finite symbolic certificate
for the whole parameter domain, not a finite parameter screen.
The compressed coefficient artifact is
[three-site-positive-certificate.json.gz](three-site-positive-certificate.json.gz).
Its SHA256 is
`5821a322d9ee946df6abbb07534063869706247b8925f93a7b9f621d13587070`.

Equivalently, (124) is the polynomial identity (121), with the scalar
certificate's coefficients and degrees. Thus it holds directly when
some k_i=0, as well as for zero internal edges. The denominators are
strictly positive: d,d_i are principal minors of the transient walk
matrix; D_one is the one-particle resolvent determinant; Zbar>0; and
det(N)>0 by the proved transience of Rcal. This proves (123)'s sign on
the whole sink-accessible kernel simplex. The endpoint-preserving
self-loop reduction then covers the original model.

There is also a useful exact factorization at the two extreme powers
of the third site's sleep odds. Write H=sum_(j=0)^4 r_2^j H_j and put

G_01=K_01+K_02 K_21, G_10=K_10+K_12 K_20,
L_1=k_1+K_12 k_2, L_2=k_2+K_20 k_0,
Lambda_2=1+r_0+r_1+d_2 r_0r_1,
F=G_01 G_10(1+r_1 L_1)Lambda_2 Z_2.

Then H_0=d_2 F and H_4=d_0 d L_2 F. Every monomial of H contains one
of the four edge products in G_01 G_10. Every positive rate monomial
does too; in odds block alpha, each of the four minimal edge products
has coefficient binom(2,alpha_0)binom(3,alpha_1)binom(4,alpha_2).
These facts are checked by `structural_audit` in the same verifier.
Neither G_01 nor G_10 divides the whole H, so the endpoint factorization
does not extend to all five layers.

In particular, all coefficients of H as a polynomial in r are
nonnegative on the physical kernel simplex, and

H(K,r)>=H(K,0)=d_2 G_01 G_10.

Consequently the off-diagonal full-minus-one increment is **strict
exactly when the added site and target lie in the same strongly
connected component**. Mutual reachability makes the displayed lower
bound positive; its failure makes every numerator monomial vanish.
This includes zero direct-sink probabilities and, after reduction,
self-loops. This strictness statement concerns the full-minus-one
comparison, not every arbitrary-input add-one comparison.

Run `python3 -B verify_three_site_certificate.py` from this directory
with SymPy available; the checked version is 1.14.0. The verifier writes
no files and uses no floating point. The fast standard-library companion
`verify_three_site_renewal.py` independently audits the physical identities.

Finally, surplus IDLA couples the binary supports of any finite all-active
input and that input plus one particle by inclusion, with fresh suffixes.
On at most three vertices, the only nontrivial new support comparison
has two common roots; one common root is N1, and full support is stationary
under an extra addition. The off-diagonal certificate and the general
added-target proof therefore establish **all finite all-active inputs on
three vertices**. Section 13's localization then proves the main
single-target inequality on any graph whenever the target's strongly
connected component has at most three vertices.

The arbitrary-graph occupation conjecture remains **OPEN**. General
activation-field positivity is false even on a three-site cyclic block,
by Section 22's counterexample; the full-field SCC theorem in Section 23
has the sharp size-two hypothesis. A short structural explanation of
the positive occupation numerator remains a separate research objective.


## 26. A pair-to-triple capacity route for general N2

**The reduction is proved; the capacity inequality on arbitrary graphs
is OPEN.** Let mu_l be the stabilization law of an all-active two-site
input A, and mu_u that of A+delta_x, x outside A. Write mu(C) for exact
terminal support mass. For every pair C put

d_C=[mu_l(C)-mu_u(C)]_+.

An inclusion allocation consists of nonnegative t_(B,C), indexed by
triples B and pairs C subset B, with

sum_(B superset C) t_(B,C)=d_C,
sum_(C subset B) t_(B,C)<=mu_u(B).                        (125)

This is feasible exactly when, for every family D of pairs,

sum_(C in D)d_C
 <=sum_(B: |B|=3, some C in D is contained in B)mu_u(B).  (126)

The equivalence follows from the finite capacitated bipartite
max-flow/min-cut theorem. In particular, (126) is independent of the
occupation target. On three sites it is precisely the stronger
capacity statement (119).

For a fixed target v, it is enough to require (126) only for families
of pairs containing v, and use only triples containing v. Such a
target-specific allocation need not be consistent with allocations
for other targets.

To prove sufficiency, use the last-arrival identity on an arbitrary
all-active input eta. With S as in Section 25, set

a_v(B)=r_v S_vv[w_v({v})-w_v(B)], v in B,
n_v(eta)=P_eta(no visit to v).

Then a_v({v})=0, a_v(B)>=0, and a_v is increasing under inclusion.
Summing (7) and using S_vv=1/[1+r_v w_v({v})] gives exactly

q_v(eta)=S_vv[1-n_v(eta)]
       +sum_(B contains v, |B|>=2)a_v(B)mu_eta(B).        (127)

The no-visit term increases when an active particle is added, by
visited-set monotonicity. The lower law has no triples, so its remaining
possible negative contribution is at most

sum_(C contains v)a_v(C)d_C
 <=sum_(B contains v, |B|=3)a_v(B)mu_u(B),

where (125) and a_v(C)<=a_v(B) prove the inequality. Thus the star
version of (126) proves N2 for this input and target, including targets
outside both input supports. Surplus IDLA and the proved N1 case
then handle repeated roots and an extra source already in A.

Equivalently, requiring all the cuts (126) says that every increasing
event generated by pair supports has larger probability under mu_u.
For the star restriction these events are
{v occupied and some site of W occupied}, W subset V\{v}.
Indeed, subtracting the probabilities gives the raw pair deficits plus
neighboring triple mass; removing pairs with nonpositive deficit
recovers the positive-part formulation. This is a stronger sufficient
statement, not a reformulation of single-target monotonicity.

The all-active creation hypothesis cannot be dropped even at the
mixed contact state with the extra particle at its active site. On
the path 0--1--2, let each endpoint jump to its neighbor or the sink
equally, let 1 jump to 0 or 2 equally, and put p_i=1/2. Starting from
(1,-1,0), where -1 denotes a sleeper, the probability that both 0 and
1 finish occupied is 45/82. From (2,-1,0) it is 45323/99876, a decrease
of 9487/99876. Thus a termwise mixed-contact pair-capacity proof fails;
the generating all-active history must remain in the comparison.
This is not a counterexample to (126) in its stated scope or to N2.

**Finite exact evidence.** The maintained verify_pair_capacity.py
checks the archived four-site kernel and Section 6's heterogeneous
conditional-order example, all twelve distinct-root additions in
each. It independently solves raw ARW absorption, checks all 63
nonempty pair-family cuts per comparison, and constructs rational
inclusion allocations satisfying every pair demand and triple supply.
All 24 comparisons pass. The smallest cut slack in the heterogeneous
case is 432635184999258001002/131921836276734018522239.
This gives no arbitrary-graph proof. The next structural question is
whether the backward deletion equations control these pair deficits
by local triple capacities without assuming the desired order.

### Joint occupation of the two original roots

A formerly open smaller inequality is now **PROVED**, for arbitrary
root sets, by Section 16's fresh-book coupling. The boundary identity
for Gamma_R implies

nu_(R,B)(R occupied)=1/m_R(B), B containing R,
m_R(R)=1/mu_R(R).

The pointwise entrance-mass bound (47a), in particular

m_R(R union {y})<=m_R(R), y outside R,

proves that adding any active input increases the probability that all
original roots remain occupied, relative to the binary-root baseline.
This is consistent with the saved negative individual-target entrance
example: it concerns one root's marginal, not their joint event. For
one root (47a) also follows from m_{ {a} }(B)=1+r_a w_a(B) and the
decrease of w_a under background inclusion.

There is an exact cone reduction for this joint event at every highest
particle layer. Write C_k for creation-cone measures supported on at
most k sleepers. Fix any k-set C
and source x, and put g_C(B)=P(C occupied after adding x to B)-1{C subset B}.
For any mu in C_k, Section 8 gives

mu=[mu(C)/mu_C(C)]mu_C+sigma,

where sigma>=0, sigma(C)=0, and sigma has at most k sleepers.
Consequently sigma g_C>=0, since C is the only initial support at this
count where the joint-event response can be negative. Equation (47c)
gives mu_C g_C>=0. Therefore

mu A_x(C occupied)>=mu(C), mu in C_k, |C|=k.              (127a)

Only nonnegativity of the remainder is used; sigma need not lie in C.
Taking k=2 proves every single-pair cut of (126), for every pair C,
not only the two original input sites. Unions of deficient pairs still
require the shared triple-capacity budget: summing individual cuts can
count the same triple more than once. The general N2 target inequality
remains open.

There is a useful quantitative form. For C subset B, |C|=2, |B|=3, put
beta_C(B)=m_C(B)/m_C(C), so 0<beta_C(B)<=1. With alpha=mu(C)/mu_C(C)
and sigma=mu-alpha mu_C as above, (21) gives the exact identity

mu(C)-(mu A_x)(C)
 =sum_(B contains C) beta_C(B) alpha (mu_C A_x)(B)
  -(sigma A_x)(C)
 <=sum_(B contains C) beta_C(B)(mu A_x)(B).              (127b)

Here the sums are over triples. Thus a sharper sufficient scalar budget
would be sum_(C subset B, v in C) a_v(C) beta_C(B)<=a_v(B). This
pointwise strengthening is false in general: on Section 22's three-site
rare-sleep kernel, p_i=1/1001, v=0, B=V, its slack after omitting the
common positive factor r_v S_vv is

-266112409230694467411733683/56456560891831241242814096384.

The quantitative pair bound remains valid. The failure concerns combining
its separate upper bounds at each triple; retaining alpha(mu_C A_x)(B)
and the nonnegative residual term in the equality may still matter.

The fixed-root confluence theorem cannot simply be stopped at a union
of pair events. For example, start active input {0,1} and stop restoration
when 0 and at least one of {1,2} are occupied. Take first records
0:(0,2), 1:(0,2), and records 2:(1,sink),(0,sink), with every other
record (1,sink). The original output is empty. Restoring root 0 first
ends at {0,2}, so the event is reached. Restoring root 1 first ends at
{1,2}; subsequently restoring 0 ends at {0,1}. Both follow permissible
missing-original-root choices, but their final emissions are respectively
(0,1,0) and (0,0,0). Thus event-stopped restoration is order-dependent.
A fixed policy might still be useful, but it needs its own image and
normalization argument: the boundary {0,2} has also forgotten original
root 1. `union_restoration_obstruction` in the pair-capacity verifier
checks the complete finite books. This is an obstruction to that direct
extension of the map, not a counterexample to the union capacity bound.

### N2 from killed sleep insertion

**Theorem.** Fix a target v. If p_z>=1/2 at every other vertex in the
strongly connected component of v, then

(mu A_x)(v occupied)>=mu(v occupied), mu in C_2,

for every source x. In particular N2 holds on arbitrary finite graphs
with all sleep probabilities at least one half. The target's own sleep
probability and probabilities outside its component may be arbitrary.
This concerns at most two common particles, including multiplicities;
it does not establish the comparison for larger common inputs.
Removing self-loops first also gives the conclusion under
r_z(1-K_zz)<=1 at those off-target sites: the effective odds are
r'_z=r_z(1-K_zz), and the active-addition and creation-cone operators
are unchanged by this reduction.

Use operators on sleeping supports avoiding v. Let Ahat_z be addition
killed on any visit to v, with Ahat_v=0, and set

Jhat_z=K_(z,s)I+sum_(y!=v)K_(z,y)Ahat_y,
Bhat_z=(Ahat_z-q_z Jhat_z)/p_z, z!=v.

The killed active additions commute. Therefore so do these operators,
including Jhat_v. On a support C avoiding z, first-step decomposition
gives delta_C Bhat_z=delta_(C union {z}). On a support containing z,
the identity Ahat_z^2=Ahat_z Jhat_z gives exactly the killed version
of (39):

delta_C Bhat_z=(1-r_z)delta_C Jhat_z
                +r_z delta_(C without z)Jhat_z^2.       (127c)

For example, expand Bhat_z^2=(1-r_z)Bhat_z Jhat_z+r_z Jhat_z^2 and
apply delta_(C without z). Thus Bhat_z is nonnegative and substochastic
when p_z>=1/2. This statement needs no restrictions on other sleep
probabilities: both Jhat_z terms are physical substochastic kernels.

Put nu=delta_empty Jhat_v, a subprobability supported on empty and
singleton supports. For D avoiding v, abbreviate w_D=w_v(D union {v}).
The deterministic insertion on vacant sites and commutation give
w_D=nu(product_(z in D)Bhat_z)1. Write d_z=1-Bhat_z 1. For distinct
a,b different from v,

w_{a}+w_{b}-w_empty-w_{a,b}
 =nu({a})delta_{a}Bhat_a d_b
  +nu({b})delta_{b}Bhat_b d_a.                           (127d)

To check it, condition nu on its at most one sleeper. At empty or a
singleton outside {a,b}, both insertions meet vacancies and the mixed
difference is zero. At {a}, the expression is Bhat_a(1-Bhat_b 1),
because Bhat_b has row mass one there. At {b}, commute the operators
and use the symmetric expression. If p_a,p_b>=1/2 both terms are
nonnegative. In the notation of (127), this proves

a_v({v,a,b})>=a_v({v,a})+a_v({v,b}).

For mu in C_2 put Delta_C=(mu A_x)(C occupied)-mu(C), which is
nonnegative for every pair C by (127a). Let N=mu L_v>=0 and
S_vv=1/[1+r_v w_v({v})]. Summing the deletion identity before and after
addition, and using A_x L_v=L_v Ahat_x, gives the exact decomposition

(mu A_x)(v occupied)-mu(v occupied)
 =S_vv N(1-Ahat_x 1)
  +sum_(C: |C|=2, v in C)a_v(C)Delta_C
  +sum_(B: |B|=3, v in B)(mu A_x)(B)
       [a_v(B)-sum_(C subset B, |C|=2, v in C)a_v(C)].   (127e)

Every term is nonnegative under the stated assumptions. For the component
restriction, note that a forced packet from v cannot return once it
leaves v's component. Sleeping sites outside that component therefore
do not change w_v; they can be replaced by the sink. If a triple contains
an outside site its bracket in (127e) is zero; otherwise (127d) applies.
This proves the theorem without a triple-allocation or source-weight
assumption.

The probability restriction enters precisely at the positivity in
(127c). The unweighted triple bracket can be negative at smaller p:
take 0->1,sink equally, 1->2 surely, 2->0 surely, and
p=(1/2,1/4,1/2), with v=0. In order D=empty,{1},{2},{1,2},

(w_D)=(13/16,9/16,5/8,1/2),
w_{1}+w_{2}-w_empty-w_{1,2}=-1/8.

This is a failure of that coefficient budget, not of N2; the three-site
theorem already covers this example. `high_sleep_budget_audit` in the
pair-capacity verifier checks the killed operator identities, commutation,
16 nonnegative two-background budgets, and 60 exact decompositions (127e)
on three fixed archived kernels, as well as this low-sleep obstruction.
The general theorem and its component restriction received a separate
AI audit. `quantitative_pair_audit` checks 198 exact identities (127b),
including repeated extra sources, and the stated rare-sleep failure.

**Bulk additions under the high-sleep hypothesis.** The same theorem
holds with A_x replaced by A_zeta for any finite nonnegative active
batch zeta, while the lower measure remains in C_2. The upper particle
count is unrestricted. This extension uses the high-sleep hypothesis
at every off-target vertex in the target component; the terminal
low-sleep extension below is not included.

First, (127a) holds with A_zeta: the response of the joint event C can
be negative only at C on supports of size at most |C|, and (47c)
applies to the all-active input C+zeta. Second, for any set D of
off-target vertices satisfying the high-sleep hypothesis,

w_empty-w_D>=sum_(z in D)(w_empty-w_{z}).                (127p)

Condition nu on its empty or singleton support. On empty or a singleton
outside D every Bhat insertion meets a vacancy, so all losses vanish.
On {a} with a in D, commute Bhat_a to the front. The remaining
nonnegative substochastic factors cannot increase the resulting row
mass. Among the individual insertions only Bhat_a has nonzero loss.
This proves (127p), hence
a_v(B)>=sum_(C subset B, |C|=2, v in C)a_v(C) for every B.
Sites outside the target component can again be discarded from w_v.
The exact decomposition (127e) now holds with A_zeta and with its last
sum over all |B|>=3; intertwining gives the first term
S_vv(mu L_v)(1-Ahat_zeta 1)>=0. Every pair-event gain and every final
bracket is nonnegative. This proves the bulk statement, which received
a separate AI audit. No induction through intermediate C_2 laws is used.

**Terminal low-sleep sites are also allowed.** The defect d_z=1-Bhat_z 1
is nonnegative at every sleep probability, even when Bhat_z has signed
entries. For z in C, write j_1=delta_C Jhat_z 1 and
j_2=delta_(C without z)Jhat_z^2 1. The released-resident comparison in
[AX Section 1](avalanche-exchange.md#1-a-comparison-for-released-resident-arrows)
gives j_2<=j_1. Explicitly, sample the first forced packet's arrow
independently of z's resident arrow. In the second system remove that
resident and hold its packet. Replay the first history, releasing the
held packet at the resident's first contact, or at the end if there is
none. Its complete visited set contains the first one. Compare before
killing, without conditioning on avoidance. Consequently

d_z(C)=1-j_1+r_z(j_1-j_2)>=0.

On supports avoiding z the defect is zero. This does not make Bhat_z
a substochastic kernel at low p: signed rows remain, and in the preceding
example d_1({1})=2.

If a cannot reach b without visiting v, both delta_{a}Jhat_a and
delta_empty Jhat_a^2 are supported on configurations avoiding b.
Equation (127c) then gives delta_{a}Bhat_a d_b=0, regardless of p_a.
The other summand in (127d) is nonnegative whenever Bhat_b>=0, since
d_a>=0 without a restriction on p_a. Hence the N2 cone theorem extends
to the following condition, with S the target's strongly connected
component:

for every a in S without v, either r_a(1-K_aa)<=1
or K(a,S without {v,a})=0.                              (127g)

Indeed every low-sleep vertex is terminal after deleting v, so each
possibly signed summand in (127d) vanishes. Self-loops and exits from S
are allowed. A terminal component containing several vertices is not
enough for this argument; the condition is on individual vertices.
In particular the central target of a bidirectional star is covered
with arbitrary sleep probabilities and leaks. The released-resident
and terminal-site extension received a separate AI audit. The verifier
retains the negative nonterminal example and a positive terminal
variant, with p=(1/2,2/3,1/4), whose bracket is 1/24.

The exact remaining low-sleep cancellation can retain the minorization
remainders instead of discarding them. For each pair C define

epsilon_C=(sigma_C A_x)(C)
 +sum_(B contains C)beta_C(B)(sigma_C A_x)(B)>=0,
sigma_C=mu-[mu(C)/mu_C(C)]mu_C.

Then (127b) gives

(mu A_x)(v occupied)-mu(v occupied)
 =S_vv N(1-Ahat_x 1)+sum_(C contains v)a_v(C)epsilon_C
  +sum_(B: |B|=3, v in B)(mu A_x)(B)
       [a_v(B)-sum_(C subset B, |C|=2, v in C)
                       a_v(C)beta_C(B)].               (127f)

The first two terms are nonnegative at every sleep probability. The
brackets need not be. A general proof would have to bound their negative
contribution by these retained terms or use a different decomposition;
their individual nonnegativity is not asserted.

Retaining only the triple part of epsilon_C is still insufficient on
the full cone C_2. For the same rare-sleep three-site kernel, target 0
and source 2, a rational cone measure makes

a_v(V)(mu A_2)(V)
 -sum_(C subset V, |C|=2, 0 in C)
       a_v(C)beta_C(V)[mu(C)/mu_C(C)](mu_C A_2)(V)<0.

After removing the common positive factor r_v S_vv, the verifier's
exact feasible witness gives
-8811470822597717846819589644771241224919/
 11046683866931891518105167894089536955986428385262336.
Its masses and all cone inequalities are reproduced by
`quantitative_pair_audit`. The full occupation response is not negative;
the exact-pair residuals or deletion-residual decrease cannot generally
be discarded. For actual two-root inputs with a distinct extra source,
this weaker source-weighted budget also fails, by the feeder example
below. In its explicit form, for lower binary roots A, x outside A and each
triple B containing v, it asks

a_v(B)mu_(A+x)(B)
 >=sum_(C subset B, |C|=2, v in C)
       a_v(C)mu_A(C)m_C(B)mu_(C+x)(B).                  (127i)

This would prove N2 by (127b), keeping the exact-pair remainder and
the no-visit gain as additional nonnegative terms. It retains both the
original input and the actual extra source; it is not a bound on entrance
masses alone. The maintained fixed-kernel audit checks 297 such budgets:
Section 12's four-cycle, the archived four-site jump kernel with all
p_i=1/1001, and the rare-sleep three-site kernel. All pass. This is finite
evidence only; (127i) is false for actual creation laws as well as for
the larger cone. The retained exact-pair remainders cannot be omitted.

**A source maximum principle.** Ordinary hitting probabilities control
some of the source weights in (127i). Let B_y=A_y+r_y(A_y-J_y) be the
full signed sleep-insertion operator, let U_y insert a sleeper at y on
supports avoiding y, and let A_y^- consist of those rows of A_y. Then

B_y-A_y=L_y(U_y-A_y^-).                                (127j)

For an all-active input eta avoiding y, compare adding a sleeping y
with holding an active y. Process eta first. At its first contact with
y the two systems agree. If it never contacts y, its residual law is
N_eta^y=mu_eta L_y and the two remaining continuations are U_y and
A_y^-. This proves (127j) on all such creation rows. When eta already
contains an active y, both sides vanish, using A_y B_y=A_y^2 and
mu_eta L_y=0. Binary creation rows form the invertible matrix P, so the
operator identity follows. No positivity of B_y is required.

Fix an exact terminal support C and a measure mu in C. Define the
source function f(y)=(mu A_y)(C), and f(s)=mu(C). For y outside C,
the insertion term cannot end at C, and (127j) gives

f(y)-sum_z K(y,z)f(z)
 =-r_y^(-1)[(mu L_y)A_y^-](C)<=0.                       (127k)

The function is also harmonic at every y with mu L_y=0. For an actual
all-active law mu_eta, put D=C without supp(eta). Stop the ordinary
jump walk on first reaching D or s. Bounded optional stopping gives

mu_(eta+x)(C)
 <=P_x(T_s<T_D)mu_eta(C)
   +sum_(j in D)P_x(X_(T_D)=j,T_D<T_s)mu_(eta+j)(C).

Here the sink-only stopping rule is used when D is empty. In particular,
exact supports contained in the original active support have decreasing
mass under an added particle. Taking C empty also proves, for every
mu in the full creation cone, mu A_x(empty)<=mu(empty): the global
nonempty event increases, without a particle-count restriction. This is
not a comparison of nonemptiness in an arbitrary proper subregion.
It has an explicit positive residual certificate. With G=(I-K_V)^(-1),

mu(empty)-(mu A_x)(empty)
 =sum_y [G(x,y)/r_y][(mu L_y)A_y^-](empty)>=0.            (127n)

This is the Green representation of (127k) when C is empty; the boundary
value at the sink is mu(empty). It uses the full cone and permits
arbitrarily many initial particles.

For binary common roots R and j outside R, set C=R union {j}. Since
mu_R(C)=0, the preceding bound becomes

mu_(R+i)(R union {j})/mu_(R+j)(R union {j})
 <=h_i(j):=P_i(T_j<T_s), for every source i.             (127l)

It permits a repeated source i in R. For a genuine replacement i
outside R union {j}, the inequality is strict whenever h_i(j)>0:
the no-i-visit history that immediately sleeps R, followed by an active
i routed to j while retaining those residents, makes the defect in
(127k) strictly positive at i. The killed walk's Green representation
then gives a strict gap.

Iterating the same comparison gives a geometric bound for an entire
replacement of roots. Let |eta|=|C| and run ordinary IDLA with parking
capacity one at C and zero outside C, killing particles at the original
sink. Write H_C(eta) for the probability that it fills all of C. Then

mu_eta(C)/mu_C(C)<=H_C(eta).                             (127o)

Hold one input particle at each currently occupied C-site. Select an
outside particle or a duplicate at a C-site, leaving common input eta'
of size |C|-1. The source comparison stops its ordinary walk at the
missing sites C without supp(eta') or the sink; the sink contribution
is zero by particle count. A successful routing increases the number
of distinct input sites in C. Iterate until C is filled or a particle
is lost. The accumulated success probability is precisely H_C(eta),
and the only nonzero terminal creation probability is mu_C(C). This
uses forcing of initial source locations through the maximum principle;
it does not assert a pathwise occupied-set coupling of ARW outputs.

There is also an order-independent lower bound for actual creation
laws. Given the target C, form a substochastic kernel Ktilde by killing
a walk with probability p_y before each ordinary jump at y outside C;
leave rows inside C unchanged. Sleep deaths are merged into the sink.
Let H_C^K and H_C^Ktilde be the corresponding target-capacity IDLA
filling probabilities. Then, for every |eta|=|C|,

H_C^Ktilde(eta)<=mu_eta(C)/mu_C(C)<=H_C^K(eta).           (127q)

For the lower inequality, let eta' be the held common input during the
routing proof of (127o), and f(y)=mu_(eta'+y)(C). The physical law with
an additional sleeper gives mu_eta' B_y>=0, so the first-step equation
implies f(y)>=q_y sum_z K(y,z)f(z). At each already occupied input
site f is harmonic. Until the walk reaches a missing target site or
the sink, every site inside C is such a held input site. Thus f is
superharmonic for Ktilde on the routing interior. The sink contribution
is zero by particle count. Iterate this lower comparison through the
same legal target-filling policy as before. Ordinary IDLA abelianness
for Ktilde identifies its success probability independently of policy.
The sharper version which omits sleep killing at other held input
sites is also valid for a specified policy; its order independence
is not asserted.

This sandwich gives a uniform limit. Put epsilon=max_(y outside C)p_y
(zero if C=V), G=(I-K_V)^(-1), and M=max_z sum_y G(z,y). Couple the
two IDLA processes with the same walks, adding independent sleep
deaths outside C. At most |C| walks are run, each stopped no later
than its original sink time. A union bound gives

0<=H_C^K(eta)-mu_eta(C)/mu_C(C)<=|C| epsilon M.

Consequently the normalized probability tends to H_C^K(eta) as all
outside-C sleep probabilities tend to zero, uniformly in the sleep
probabilities on C and the input eta of this fixed particle count.
The first-contact, maximum-principle, Green, replacement, and IDLA
arguments (127j--q) received an independent AI proof audit.

For example, let A={v,a}, x=b, C={v,b}, B={v,a,b}. Applying (127l)
first with common root {v}, then with common roots {v,b}, yields

[mu_A(C)/mu_C(C)] [mu_(C+b)(B)/mu_B(B)]
 <=h_a(b)h_b(a).                                       (127m)

Thus the other pair's source weight is suppressed by an ordinary
round-trip probability when the upper input is this triple. This
controls that triple's term; it does not by itself handle all other
possible terminal triples. The first-contact identity has an exact
1,024-entry check on the archived four-site kernel. The maintained
`source_maximum_principle_audit` also checks 64 empty-event Green entries,
128 source bounds, 48 strict replacements, 324 exact-support decreases,
and 70 target-capacity IDLA bounds, including reversed routing orders.

**Two uniform-parameter source repairs.** The ordinary hitting factor
can resolve examples in which the pointwise beta budget is negative.
On three sites V={v,a,b}, take A={v,a}, extra source b, and C={v,b}.
Write u_D=w_v({v})-w_v(D), tau=mu_V(V), alpha=mu_A(C)/mu_C(C),
kappa=h_b(a), and lambda=h_a(b)h_b(a). Surplus IDLA gives
mu_(C+b)=(1-kappa)mu_C+kappa mu_V. Thus the actual cross-source weight
gamma=alpha kappa is at most lambda. The original pair's minorization
remainder is zero. With rho_C=(sigma_C A_b)(C)>=0, (127b) gives

[q_v(V)-q_v(A)]/(r_v S_vv)=u_C rho_C+tau T,
T=u_V-u_A beta_A(V)-u_C beta_C(V)alpha kappa.

For bar_alpha=h_a(b), this has the positive decomposition

T=[u_V-u_A-lambda u_C]
  +u_A[1-beta_A(V)]+u_C kappa[bar_alpha-alpha]
  +u_C alpha kappa[1-beta_C(V)].

All terms except the first are nonnegative by the general theorems.
The following two fixed kernels have positive first terms throughout
their uniform-sleep parameter ranges. These are structural certificates
for these comparisons, not a proof of (127i) on arbitrary graphs.

First use Section 22's rare-sleep kernel, now with common odds r>0 and
target 0. Its ordinary hitting probabilities are h_1(2)=1/4 and
h_2(1)=3/34, so lambda=3/136. Direct processing of the killed domain
gives

u_{01}=3r(r+2)/[512(1+r)^2],
u_{02}=3(205r^2+5574r+5120)/[32768(1+r)^2],
u_V=u_{02}+u_{01}(1/r+7/128).

For clarity these are exact two-buffer identities. Put p=1/(1+r),
q=1-p, c=p+q/2, s=1/64 and t=3/64. Starting active 2 with initial
sleep indicators i at 1 and j at 2, its no-0-visit probability is

F_(i,0)=p+q(s+t c/2^i),
F_(i,1)=p(s+t c/2^i)+q[s^2+2st c/2^i+t^2 c/2^(i+1)].

The initial forced 0-jump enters 2 with probability one half, so
u_{01}=(F_00-F_10)/2, u_{02}=(F_00-F_01)/2 and
u_V=(F_00-F_11)/2. The two positive first terms are

u_V-u_{02}-lambda u_{01}
 =3(r+2)(71r+2176)/[1114112(1+r)^2],
u_V-u_{01}-lambda u_{02}
 =3(19037r^2+733590r+698368)/[4456448(1+r)^2].

Hence both actual root comparisons with target 0 hold by this source
budget for every r>0, even though their source-free beta triangle is
negative at r=1000. The verifier checks these small polynomial identities
and independently audits the two residual components and occupation
gains with a 28-state raw solve at r=1000; no large beta numerator is
needed for the positive proof.

Second, even a directed cycle can have a negative source-free beta
triangle. Let 0 jump to 1 with probability 1/4 and to the sink otherwise,
1 jump surely to 2, and 2 jump to 0 or the sink equally. Put all sleep
probabilities equal to p. At p=1/10, target 0,

(w_empty,w_{1},w_{2},w_{1,2})
 =(719/800,1321/1600,1339/1600,2521/3200),
beta_{01}(V)=7618127/8649734,
beta_{02}(V)=50085190/56223271.

The weighted slack is -772388791/89957233600. An exact univariate
certificate shows it remains negative for 0<p<=1/10 and tends to
-1/64 as p decreases to zero; both beta ratios tend to 7/8. Thus the
failed entrywise bound cannot be rescued by assuming a directed cycle.
For the actual source comparison, however, lambda=1/8 and q=1-p give

u_{01}=q(1+3p)/16, u_{02}=q(1+p)/16,
u_V=(3+6p-5p^2)/32,
u_V-u_{01}-lambda u_{02}=(3+8p+5p^2)/128>0,
u_V-u_{02}-lambda u_{01}=(3+13p+9pq)/128>0.

The same decomposition therefore proves both specified actual-source
budgets for every 0<p<1. The maintained cycle audit checks the exact
raw and pair-renewal fractions, the negative interval certificate, and
these positive quadratic identities. The arbitrary-cycle N2 extension
is not inferred from this one kernel.

**Actual-source obstructions to (127i).** A single-leak three-cycle
already suffices. Take v->a,sink equally, a->b surely, b->v surely,
with (p_v,p_a,p_b)=(1/2,t,1/2). Put A={v,b}, extra source a,
C={v,a}, B={v,a,b}, and D=t^2+15t+40. Exact pair renewal and
the killed two-buffer calculation give

u_A=(1-t)/4, u_C=1/4, u_B=(1+t)/4,
beta_A(B)=(36-t)/D,
beta_C(B)=2(t^2+13t+10)/D,
gamma=[mu_A(C)/mu_C(C)][mu_(C+a)(B)/mu_B(B)]=2/(t+7).

The duplicate at a goes directly to vacant b, so mu_(C+a)=mu_B as
whole laws; gamma is the actual source factor. The slack normalized
by r_v S_vv mu_B(B) is

T(t)=u_B-u_A beta_A-gamma u_C beta_C
 =[t^4+22t^3+193t^2+596t-12]/[4(t+7)D]<0,
 0<t<=1/50.

The numerator is increasing on t>0 and equals -16399/6250000 at
t=1/50. At t=1/100 the normalized slack is
-602067799/112580880400; its boundary limit is -3/280.
The exact-pair remainder rho_C=(sigma_C A_a)(C) restores the budget:

u_C rho_C/mu_B(B)=(9-t)(t+5)/[(t+7)D],
T(t)+u_C rho_C/mu_B(B)=(t^3+15t^2+84t+24)/(4D)>0.

The maintained cycle audit verifies the polynomial identities and the
failure interval, and independently checks the renewal, inverse boundary
ratios, and duplicate-source identity by 28-state raw absorption.
The occupation gain at t=1/100 is 49683002/112862332601>0.

A second example shows how a feeder enforces the same constraint.
Let the vertices be
v,a,b,z, with v->a,sink equally, a->b surely, b->v,sink equally,
and z->a surely. Put p_v=p_b=p_z=1/2 and p_a=t. Take lower input
A={v,z}, extra source a, and terminal triple B={v,a,b}. Thus the
extra source is distinct from the original roots and is the only
low off-target site when t<1/2.

Set R={v,a}, C={v,b}, pi=mu_B and Q=t^2+31t+208. On exact supports
omitting z, mu_A=mu_R/2: sleeping z remains forever, while its jump
branch enters a. Surplus IDLA gives the whole-law identities
mu_(R+a)=mu_(C+a)=pi. The three-site renewal and killed two-buffer
calculation give

pi(B)=4t/(12+t),
rho=mu_R(C)/mu_C(C)=(1-t)(t+16)/(16-t),
beta_R(B)=(5t^2+115t+312)/(2Q),
beta_C(B)=(7t+188)/Q,
u_R=(3-t)/16, u_C=3(1-t)/16, u_B=(9+t)/32,

where a_v(D)=[8/(15+t)]u_D. Consequently the slack in (127i),
divided by mu_(A+a)(B)=pi(B)/2 and the common factor 8/(15+t), is

E(t)=u_B-u_R beta_R-u_C rho beta_C
 =-(3t^4+110t^3+794t^2-2539t+192)
   /[2(16-t)(t^2+31t+208)]<0, 0<t<=1/16.

The numerator is at least 192-2539/16=533/16 on this interval.
At t=1/16, E=-53071/9745760 and the unscaled inequality has slack
-424568/14165766715. The target gain itself is positive.

The omitted exact-pair remainder supplies the missing compensation.
For sigma=mu_R-rho mu_C, the same whole-law identities give
sigma A_a=(1-rho)pi. In the feeder example the pair-C remainder is
therefore (1-rho)pi(C)/2. Including it changes the normalized budget to

E(t)+u_C(1-rho)pi(C)/pi(B)
 =(2t^2+43t+51)/[2(t^2+31t+208)]>0.

The two feeder branches also give the full target gain directly as
[mu_B(v occupied)-mu_{ {v} }(v occupied)]/2>=0 by root retention;
at t=1/16 it equals 7922477/499968237. This is a failure of dropping
pair remainders, not of N2. More generally, a new vertex with no
incoming arrows and a jump to an old root converts a repeated-source
comparison on the old graph into a distinct-source comparison on
terminal supports omitting the new vertex. Distinct-source budgets
therefore cannot ignore these repeated-source constraints.

**Two local lower bounds for the retained pair credit.** Let the added
source be a, and for each pair C put tau_C=sigma_C A_a. The exact
operator identity A_a L_a=0 gives tau_C L_a=0, even though sigma_C
need not belong to the creation cone. For C avoiding a and B=C union {a},

tau_C(C)=r_a sum_(T contains a)tau_C(T)W_a(T,C)
        >=r_a W_a(B,C)tau_C(B).                         (127x)

Different triples containing v and a use different pairs C=B without a,
so these lower bounds can be summed without reusing pair credit.
For the pair D={v,a} containing the source, the immediate-sleep event
instead gives tau_D(D)>=p_a sigma_D({v}). This one pair's credit is
shared by all relevant triples and must be budgeted only once.

For example, the feeder kernel above has W_a(B,{v,b})=5/16, independent
of t. At t=1/16, the smaller credit from (127x), normalized by the
upper triple mass, is 50625/1114112; its sum with E(t) is
798257467/19959316480>0. The maintained audit proves positivity of
this smaller repair throughout 0<t<=1/2 by a positive Bernstein
expansion. In the single-leak three-cycle at t=1/100, the other pair
D={v,a} supplies the repair: sigma_D({v})=80000/491401, so its direct
sleep contributes normalized credit 40100/491401, exceeding -T(t).
These local lower bounds do not suffice to repair every triple. Change
the feeder example to b->v surely, leaving all other routes and
p_v=p_b=p_z=1/2, p_a=1/16 unchanged. For the same A,R,C,B, the
normalized source budget is -83033/1331087. The containing-source
remainder tau_R(R), including its direct-sleep lower bound, is zero.
The self-triple credit is 21825/1040384, leaving
-451458511/10904264704<0. The full pair-C credit is 174600/1331087
and gives repaired slack 721/10481>0. Here the complete deletion is

tau_C(C)=15[tau_C({v,a})/2+tau_C({a,b})/4+tau_C(B)/8].

The two source-containing pair terms cannot be dropped; even
sigma_C({v})=sigma_C({b})=0, so direct-sleep bounds on their singleton
predecessors give no replacement. A maintained 51-state raw audit
reproduces these fractions and the positive target gain
2167809/76982945. This obstructs that local repair, not the full
residual budget or transfers between triples.

**One-site continuation by repeated resident releases.** There is an exact
positive representation of changing only one sleep probability. Fix base
odds r^0 and a vertex a, keep all other odds fixed, and increase r_a to
r_a^1>=r_a^0. All operators with superscript 0 use the base parameters.
Let R_a^0 leave a support avoiding a fixed; from C containing a, remove
its resident, force one K-step, and stabilize at the base parameters:

delta_C R_a^0=delta_(C without a)J_a^0, a in C.

This is a Markov kernel which never increases particle count. Put
t=p_a^0(r_a^1-r_a^0). Then the parameter-transfer kernel is

E_0 P_1=[I+t(I-R_a^0)]^(-1)
 =(1/(1+t))sum_(k>=0)(t/(1+t))^k (R_a^0)^k.             (127h)

To prove it, use that each row of the complete-block inverse E is affine
in r_a in the surplus-IDLA coordinates. Rows avoiding a have derivative
zero. For a row C containing a, multiplication by P_0 converts its
derivative to

delta_(C without a)(A_a^0-J_a^0)
 =p_a^0[delta_C-delta_(C without a)J_a^0].

The equality is the vacant-site first-step identity. Hence
E_1 P_0=I+t(I-R_a^0), and inversion proves (127h); the geometric series
converges because R_a^0 is stochastic. Thus one can sample an independent
geometric number of forced resident releases at a, each time using the
current support and fresh base-parameter books. The kernel is the
previously proved parameter transfer, now expressed entirely at the
base probability. This does not transfer the occupation sign by itself:
the resulting target payoff need not be an increasing function of the
initial sleeping support. The maintained verifier checks all 256 entries
on the archived four-site kernel, changing p_a from 1/2 to 1/4.

The source constraint is preserved by this parameter flow, but does
not close a scalar sign inequality. Fix a!=v, vary r=r_a, write
p=1/(1+r), q=1-p, delta_r=mu_(A+a)^r-mu_A^r, and D(r)=delta_r f_v.
The lower no-a-visit law N=N_A^a and the forced no-return kernel W_a
are independent of r. Thus L_a=P_avoid-rW_a and

delta_r L_a=-N, delta'_r=p delta_r(R_a-I),
delta'_r L_a=delta_r W_a.

In fact p(R_a-I)L_a=W_a, so the differentiated constraint is already
automatic. Let u_r be delta_r restricted to a-containing supports,
J_r=delta_r 1{a and v occupied}, and k_r(C) the probability that a
is visited during the forced-resident-release continuation from C and
v is finally occupied. Put C_r=u_r k_r. Splitting that continuation
according to whether it visits a gives the exact scalar identity

rD'(r)=p[D(r)+Nf_v]-J_r+q C_r.                         (127h')

The missing comparison at a hypothetical zero of D is q C_r+pNf_v>=J_r.
It concerns the signed terminal background weighted by the target
reward after a return, and is not implied by the no-return constraint.
The stronger sufficient inequality giving [(1+r)D]' >=0 is false:
on the single-leak three-cycle above, with t=p_a,

D/t=2(t^3+15t^2+84t+24)/[(t+7)(t+4)(t^2+15t+40)],
partial_t(D/t)|_(t=0)=459/3920>0.

Thus D/t decreases with r near the zero-sleep boundary despite D>0.
The exact parameter and return identities received a separate AI audit;
they do not prove continuation of the high-sleep sign.

**A positive source-switch reduction with a common target root.**
Remove self-loops, fix v, and write mu_z=mu_{ {v,z} }, z!=v,
f_v(C)=1{v in C}, and D_z(x)=(mu_z A_x-mu_z)f_v, with D_z(s)=0.
For x outside {v,z} with p_x>=1/2, set

alpha_(z,x)=mu_z({v,x})/mu_x({v,x}),
sigma_(z,x)=mu_z-alpha_(z,x)mu_x.

Minorization gives sigma>=0, sigma({v,x})=0 and 0<=alpha<=1.
The ordinary sleep insertion B_x is stochastic at these parameters.
On supports of size at most two, (B_x-I)f_v can be negative only at
{v,x}: inserting at a vacant x preserves the target bit, and a
target-free support has nonnegative response by positivity of B_x.
Hence h_(z,x)=sigma_(z,x)(B_x-I)f_v>=0. Since x is a common root
of mu_x, mu_x B_x=mu_x A_x=mu_x J_x. First-step decomposition gives

D_z(x)=q_x sum_y K(x,y)D_z(y)
        +p_x alpha_(z,x)D_x(x)+p_x h_(z,x).             (127r)

At either common root x in {v,z}, D_z(x)=sum_y K(x,y)D_z(y),
with no restriction on p_x. These equations are a positive renewal
on states (z,x). Treat (z,a) with a outside {v,z} and p_a<1/2 as
boundary states. Every interior row has total transition weight at
most one. Following a positive K-path reaches either such a boundary
or the sink with positive probability; every q_x is positive. Finiteness
therefore makes the interior kernel Q transient, and

D_interior=(I-Q)^(-1)(h+H D_boundary),
(I-Q)^(-1), h, H>=0.

Thus, with exactly one low off-target site a, it suffices for all
common-target two-root comparisons to prove those with a as the extra
source and a outside the common pair. This is a reduction, not a proof
of that remaining boundary sign.

The source-switch identity holds more generally for every normalized
mu in C_2, using rho_x=mu_{ {v,x} }, alpha=mu({v,x})/rho_x({v,x})
and sigma=mu-alpha rho_x. At x!=v with p_x>=1/2 the only dangerous
support is still {v,x}, so the remainder is scored nonnegatively at
once; it need not lie in the cone. For an actual pair input A, common
root sources remain harmonic. The target source has its own positive
row at every sleep probability and on the full creation cone:

D_mu(v)=sum_y K(v,y)D_mu(y)+c_mu(v),
c_mu(v)=r_v^(-1)(mu L_v)(1-A_v^- f_v)>=0.               (127s)

Indeed divide (127j) at y=v by r_v and use U_v f_v=1. Combining
these rows gives a finite positive renewal on labels (A,x), or on
the original cone measure together with the pair laws rho_x. The
only remaining boundary states have a low-sleep off-target extra
source that is not a common root. The target-source boundary is
eliminated without assuming its response already nonnegative. Every
switch enters a pair containing v; subsequent switches retain v.
Thus the common-target reduction extends to all two-root inputs.
With no low off-target sites this is a second proof of the high-sleep
N2 theorem, though it does not give its bulk-addition extension.

Intermediate resident-release laws
cannot be presumed to lie in the reference creation cone: on
v->a,sink equally, a->b surely, b->v,sink equally, all reference
probabilities 1/2, lambda=mu_{ {v,b} } R_a satisfies
(lambda L_v)({a,b})=-66/895. The geometric representation (127h)
therefore does not directly permit repeated application of the
high-sleep cone theorem.

### A small-sleep theorem under a hitting hypothesis

**Proved; two independent AI proof audits passed.** Let K be a finite,
sink-accessible jump kernel whose internal graph is strongly connected.
Suppose

h_i(j)=P_i(T_j<T_s)<1 for every distinct i,j.            (127t)

Equivalently, after deleting any one vertex, every remaining vertex
still has a path to the sink. There is delta(K)>0 such that the
single-target monotonicity conjecture holds for every finite all-active
input whenever 0<p_i<delta(K) for all i. The sleep probabilities need
not be comparable. The threshold is uniform over inputs and targets
on this fixed graph; no bound uniform in graph size is asserted.
Every binary comparison A to A union {x}, x outside A, is strict.

Here is the proof. Write G=(I-K_V)^(-1) and p_C=product_(c in C)p_c.
For any fixed finite input eta and terminal support C, the raw finite
absorption problem is analytic in all p_i in a neighborhood of zero.
At zero sleep no new sleeper is created. The potential
sum_i n_i E_i T_s, counting both active and sleeping particles, decreases
in conditional expectation by one per jump. Its finite initial value
bounds the expected number of jumps until no active packet remains.
Thus the transient matrix is invertible at zero. Since occupation of C requires
a sleep at each of its vertices, mu_eta(C) vanishes identically when
any p_c=0 and its analytic extension is divisible by p_C.

Set the outside-C sleep probabilities to zero. Trace every active
excursion onto C, merging excursions lost to the sink into killing.
The internal trace kernel K^C satisfies

(I-K^C)^(-1)=G_CC.

The remaining ARW has unchanged sleep probabilities on C. Surplus IDLA
with capacity one on C and zero outside C first produces an active
binary subset of C, with fresh suffixes. Only the full subset can end
with all |C| sites occupied. Consequently, for any particle count,

mu_eta(C)=H_C^K(eta)/Z_C^trace, when p_(V without C)=0.

Here H_C^K is the probability that ordinary target-capacity IDLA fills
C, processing all particles even after losses, and Z_C^trace is the full
inverse diagonal (44) for the trace graph. In particular the full
active input C is stationary for subsequent additions on that graph.
The inverse-diagonal formula gives
p_C Z_C^trace -> det(I-K^C)=1/det(G_CC). Analytic divisibility therefore
identifies the leading coefficient in all parameter directions:

mu_eta(C)=p_C[det(G_CC)H_C^K(eta)+O(||p||)].             (127u)

This is an asymptotic for each fixed input; uniformity over all input
sizes will follow from a separate finite binary reduction below.

Fix binary A, x outside A and target v, and put
F_(A,x,v)(p)=[q_v(A union {x})-q_v(A)]/p_v. Analyticity and the
necessary terminal sleep at v extend F analytically through p_v=0.
If v is outside A, the singleton case of (127u) gives

F_(A,x,v)(0)=G_vv h_x(v) product_(a in A)[1-h_a(v)]>0.   (127v)

Indeed target-capacity IDLA on {v} fills it iff at least one of the
independent source walks visits v. Strong connectivity makes h_x(v)>0,
and (127t) makes every escape factor strictly positive.

Now suppose v is in A. In (127) both no-visit probabilities vanish,
so the terms S_vv[1-n_v] cancel between the inputs. The factors
a_v(B)=(1-p_v)[w_v({v})-w_v(B)]/
[p_v+(1-p_v)w_v({v})] are analytic near
zero and their common prefactor tends to G_vv. For c!=v let
rho_(v,c) be the probability that a forced departure from v first
reaches c before returning to v or hitting the sink. At zero sleep,
adding one resident at c changes the forced packet's target-hitting
probability by

u_(v,c)^0=rho_(v,c) h_c(v)[1-h_c(v)].

On first reaching c there are two ordinary active walks instead of
one; their probabilities of both avoiding v multiply. The two-site
trace inverse also gives rho_(v,c)det(G_{ {v,c} })=G_vc. Combining
these identities with the pair case of (127u), while supports of size
at least three contribute only higher total degree, yields

partial_(p_c) F_(A,x,v)(0)
 =G_vc G_cv[1-h_c(v)]h_x(c)
     product_(a in A without {v})[1-h_a(c)]>0,
     c outside A.                                    (127w)

For this IDLA difference, park the common particle at v first. The
failure probability to fill c is the product of its avoidance
probabilities from the other sources; the extra walk multiplies that
failure probability by 1-h_x(c). The expression in (127w) is strictly
positive by the same two graph assumptions.

The whole small parameter cube, rather than only fixed rays, follows
from an exact vanishing face. Keep arbitrary sufficiently small
probabilities on A and set all outside-A probabilities to zero. Trace
onto A. Its full all-active law is stationary, so the extra particle,
whether killed before reaching A or added at an A-site, does not change
the output law. Thus F vanishes identically on this face. By continuity,
all derivatives in (127w) remain positive on a sufficiently small
parameter cube. Holding p_A fixed, turn on the outside coordinates
one at a time and integrate these derivatives from the vanishing face.
This proves F>0 throughout the positive cube. In the case v outside A,
continuity of (127v) directly proves the same conclusion.

There are finitely many binary triples (A,x,v). Take the minimum of
their positive cube sizes. Every binary inclusion then increases each
target probability, by successive one-site additions. Finally surplus
IDLA represents every finite all-active input as a mixture of binary
creation laws with fresh suffixes. Ordinary IDLA occupied supports can
be coupled increasingly after an extra particle. The already proved
order of their binary continuation probabilities proves the theorem
for all finite inputs with this same delta(K).

Section 13 also localizes this theorem. For a fixed target component S,
kill jumps on exit from S and assume (127t) for this local kernel.
Small enough sleep probabilities on S then suffice for monotonicity
of its target in the full graph, with arbitrary sleep probabilities
and component sizes outside S. The local theorem handles all finite
conditional arrival vectors, so its threshold does not depend on the
upstream input or emission counts.

There is a weaker hypothesis for a single specified target v:

h_i(d)<1 for i in S without {v} and d in S without {i}. (127y)

Only the target itself may be forced through another transient vertex
before reaching the sink. Equivalently, deleting one vertex can
disconnect at most v from the sink. Equations (127v,w) show directly
that this hypothesis makes all coefficients required for this target
strictly positive. The same proof therefore gives its all-input
monotonicity on a small heterogeneous cube, with arbitrary parameters
outside S. This target-specific extension received a separate AI audit.

The vanishing coefficients identify the remaining obstruction exactly.
Say d dominates i toward the sink if h_i(d)=1. For a common-target
binary A, the quadratic coefficient at a missing site c is positive
iff c is not dominated by v and c dominates no off-target root in A.
Thus every quadratic coefficient vanishes precisely when each missing
site either is dominated by v or dominates an off-target common root.
If an off-target i has a proper dominator d, a full-minus-one input
already exposes such a zero: take A=S without {i} when d=v, and
A=S without {d} otherwise. This is a limitation of the first/second
degree proof, not a zero-response theorem or a counterexample.
One positive quadratic coefficient suffices along each fixed positive
scaling ray, but does not alone control a heterogeneous cube when
other outside-coordinate derivatives vanish.

The threshold can be made explicit; its numerical derivative bound
received an independent AI audit. For n=|V|>=2 under (127t), put

M=max_i E_i T_s, h=min_(i!=j)h_i(j),
e=min_(i!=j)[1-h_i(j)], kappa=h^3 e^(n-1).

Then one valid choice is

delta(K)=kappa/[768 n^4 M^3].                           (127z)

For the fixed-target version (127y), replace e by
e_v=min_(i!=v, d!=i)[1-h_i(d)]. For n=1 the assertion holds at every
sleep probability. This bound is small and is not uniform over graphs.

To verify the bound, include every transient raw state with at most n
particles. At parameters p_i<=1/2, the same potential drops by at
least one half per step, so the fundamental matrix has row-sum norm
at most L=2nM. For a terminal occupation payoff, its first parameter
derivative solves a Poisson equation whose forcing is the difference
between a sleep continuation and a jump continuation. Those continuation
values lie in [0,1], so the forcing has norm at most one. Each derivative
Q_i of the transient matrix has row-sum absolute norm at most two.
Differentiating the Poisson equation again gives bounds L, 4L^2 and
24L^3 for first, second and third mixed derivatives of the payoff.

The difference D of two initial rows doubles these bounds. Since
D vanishes on p_v=0,
F=D/p_v=int_0^1 partial_v D(t p_v,p_(V without v))dt.
For a cube of side delta<=1/2, integrating parameter derivatives gives

|F(p)-F(0)|<=32 n^3 M^2 delta,
|partial_c F(p)-partial_c F(0)|<=384 n^4 M^3 delta,
c!=v.

Every positive coefficient required in (127v,w) is at least kappa:
use G_ii>=1, G_ij=h_i(j)G_jj, and at most n-1 escape factors.
The choice (127z) keeps the second displayed variation at most
kappa/2 and the first at most kappa/(24nM). Thus it supplies all the
strict signs used in the cube proof. The derivative estimates also
hold on boundary faces by the same finite-matrix formulas.

The maintained `rare_sleep_expansion_audit` independently differentiates
the raw transient absorption matrix at zero, without assuming (127u).
On a four-site kernel with positive direct leaks and on a four-site
nearest-neighbor path killed beyond its two endpoints, it verifies
256 first derivatives and 384 second derivatives against (127v,w).
The two systems have 118 and 101 transient states. Exact zero-outside
calculations also check the trace diagonal, 32 target-IDLA identities
including unequal particle counts, and stationarity on the vanishing
face. These are checks of the saved proof, not a substitute for it.

### Fixed-ray extensions with outstanding audit

The following derivation retains its original reported status. Its general statements await an independent proof audit, and concern fixed rays rather than a heterogeneous cube.

#### Single-leak cycles along fixed positive sleep-rate rays

Scope. Sites 0,...,n-1 follow this order, with i->i+1 surely for i<n-1,
and n-1->0 with probability c, to the sink with probability d=1-c.
Assume 0<c<1. Put p_i=epsilon lambda_i with every lambda_i>0 fixed.
All assertions below are as epsilon decreases to zero on this fixed ray.
They do not assert a uniform small heterogeneous parameter cube.

Let A be binary, x outside A, k=|A|, and Delta_v=q_v(A+x)-q_v(A).
Write r=min{i:i not in A}; thus 0,...,r-1 are precisely the initial
consecutive roots before the first vacancy. Put h_x(j)=1 if x<=j,
and h_x(j)=c if x>j. This is the ordinary hitting probability of j.

1. A general-target prefix case.

Let b=min(A\{0,...,r-1}), or b=n when the set is empty. If v<b, then

Delta_v = epsilon^(r+1) C + O(epsilon^(r+2)),
C = d^(k-r-1) sum_{0<=i_0<...<i_r<b, v in {i_0,...,i_r}}
       h_x(i_r) product_t lambda_(i_t)
       [1{v=i_r}+c 1{v<i_r}] > 0.

For v>=b this statement makes no claim about the leading term.

2. Every binary input when the target is the leaky vertex v=n-1.

Set L_i=sum_{j=i}^{v-1}lambda_j. When v is not in A, with sorted roots
a_1<...<a_k, define the positive no-visit coefficient

E_A = (lambda_v/d)
      sum_{0<=j_1<...<j_k<v, a_t<=j_t for all t}
             product_t lambda_(j_t).

For k=0 the empty product and sum are one.

If v is not in A and k-r<=1, then

Delta_v=epsilon^(k+1) E_A+O(epsilon^(k+2)).

In every other case, set
b=min(A\({0,...,r-1} union {v})), or b=v if this set is empty.
Then b>r. Define

J_A = c d^(k-3) lambda_v sum_{j<b} h_x(j)lambda_j L_j,             r=0;

J_A = c^2 d^(k-r-3) lambda_v
      sum_{0<=i_0<...<i_r<b} h_x(i_r)
             product_t lambda_(i_t) L_(i_0),                    r>=1.

Then J_A>0 and

Delta_v=epsilon^(r+3)
        [J_A+1{v not in A, k-r=2} E_A]+O(epsilon^(r+4)).

In particular, for two original roots and leaky target:
- A={v,y}, y>0: positive cubic, with
  C=(c/d)lambda_v sum_{j<y}h_x(j)lambda_j L_j.
- A={v,0}: positive quartic, with
  C=(c^2/d^2)lambda_v sum_{i<j<v}h_x(j)lambda_i lambda_j L_i.
- A={a,b}, a<b<v: positive cubic, with
  C=E_A+(c/d)lambda_v sum_{j<a}h_x(j)lambda_j L_j.

One arbitrarily degenerate example: for n>=3, A=V\{n-2}, x=n-2,
v=n-1, and lambda_i=1, the first gain is

Delta_v=[(n-1)c^2/d^2] epsilon^(n+1)+O(epsilon^(n+2)).

More generally its leading polynomial is (c^2/d^2)p_V sum_{i<v}p_i.
The n=2 case belongs to r=0 and instead has leading c p_v p_0^2/d^2.

Thus, for any fixed c,n and positive rate vector lambda, sufficiently
small epsilon gives all binary add-one comparisons at the leaky target.
Finitely many comparisons and surplus IDLA then give the same threshold
for every finite all-active input at that target. General nonleaky
targets beyond case 1 remain outside this proof.

Proof ingredients and cancellation.

(a) Analyticity and trace IDLA are exactly RN (127u). On the face where
only S may sleep, ordinary target-capacity IDLA gives a mixture of binary
inputs on the trace cycle S, with fresh suffixes. Its cycle product is c.
The full active input on S is stationary for further additions.

(b) Two elementary trace-cycle coefficient vectors. Let pi_S be the full
active law. If j is the last site of S (the trace leaky vertex), then

[p_S](pi_S-mu_(S\{j}))
  = (delta_S-d delta_(S\{j})-c delta_empty)/d.                  (I)

Every lower-degree coefficient is zero. Here delta_B is point mass at
the exact occupied support B, and p_S=product_{i in S}p_i.

If v is the last site and j the penultimate site of S, then, restricted
to output supports containing v,

[p_S](pi_S-mu_(S\{j}))
  = delta_S/d-delta_(S\{j})-(c/d)delta_{ {v} }.                 (II)

To prove (I), RN (113),(113a) give
pi_S-mu_(S\{j})=Theta_S(S)/(d r_S), with r_i=(1-p_i)/p_i.
At zero sleep, for a mixed full input whose first active vertex is i,
only its initially sleeping prefix can remain: with probability
d^(|S|-i) every active packet leaves without a circulation; otherwise
the output is empty. In the full finite difference, all first-active
classes cancel except the last-active class and the all-sleeping input.
This gives delta_S-d delta_(S\{j})-c delta_empty.

For (II), mark C=S\{v}; (113) instead gives
pi_S-mu_(S\{j})=Theta_S(C)/(d r_C). Set all off-v sleep probabilities
to zero. If the first active nonleaky vertex has index i, the surviving
prefix plus v has probability p_v d^(|S|-i-1). Total v occupation is
p_v/(d+c p_v), independently of the initial flags. Taking the C finite
difference cancels all classes except all C asleep and only j active.
Its coefficient of p_v is delta_S-d delta_(S\{j})-c delta_{ {v} },
which gives (II). The statement concerns supports containing v only.

(c) Proof of case 1. The first r original particles deterministically
fill the first r sites of any trace set S of size r+1. Its last site j
remains unfilled with probability d^(k-r) when j<b; if j>=b it is filled
surely. The extra particle increases its filling probability by
h_x(j)d^(k-r). Every face of size at most r is already stationary.
Consequently the first possible coefficients are squarefree, on these
sets S. Substitution of (I), followed by taking the v marginal, gives C.

(d) For the leaky target, let w_v(B) be the forced-jump no-return
probability and u_B=w_v({v})-w_v(B). Put u_{ {v} }=0. If B contains v
and has first off-target resident i, the deterministic open queue gives

u_B=c sum_{l=i}^{v-1}p_l+O(||p||^2).                            (III)

A single successful sleep must occur before this first resident to stop
the forced packet: after that resident is woken at least two particles
must be absorbed. This proves (III). Also r_v S_vv=1/d+O(||p||).
The exact occupation identity is

Delta_v=S_vv(N_A^v(1)-N_(A+x)^v(1))
        +r_v S_vv sum_{B contains v, |B|>=2}u_B[mu_(A+x)(B)-mu_A(B)].

(e) Suppose v is an original root, or v is absent and k-r>=2. Every
trace face S containing v with |S|<=r+1 is already full and stationary.
For |S|=r+2, write its off-v sites as i_0<...<i_r=j. If j<b, original
IDLA fills S\{j}, and the probability that j remains vacant is
d^(k-r-1). The extra fills that remaining site with probability h_x(j).
If j>=b the face is already full. These assertions hold also when v is
not an original root: one of the at least two post-prefix roots fills v,
and a root at or before j would fill j deterministically.

Thus every leading coefficient of the exact-support differences with v
occupied is squarefree of degree r+2, and is (II) times
h_x(j)d^(k-r-1). Multiply by (III). For r>=1, S and S\{j} have the
same first resident i_0, so their coefficients combine as
(1/d)-1=c/d. This gives J_A with its factor c^2. For r=0 the lower
support is {v}, whose u is zero, giving the stated factor c instead.
The tail L_(i_0) runs through ALL original nonleaky sites, including
sites outside S; omitting those sites would give a wrong coefficient.

(f) If v is absent from A, avoiding v requires retaining all k particles
on the deterministic open line. Its leading no-visit probability is
epsilon^k times the parking sum in E_A without lambda_v/d. The upper
input requires at least k+1 sleeps (or has no no-visit branch if x=v).
Hence its contribution to Delta_v is epsilon^(k+1)E_A.
When k-r<=1, every face containing v of size at most k is already full;
the u-weighted remainder starts at degree k+2. Otherwise compare k+1
with r+3 in (e). This proves all cases in section 2.

**Agent-reported exact audits (not rerun at the stop).**
The maintained `verify_cycle_rare_sleep.py` differentiates the raw
112-state four-cycle absorption problem at zero sleep through degree 5.
For c=1/2, lambda=(1,1,1,1), all 128 binary/target comparisons have a
positive first nonzero coefficient; the maximum degree is 5. All 48
two-root comparisons have degree at most 4. The leaky-target formulas
and case-1 formulas match 32 and 82 entries respectively. A second fixed
identity audit c=2/3, lambda=(2,3,5,7) checks the same formulas exactly.
These two exact finite systems do not establish the remaining targets
on arbitrary cycle lengths, or uniformity over heterogeneous rates.


#### Fixed-root boundary extension: proof sketch awaiting audit

Erdos proposed keeping A={v,a} and its two sleep probabilities fixed in
(0,1), while taking every nonroot sleep probability small. On the face
where only A union {c} may sleep, trace onto these three vertices; the
claimed exact reduction is

Delta_v(A,x)=h_x(c)[q_v^trace(A union {c})-q_v^trace(A)].

Here x is the extra source and c is the sole nonroot allowed to sleep.
For a strongly connected kernel and h_c(v)<1, the proposed positive
p_c derivative would give a small nonroot cube for this fixed A by
coordinate integration from the exact zero-outside-A face. The threshold
would depend on the fixed root probabilities (uniform on compact subsets
of their open parameter square), not permit those roots to tend to zero
without further work.

After tracing and removing self-loops, index 0=v, 1=c, 2=a. In the
notation of (122), let H_3=[r_1^3]H, delta=det(I-K), and

D_*=1+d_2 r_0+d_0 r_2+delta r_0 r_2,
Lambda=1+r_0+r_2+d_1 r_0 r_2,
Z_trace=1+delta(r_0/d_0+r_2/d_2+r_0 r_2).

The proposed derivative is
r_0 d_1 H_3/[D_* d_0^2 d_2^2 delta Z_1 Lambda Z_trace].
The boundary coefficient claimed to ensure strictness is

[r_2^0]H_3=d_2^2 g_01 g_10 ell_1(1+r_0)(1+d_2 r_0),
g_01=K_01+K_02 K_21, g_10=K_10+K_12 K_20,
ell_1=k_1+K_12 k_2=d_0[1-h_1(0)].

The positive polynomial certificate supplies nonnegativity of the other
terms if the indexing and boundary identities are verified. Restoring a
trace self-loop multiplies the p_c derivative by 1/(1-K^trace_cc).
The standalone `verify_three_site_boundary.py` reconstructs four
coefficient identities from (122). The final agent report confirms only
checks 1-2 were collected in its fresh run; checks 3-4 were uncollected
when interrupted. The extension and strictness when h_c(v)=1 remain
**unproved here**. Do not confuse [r_1^3]H with (122)'s endpoint layers
in powers of r_2.

#### One pendant leaf: candidate extension, unaudited

Poincare proposed extending (127t) to a core satisfying that condition
with one bidirectional pendant leaf. The exceptional comparisons were
reported to reduce to traced full-minus-attachment and full-minus-leaf
comparisons. General cubic formulas and the heterogeneous-cube assembly
were not audited or saved as a complete proof.

One agent-reported raw calculation used c->j,ell,s with probabilities
1/4,1/4,1/2, j->c,s equally, and ell->c surely (29 transient states).
Writing D_(missing z;v)=q_v(V)-q_v(V without {z}), it gave

[p_c p_j p_ell]D_(missing c;j)=76/315,
[p_c p_j p_ell]D_(missing ell;c)=4/525,
[p_c p_ell^2]D_(missing ell;c)=96/175.

These are reported finite checks only; no independent artifact was
saved for them, and the proposed extension remains open.

#### Low-source insertion: the pointwise pair shortcut fails

The proposed shortcut was (B_a f_v)(C)>=0 for v outside C, |C|<=2.
Vacant-a rows are immediate. At C={a}, N1 gives
S_aa(B_a f_v)({a})=q_v(2a)-q_v(a)>=0. The remaining two-resident
case is false. Russell reported an exact calculation under both raw
priorities with v=0, b=1, a=2, C={1,2},

K=[[0,1/2,0,1/2], [1/2,0,0,1/2], [1,0,0,0]],
p=(1/2,1/2,1/64),
(B_2 f_0)({1,2})=-109/780.

Columns are 0,1,2,sink. No script was saved for this inline calculation;
main preserved its parameters and result without rerunning it at the
stop. This is an agent-reported operator-row obstruction, not a
counterexample to all-active monotonicity.

The narrower remainder question remains open: do sigma>=0,
sigma({v,a})=0, sigma L_v=0 and sigma L_a>=0 force
sigma(B_a-I)f_v>=0 on supports of size at most two? One exact LP on the
rare three-site kernel from section 22 had minimum zero, but all its
allowed insertion rows were already nonnegative, so the deletion
constraints played no role. It supplies no general certificate.

**Verification and continuation status.** Before the stop, the maintained
pair-capacity suite passed; the subsequently integrated three-cycle,
feeder and cheap-credit obstruction functions passed individual exact
runs. The small-sleep raw derivative, zero-outside trace, source-switch
and bulk-addition checks also passed. The two newly preserved standalone
verifiers above retain their distinct agent-report/uncollected statuses.
### A general positive renewal for the pair layer

The pair renewal underlying Section 25 extends to arbitrary graphs.
Let ell_C=mu_A(C), u_C=mu_(A+x)(C) for pairs C, and
tau_B=mu_(A+x)(B) for triples B. For each pair choose

i(C)=x if C=A; otherwise any element of A without C.

Define nonnegative matrices, with pair rows and pair/triple columns,

P_(C,D)=r_i 1{i in D}W_i(D,C),
Q_(C,B)=r_i 1{i in B}W_i(B,C), i=i(C).

Last-arrival deletion at i(C) gives exactly

(I-P)ell=kappa e_A, (I-P)u=Q tau,
kappa=N_A^x(A)>0.

Here kappa is the exact-support no-x-visit probability from active A,
equivalently its stabilization diagonal on the domain killed at x.
No reciprocal normalizer is asserted on a larger graph. For n=3 the
killed domain equals A, recovering kappa=1/Z_x. The upper source
vanishes because every selected i(C) is an upper root. The lower
source vanishes except at C=A because otherwise it is a lower root.
Singleton predecessors cannot end at a pair by particle conservation.

The inverse G=(I-P)^(-1) is nonnegative. To prove transience, let
M=binom(n,2) and average all two-root creation laws:

h_C=M^(-1) sum_(D: |D|=2)mu_D(C).

Applying the selected deletion rows gives h-Ph=b with
b_C=M^(-1)sum_D N_D^(i(C))(C)>=p_C/M>0. The term D=C has the
positive immediate-sleep history and i(C) outside C. Also 0<h_C<=1.
Thus Ph<=(1-epsilon)h for epsilon=min_C p_C/M>0, and
rho(P)<=1-epsilon<1. This does not bound the ordinary row sums of P.
Consequently ell-u=kappa Ge_A-GQ tau.

For a family E of pairs put phi_E=G^transpose 1_E and let Gamma(E)
be the triples containing at least one pair of E. The capacity cuts
are now equivalent to the exact potential inequalities

kappa phi_E(A)
 <=sum_B tau_B[(Q^transpose phi_E)_B+1{B in Gamma(E)}].

Indeed this is the signed cut ell(E)-u(E)<=tau(Gamma(E)). Applying it
to the positive-deficit subfamily of any D recovers (126); the converse
is immediate. Restriction to target stars preserves this equivalence.

This gives a positive, transient renewal, **not the missing capacity
proof**. Neither Q nor GQ has containment support: deleting a sleeper
from a triple can move another resident to a site outside that triple.
Nor do their columns have the required mass bound; r_i W_i can be
arbitrarily large as p_i decreases to zero. The remaining inequalities
must use the actual triple law or additional deletion constraints. The
maintained pair-capacity verifier checks both missing-root policies on
the two four-site kernels, including direct killed-domain source laws,
48 pair renewals and their strictly positive contraction defects.

## 27. Product backgrounds and directed-cycle queue reversal

**Proved:** on every directed cycle, the full all-active input dominates
every full-minus-one input for each joint-occupation event. On a cycle
with only one leaking edge, every full-wake field coefficient is
nonnegative. Neither assertion alone proves arbitrary-input N2 on a
cycle of more than three vertices.

Let K(i,i+1 mod n)=a_i, K(i,sink)=1-a_i, and product_i a_i<1.
Fix the missing site j. Order the remaining sites 1,...,m clockwise
after j, with m=n-1, relabelling their parameters accordingly. Initially
assume 0<a_i<1. Define

gamma_i=1-product_(k=i)^m a_k,
beta_i=p_i/(p_i+q_i gamma_i),
z_j=product_(i!=j)(p_i+q_i gamma_i),
mu_beta=product_i Ber(beta_i), nu_p=product_i Ber(p_i).

Here gamma_i is the probability that a jump from i reaches the sink
before j. Holding the particle at j, stabilize the other full active
sites while killing on a visit to j. Their no-visit terminal measure is

z_j mu_beta.                                            (128)

One way to check (128) is to process the cut path in order. At each
site, the initial active resident guarantees a positive incoming total.
The final resident is a fresh Ber(p_i) coin; every emitted particle
contributes the downstream sink-before-j probability gamma_i.
The weights telescope, giving independent unnormalized weights p_i for
an occupied site and q_i gamma_i for an empty one.

Let U_j(T) be the target-T probability from one active at j and an
independent mu_beta sleeping background on the other sites, retaining
all returns. First-contact cancellation at j gives

Theta_T({j})=z_j[beta_(T without j)-U_j(T)].              (129)

This is a finite-histories identity: histories which hit held j erase
its initial active/sleeping distinction; on the no-contact branch,
sleeping j gives target probability beta_(T without j).

### A queue lemma which retains every return

Consider the ordered suffix 1,...,m with incoming batch M>=0. Its
initial binary residents B have density f>=0 relative to mu_beta.
At a visited buffer i, draw Z_i~Ber(p_i), emit
L_i=M_i+B_i-Z_i, and pass Binomial(L_i,a_i) particles to the next
buffer. If M_i=0, leave Z_i=B_i and all later arrivals zero.
The final passed batch H returns to j. This is a legal ARW batch
schedule, using fresh coins and arrows.

For i<m put

a_i^*=a_i gamma_(i+1)/gamma_i in [0,1].

For h>=1 the following exact reversal holds:

P_f(Z=z,H=h)
 =mu_beta(z) gamma_1^M (a_m/gamma_m)^h
   E_[b~nu_p, tilted routes][
       f(b) 1{every buffer visited} binom(L_m,h)].        (130)

Inside this expectation start M_1=M, set L_i=M_i+b_i-z_i, and route
M_(i+1)~Binomial(L_i,a_i^*) until extinction. The binomial coefficient
is zero for L_m<h. To derive the formula, positive output requires
every buffer to be visited. The local change of measure is

mu_beta_i(b_i) nu_p_i(z_i) gamma_i^(M_i+b_i-z_i)
 =gamma_i^M_i mu_beta_i(z_i) nu_p_i(b_i).

Multiplying and tilting each internal binomial route telescopes the
gamma factors. The final route contributes binom(L_m,h)(a_m/gamma_m)^h.

For fixed b and common Bernoulli routing trials, increasing z decreases
every subsequent queue length, the full-visit indicator and the final
binomial coefficient. Thus, for every starting f, conditioning on any
positive return H=h leaves a **coordinatewise decreasing density
relative to mu_beta**. This statement is about the density, and is
stronger than stochastic domination.

For zero returns, draw Z~mu_beta, independent P_i~Ber(p_i), and construct
an auxiliary queue. Start Q_1=M. If Q_i>0, set

Btilde_i=P_i, R_i=Q_i+P_i-Z_i;
Q_(i+1)~Binomial(R_i,a_i^*) for i<m.

If Q_i=0, put Btilde_i=Z_i and keep the remaining queue zero. The last
queue is absorbed. The same change of measure gives, for every phi,

E_f[phi(Z);H=0]
 =gamma_1^M E_aux[phi(Z)f(Btilde)].                       (131)

At M=0 this includes Btilde=Z. In particular, f=1 gives an exact reset:
P(H=0)=gamma_1^M and Z conditional on H=0 has law mu_beta.

If f is decreasing, (131) yields the terminating target budget

P_f(Z_T=1,H=0)<=beta_T P_f(H=0).                         (132)

Here is the required coupling, rather than an appeal to ARW
attractiveness. Conditioning Z_T=1 increases the auxiliary service
probabilities from beta_i to one at T. Couple the resulting queue below
the original queue. When both queues are positive, use the same P_i
and ordered service coins, then monotonically coupled binomial routes;
their reconstructed Btilde bits agree. When the smaller queue is empty
and the other positive, couple its service/output coin above the other's
P_i, using beta_i>=p_i. When both are empty, order their service coins.
Queue order and coordinatewise order of the reconstructed vectors
persist. Thus Btilde conditional on Z_T=1 stochastically dominates
unconditioned Btilde. Decreasing f and (131) prove (132).

### Closing the circulation renewal

Initially f=1. Condition on the untimed local history at j and the
successive boundary batch counts M_k,H_k, leaving suffix residents,
internal outflows and routing details unobserved. Every positive-return
pass makes the conditional suffix
density decreasing by (130), irrespective of its previous density.
Processing j uses fresh records and does not inspect the suffix.
Consequently conditioning on the next outgoing batch preserves that
property; mixtures also preserve it. If j absorbs or loses the packet,
treat the next suffix input as M=0.

At the terminating zero-return pass, (132) bounds each suffix joint
occupation probability by beta_T. Termination is almost sure, with
finite mean: there are at most n particles, and a uniformly positive
probability of taking no sleep until a leaking edge and losing the
whole packet there. Averaging the terminal bound gives

U_j(T)<=beta_T for T avoiding j,
U_j(T)<=U_j(T without j)<=beta_(T without j) in general.

Together with (129) and (113a),

Theta_T({j})>=0,
mu_(V,T)>=mu_(V without j,T).                            (133)

Zero or unit a_i follow by continuity at sink-accessible kernels;
the finite transient ARW resolvent stays invertible. Self-loops can
be removed by the endpoint-preserving reduction. A one-vertex cycle
has occupation probability p/[p+(1-a)q] and satisfies the assertion
directly.

If the unique leaking edge is n-1 to 0, take the usual order
0,...,n-1. For C nonempty and j=max C, repeated use of the leak-free
case of (113) gives

Theta_T(C)=r_(C without j)Theta_T({j})>=0.                (134)

This proves the whole joint-target field on a single-leak cycle,
including arbitrary extra particles through (85a). The argument for
several leaking edges proves (133); it does not assert (134) there.

A separate AI audit checked the change of measure, zero-return coupling,
boundary-history conditioning and endpoint limits. The maintained
`verify_cycle_queue.py` checks (130)--(131) on complete point-mass bases
for three buffers with both single and varied leaks, and (132) on all
19 nonempty Boolean down-sets plus a positive nonproduct decreasing
density. Incoming batches range from zero through four. Independent
raw ARW absorption on heterogeneous four-cycles checks the no-contact
product laws, target budgets and singleton/full-minus-one identities;
for the single-leak case it also checks every field coefficient and
marked-set reduction. The proof above establishes the arbitrary-size
claim; these finite calculations audit its identities.



### What survives when the cut domain branches

The zero-return identity has a useful DAG extension. Let R be an
acyclic directed domain, with exits to a portal j or to the sink, and
assume h_i=P_i(sink before j)>0. Put

beta_i=p_i/(p_i+q_i h_i),
K^h_(i,k)=K_(i,k)h_k/h_i, K^h_(i,sink)=K_(i,sink)/h_i.

These transformed rows sum to one after the portal is removed.
Let M be a fixed vector of incoming active batches, and start the
sleeping residents with density f relative to the product mu_beta.
Process vertices in topological order, and let H count all exits to j.
There is an exact auxiliary representation

E_f[phi(Z);H=0]
 =product_i h_i^M_i E_aux[phi(Z)f(Btilde)].               (135)

To generate the auxiliary law, draw Z~mu_beta. At a vertex with
positive total incoming queue Q_i, take Btilde_i~Ber(p_i), independently
of Z_i and the past, and route Q_i+Btilde_i-Z_i particles using K^h.
At an inactive vertex set Btilde_i=Z_i. Sink exits disappear.
The Bernoulli and routing factors telescope exactly as in (131).

In particular, the product background resets conditional on no returns.
For decreasing f, the joint target budget (132) remains valid. The
auxiliary comparison uses the same local coupling of larger service
probabilities and reconstructed bits, now with monotonically coupled
multinomial routes along the DAG. Incoming queue vectors stay ordered
in topological order. This establishes the terminating-pass bound;
it does not by itself show that the cone survives positive returns.



Positive-return closure of the decreasing-density cone is false on a
branching DAG. Take two buffers a,b, with

a->(b,j,sink) probabilities (9/10,1/20,1/20),
b->(j,sink) probabilities (1/2,1/2),
p_a=p_b=99/100.

Both h-values are 1/2 and beta_a=beta_b=198/199. Send M=1 to a and
start with density f(B)=1{B_b=0}/(1-beta_b), a decreasing density
relative to mu_beta. There are at most two particles. For the event of
exactly one portal return, write
D(z)=P_f(Z=z,H=1)/mu_beta(z). Exact calculation gives

D(0,0)=2014477/500000, D(1,0)=2147409/200000,
D(1,0)-D(0,0)=6708091/1000000>0.

The same failure persists after mixing f with a small positive constant.
For a compact derivation, allow p_a=u, p_b=v. Conditional on L departures
from a,

E[1{Z_b=0}t^H | L]
 =v[(1+t)/20]^L+(1-v)[(1+t)/2]^L.

Taking the coefficient of t and averaging the two possible initial
a-residents gives

D(1,0)-D(0,0)
 =[(1+v)/(1-v)] [9v(u+5)-50]/100.

In a chain, a positive return forces every buffer to be visited. Here
the shortcut to the portal permits avoiding b, which defeats that
pointwise density argument. This example alone does not show that the
chosen starting density is generated by the original renewal.

There is also a reachable failure. Start the same suffix from mu_beta,
with one incoming particle, and condition on H_1=3. This has probability
970299/79202000000 and leaves the suffix empty. Let the portal send
the three returning particles back to the entrance, and condition on
H_2=1. The second-pass joint masses at (00,01,10,11) are
(3297,395604,788040,38811960)/80000000. The resulting conditional density
relative to mu_beta increases from 00 to 10 by
27047483/39998901>0. A portal with p_j=1/2 and a deterministic jump to
the entrance realizes the requisite emptying decisions with positive
probability. Thus conditioning on reachable return histories does not
restore the density invariant. Section 28 replaces it with deletion
inequalities which do survive these histories.

Weakening the hypothesis to ordinary stochastic domination does not
repair the terminating step. For the same parameters and M=1, the law

nu(0,0)=1/39601, nu(0,1)=198/199,
nu(1,0)=198/39601, nu(1,1)=0

is stochastically dominated by mu_beta, as all upper-set inequalities
verify. Nevertheless,

P_nu(Z_b=1,H=0)-beta_b P_nu(H=0)
 =368899839/78805990000>0.

Thus even a single-target zero-return budget can fail under this weaker
hypothesis. The maintained cycle-queue verifier checks both obstructions
against raw ARW absorption, and checks (135) on all input/output point
bases and its budget on all five nonempty two-buffer down-sets. These are
failures of proposed invariant classes, not all-active monotonicity
counterexamples. The unresolved extension must retain more about how
the background and the returning packet were created.




### Partial root sets and the unresolved reference comparison

Cut a directed cycle at j and order the remaining path 1,...,m.
For a finite all-active path input eta let nu_eta(B) be its terminal
measure restricted to no visits to j. With h_i=1-product_(k=i)^m a_k,
ahat_i=a_i h_(i+1)/h_i, beta_i=p_i/(p_i+q_i h_i), and
t_i=(p_i+q_i h_i)/h_i, there are two positive representations:

nu_eta(B)=(h^eta/h_B) P_eta^(ahat,p)(Z=B)
 =h^eta E_eta^(ahat,beta)[1{Z=B} product_(i visited)t_i].

The second formula uses the forward path order, so each visited site's
last coin contributes once. In particular sum_B h_B nu_eta(B)=h^eta.
With at most two roots the queue takes only values 0,1,2, and the
visited set is a union of at most two intervals.

Fix the future active inputs and let F_i(k) be their no-j-visit
probability with k additional incoming particles at i. The positive
representation shows F_i(k+1)>=h_i F_i(k): after dividing by h_i^k,
common beta-coins and monotone binomial routes make the visited set
increase with k, while every t_l>=1.

For Psi_i(s)=E[F_(i+1)(Binomial(s,a_i))] this gives
Psi_i(s+1)>=h_i Psi_i(s). Hence, conditional on any incoming total n>0
and no portal exit, the site's occupation probability is

p_i Psi_i(n-1)/[p_i Psi_i(n-1)+q_i Psi_i(n)]<=beta_i.

Unvisited sites are empty. Averaging over incoming counts proves the
same bound conditional on any clockwise output prefix. Thus the
normalized partial creation law is stochastically dominated by
mu_beta. If K_M^0 is the zero-return suffix kernel, Abelianity gives

nu_eta K_M^0=nu_(eta+M delta_1).

Consequently the terminating product budget extends to these partial
creation backgrounds. Section 28 proves it for the larger no-portal
deletion cone.

This class need not have decreasing density. On the four-cycle with
all a_i=p_i=1/2, cut at j=0 and use roots {1,3}. Then
h=(7/8,3/4,1/2), beta=(8/15,4/7,2/3), and

nu(1)=93/128, nu(empty)=11/128, nu({1})=1/8.

The normalized density relative to mu_beta is 55/31 at empty but 70/31
at {1}. Moreover its actual occupation probability at 1 is 16/31,
strictly below beta_1=8/15. Bounding the continuation by the latter does
not compare it with the former.

The stronger pointwise comparison against the actual baseline is
false even for a target at an original root. Take
K_01=K_12=K_23=1, K_30=2/3, K_3s=1/3,
p=(1/3,2/5,3/7,4/9), roots A={0,2}, extra site j=1 and target 2.
The mixed-minus-active coefficient is

q_2(1,-1,1,0)-q_2(1,1,1,0)
 =-15889009571721264/1490750968359917015<0.

Thus the product budget cannot simply be replaced by the partial-root
marginal, even on a single-leak cycle. These are not add-one occupation
counterexamples.

The exact source identity clarifies what remains. For fixed all-active
eta set D_x=q_v(eta+delta_x)-q_v(eta), D_sink=0. Let theta_x be the
mixed-minus-active coefficient for an extra initially sleeping particle
at x when x is outside the input support, and set theta_x=0 on that
support. The addition algebra gives

(I-K)D=diag(1/r_x)theta, D=(I-K)^(-1)diag(1/r_x)theta.

The negative example refutes pointwise nonnegativity of this forcing.
A successful use of the partial-root representation must combine the
source terms with their Green weights. The maintained cycle-queue
verifier checks the positive transfer formula, three zero-pass
compositions, prefix bounds, both exact obstructions and this source
identity.


## 28. Exit counts and backward creation cones

The failed density invariant has a useful replacement. The results in
this section concern an open domain R with two absorbing boundaries:
the original sink and a portal j. Let h_i be the ordinary probability
of reaching the original sink before j, and initially assume h_i>0.
Put h_B=product_(i in B)h_i and beta_i=p_i/(p_i+q_i h_i).

Write W_i^k(B,C) for the forced-deletion kernel which starts from
sleeping B containing i, removes that sleeper, forces its first jump,
stabilizes without revisiting i, and records exactly k portal exits.
Define L_i^0 by the identity row on supports avoiding i and
-r_i W_i^0 on supports containing i. The no-portal creation cone is

C^0={mu>=0: mu L_i^0>=0 for every i}.                    (136)

### Exact exit-count slices and the Doob transformation

For any all-active input eta, let mu_H(C) record terminal support C
and exactly H portal exits, and let N_H^i be its no-i-visit submeasure.
The last-arrival bijection preserves portal counts additively, giving

mu_H(C)=N_H^i(C)
 +r_i sum_(B contains i) sum_(k=0)^H mu_(H-k)(B)W_i^k(B,C),
 C avoiding i.                                         (137)

Consequently every exact total-exit slice mu_H lies in C^0. This is
valid on arbitrary finite sink-accessible domains. It does not assert
the same thing after imposing additional intermediate exit-count
conditions.

Let Khat_(i,l)=K_(i,l)h_l/h_i and
Khat_(i,sink)=K_(i,sink)/h_i, removing portal jumps. This is a stochastic
kernel on R plus the original sink, with the same sleep probabilities.
The pathwise change of measure gives

W_i^0(B,C)=(h_B/h_C) What_i(B,C).

Thus mu belongs to C^0 exactly when lambda(B)=h_B mu(B) belongs to
the ordinary creation cone for Khat. In particular, adding an active
particle and retaining only zero-portal-exit histories preserves C^0:
if A_y^0 is that subprobability addition kernel, then

(mu A_y^0)(C)h_C=h_y(lambda Ahat_y)(C),

and Section 7 gives ordinary-cone invariance.

There is also an arbitrary-graph survival bound for genuine all-active
inputs. Let F(eta) be their probability of avoiding the portal. Then

F(eta+delta_y)>=h_y F(eta),
F(eta+zeta)>=h^zeta F(eta).

Freeze eta and force the marked extra particle along an ordinary walk
until sink or portal: after each lone sleep, reactivate it as part of
the next acceptable toppling. At frozen roots, jumps leave exactly the
original active population. On reaching the original sink, the remaining
configuration is precisely eta, with independent fresh stack tails.
The probability that this walk and the subsequent stabilization both
avoid the portal is h_y F(eta). Their concatenation is an acceptable
stabilizing sequence. Least action bounds the legal odometer by its
consumed prefixes, so legal stabilization cannot introduce a portal jump
absent from those prefixes. This proves the inequality; iteration gives
the second one. Multiplicities and self-loops cause no change. Initially
sleeping residents cannot be frozen this way, so the all-active
hypothesis is essential. This survival ratio does not itself give an
occupation budget.

### Product domination for the whole cone on a DAG

**Theorem.** On an acyclic domain, every normalized mu in C^0 is
stochastically dominated by the product law mu_beta. More precisely,
each occupation probability conditional on an earlier topological
output prefix is at most beta_i. Also, for every y,

mu A_y^0 1>=h_y mu 1.                                   (138)

Here is an induction proving both assertions. Let i be the first
vertex in a topological order, and split mu into measures mu_0,mu_1
on the remaining domain according to its occupation bit. Deletions
at later sites cannot alter this bit, so both measures belong to the
remaining domain's cone. The deletion constraint at i is exactly

mu_0>=r_i mu_1 J,

where J=K_(i,sink)I+sum_(l in R without i)K_(i,l)A_l^0
is the no-portal forced-jump kernel on the remaining domain.

By induction, nu J1>=h_i nu1 for every nu in its cone, and J preserves
that cone. Summing the displayed constraint gives
mu_0 1>=r_i h_i mu_1 1, hence the first-bit bound beta_i.
Restricting successively to any earlier output prefix gives all the
claimed conditional bounds. Coupling each conditional Bernoulli with
an independent Ber(beta_i) coin proves full stochastic domination.

For (138), sources later than i are handled by induction separately
on mu_0 and mu_1. At i, the one-particle addition has no returns from
the later domain. Its no-portal survival mass is

p_i(mu_0+mu_1 J)1+q_i(mu_0 J+mu_1 J^2)1
 >=(p_i+q_i h_i)(mu_0+h_i mu_1)1
 >=h_i(mu_0+mu_1)1.

The last inequality is
(1-h_i)[p_i mu_0 1-q_i h_i mu_1 1]>=0, already supplied by the
deletion constraint. This completes the induction. The empty domain
is immediate. Repeated zero-exit additions give the analogous
product h^eta bound for any fixed active input.

This result is stronger than a bound just for genuine partial-root
laws. It does not follow from ordinary stochastic domination alone,
as Section 27's zero-return obstruction demonstrates.

### Positive return counts along a path with portal shortcuts

**Theorem.** Order R={1,...,m}. From i<m allow jumps only to i+1,
the portal j, or the sink; from m allow only portal or sink jumps.
For every incoming batch M at site 1 and every exact portal return
count H, the corresponding coefficient kernel preserves C^0.

Prove this by induction on the path length. On the suffix after its
first site i, let A_(i+1)(z) be active addition with each portal exit
weighted by z, and put

J_i(z)=K_(i,sink)I+K_(i,j)z I+K_(i,i+1)A_(i+1)(z).

At the last site omit the final term. All coefficients are nonnegative.
By induction each coefficient of J_i(z), and of every nonnegative
power, preserves the remaining cone.

Split an initial cone measure into mu_0,mu_1 at i. For M>=1, the
output measures at its two occupation bits are exactly

nu_1(z)=p_i[mu_0 J_i(z)^(M-1)+mu_1 J_i(z)^M],
nu_0(z)=q_i[mu_0 J_i(z)^M+mu_1 J_i(z)^(M+1)].

Both coefficient families satisfy all downstream deletion constraints.
The remaining constraint is positive coefficientwise because

nu_0(z)-r_i nu_1(z)J_i(0)
 =q_i[mu_0 J_i(z)^(M-1)+mu_1 J_i(z)^M]
       [J_i(z)-J_i(0)]>=0.                              (139)

No commutation between J_i(z) and J_i(0) is used; the displayed factor
order is essential. The M=0 kernel is the identity. This proves the
claim for every batch and return count.

### Disjoint branches and outward forests

The entrance result extends to an outward directed forest: internal
edges run from parents to children, and each vertex may also jump to
the portal or sink. For every root, every incoming batch, and every
exact portal-exit count, the addition coefficient preserves C^0.

First note a tensor property that does not require independence of the
input measure. For a disjoint union of domains, fixing all other
component configurations leaves a cone measure on the selected
component. A positive kernel on that component which preserves its
cone therefore preserves its own constraints on the whole union.
Deletion constraints in other components commute with the kernel,
so their nonnegativity is preserved too.

Induct from the leaves. After removing a root i, its child subtrees
are disjoint. By the tensor property, each coefficient of the
child-entry kernel A_c(z) preserves the remaining forest's cone.
Replace the path's jump operator by

J_i(z)=K_(i,sink)I+K_(i,j)z I+sum_(c child of i)K_(i,c)A_c(z).

Its coefficients and all its powers preserve that cone. The same
two output formulas and the ordered factorization (139) give the
root deletion constraint. This completes the induction. Mixtures
of root entries, and any sequence of such entries with specified
total or successive exit counts, preserve C^0 as well.

### Closing feedback through forest roots

Now adjoin j to an outward forest R. Its jumps may go to any forest
root or the original sink; forest vertices may return directly to j.
Assume the sink is accessible from every vertex. For every increasing
terminal payoff f on occupied sets,

mu_V(f)>=mu_(V without j)(f).                           (140)

Thus this full-minus-one comparison has **occupied-set stochastic
domination**, stronger than just the specified single-target order.

Indeed, the full-root no-j-visit law on R is z_j mu_beta, by the DAG
product formula. It belongs to C^0 by (137). Start an active particle
at j against mu_beta. Condition only on its untimed local history and
successive boundary batch counts. By (139) and the tensor argument,
every positive-return conditional forest law stays in C^0. Fresh processing at j and
conditioning on its next outgoing batch preserve this property.
At termination, (138)'s domination theorem bounds the suffix law
stochastically by mu_beta; the same applies if the next batch is zero.
Absorption has finite mean by sink accessibility and the bounded
particle count.

Let U_j(f) be the final expectation in this continuation. Since f is
increasing,

U_j(f)<=E_mu_beta[f(B union {j})].

First-contact cancellation and the full-minus-one addition identity give

z_j[E_mu_beta f(B union {j})-U_j(f)]
 =r_j sigma_j[mu_V(f)-mu_(V without j)(f)],

where sigma_j is the strictly positive ordinary sink-before-first-return
probability at j. This proves (140). The identity holds for arbitrary
payoffs by linearity of the joint-occupation basis.

Zero h_i are handled by adding vanishing direct-sink leakage within the
same graph class and taking the finite-resolvent limit. Self-loops are
removed by the established endpoint-preserving reduction. Every directed
cycle belongs to this class for every choice of j, so (140) strengthens
Section 27's cycle comparison to full occupied-set order. Both disjoint
branches and shortcut returns to j are allowed. This does not assert
higher marked-field positivity, or arbitrary-input N2 on a larger graph.

### Exact limits of return-count closure

Entry at a later site can fail even on a two-vertex path. Use Section
27's two-buffer example: 0 jumps to 1, portal, sink with probabilities
9/10,1/20,1/20; 1 jumps to portal or sink equally; p_0=p_1=99/100.
Then h_i=1/2 and beta_i=198/199. Start from mu_beta and add one active
particle at site 1. The event of two portal exits has probability
99/39800 and conditional law

nu(empty)=1/199, nu({0})=198/199.

Its deletion constraint is
(nu L_0^0)({1})=-891/99500. Thus arbitrary entry locations are excluded.

Merging branches also fails, even with entry at the unique source and
the correct full-root no-portal background. Take R={0,1,2,3}, with
0 jumping to 1 or 2 equally, 1 and 2 jumping surely to 3, and 3 jumping
to the portal or sink equally. Put p_i=1/2, so h_i=1/2 and beta_i=2/3.
Let nu_H be the unnormalized output law from mu_beta plus one active
particle at 0, with exactly H portal exits. Then

nu_3({0,3})=0, nu_3({0,1})=nu_3({0,2})=1/648,
(nu_3 L_1^0)({0,3})=-1/1296.

Here is a direct count. Three exits and two terminal sleepers require
all four initial sleepers. If 0 remains occupied, only one particle
leaves it, so only one branch is visited and the other keeps its sleeper.
Consequently {0,3} is impossible. The event ending at {0,1} has initial
mass (2/3)^4, then root sleep, routing to 2, site-2 emptying, and site-3
emptying, each with probability 1/2; all three departures from 3 must
use the portal, with probability 1/8. Their product is 1/648. Finally
W_1^0({0,1},{0,3})=1/2 and r_1=1. The total H=3 mass is 175/5184,
so the normalized residual is -4/175.

For this same constraint the unnormalized residuals at H=0,...,5 are
(0,0,-1/432,-1/1296,0,0). Summing all return counts still gives -1/324.
Thus forgetting the return count does not repair this cone invariant.

This is a counterexample to the proposed return-conditioned invariant,
not to occupied-set order or the main all-active conjecture. The tensor
argument uses disjoint subtrees and cannot cover this merge.

### Static domination on a two-site cyclic domain

The static theorem has a cyclic extension, separate from positive-count
closure. On a two-site domain with h_1,h_2>0, let pi_0 be the normalized
no-portal full-root law. Every normalized mu in C^0 is stochastically
dominated by pi_0.

Remove self-loops and write K_12=a, K_21=b, K_1s=u, K_2s=v, with the
remaining row masses going to the portal. Set q_i=1-p_i, r_i=q_i/p_i.
For mu=(e,x,y,z), indexed by empty,{1},{2},{1,2}, direct forced-deletion
enumeration gives exactly

e>=A x+B z, e>=C y+D z,
y>=alpha x+beta z, x>=gamma y+delta z,

where

A=r_1(u+a q_2 v), B=r_1 a q_2 v^2,
C=r_2(v+b q_1 u), D=r_2 b q_1 u^2,
alpha=r_1 a p_2, beta=r_1(u+a p_2 v),
gamma=r_2 b p_1, delta=r_2(v+b p_1 u).

These beta,delta are local scalar coefficients, unrelated to the DAG
product parameters above. Both A and C are positive because h_i>0.
Write pi_0=(E,X,Y,1)/Z, Z=E+X+Y+1. Its four deletion constraints are
equalities, so

E=A X+B=C Y+D, Y=alpha X+beta, X=gamma Y+delta.

The cone is generated by the following five rays:

(E,X,Y,1), (1,0,0,0), (A,1,alpha,0),
(C,gamma,1,0), (AC,C,A,0).

Indeed, subtract z(E,X,Y,1). Since alpha gamma=ab q_1q_2<1, the
remaining cross inequalities force x,y>=0 and then e>=0. In the z=0
section, set e=1. The feasible polygon is
0<=x<=1/A, 0<=y<=1/C, alpha x<=y, gamma y<=x.
Its vertices give the other four rays; the corner ordering follows from

A-C alpha=r_1u(1-ab q_1q_2)>=0,
C-A gamma=r_2v(1-ab q_1q_2)>=0.

It suffices to check normalized rays against the three nontrivial
upsets other than full occupation. Put L=A beta-alpha B=r_1^2u(u+av).
The table gives the numerator of pi_0(U)-nu(U), over Z times the
unnormalized ray's total mass:

| Ray | 1 occupied | 2 occupied | Nonempty |
|---|---:|---:|---:|
| (A,1,alpha,0) | r_1a(1-v)(p_2+q_2v) | 1+A+beta+L | A-B+L |
| (AC,C,A,0) | (C+1)(A-B)+D | (A+1)(C-D)+B | AC-BC-AD |

The other singleton ray follows by exchanging sites; the empty ray is
immediate. Every non-full ray has zero mass at {1,2}. All table entries
are nonnegative: A>=B, C>=D, L>=0, and

(AC-BC-AD)/(r_1r_2)
 =uv[1+ab q_1q_2(1-u-v)]
  +a q_2v^2(1-v)+b q_1u^2(1-u)>=0.

The bracket is at least 1-ab q_1q_2>0. This proves domination for every
cone measure by convexity, including zero direct-sink probabilities.
The self-loop reduction preserves h and endpoint laws; it has
W_i^0=(1-K_ii)(W_i^0)' and r_i'=r_i(1-K_ii), so it preserves C^0 too.

### Static domination assembles across strongly connected components

The preceding static result extends to every domain whose strongly
connected components have at most two vertices. Let pi_R^0 be its
normalized no-portal full-root law. Then, for every normalized mu in C^0,

mu <=_st pi_R^0,
mu A_zeta^0 1>=h^zeta.                                 (141)

Here A_zeta^0 adds a fixed active batch zeta and keeps only no-portal
histories. This theorem allows arbitrary edges between the components;
no disjoint-branch hypothesis or positive-count closure is asserted.

First assume every h_i>0. Write pihat_R for the full-input stationary
law under the Doob kernel, and Z for the normalization of h_B^(-1)pihat_R.
The following eigenmeasure identity holds on every finite domain,
without any component-size assumption:

pi_R^0(B)=Z^(-1)h_B^(-1)pihat_R(B),
pi_R^0 A_zeta^0=h^zeta pi_R^0.

Indeed, multiply the left-hand side at C by h_C and use the Doob
change of measure, followed by pihat_R Ahat_zeta=pihat_R.

The reference law factorizes across the components. Process them in
topological order under the Doob kernel. Given the previous components'
histories and emitted batches, the next component starts with full
active support plus a fixed extra batch. Full-input stationarity makes
its terminal support law independent of that upstream history. Thus
pihat_R is the product of the component stationary laws, and the h_B
tilt shows pi_R^0 is the product of corresponding laws pi_C^0. Do not
condition additionally on the current component's own emitted counts.

The effective two-boundary kernel for a component C retains internal
K entries and has boundary weights

K^C_(i,sink)=K_(i,sink)+sum_(l outside C)K_(i,l)h_l,
K^C_(i,j)=K_(i,j)+sum_(l outside C)K_(i,l)(1-h_l).

Its avoidance potential is h restricted to C. Its Doob transform is
the ordinary component kernel with outgoing transitions killed, so
its normalized no-portal full-input law is exactly pi_C^0.

Induct on the number of components. Take the first component C and
put D=R without C. For a cone measure mu, fixing its configuration
B in C leaves a cone submeasure on D: downstream deletion cannot
alter C. By induction every normalized such measure nu_B is dominated
by pi_D^0. The survival payoff

F_zeta(E)=P(no portal exit from sleeping E plus active zeta)

is decreasing in E. Under common stacks, a legal sequence for the
smaller configuration remains legal for the larger one; consumed portal
jumps remain consumed. Hence the reference eigenmeasure identity gives
nu_B A_zeta^0 1>=h^zeta.

For a deletion at i in C, first stabilize inside C, avoiding a revisit
to i and any direct portal exit, and hold the emitted downstream batch
zeta. Conditional on the starting B, this internal history is independent
of the downstream background. The preceding survival bound lower-bounds
its downstream continuation by h^zeta. Therefore the downstream-averaged
deletion kernel is at least the effective component deletion kernel.
Summing the global i-deletion inequality over downstream outputs shows
that the C marginal belongs to the effective component's C^0. The
one- or two-site static theorem bounds it by pi_C^0.

Conditional downstream laws satisfy the induction hypothesis. Successive
conditional Strassen couplings, with the same fixed reference marginal
at each component, give domination by the product pi_R^0. Applying the
decreasing payoff F_zeta proves the survival inequality in (141), so the
induction is noncircular.

Zero h causes no fixed-measure perturbation problem. Within a component,
either all h_i are positive or all are zero. In the zero case the effective
sink weights vanish. Its no-portal full-input law is concentrated on full
occupation, by particle conservation, and domination is automatic. The
reference factorization and eigenmeasure identity extend by continuity;
the same conditional-component induction then applies. Self-loop removal
again preserves the cone.

### The remaining three-site static budget

There is a small exact reduction, still not a proof of three-site static
domination. Assume h_i>0 and write pi=pi_0. Section 8's full-support
minorization, after the Doob tilt, gives
mu>=[mu(V)/pi(V)]pi. Subtracting this component leaves a nonnegative
measure in C^0 supported on at most two sleepers: unlike a proper-root
peel, all deletion residuals of the subtracted full law vanish.

For such an unnormalized measure on three vertices, write
b_i=mu(R_i), s_i=mu({i}), e=mu(empty), R_i=V without i. Its pair
constraints are b>=Qb, where

Q_ij=r_i W_i^0(R_j,R_i), i!=j, Q_ii=0.

Put G=(I-Q)^(-1)>=0 and t=(I-Q)b>=0, so b=Gt. Positivity and
transience follow from the full reference law's strict pair forcing,
as in Section 25. The remaining constraints are

s_j>=a_(j,i)s_i+B_(j,i)b, i!=j,
e>=c_i s_i+d_i b,

where a_(j,i)=r_iW_i^0({i},{j}),
(B_(j,i))_k=r_iW_i^0(R_k,{j}),
c_i=r_iW_i^0({i},empty), and
(d_i)_k=r_iW_i^0(R_k,empty). Terms whose starting support avoids i
are zero. The least singleton vector satisfies

s_j^*=max_(i!=j)[a_(j,i)s_i^*+B_(j,i)Gt],
e^*=max_i[c_i s_i^*+d_iGt].

There are eight policies sigma, choosing one predecessor sigma(j)!=j
in each singleton row. Let A^sigma,B^sigma contain those selected rows.
Each policy has a positive inverse I-A^sigma; scaling by pi({i}) gives
a common strict contraction. Monotone value iteration and a greedy
policy at its fixed point show that one policy attains s^* in every
coordinate. Thus the minimum required total mass is

M(t)=1^transpose Gt+1^transpose s^*+e^*
    =max_(sigma,a) ell_(sigma,a)t,
ell_(sigma,a)=(1^transpose+d_a)G
 +(1^transpose+c_a e_a^transpose)(I-A^sigma)^(-1)B^sigma G.

There are only 24 row vectors ell. For a nonempty pair family J, put
q_J=pi(V)+sum_(i in J)pi(R_i). Domination for the upset generated by
these pairs is equivalent to

1_J^transpose Gt<=q_J M(t) for every t>=0.

Finite linear-programming duality makes this equivalent to existence
of probabilities omega_(sigma,a) such that

sum_(sigma,a)omega_(sigma,a)ell_(sigma,a)
 >=(1_J^transpose G)/q_J coordinatewise.

This is a three-coordinate convex-hull budget. It addresses the seven
pair-generated upsets, not all three-site upsets. Even when h=1,
Section 25 supplies only its coordinate-axis tests: M(e_x)=Z_x and
(1_J^transpose G)_x<=q_J Z_x. A single mixture of policies must meet
all three coordinate requirements; that step remains open.

The normalization issue is real. On the symmetric three-site kernel
K_ij=1/4 for i!=j, K_is=1/2, with p_i=1/2 and no portal edges, the law

mu(empty)=1161/11571,
mu({i})=1518/11571 for each i,
mu(R_i)=1952/11571 for each i, mu(V)=0

belongs to C_2. All empty and singleton deletion residuals vanish;
every pair residual is 16/133. Matching its pair vector by two-root
creation laws uniquely forces coefficient 61/133 on each law, totaling
183/133>1. Thus those laws plus a nonnegative one-sleeper remainder
cannot give a convex decomposition. Equivalently,
M(1,1,1)=133/16<sum_x M(e_x)=183/16. The sources share lower-layer
mass requirements. This obstructs that mixture proof, not the proposed
static domination. `source_sharing_audit` in the cone verifier checks
all residuals and the forced mixture coefficients exactly.

### Weighted intertwining and a dominator criterion

There is an exact algebraic reason that some positive-count kernels
preserve the cone. Let A_y(z) weight every portal exit by z, and write
W_i(z)=sum_k z^k W_i^k and L_i(z)=I_i-r_iW_i(z). Here I_i has identity
rows on supports avoiding i and zero rows on supports containing i.
Let B_i(z) be addition at y killed on any visit to i, on supports
avoiding i, with the same exit weights. Then

A_y(z)L_i(z)=L_i(z)B_i(z).                              (142)

Apply the count-preserving last-arrival identity to each binary
all-active input, before and after adding y. Sequential stabilization
preserves both the no-i-visit event and total portal counts, giving
P(z) times the two sides of (142). The polynomial matrix P(z) is
invertible over the rational-function field because P(1), with portal
and sink merged, is invertible. Cancel P(z). No positive h assumption
or constant-term invertibility is needed.

Define D_i(B,C)=1{i in B, C=B without i}. Its rows avoiding i are zero.
On the domain killed at i, put

J_i(z)=K_(i,sink)I+K_(i,j)zI
       +sum_(x!=i)K_(i,x)A_x^[i](z),

where A_x^[i] is killed addition. A forced self-loop at i is an immediate
forbidden revisit and contributes zero. Then W_i(z)=D_iJ_i(z).
All weighted killed additions commute at a common z, so J_i(z) commutes
with B_i(z). Write J_i^+(z)=J_i(z)-J_i(0) and

T_i(z)=A_y(z)D_i-D_iB_i(z)>=0 coefficientwise.

The inequality is a literal inclusion of histories: the second term
requires an initially occupied i never to be visited; these histories
are among those in the first term, which requires only terminal
occupation at i. The portal counts agree.

The general decomposition is

A_y(z)L_i^0=L_i^0B_i(z)+r_i T_i(z)J_i^+(z)
             +r_iD_i[J_i(0),B_i(z)].                  (143)

Thus the final commutator is the only potentially signed term.

**Dominator lemma.** Suppose every directed path from y to the portal
passes through i. Then every coefficient A_(y,H) preserves the single
i-deletion inequality. Indeed, B_i(z) is independent of z: without
visiting i, the active population cannot reach a portal exit. Awakening
sleepers causes no exception, since each newly active site is reachable
from y without i. Commutation at z and at 0 separately gives
[J_i(0),B_i]=0. Multiplying (143) by mu>=0 with mu L_i^0>=0 proves
the claim. When y=i, B_i=0 and the output inequality is automatic.

If every vertex dominates the portal from y, each coefficient preserves
the whole cone C^0. Products and coefficient convolution give every
incoming batch and exit count. In particular, order the vertices
1,...,m, allow jumps i->k only when k<=i+1, and allow portal exits only
from m. With entry at 1, arbitrary backward jumps, sink leaks, and
self-loops are permitted. Every vertex is a dominator, so countwise
closure holds even when the domain is cyclic. This includes a directed
cycle with its portal exits confined to the predecessor of the entry.

For two vertices, allowing portal exits also from the entry can destroy
closure. In columns 1,2,portal,sink, take

K_1=(0,1/2,1/3,1/6), K_2=(1/3,0,1/6,1/2), p=(1/3,2/5).

Here h=(1/2,2/3) and pi_0=(31,39,38,42)/150 in support order
empty,{1},{2},{1,2}. After one active entry at 1, the event of two
portal exits has probability 428/4725 and conditional law
(260,21,147,0)/428. Its site-2 deletion residual at {1} is

21/428-(1/6)(147/428)=-7/856,

or -1/1350 before normalization. The static domination theorem remains
valid on this domain. The failure concerns the return-conditioned law.

### Complete feedback as the remaining cone question

The feedback proof has a general sufficient form. For any open domain
R, let pi_R^0 be its normalized no-portal full-root law. Adjoin portal
j and start it with one active particle, independently of a background
distributed as pi_R^0. Let F_j be the kernel on R obtained after the
entire enlarged graph has stabilized, forgetting j's terminal bit.
If

pi_R^0 F_j <=_st pi_R^0,                                (144)

then the full-minus-j occupied-set comparison follows by the same
first-contact identity used in (140). The two sides of that identity
use z_j pi_R^0 for the no-j-visit law; no product structure is needed.
For increasing f, dominate the terminal R configuration and replace
the terminal j bit by 1.

Static domination plus pi_R^0 F_j in C^0 would suffice for (144).
The stronger proposed statement is

mu in C^0 implies mu F_j in C^0.                        (145)

**General complete-feedback preservation (145) is FALSE**, even on a
two-vertex DAG if entry is allowed at a later vertex. The remaining DAG
question restricts portal departures to sources of R. In that class
(145) would combine with static domination to prove full-minus-portal
occupied-set order with arbitrary merges. General static domination and
the occupation budget (144) remain open; failure of this sufficient cone
invariant proves neither false.

There is bounded exact evidence for (145). Close the two-portal cyclic
pair above by giving j sleep probability 1/2 and a deterministic jump
to its entry. Also close the merging diamond by the same portal rule,
with all five sleep probabilities 1/2. Both entire feedback kernels
preserve their full cones: 4 and 32 rational primal/dual LP certificates,
respectively. In particular, the diamond's intermediate residual
-1/324 is not a terminal obstruction. Thus complete feedback genuinely
repairs some intermediate cone failures, but it does not do so universally.

Here is the simplest terminal obstruction. Let R={a,b}, with a jumping
surely to b and b jumping to the portal or sink equally. Put all sleep
probabilities equal to 1/2 and let j jump surely to b, the later vertex.
Then pi_R^0 is product Ber(2/3). Complete feedback leaves a untouched
and changes b to Ber(22/35), independently. In bit order 00,01,10,11,
the terminal law is (13,22,26,44)/105. Since
W_a^0({a},{b})=1/2 and W_a^0({a,b},{b})=1/4,

(pi_R^0 F_j L_a^0)({b})=22/105-26/210-44/420=-2/105.

Thus (145) fails even for one-site components. The terminal law is
nevertheless stochastically dominated by pi_R^0. This example explains
why the initial entrance restriction cannot be dropped merely by
waiting for all feedback to finish.

For the counterexample use R={0,1,2}, with rows in columns
0,1,2,portal,sink,

K_0=(0,0,1/2,1/4,1/4),
K_1=(1/2,0,0,1/4,1/4),
K_2=(15/16,3/64,0,1/128,1/128).

Give all four physical vertices sleep probability 1/2 and let j jump
surely to site 2. All h_i=1/2. The actual normalized no-portal full-root
law, in bit order 000,001,010,011,100,101,110,111, is

pi_R^0=(2779822848587,12136565503936,5786344507216,
        24809449022720,10925713763524,37005039591680,
        22540864143680,74806876072960)/190790675454303.

It has every deletion residual zero. Nevertheless its entire stabilized
feedback law nu=pi_R^0 F_j has

(nu L_1^0)({0,2})
 =-429949477180782443029835986988779712426332577408224
   /539009793644576310168618312663274330682808973998650717<0.

This is a finite positive-parameter counterexample, with no small-sleep
limit. All increasing-event budgets pi_R^0(U)-nu(U) still pass in this
example, so it specifically separates cone preservation from the
occupation budget. It does not refute the main conjecture.

The initial hypothesis cannot simply be weakened to stochastic
domination by pi_R^0, even for complete feedback on two vertices.
In the two-portal pair above, with p_j=1/2 and j jumping surely to
site 1, delete site 1 from a pi_0 sample whenever both sites are occupied.
The resulting law nu=(31,39,80,0)/150 in support order
empty,{1},{2},{1,2} satisfies nu<=_st pi_0. Nevertheless

(nu F_j)(2 occupied)-pi_0(2 occupied)=294469/31635375>0.

The reference probability is 8/15. Thus the successful cone certificates
in this example use more than an initial occupied-set domination bound.

**Entry at a DAG source also fails.** Prepend a vertex d to the merging
diamond: j jumps surely to d; d jumps to 0 with probability epsilon and
to the original sink otherwise; 0 jumps to 1 or 2 equally; 1 and 2 jump
surely to 3; and 3 jumps to j or the sink equally. Put every sleep
probability, including that at j, equal to one half. The portal enters
the unique source of the open DAG. Its reference is

pi_epsilon^0=Ber(2/(4-epsilon)) at d times Ber(2/3)^4.

For the deletion-1 column at exact support {0,3}, write

phi(B)=1{B={0,3}}-(1/2)1{B={0,1}}-(1/4)1{B={0,1,3}}.

Thus phi requires d and 2 absent, and pi_epsilon^0 phi=0. Nevertheless

pi_epsilon^0 F_j phi=-epsilon/1296+O(epsilon^2).          (145)

To see the linear term, use the open diamond residuals resolved by
portal count H=0,...,5:
(a_H)=(0,0,-1/432,-1/1296,0,0). Conditional on j's first jump, at
epsilon=0 the initial and final d-bits B,Z are independent Ber(1/2),
so its first batch has N=1+B-Z departures and E N=1. Marking one
departure to the diamond contributes N epsilon. For H>=2 the portal
must send at least H-1 particles back to d. Its final bit is therefore
a fresh Ber(1/2), independent of its previous bit and the return count.
The initial jump and final absence at d each contribute one half:
the derivative is (1/2)(1)(1/2)(a_2+a_3)=-1/1296. The H=0,1
residuals vanish separately. Histories with no diamond entrance have
zero residual for every initial d-law, so differentiating its reference
parameter adds nothing. Finite absorption is analytic near epsilon=0.

At epsilon=1/4860, the maintained verifier stores and checks the exact
rational residual, which is less than -1/7000000. An independent
112-state return renewal matches all 1,024 entries of the complete
feedback kernel, retaining every portal return. All 80 adjacent-coordinate
comparisons of the final density relative to pi_epsilon^0 still decrease.
Thus even source-only DAG entry cannot preserve the deletion cone;
the weaker reference domination or density question survives this example.

The maintained `verify_feedback_cone.py` reproduces these checks.
It constructs binary stabilization P by the complete-block inverse in
surplus-IDLA coordinates and computes active addition as P^(-1)D_jP.
It checks every entry against the archived four-site P and inverse,
checks the three- and four-site feedback kernels independently by raw
ARW absorption, and checks the first-contact vector identity in all four
original cases. All three terminal counterexamples use pi_R^0 itself;
it is an exact LP maximizer for the earlier cyclic and later-entry
counterexamples. The source-entry case has its separate complete-return
check above. These are fixed-kernel
certificates, not an arbitrary-graph theorem or a parameter search. The
remaining feedback task is to bound increasing payoffs without requiring
cone preservation; source entry alone does not control (143)'s commutator.

The standard-library verifier `verify_backward_queue_cone.py` checks
the path, fork, and cyclic identities against raw rational ARW absorption,
including 706 primal/dual cone certificates, 336 upper-event comparisons
on two four-vertex feedback graphs, 2,304 weighted-intertwining coefficient
identities, and 768 dominator-identity coefficients. It also verifies the
scope counterexamples and that Section 27's reachable density obstruction
remains in C^0. These finite checks audit the proofs and delimit their
scope; they are not the general arguments.


## 29. Parameter transport, exclusion, and fractional stabilization

**Proved identities and exact finite certificates, 2026-10-19; independent
AI audits passed for the flow, (150a)--(150d), and the logarithm obstruction.**
This exploration starts from the saved stochastic
parameter transfer, rather than from another proposed occupation inequality.
Write r_i=(1-p_i)/p_i, P_r for binary all-active stabilization, E_r=P_r^{-1},
A_z(r) for active addition on stable supports, and
J_x(r)=K(x,s)I+sum_y K(x,y)A_y(r). Odds r_i=0 are allowed here; P_0=I.

The familiar baseline is the IDLA intertwining P_r A_z=D_z P_r, where
D_z adds a walk until the first vacancy or the sink. Its stationary and
strong-stationary-time interpretation is established in
[Levine--Liang, §§1.4,2](https://arxiv.org/html/2110.14008v3#S2).
The new viewpoint below is a consequence of our complete-block and
parameter-transfer identities, not a claim of literature novelty.

### A positive flat connection in sleep parameters

For r<=s the kernel T(r,s)=E_r P_s is stochastic: independently wake each
current sleeper i with probability 1-r_i/s_i, then stabilize at s. A
coordinate with r_i=s_i=0 is left unchanged. The cocycle identity is
T(r,s)T(s,t)=T(r,t). Define G_x(r)=E_r partial_(r_x)P_r. Then

G_x(r)(C,.) = p_x[delta_(C without {x})J_x(r)-delta_C], x in C;
G_x(r)(C,.) = 0,                                      x notin C.  (148)

This formula gives a genuine Markov generator. Its event removes the
resident at x, gives it one K-jump, and stabilizes the remaining sleeping
background plus that particle. The original sink deletes the particle.
The event cannot increase particle count, although it can relocate or
remove several residents. Its clock rate is p_x per unit jump odds.

To derive (148), write the row of the inverse in the fixed IDLA algebra as
E_r(C,.)=delta_empty product_(i in C)[D_i+r_i(D_i-J_i^D)],
where J_i^D=K(i,s)I+sum_y K(i,y)D_y. This is multiaffine in r.
Differentiating at x and right-multiplying by P_r gives
(partial_x E_r)(C,.)=delta_(C without {x})(A_x-J_x)E_r.
At a vacant x, A_x=p_x B_x+(1-p_x)J_x and B_x inserts its sleeper.
Thus -partial_x E_r P_r is exactly (148). Differentiating E_r P_r=I
identifies this matrix with G_x.

Along any coordinatewise nondecreasing differentiable path r(t),

partial_t T(r(0),r(t))=T(r(0),r(t)) sum_x r'_x(t)G_x(r(t)).       (149)

Hence binary ARW stabilization P_r is the terminal law of an inhomogeneous
Markov process on stable supports, starting at the binary input at odds
zero. In coordinates u_x=log(1+r_x)=-log p_x, each site's potential event
clock has rate one per unit u_x; absent residents suppress the event.
The total number of potential events on a finite monotone path is Poisson
with mean sum_x[u_x(final)-u_x(initial)]. The event kernels still involve
ARW avalanches, so this is a structural representation, not a cheaper
independent-particle sampling algorithm. It gives a coupling under which
particle count decreases when sleep probabilities decrease.

The compatibility identities are

partial_x G_y-partial_y G_x+[G_x,G_y]=0,
partial_x A_z=A_z G_x-G_x A_z.                                  (150)

They follow directly by differentiating E_r=P_r^{-1}. The connection is
flat but its generators at one parameter vector need not commute. For the
two-site kernel K(0,1)=K(1,0)=K(0,s)=K(1,s)=1/2 at r=(1,1),
[G_0,G_1]({0},empty)=2/45. The parameter dependence in (150) cancels this
commutator. Ordinary inclusion order is not preserved by the flow: in
that same example,

(delta_{ {0,1} }-delta_{ {0} })G_1 f_0=-22/195,                  (151)

where f_v(C)=1{v in C}; the two initial target bits are equal.

The flow also gives a conservation identity for potential-weighted
occupations. Let M=I-K_V and G=M^{-1}. For a stable C containing x, let
U(C,x;z) be the expected number of departures from z during the forced
resident-release avalanche, including its compulsory first jump from x.
These expectations are finite at positive sleep and jump probabilities.
Particle conservation, averaged over the independent jump destinations,
gives

E[1_terminal]-1_C=-U(C,x;.)M.

For any fixed all-active input eta, surplus IDLA writes mu_eta^r=nu_eta P_r
with nu_eta independent of r. Consequently partial_x mu_eta^r=mu_eta^r G_x,
and, writing q(eta;r) for its row of occupation probabilities,

partial_(r_x) q(eta;r) G
 =-p_x E_(C~mu_eta^r)[1{x in C} U(C,x;.)] <= 0.              (150a)

Thus every weighted occupation sum with weight h=G f, f>=0, decreases
as jump odds increase. Equivalently these are the nonnegative excessive
potentials h with (I-K_V)h>=0. For h=1, the right side counts expected
sink losses; for h=G 1, it counts all departures in the release avalanche.
In logarithmic coordinates u_x=-log p_x, the factor p_x disappears.
This parameter-response identity compares weighted means, not individual
occupation probabilities or inputs differing by one particle.

### The same transport orders statistical information

Regard the binary all-active input B as an unknown parameter and its
terminal sleeping support as an observation. The factorization
P_s=P_r T(r,s) says that the experiment at larger jump odds is obtained
by applying one stochastic postprocessing rule to the experiment at r,
independently of B. Thus P_s is a garbling of P_r in the statistical
comparison of experiments. Surplus IDLA extends the same statement to
any family of finite all-active inputs, with the same postprocessing T.

For the full binary-input experiment the parameter order is exact:
P_s is a stochastic postprocessing of P_r if and only if r<=s,
for finite strictly positive odds. Invertibility makes the candidate
postprocessing unique, namely E_r P_s. Its singleton row satisfies

    T(r,s)({i},empty)=(1-r_i/s_i)P_s({i},empty).

This follows by multiplying the singleton inverse row
(1+r_i)delta_{i}-r_i sum_j K(i,j)delta_j-r_i K(i,sink)delta_empty
by P_s and using the single-particle first-instruction equation at s.
The final factor is strictly positive by a finite all-jump route to the
sink. Positivity therefore forces r_i<=s_i at every i. This necessity
uses the entire binary-input experiment; a restricted input family can
have additional coincidences.

In particular, for any such eta,xi and any convex f,

    D_f(mu_eta^s || mu_xi^s) <= D_f(mu_eta^r || mu_xi^r),       (150b)

with the usual extended-value convention. To see this directly when
the reference law is positive, its output likelihood ratio is the
conditional expectation of the input likelihood ratio under the joint
law nu(C)T(C,D). Conditional Jensen gives (150b); approximation handles
zeros. Total variation and relative entropy are special cases. For every
prior on the input label, mutual information with the observed support
also decreases, by the same joint Markov construction.

The full-input stationary laws pi_r=delta_V P_r satisfy pi_s=pi_r T too.
Hence (150b) also applies with mu_xi replaced by pi: distance of a fixed
all-active input from the stationary law decreases with jump odds.
This applies to a fixed mixture of inputs independent of the instruction
books, including a fixed driven law started empty after a fixed number of
additions. The mixing law must be the same at r and s.
It does not compare arbitrary fixed sleeping initial backgrounds.

If s_i>r_i for every site, T(C,empty)>0 for every C: wake every resident
and follow positive-probability routes that lose all of them. Consequently
beta=min_C T(C,empty)>0, and subtracting the common atom gives

    ||mu T-nu T||_TV <= (1-beta)||mu-nu||_TV.                 (150c)

This strict contraction statement concerns finite positive endpoint odds;
it need not hold for a change at only one site, since supports avoiding
that site are fixed. No individual target sign follows from information
contraction.

For completeness, strictly positive laws evolving with a coordinate
generator G_x obey an exact entropy dissipation formula. Put u_C=mu_C/nu_C.
Differentiating sum_C mu_C log u_C and using zero row sums gives

    partial_x D(mu||nu)
      =-sum_(C!=D) nu_C G_x(C,D)
           [u_C log(u_C/u_D)-u_C+u_D] <=0.                  (150d)

The nonnegative brackets describe information lost during the same physical
resident-release events as (148). These are consequences of the already
proved parameter transport, not a new occupation comparison.

There is a simple boundary to this information loss. If one retains the
entire directed jump-count vector J as well as the terminal support C,
particle conservation reconstructs the all-active input exactly:

    eta_i=1{i in C}+sum_y J_(i,y)-sum_(y in V)J_(y,i),

where outgoing sums include sink jumps. Thus discarding the currents is
essential to the statistical degradation of the support observation.
This reconstructs initial particle counts, not instruction books; it does
not remove the record-matching obstruction F3.

### The all-sleep tangent is long-range exclusion

At r=0, all lone particles sleep immediately. A release from occupied x
moves one particle along an ordinary K-walk until it finds a vacancy in
C without {x}, including a return to x, or reaches the sink. Other
residents are retained. This is the finite killed version of long-range
exclusion; compare [Liggett, *Long Range Exclusion Processes*](https://doi.org/10.1214/aop/1176994618).
Its attractiveness has a direct proof here. For C subset D use the same
walk from x. If x belongs to both, the first vacancy for C without {x}
comes no later than that for D without {x}. The former landing site is
either already present in the upper configuration or is also its first
vacancy. If x belongs only to D, removing it cannot remove a lower
resident. Thus the single-event map preserves inclusion. Any nonnegative
combination of these generators is attractive. Equation (151) shows why
this boundary structure does not persist throughout the parameter flow.

### A stochastic fractional power can lose target monotonicity

There is also a distinct homogeneous-embedding question: is log P_r a
Markov generator? The general answer is **not established here**. Positive
inhomogeneous transport alone does not answer it. This is the usual
[Markov embedding problem](https://arxiv.org/abs/1903.08736), applied to the
binary stabilization matrix rather than to the driven addition chain.

A fixed three-site example gives a precise distinction. Let
0->1->2 deterministically, and 2->0 or s with probabilities 1/2 each.
Put r_i=epsilon, so p_i=1/(1+epsilon), and let
L(epsilon)=log P_epsilon=sum_(j>=1)epsilon^j L_j near epsilon=0.
The complete-block inverse E_epsilon is a polynomial of degree three.
Writing H=E_epsilon-I determines the Taylor coefficients exactly from
log P_epsilon=-log(I+H).

There are 34 off-diagonal entries allowed by nonincrease of particle
count. Exact rational coefficient calculation gives:

| First nonzero degree | Number of entries | Smallest coefficient |
|---|---:|---:|
| 1 | 16 | 1/2 |
| 2 | 13 | 1/8 |
| 3 | 5 | 1/24 |

All other off-diagonal entries vanish identically by the count-block
structure. Since log P has zero row sums and is analytic near zero,
L(epsilon) is a Markov generator for all sufficiently small positive
epsilon. Thus every P_epsilon^t=exp(t L(epsilon)), t>=0, is stochastic.
No quantitative neighborhood size is asserted.

Nevertheless take A={2}, B={1,2}, and target 2. Define
ell(M)=(delta_B-delta_A)M f_2. The exact coefficients are

ell(L_1)=0, ell(L_2)=-1/4, ell(L_3)=1/6, ell(L_4)=-1/24,
ell(L_1^2)/2=3/4.

There is a short direct check of the decisive quadratic coefficient.
Write E_epsilon=I+epsilon B_1+epsilon^2 B_2+epsilon^3 B_3. If R_S(C,.)
releases all residents in S with compulsory first jumps, then performs
IDLA, the inverse gives

    B_1(C,.)=|C|delta_C-sum_(x in C)R_{x}(C,.),
    B_2(C,.)=binom(|C|,2)delta_C
              -(|C|-1)sum_(x in C)R_{x}(C,.)
              +sum_({x,y} subset C)R_{x,y}(C,.).

For A={2} and B={1,2}, direct IDLA releases give

| Support | B_1 f_2 | B_2 f_2 | B_1^2 f_2 |
|---|---:|---:|---:|
| A | 1 | 0 | 1 |
| B | 1 | 1 | 5/2 |

Since L_1=-B_1 and L_2=B_1^2/2-B_2, this proves ell(L_2)=-1/4
and ell(L_1^2)/2=3/4 without using the logarithm-series verifier.

Consequently, for each fixed t,

(delta_B-delta_A)P_epsilon^t f_2
 =t(3t-1)epsilon^2/4+O(epsilon^3).                            (152)

In particular its stochastic quarter-power has target increment
-epsilon^2/64+O(epsilon^3)<0. At time one the leading coefficient is
+1/2, consistently with the proved three-site ARW theorem. The fractional
kernel is a different process, not an ARW counterexample. A monotone
complete stabilization need not admit monotone fractional iterations.

`verify_sleep_transport.py` checks (148) against independent raw resident
release, and checks the transfer law by independently waking and stabilizing
under both coordinate orders. The two fixed motifs give 224 generator
entries, 208 curvature entries, four physical path comparisons and 80
zero-odds stabilization entries. The logarithm certificate evaluates the
whole inverse polynomial against a separate 28-state raw absorption system;
it also verifies all 34 leading signs and (152) in rational arithmetic.
The Markov-flow and analyticity arguments supply the general/existence
statements; the finite computations certify the explicitly scoped example.
An independent 12-state IDLA calculation, importing no project code,
reproduced the leading-sign table and the coefficients through degree three.
The displayed degree-four coefficient was checked by the main verifier only.


## 30. Loss markers as mixtures of forest laws

**Derived identities, with exact checks and main-thread review, 2026-10-19.**
These describe auxiliary markers, not an occupation comparison.

Fix an initial root set Q. Write k_i=K(i,s), p_i for the sleep probability,
q_i=1-p_i, r_i=q_i/p_i, d_R=det(I-K_RR), and

\[
Z_R=d_R\sum_{C\subseteq R}\frac{r_{R\setminus C}}{d_C},\qquad
g_R=p_R Z_R/d_R.
\]

Subscripts on scalar vectors denote products over the indicated set. Every
exit outside Q retains its label; s denotes the original sink. All d_R are
positive. In the creator sampler, C has law

\[
\nu_Q(C)=\frac{d_Q}{Z_Q}\frac{r_{Q\setminus C}}{d_C}.
\]

Independently, F is the exit forest on Q with probability w(F)/d_Q. Conditional
on C, the balanced prefix alpha has law d_C w(alpha), supported in C. Let
L be its vertex range. The released-input loss set and its companion are

\[
E=\{i\notin C:F_i=s\},\qquad D=\{i:F_i=s\}.
\]

D records reserved forest arrows, not actual physical particle losses.

**An exact determinantal mixture.** Put M=I-K_QQ, G=M^{-1}, and
H=G diag(k). Fixing the sink arrows at a set S and deleting those vertices
leaves the forest on Q\S, with arrows into S treated as labelled exits.
Consequently

\[
\mathbb P(S\subseteq D)=k_S\frac{d_{Q\setminus S}}{d_Q}
 =\det H_{SS}.                                                   \tag{M1}
\]

This is the directed forest-root determinantal law; compare
[Avena--Gaudilliere, Theorem 3.4](https://arxiv.org/pdf/1310.1723).
The displayed identity also proves it directly in the present notation.

The joint weight of C and alpha is

\[
\mathbb P(C,\alpha)=\frac{d_Q}{Z_Q}
 r_{Q\setminus C}w(\alpha)\,\mathbf1\{L\subseteq C\}.
\]

Thus, conditional on alpha, C contains L and includes each site outside L
independently with probability p_i. This conditional law depends only on L.
The forest is independent of both variables. It follows that

\[
 E\mid L \quad\text{is D independently thinned with probabilities}\quad
 t_i(L)=q_i\mathbf1\{i\notin L\}.                                \tag{M2}
\]

In particular its conditional inclusion kernel can be taken as
H diag(t(L)). Its probability generating polynomial is

\[
\mathbb E[z^E\mid L]
 =\det\bigl(I+H\operatorname{diag}(t_i(L)(z_i-1))\bigr).           \tag{M3}
\]

Conditioning on C instead gives the analogous restriction of H to Q\C.
There is no averaging of the kernel in (3): the determinant itself must be
averaged over L.

If K_QQ is reversible with respect to positive weights m_i, then
B=diag(sqrt(m)) M diag(1/sqrt(m)) is symmetric positive definite.
Its quadratic form is the internal Dirichlet form plus exit killing, so
B >= diag(k). Hence

\[
\widehat H=\operatorname{diag}(\sqrt{k})B^{-1}
                 \operatorname{diag}(\sqrt{k})
\]

is a positive contraction with the same principal minors as H. This includes
zero k_i by continuity. Thinning gives the positive contraction
diag(sqrt(t)) \widehat H diag(sqrt(t)). Thus each conditional law (2) is
strongly Rayleigh in the reversible case, by
[Borcea--Branden--Liggett, Proposition 3.5](https://arxiv.org/pdf/0707.2340).
For directed kernels, (1)--(3) remain determinantal identities without a
Hermitian or strong Rayleigh assertion.

The range mixing law is finite and exact, despite unbounded prefix lengths:

\[
\mathbb P(L=T)=\frac1{g_Q}\sum_{A\subseteq T}
                   (-1)^{|T\setminus A|}g_A.                    \tag{M4}
\]

Indeed, the marginal weight of alpha is p_L w(alpha)/g_Q, and g_A sums these
unnormalized weights over words supported in A. Equivalently L is the range
of the retained Poisson soup from §19. Its atoms need not be the raw popped
cycles, and nothing here erases the balanced prefixes from physical ARW.

**A reversible three-site counterexample.** Take Q={0,1,2}, p_i=1/2,
K(i,j)=1/4 for distinct internal vertices, K(i,i)=0, and K(i,s)=1/2.
Then

\[
d_1=1,\ d_2=15/16,\ d_3=25/32;\qquad
Z_1=2,\ Z_2=61/16,\ Z_3=53/8.
\]

Here d_j and Z_j depend only on the set cardinality. Formula (64) gives

\[
\mathbb P(i\in E)=61/212,\quad
\mathbb P(i,j\in E)=4/53,\quad
\mathbb P(Q\subseteq E)=1/53.
\]

For e_j the elementary symmetric polynomials in z_0,z_1,z_2, the full law is

\[
\boxed{\mathbb E[z^E]=\frac{73+33e_1+12e_2+4e_3}{212}.}          \tag{M5}
\]

Thus each individual empty/singleton/pair/triple subset has probability
73/212, 33/212, 12/212, 4/212, respectively. Its count polynomial is

\[
\mathbb E[z^{|E|}]=\frac{4z^3+36z^2+99z+73}{212}.
\]

The numerator's discriminant is -15984, so it has nonreal conjugate zeros.
A real stable generating polynomial must have a real-rooted diagonal
specialization. Therefore E is **not strongly Rayleigh**, even with a
reversible K. It is not a Hermitian positive-contraction DPP either. In fact
no Hermitian correlation kernel is possible, since its count polynomial
would then factor into linear factors with real zeros.

The same example makes (2) explicit. Let J denote the all-ones 3x3 matrix.
Then H=(2I+J)/5, and

\[
\mathbb P(L=\varnothing)=50/53,\quad
\mathbb P(L=\{i,j\})=5/318\ \text{for each pair},\quad
\mathbb P(L=Q)=1/106.
\]

There are no singleton ranges. For L empty, E has kernel H/2; for a pair
range, E is a Bernoulli(3/10) point at the remaining vertex; for the full
range, E is empty. Thus (5) is a five-component mixture of strongly Rayleigh
laws. On the diagonal the mixture reads

\[
\frac{50}{53}\frac{(1+z)(4+z)^2}{50}
 +\frac{15}{318}\frac{7+3z}{10}+\frac1{106}.
\]

The marker covariance is still -329/44944. This example is consistent with
§19's conditional pairwise negative correlations. Stability fails at real
negative fields: for the unnormalized numerator f in (5),

\[
\partial_0 f\partial_1 f-f\partial_{01}f
 =213+104z_2+12z_2^2,
\]

which is -37/3 at z_2=-13/3. The failure does not itself refute any
negative-dependence property restricted to positive fields.

**A companion with a surviving sign.** For disjoint A,B subset Q, independence
of C and F gives the mixed identity

\[
\boxed{\mathbb P(A\subseteq E,\ B\subseteq D)
 =\frac{a_A Z_{Q\setminus A}}{Z_Q}\,
    k_B\frac{d_{Q\setminus(A\cup B)}}{d_{Q\setminus A}},
 \qquad a_i=r_i k_i.}                                         \tag{M6}
\]

In particular for i != j,

\[
\operatorname{Cov}(\mathbf1_{i\in E},\mathbf1_{j\in D})
 =-\mathbb P(i\notin C)\,k_i k_j G_{ij}G_{ji}\leq0.             \tag{M7}
\]

Conditional on C the same assertion holds, with the prefactor replaced by
the fixed indicator i notin C; conditional on L it is replaced by t_i(L).

In the saved two-site half-sink example, conditional on C={1}, E_0=D_0.
The forests (s,s), (s,0), (1,s) are equiprobable. Hence

\[
\operatorname{Cov}(E_0,D_1\mid C=\{1\})=-1/9,
\qquad
\operatorname{Cov}(E_0,X_1\mid C=\{1\})=14/135.
\]

For the second value, X_1 has conditional means 1 and 8/15 according as
E_0 occurs or not, and P(E_0|C)=2/3. Thus the raw forest-root companion has
an exact negative correlation while terminal occupation does not. The
reserved arrow F_1=s may never be emitted by the physical particle.

**A connected-loop identity for comparison with the transfer-flow branch.**
Let w_T be §19's rejected-loop intensity with vertex range exactly T, and
W(R)=sum_{nonempty T subset R} w_T. That section proves
Z_R=p_R^{-1} exp(-W(R)). Substitution into (64) gives

\[
\boxed{\mathbb P(S\subseteq E)
 =q_S k_S\exp\left(\sum_{T\subseteq Q:T\cap S\ne\varnothing}w_T\right).}
                                                                    \tag{M8}
\]

When all k_i in S are positive, its logarithm is a singleton sum plus a
nonnegative loop-coverage sum. In particular

\[
\frac{\mathbb P(i,j\in E)}{\mathbb P(i\in E)\mathbb P(j\in E)}
 =\exp\left(-\sum_{T\supseteq\{i,j\}}w_T\right).                \tag{M9}
\]

More generally, for disjoint A,B with |B|>=2, if
h(S)=log P(S subset E),

\[
\Delta_B h(A)=(-1)^{|B|+1}
 \sum_{T\supseteq B:T\cap A=\varnothing}w_T.                    \tag{M10}
\]

Only ranges of closed walks contribute. These are logarithmic inclusion
interactions, not cumulants of the probability generating polynomial and
not entries of log(P_r). In the construction (2), ordinary forest-root
exclusion is partially cancelled by positive correlations of the vacant
retained-loop range; precisely the rejected-loop intensity remains in (9).
No statement identifying this expansion with a homogeneous Markov
generator has been proved here.

**Distinguish nonsymmetric determinant kernels.** The non-stability result
does not rule out arbitrary nonsymmetric determinantal representations.
For the directed three-cycle K(i,i+1)=K(i,s)=1/2, p_i=1/2, the marker law is

\[
\mathbb E[z^E]=\frac{\prod_i(3+z_i)-7}{57}.
\]

It is not strongly Rayleigh, but it has a real nonsymmetric inclusion kernel
J=A^{-1}, where

\[
A=\begin{pmatrix}4&-1&0\\0&4&-1\\-7&0&4\end{pmatrix}.
\]

Indeed det(A+diag(z_i-1))/det(A) is the displayed polynomial. Existence of
a single nonsymmetric kernel for general loss-marker laws is not settled
by this calculation. The conditional mixture (2) is general.

**Scope and verification.** The mixture and the mixed forest-root identities
give a precise nearby object; they do not compare physical ARW occupations.
An extension having E as a coordinate marginal cannot be strongly Rayleigh
in general, since marginalization preserves stability. This does not rule
out broader conditioning followed by nonlinear maps.
Known loop- and forest-conditioned occupation obstructions still apply.
Nothing here upgrades §22's failed proper-target coefficient signs.

The new project file
`verify_marker_geometry.py` verifies (1)--(6) on
exactly three fixed examples, using rational arithmetic and finite forest
enumeration; the infinite prefixes are integrated using (4). It also checks
the two non-stability examples, the directed nonsymmetric kernel, and the
two-site covariance calculation. Result: 63 mixed marker/root identities
and 19 conditional marker-pair checks pass. No sampling or parameter sweep.


### Positive external fields: a scoped Rayleigh theorem

**Proved, with an exact generic polynomial certificate and main review.**
For F(z)=E[z^E], write Delta_ij F=F_i F_j-F F_ij. Under positive
external fields, Cov_z(E_i,E_j)=-z_i z_j Delta_ij F/F^2.
The marker law is Rayleigh, meaning Delta_ij F>=0 for every z>0,
whenever every SCC of K_QQ has at most three vertices. On arbitrary
root domains the inequalities hold at least for 0<z_i<=1. These are
auxiliary-marker statements, not occupation comparisons.

Here is the three-site coefficient that is not supplied by conditional
negative correlation alone. Remove self-loops by c_i=1-K(i,i),
K'(i,j)=K(i,j)/c_i for j!=i, and r'_i=c_i r_i. Both Z_R and
a_i=r_i K(i,s) are unchanged, so assume K(i,i)=0. Put Q={i,j,k},
b_l=1+r_l, kappa=K(k,s),

    u=K(i,j)K(j,i), v=K(i,k)K(k,i), w=K(j,k)K(k,j),
    tau=K(i,j)K(j,k)K(k,i)+K(i,k)K(k,j)K(j,i),
    d=1-u-v-w-tau>0, H=b_i b_j-1.

Then Z_i=b_i, Z_ij=b_i b_j-u(b_i b_j-1), and

    Z=1+d[b_i b_j b_k-1+r_i w/(1-w)+r_j v/(1-v)+r_k u/(1-u)].

Since Delta_ij F is independent of z_i,z_j, write it as alpha+beta t+
gamma t^2 for t=z_k. Conditional NC given E_k=0, proved in (64)--(65),
gives alpha>=0. Substituting m_S=a_S Z_(Q\S)/Z gives

    gamma=(a_i a_j a_k^2/Z^2)uH,
    beta=(a_i a_j a_k/Z^2) B,
    B=[2u+tau+r_k{2u(1-kappa)+tau}]H
      +r_k[tau+u(v+w+tau)/(1-u)]
      +r_i w(u+v+tau)/(1-w)+r_j v(u+w+tau)/(1-v)>=0.         (M11)

Equivalently B=Z_jk b_i+Z_ik b_j-Z-Z_ij b_k-2a_k uH. All
denominators are positive and 0<=kappa<=1. Clearing the three
denominators gives the exact sparse-polynomial identity checked by
`symbolic_rayleigh_middle_coefficient()` in `verify_marker_geometry.py`.
The previous triangle and directed-cycle examples have respectively
Delta_ij=(213+104t+12t^2)/212^2 and 7(3+t)/57^2. They are Rayleigh
although their strong Rayleigh property fails as already recorded.

Principal determinants factor across the SCCs of K_QQ, hence Z_R does
too. The marker inclusion moments factor, so different SCC marker
vectors are independent. This proves the stated component theorem;
it is not a theorem about selecting three coordinates of a larger SCC.

For the one-sided field claim on arbitrary Q, fix fields 0<z_b<=1
at sites B and put c_b=a_b(1-z_b). Replace Phi_B from Section 19 by

    Phi_(B,z)(U)=sum_(D subset B)(-1)^|D| c_D Z_(U union (B\D)).

This is positive by its field-weight interpretation and obeys
Phi_(B+b,z)(U)=Phi_(B,z)(U+b)-c_b Phi_(B,z)(U). The proof following
(65) applies unchanged for 0<=c_b<=a_b: log of the ratio being
subtracted has nonnegative nonconstant Boolean coefficients, so the
convergent expansion of log(1-ratio) preserves nonpositive higher
log-Mobius coefficients. For R=Q\B, the tilted inclusion moments are
a_S Phi_(B,z)(R\S)/Phi_(B,z)(R). Log-submodularity gives their pairwise
NC. A tested pair's own field values do not affect its Rayleigh sign.
For z_b>1 the subtraction becomes addition and this proof gives no sign.
Unrestricted fields in larger SCCs remain open here.

The nearby principal-minor law has Lorentzian homogenization
det(w_0 I+diag(w_i)(I-K_QQ)), by
[Brändén--Huh, Theorem 4.15](https://web.math.princeton.edu/~huh/Lorentzian.pdf).
This does not classify the marker law: its polynomial is instead
Z_Q^(-1) sum_S a_S Z_(Q\S) product_(i in S)(z_i-1). A transformation
preserving the Lorentzian property between these formulas is missing.

For a directed simple cycle of any length n>=2 there is a stronger
closed formula. Use effective parameters after the preceding self-loop
reduction. If t_i=K(i,i+1) are the cycle-edge probabilities, put
theta=product_i t_i<1, beta=theta(1-product_i p_i), and
u_i=(1-p_i)K(i,s). Every proper induced subgraph is acyclic, so
Z_R=product_(i in R)(1+r_i) for R proper, whereas
Z_Q=product_i(1+r_i)(1-beta). The inclusion moments then give

    F_E(z)=[product_i(1-u_i+u_i z_i)-beta]/(1-beta),
    Delta_ij F_E=beta u_i u_j/(1-beta)^2
                      product_(h distinct from i,j)(1-u_h+u_h z_h).
                                                                  (M12)

Thus this law is Rayleigh under every positive external field. It is
a Bernoulli product law with mass beta subtracted from its empty atom,
then renormalized. This subtraction is valid because
product_i(1-u_i)>=theta>beta. It also survives arbitrary positive
fields: replace u_i by u_i z_i/(1-u_i+u_i z_i) and beta by
beta/product_i(1-u_i+u_i z_i). Together with SCC factorization, the
Rayleigh theorem therefore allows each component to be either of size
at most three or a directed simple cycle of arbitrary length.

Distinct exit labels can be retained as auxiliary markers indexed by
(root,label). Set u_(i,l)=(1-p_i)K(i,l) and replace the i-th factor
in (M12) by 1-sum_l u_(i,l)+sum_l u_(i,l)z_(i,l). Different roots
give the same factorized Rayleigh difference; two different labels at
one root are mutually exclusive, so their mixed derivative is zero.
For exits to outside vertices this is only a marker identity, not the
original-sink endpoint deletion identity (64). The focused
`labelled_cycle_audit()` compares all 81 marker masses on one heterogeneous
four-cycle with self-loops and two exit labels, by direct latent-subset
and forest enumeration. Main reran it successfully. No larger-SCC
parameter sweep was performed.

The cycle law in fact has negative association after arbitrary positive
fields and positive-probability partial assignments (CNA+), including
all root-and-exit-label coordinates. The following elementary lemma is
enough. For a measure mu on binary vectors with minimal atom 0, let
nu=(mu-beta delta_0)/(1-beta), where 0<=beta<=mu(0) and beta<1. Then

    Cov_nu(f,g)=Cov_mu(f,g)/(1-beta)
      -beta[mu(f)-f(0)][mu(g)-g(0)]/(1-beta)^2.               (M13)

Thus feasible subtraction at zero preserves negative association.
The categorical product base is NA: within a root, disjoint increasing
functions have covariance
-(sum_(a in I)u_a[f(e_a)-f(0)])(sum_(b in J)u_b[g(e_b)-g(0)]).
Induction over independent roots using conditional covariance proves
the assertion for disjoint coordinate sets spanning several roots.
Fields and assignments preserve this categorical product structure.

For the subtracted law, a positive field z replaces mu by mu^z and
beta by beta/mu(z^X), with atom feasibility preserved. Conditioning
on any assignment containing a one removes the subtracted atom.
Conditioning on only zeros gives the same form, with base mu(.|H)
and subtraction beta/mu(H); positive probability ensures mu(H)>beta.
These identities and (M13) prove CNA+. Independent cycle components
retain the property. This conclusion is stronger than the pairwise
Rayleigh calculation for this class; it does not establish the
marker--occupation covariance sign in (65).


## 31. Finite boundary responses and associative gluing

**Derived identities, with exact checks and main-thread review, 2026-10-19.**
The full state includes local sleeping supports and pending arrivals; no
occupation comparison or cone invariance is asserted.

### The response object

Let D be a finite ARW domain, I a set of input sites in D, and O a finite
set of labeled absorbing exits. Include the original sink as an exit if
its count is wanted; otherwise set its generating variable to 1. Jumps
within D retain their usual kernel. Assume a killed exit is accessible
from every vertex and 0<p_x<1. All statements below concern the law
averaged over independent instruction stacks.

For a sleeping support B subset D and one active addition at i in I,
stabilize completely inside D. Let C be the terminal sleeping support
and k=(k_o:o in O) the exit-count vector. Define the matrix polynomial

    A_i^D(z)[B,C] = sum_k P_B(C,k | active addition at i) z^k.       (B1)

The matrices have 2^|D| rows and columns. If the sink is counted, mass
conservation gives |C|+sum_o k_o=|B|+1, so their degree is at most |D|+1.
Thus a finite list of finite polynomial matrices describes responses
to arbitrarily many successive inputs; it is not necessary to tabulate
each possible incoming batch separately. Coefficients can be rational
functions of the underlying kernel and sleep parameters; for fixed
rational parameters they are rational numbers.

**Exact identity.** For any input word i_1,...,i_m,

    A_(i_1)^D(z)...A_(i_m)^D(z)

is the joint generating matrix of terminal support and total labeled
exits. In particular,

    A_i^D(z) A_j^D(z) = A_j^D(z) A_i^D(z).                       (B2)

Proof: sequential stabilization and simultaneous addition have the same
terminal support and boundary instruction counts by the Abelian property.
Unused stack tails after each complete stabilization are independent
fresh stacks. Conditioning on previously revealed supports and exit
counts therefore gives precisely matrix multiplication. The equality
holds for arbitrary initial sleeping B; it is an addition/stabilization
identity, not sleeping-input attractiveness.

The necessary retained information is the joint law of C and k, not
their separate marginals. For a fixed routing network and terminal
observables, exit chronology inside a completed block is unnecessary.
If one conditions on intermediate exit counts, that conditioning must
be applied to the corresponding matrix coefficients in their actual
order. Its effect on the next pass is retained by the posterior over C.

Activity flags are unnecessary at these particular stopping times,
because every module is locally stable. This does not contradict RN §24:
stopping before local stabilization requires activity information.
For fixed, already exposed instruction stacks, (B1) is not a sufficient
deterministic state: consumed stack indices or unread tails must also
be retained. The finite reduction here uses iid tails at stopping lines.

### Exact gluing, including feedback

Take disjoint domains D_1,...,D_r and route specified output ports to
specified input sites of other domains. Remaining outputs stay absorbing.
Assume the resulting physical graph can reach its absorbing exterior
from every vertex (the original sink after all feedback is closed). Ports here are
directed exit/entry interfaces, not duplicate copies of a shared physical
vertex. A shared vertex can be assigned to one component, or represented
by a separate singleton component.

At a macro time retain

    (B_1,...,B_r; m; ell),

where B_a is the sleeping support stored in component a, m_i counts
pending active arrivals at input i, and ell counts external exits already
produced. Pick any pending input i, remove one message, sample (C,k)
jointly using (B1), replace the selected B_a by C, and route all components
of k. Internal exits become pending arrivals; external exits increment ell.
Do not replace the joint distribution of these states by a product law.

With initial particle mass N, the sum of stored sleepers, pending messages
and recorded external exits is N. Hence the state space is finite.
For an elementary response of the combined domain, N<=sum_a |D_a|+1.

**Why delayed arrivals are legitimate.** A stored sleeper together with
an undelivered active message at its site represents two active particles
in the physical configuration. The representation need not itself look
like the physical active/sleeping configuration. Nevertheless an active
addition commutes with each toppling that was legal before that addition.
For example, a sleep instruction on one active particle followed by an
addition produces two active particles; adding first and then consuming
that sleep instruction produces the same two active particles. Jump
instructions satisfy the same identity. Thus the local topplings made
while other messages are delayed remain legal after those messages are
inserted, and give the same reconstructed physical configuration.
Induction over blocks embeds the queue procedure into a legal ARW
toppling sequence. Empty queues and locally stable supports mean the
physical configuration is stable.

**Termination.** Let h_x be expected ordinary-walk time from x to the
absorbing exterior on the glued graph, and let Phi be the sum of h over all
particles, including pending ones at their destination sites. Each raw
instruction at an active x has conditional expected change

    E[Delta Phi | current information, topple x] = -(1-p_x).

Sleep changes no particle positions, while the jump equation is
h_x=1+sum_y K(x,y)h_y. Put q_min=min_x(1-p_x)>0. Stopping after any
fixed number of instructions bounds their expected number before
stabilization by Phi_initial/q_min. Every response call consumes at
least one instruction. Consequently the queue chain is transient and
has finite mean absorption time, for any causal fixed priority rule.
This argument also justifies the unbounded number of possible returns.

Write Q for the transient-to-transient queue matrix, S for its
transient-to-terminal matrix, and J for the selector of initial states
(B, one pending input, zero external counts). The exact composition is

    R_glued = J (I-Q)^(-1) S = J sum_(t>=0) Q^t S.               (B3)

This is the full joint response, including terminal configurations and
external exit counts. It has nonnegative coefficients without requiring
any creation-cone invariant.

**Associativity in law.** Every grouping of components computes the
same killed ARW response on their union. Thus, with consistent labels,

    (R_1 glue R_2) glue R_3 = R_1 glue (R_2 glue R_3).            (B4)

Grouping may change the queue chain and scheduler, so the two Q matrices
need not be identical. Abelianity identifies their terminal responses.
Within any one such chain, feedback elimination is ordinary associative
matrix elimination. For a partition of transient states into U,W,

    H_U = [I-Q_UU-Q_UW(I-Q_WW)^(-1)Q_WU]^(-1)
          [S_U+Q_UW(I-Q_WW)^(-1)S_W].                          (B5)

This follows by substituting the W equation into H=QH+S. Transience
makes all displayed inverses legitimate. No signed-renewal convergence
assumption is used.

The stochastic-processor viewpoint has a primary-source antecedent in
[Bond--Levine, Abelian networks I, §3.9 and Lemma 4.7](https://lionellevine.github.io/abelian-networks-I-foundations.pdf).
Section 3.9 lists ARW among stochastic abelian networks; Lemma 4.7 treats
a halting subnetwork as one processor. Their stochastic theory is deferred
there. The iid-tail compression and finite queue proof above specify what
is used here; no novelty claim or general theorem about arbitrary
stochastic transducers is intended.

### A minimal noncommuting response

Consider one site with sleep probability p, no self-loop, and jump
generating polynomial w(z)=s+sum_o a_o z_o, where s is the probability
of an uncounted original-sink jump and s+sum_o a_o=1. Put q=1-p. In
state order empty, sleeping, its elementary response is

              [ q w      p   ]   [1] [q w, p].
    A(z) =    [ q w^2    p w ] = [w]

An input at an occupied site forces one departure before the remaining
particle has its fresh sleep/jump decision. The outer-product expression
gives the exact local identity

    A(z)^2 = w(z) A(z),  A(z)^m = w(z)^(m-1) A(z), m>=1.       (B6)

Now set p=1/2 and let each departure go to one portal or the original
sink with equal probability. Then w(z)=(1+z)/2 and A(z)=A_0+z A_1+z^2 A_2,
where

    A_0 = [1/4  1/2],   A_1 = [1/4  0  ],   A_2 = [0    0].
          [1/8  1/4]         [1/4  1/4]         [1/8  0]

Their coefficient commutator is

    A_0 A_1 - A_1 A_0 = diag(1/8,-1/8).                       (B7)

Starting empty and making two additions, the unnormalized terminal
rows for successive portal counts (0,1) and (1,0) are respectively

    e_empty A_0 A_1 = (3/16, 1/8),
    e_empty A_1 A_0 = (1/16, 1/8).

Thus terminal occupation conditioned on those histories is respectively
2/5 and 2/3, although both histories have total portal count one. Summing
the two rows gives the correct total-count-one row (1/4,1/4).

This distinguishes the commutative addition algebra at a common exit
variable from the noncommutative algebra of conditioned pass operators.
It is a minimal warning relevant to RN §§27--28, not a proof of any new
cone failure. Associative gluing does not require commuting coefficients.

### The smallest feedback check

Take two physical sites a,b, with a->b surely, b->a or sink equally,
and both sleep probabilities 1/2. Start with one active particle at each.
Let L_a,L_b be the terminal laws of a single active particle at the named
site, and F the two-particle law. Writing delta_uv for the sleeping-bit
state (a=u,b=v), direct one-particle renewal gives

    L_a = (1/2)delta_10+(1/2)L_b,
    L_b = (1/2)delta_01+(1/4)delta_00+(1/4)L_a
        = (2 delta_00+4 delta_01+delta_10)/7.

Process a first in the two-particle system. The complete-cycle branches
return to two active particles with fresh tails. The exact renewal is

    F = (1/4)delta_11+(1/8)delta_10+(3/8)F+(1/4)L_b.

Consequently, in bit order 00,01,10,11,

    F = (4,8,9,14)/35.                                        (B8)

The occupation covariance is -16/1225, so even this tiny feedback law
is not the product of its one-site marginals. The glued two-state
responses recover (B8) and every stable-background single-addition row.

### A limitation on boundary-only memory

There is no uniform finite number of hidden memory states depending
only on the number of ports, if the interface must represent arbitrary
sleeping backgrounds and exact output counts.

For a directed n-site path, use one inlet at its first site and one
counted outlet after its last site. Let 0<p_i<1. Start from any sleeping
support B of size k, and inject M=n+1 active particles at the inlet.
Processing the path in order necessarily visits every site. Its final
bits Z_i are independent Ber(p_i), and its outlet count is

    H=M+k-sum_i Z_i.

Choose one initial background for each k=0,...,n. Their outlet generating
polynomials are

    G_k(z)=z^(k+1) product_(i=1)^n [p_i+(1-p_i)z].              (B9)

These n+1 polynomials are linearly independent. Any exact r-state
stochastic transducer, even allowing each B to be encoded by a distribution
over its hidden states, expresses these fixed-input output laws in the
span of at most r state-conditioned laws. Hence r>=n+1. The claim assumes
no additional B-dependent parameter or externally retained mass register;
such a register would itself be memory.

This lower bound does not show that all 2^n sleeping supports are needed,
or identify a minimal realization on a fixed graph. It does show why a
fixed-dimensional linear response determined only by port number cannot
represent this exact interface.

### Checks and remaining questions

`verify_boundary_response.py` is the sole new
project file. Run from the worktree with

    python3 -B verify_boundary_response.py

It checks (B6)--(B8), all 30 nonzero elementary-response outcomes of the
two-site feedback motif against raw ARW, and all 171 nonzero outcomes
for two parenthesizations of a three-site cycle against raw ARW. The
three-site graph is a->b->c, with c->a/sink equally and all p_i=1/2.
The two queue chains have 97 and 105 transient states; the direct raw
chain has 98. Zero entries are included in the dictionary equalities.
It also checks input-order independence of the open two-site block and
the all-active two- and three-site laws against simultaneous raw starts.
These are fixed exact audits, not a parameter search.

The main object (B1) and gluing identities (B3)--(B5) are proved above;
the computations audit them on small motifs. The open structural question
is useful compression: when can two sleeping supports be identified
without changing their joint exit/future-response laws? Occupation
monotonicity, cone preservation, and any useful order on these response
matrices remain unproved. A category or stochastic-connection formalism
would require separate definitions and compatibility statements; the
associativity and coefficient noncommutativity here do not establish one.

A relation to the parameter connection of section 29 remains to be developed.


### Raw stacks and reset books: the precise pathwise dictionary

The iid reset processor itself is already in
[BS §2](evidence/notes/buffer-switching.txt), R §1, and RN §11 (17).
The present follow-up supplies an explicit parser and its quantifiers.
Fix the initial sleeping background b. At each site use independent
B_x(m)~Ber(p_x), m>=1, with B_x(0)=b_x, and independent arrows Y_x(k)~K(x,.).
The cumulative output after m processed arrivals is

    F_x(m)=m+b_x-B_x(m),
    F_x(m)-F_x(m-1)=1+B_x(m-1)-B_x(m) in {0,1,2}.              (B10)

If N_yx(k) counts y->x arrows in the first k arrows at y, the terminal
processed-input counts solve

    u_x=alpha_x+sum_y N_yx(F_y(u_y)),
    terminal bit at x=B_x(u_x).                               (B11)

Legal processing selects the least nonnegative integer solution. Indeed,
before a first attempted processing beyond any proposed solution v,
every other count is at most v_y. Monotonicity of all arrow-prefix counts
then bounds the available arrivals at x by v_x, a contradiction.
Frozen books therefore define an abelian processor network with tape
cursors; the two-state processor description averages over fresh tails.

To parse raw instruction stacks, retain a virtual bit c, initially b_x.
If c=1, consume sleeps through the next jump and append that jump's arrow
to Y_x. Then read one fresh instruction: a sleep sets the next B_x to one;
a jump sets it to zero and appends its arrow. When c=0, omit the compulsory
jump block. Replace c by the new bit and repeat. This parser depends on b
but not on alpha or the execution schedule. It consumes every raw jump
in its original order. The selected decision instructions are fresh;
their sleep indicators are iid, and the entire arrow sequence is independent
of the indicator sequence. Thus the parsed B and Y have the required joint law.

For legality, a processor boundary with q pending messages and virtual
bit c represents q+c active particles when q>0, or the stable bit c
when q=0. If c=1 and a message is processed, there are at least two active
particles, so the sleeps preceding the compulsory jump fail legally.
The decision instruction is then legal too. Its sleep succeeds only if
the physical population is one; otherwise c'=1 and the remaining pending
messages represent the unchanged active population. Self-jumps are delivered
back to the queue immediately in this accounting. Every processor block
is consequently a legal raw instruction block, including crowded decisions.

Conversely, independent geometric counts G_x(m), with
P(G_x(m)=g)=(1-p_x)p_x^g, restore the omitted sleeps before each compulsory
jump. The raw prefix length after m records is

    t_x(m)=m+sum_(j=1)^m B_x(j-1)[G_x(j)+1].                  (B12)

The reconstructed raw stacks are iid. On a common full-probability event,
for fixed b and every finite active input alpha, the two constructions
have identical terminal states, directed jump counts, and raw odometers
t_x(u_x). Termination follows from the sink-hitting-time potential; take
the countable intersection over alpha. The qualification "fixed b" matters:
on a one-site sink graph, raw prefix S,J_sink,J_sink gives different first
reset bits when parsed from b=0 and b=1. A common reset tape across changed
backgrounds is generally a different coupling from common raw stacks.

Terminal reset bits are not fresh iid coins. On one site with self-loop
probability r<1, a lone particle has

    P(B(N)=1)=p/[1-(1-p)r],                                  (B13)

because P(N=m,B(N)=1)=[(1-p)r]^(m-1)p and the corresponding final-zero
probability has last factor (1-p)(1-r). It equals 2/3 at p=r=1/2.
Conditioning on N=m leaves the same bias. Eliminating self-loops changes
the reset parameter to p'=p/[1-(1-p)r], as already proved in RN §21;
retaining self-messages instead keeps the original p. Fresh unused tails
do not make the selected terminal bit fresh.

`verify_boundary_response.py --reset-only` checks one fixed two-site raw-book
cylinder with failed sleeps and a self-loop, for three active inputs under
both opposite processor schedules and raw priorities. It verifies (B12)
and the exact rational case of (B13). This focused check passed separately
from the earlier gluing audits. No occupation order follows from (B10):
the final bit is sampled at the endogenous index u_x.


## 32. Dominator matroids and multivariate cancellation faces

**Derived statements, reviewed by main; one exact symbolic example,
2026-10-19.** Analytic assertions concern each fixed finite input.
The general occupation conjecture remains open.

### The sets that target-capacity IDLA fills surely form a laminar matroid

Let K be a finite sink-accessible jump kernel on V with absorbing sink s.
Use only positive-probability directed edges when discussing paths. For
d in V define

    B_d = {a in V : every a-to-s path visits d},
    m_d = sum_{a in B_d} eta(a).

In particular d belongs to B_d. Let H_C(eta) be the probability that
ordinary capacity-one IDLA, with parking allowed only at C, fills every
site of C from the finite all-active input eta. Process every particle,
including those that are killed before parking.

**Proposition.**

    H_C(eta)=1  iff  |C intersect B_d| <= m_d for every d in V.       (V1)

This is a statement about certainty, not a formula for H_C when it is
strictly between zero and one.

**Proof of necessity.** Suppose |C intersect B_d|>m_d. From every source
a outside B_d choose a finite positive-probability path to s avoiding d.
Such a path avoids all of B_d: if it visited b in B_d, its suffix would
be a b-to-s path avoiding d. With positive probability all particles
initially outside B_d follow these paths until parking or killing.
None parks in B_d. The m_d particles initially inside B_d can occupy
at most m_d of its requested sites. Thus C is not filled on this event.

**Proof of sufficiency.** Suppose c in C is left vacant. A particle
starting in B_c cannot leave B_c without visiting c: an exit to a vertex
outside B_c would furnish a continuation to s avoiding c. Since c is
vacant throughout IDLA, such a particle cannot visit it, either. Every
one of the m_c particles must therefore park at a distinct site of
(C intersect B_c) minus {c}. Sink accessibility ensures their walks
terminate almost surely. Hence m_c <= |C intersect B_c|-1, contradicting
the assumed inequality for d=c. This proves (V1).

The family {B_d} is laminar. Indeed, if B_c and B_d intersect at a, take
a simple a-to-s path containing c before d. Its a-to-c prefix avoids d.
Any c-to-s path avoiding d would contradict domination of a by d.
Therefore d dominates c and B_c is contained in B_d. The opposite order
gives the opposite inclusion.

Consequently

    I_eta = {C subset V : |C intersect B_d| <= m_d for all d}

is the family of independent sets of a laminar matroid. This is the
standard terminology in [Fife and Oxley, *Laminar Matroids*, Introduction](https://www.math.lsu.edu/~oxley/LaminarMatroids_revised.pdf).
For completeness, the exchange proof is short. If I,J satisfy the
constraints, |I|<|J|, and no j in J minus I can be added to I, each such
j lies in an I-tight B_d. Inclusion-maximal tight sets are disjoint.
On their union U, |J intersect U|<=|I intersect U|; outside U,
J is contained in I. This contradicts |I|<|J|.

Thus this matroid depends on the support of K and the input counts,
not on the positive values of the transition probabilities. Adding a
particle at x raises the capacities precisely along the dominators of x.

Two specializations connect it directly to RN Section 26. Under (127t),
B_d={d}; for a binary input A the independent sets are exactly the
subsets of A, and the circuits are the singletons outside A. On a
single-leak directed cycle, ordered just after the leak, B_j={0,...,j}.
The independent sets instead satisfy the prefix inequalities

    |C intersect {0,...,j}| <= |A intersect {0,...,j}|.

If A has an initial consecutive prefix of r roots and its next root is
b, the circuits of smallest size are exactly the (r+1)-subsets of
{0,...,b-1}. These are precisely the index sets occurring in the first
prefix leading polynomial in the old cycle handoff. This observation
concerns support and vanishing faces; it does not prove its coefficients.

### Source insertion, circuit roots, and trace operations

Write M_eta for this matroid and D(a)={d:a in B_d} for the dominator
chain of a. Give each source a exactly eta(a) distinguishable tokens,
each eligible to represent one vertex of D(a). Then C is independent
exactly when its vertices match injectively to eligible source tokens.
Necessity is the B_d capacity constraint. For sufficiency, the token
neighbors of T subset C are the union of B_c, c in T, counted with source
multiplicities. Inclusion-maximal such regions are disjoint, and their
capacity inequalities imply Hall's condition. Thus M_eta has a concrete
transversal presentation as well as its laminar presentation.

Adding one source is the matroid union with a rank-one matroid supported
on its dominator chain:

    M_(eta+delta_x)=M_eta vee U_(1,D(x)),
    I_(eta+delta_x)=I_eta union {I union {d}: I in I_eta,d in D(x)}. (V1a)

This follows by assigning or removing the new token in a matching. Every
restricted rank rises by at most one, but this is not generally global
rank truncation. Let J(M) be the squarefree circuit ideal, extending
deletion ideals to the common polynomial ring. Then

    J(M_(eta+delta_x))=intersection_(d in D(x)) J(M_eta\d).      (V1b)

A set S remains dependent exactly when S without {d} was dependent for
every eligible d. This proves (V1b), including d not in S. On the path
0->1->2->sink with eta=delta_0, adding at 1 changes U_(1,3) to U_(2,3):
the three pair generators become the single generator p0 p1 p2.
New circuits can therefore appear as old ones are destroyed.

Every circuit C has a distinguished root d. Choose an inclusion-minimal
violated B_d. Minimal dependence gives C subset B_d and |C|=eta(B_d)+1.
Every immediate dominator-child region B_e is unviolated, and these
regions partition B_d without {d}. Hence

    eta(B_d)+1=|C| <= 1+sum_children eta(B_e)
                       =1+eta(B_d)-eta(d).

Equality is forced, proving

    d in C, eta(d)=0,
    |C intersect B_e|=eta(B_e) for each child e.              (V1c)

Conversely, choose in each child region an independent set I_e of size
eta(B_e), if such choices exist. Then {d} union union_e I_e is a circuit
when eta(d)=0: deleting any one element repairs the B_d constraint and
all ancestor constraints; all child constraints already hold. Thus d is
the unique vertex of C nearest the sink in the dominator tree. Inserting
at x destroys this circuit exactly when x in B_d. If x is inside, all
violated ancestor capacities rise by one; if outside, B_d still fails.

Tracing K onto C and the sink gives exactly one leaking vertex, d.
Every other c in C must visit d before the sink, while a simple d-to-sink
path avoids all such c. Adjoining a target v yields the exact classification

| Dominator relation | Directly leaking vertices in the trace onto C union {v} |
|---|---|
| v in B_d | d |
| v not in B_d and d in B_v | v |
| d and v incomparable | d and v |

For incomparable d,v, simplify a d-to-sink path avoiding v; it also avoids
C without {d}. A v-to-sink path avoiding d cannot visit any other C-site.
Both leaks are positive. These one/two-leak kernels describe the minimal
circuit-supported faces below, without supplying their coefficient signs.
For arbitrary retained sets, collective cutsets still matter.

Two further operations have precise IDLA meanings. If eta is supported
on W, tracing excursions outside W gives

    M(K^W,eta)=M(K,eta)|W.                                    (V1d)

If sources lie outside W, the independent sets on the right are instead
the intersection over every positive-probability entrance configuration
zeta of I(M(K^W,zeta)). One cannot replace that random entrance law by its
mean. For supplied roots R, contraction gives

    M_(eta+1_R)/R=M_eta\R.                                   (V1e)

Indeed |R intersect B_d| cancels from both sides of each capacity constraint.
Physically, the supplied particles can park immediately at R. Contracting
an arbitrary independent R in M_eta need not correspond to spending a
fixed selection of the original sources: with a->r, b->r, r->sink and
eta=delta_a+delta_b, contraction gives U_(1,2) on {a,b}, whereas spending
either source to fill r makes its own source site a loop in the remaining
sure-filling matroid.

### Circuits give an exact analytic support restriction for ARW

Write D_v(p)=q_v(eta+delta_x;p)-q_v(eta;p), and let C_eta be the circuits
(minimal dependent sets) of the matroid above. RN Section 26 (127u) gives
analyticity at zero and the trace-IDLA reduction. On the parameter face
p_(V minus S)=0, with S in I_eta, the target-capacity exploration fills
all of S surely. Its fresh continuation is the full active law on the
trace graph S. That law is stationary under any further active addition.
Hence the *entire output-law difference* is zero on this face.

For a Taylor coefficient [p^m]D_v to survive, its support must therefore
be dependent and contain v. Equivalently, in the ring of analytic germs,

    D_v belongs to < product_{i in C union {v}} p_i : C in C_eta >. (V2)

The proof is just coefficient restriction to each coordinate face,
together with divisibility by p_v. A dependent support contains a
circuit; grouping the convergent series by one of the finitely many
resulting monomial divisors gives (V2).

In particular, for fixed positive alpha_i and p_i=epsilon^{alpha_i},

    D_v = O(epsilon^b),
    b = min_{C in C_eta} sum_{i in C union {v}} alpha_i.           (V3)

If there are no circuits, D_v is identically zero. Neither equality in
(V3) nor the sign of its first surviving coefficient is asserted. Scalar
target projection can cancel further terms. Also (V1) identifies exactly
the *sure-saturation* faces, not every face on which a specified scalar
response happens to vanish.

### A complete anisotropic calculation on one three-cycle

Take 0->1->2 deterministically and 2->0 or s with probability 1/2 each.
The lower input is A={0,2}, the added particle starts at 1, and v=2.
All inputs are active. Put

    P=p0 p1 p2,  q_i=1-p_i,  s0=p0+p1,  t0=p0 p1,
    H=1-(q0 q1 q2)/2,

        [ 1          -q0           -q0 p1 ]
    M = [ -q1 p2/2    1            -q1    ],
        [ -q2/2      -q2 p0/2       1     ]

    B=s0-s0^2+(1+p2)t0+(4-p2)s0 t0-(3-p2)t0^2.

An exact rational identity is

    D_2 = P q2 B / [4 H (1+P) det(M)].                          (V4)

Here det(M)(0)=H(0)=1/2, so the denominator in (V4) equals one at zero.
For clarity,

    det(M)=1-(q0 q1 p2+q0 q2 p1+q1 q2 p0+q0 q1 q2)/2
             -P q0 q1 q2/4.

**Derivation of (V4).** Let beta_i and gamma_i be the full-input and
lower-input probabilities of the exact terminal pair V minus {i}.
The killed two-site normalizer is Z_i=1/(p_j p_k), and the full one is
Z_V=(1+P)/(2P). Specializing RN (115), and writing L=(1+P)/2 and
t=(q0 p1 p2,q1 p2,q2)^T, gives

    M diag(p)(beta-gamma) = P[(t/2)/L-e_1].                     (V5)

The forced no-return probabilities from target 2 are particularly
elementary on the open path 0->1->2:

    w({2})     = 1/2 + (p0+q0 p1)/2,
    w({0,2})   = 1/2 + p0 p1/2,
    w({1,2})   = 1/2 + p0/2,
    w(V)       = 1/2.

For example, after forcing 2->0 with a resident at 0, avoidance of 2
requires two particles to sleep at 0 and 1, with probability p0 p1.
With residents at both 0 and 1, three particles cannot be absorbed on
the two-site open path, so a return to 2 is forced.

The occupation identity (127) now becomes

    D_2 = q2/(2H) [
        (p0+p1-2p0 p1)(beta_1-gamma_1)
        +q0 p1(beta_0-gamma_0)
        +(p0+p1-p0 p1)P/L ].                                  (V6)

Solving the three-by-three system (V5) and substituting in (V6) yields
(V4). The new verifier constructs its adjugate and verifies the polynomial
identity by exact rational arithmetic. It also checks (V4) independently
against the raw absorption equations, as detailed below.

Since B=s0+O(||p|| s0), (V4) proves the uniform multivariate asymptotic

    D_2=P[(p0+p1)+O(||p||(p0+p1))].                            (V7)

The norm can be the maximum norm and the error bound holds throughout
a sufficiently small nonnegative cube. It does not require comparability
of the coordinates. This stronger remainder is essential: an undivided
total-degree error O(||p||^5) would not prove the following assertion.

For *every fixed positive* alpha=(alpha0,alpha1,alpha2), setting
p_i=epsilon^{alpha_i} gives

    D_2 ~ a(alpha) epsilon^{alpha0+alpha1+alpha2+min(alpha0,alpha1)},
    a(alpha)=1 if alpha0!=alpha1, and 2 if alpha0=alpha1.         (V8)

With additional fixed multipliers lambda_i>0, multiply the leading
coefficient by lambda0 lambda1 lambda2 and by lambda0, lambda1, or
lambda0+lambda1 according to which of alpha0,alpha1 is smaller or equal.

The Newton polyhedron (convex hull of Taylor exponents plus the
nonnegative orthant) is exactly

    conv{(2,1,1),(1,2,1)} + R_{≥0}^3.                          (V9)

Indeed the exact numerator has precisely those componentwise-minimal
exponents, both with positive coefficient, and its denominator is an
analytic unit. There is one transition wall, alpha0=alpha1.

#### Where the cancellations occur

The dominator regions are B_0={0}, B_1={0,1}, B_2=V, with input capacities
1,1,2. The only circuit is {0,1}. Formula (V2) thus forces the factor
p0 p1 p2 for target 2. It does not force the extra p0 or p1 in (V7).

At the squarefree cubic, the output-law difference restricted to supports
containing 2 is

    [p0 p1 p2](mu_upper-mu_lower)|_{2 occupied}
       = 2 delta_{012} - delta_{02} - delta_{2}.                (V10)

Its target marginal is 2-1-1=0. Thus the matroid detects the required
support, but the first marginal has a second power at 0 or 1. These
are the two exponents in (V9). In (V6) that extra power comes from one
sleep in an auxiliary forced-return avoidance probability. This is
not a proved positive pairing of complete upper and lower histories.
There is no positive coefficient theorem here: in fact

    [p0^3 p1 p2]D_2 = -2.                                     (V11)

**Exact verification.** The standard-library script
`verify_sleep_valuations.py` passed:

- the adjugate identities and the compact polynomial numerator identity;
- all 56 occupation Taylor coefficients of total degree at most five,
  from a 28-transient-state raw ARW absorption system;
- the four cubic exact-support coefficients in (V10);
- one parameter point p=(1/3,2/5,3/7), under both opposite raw toppling
  priorities, with gain 12670/1398193.

These are one symbolic example and one fixed rational check, not a
parameter search. The leading homogeneous polynomial in (V7) agrees
with the old fixed-ray handoff. The additional content is the exact
support restriction, the anisotropic conclusion, and the explicit
cancellation interpretation. This does not audit the arbitrary-length
cycle or pendant extensions in the old handoff.

### Minimal sleep histories: a positive statement with a precise boundary

Fix a legal deterministic state-dependent scheduler and a terminal event E for one
finite input. Discard ineffective sleep instructions at crowded sites.
For each admissible complete history ending in E, count its successful
sleeps by site, obtaining n in N^V. Let H_E be the componentwise-minimal
vectors so obtained. There are finitely many: erasing a repeated raw
state from a history removes a cycle, preserves admissibility and its
terminal event, and can only decrease every sleep count. Simple paths
in the finite raw state graph suffice to witness all minimal vectors.

For the nonnegative probability f_E(p), its Taylor support has exactly
the same componentwise-minimal exponents H_E, and their coefficients
are positive. To see this, group histories by their sequence of sleep
transitions. The jump-only fundamental matrix is analytic at zero
and has nonnegative entries there; its transience follows from the
finite expected sink-hitting-time potential. A sleep sequence with
count vector n contributes p^n times an analytic function with positive
constant term when the sequence is feasible. Every further Taylor
monomial is componentwise at least n. The expansion over sleep sequences
converges near zero, since each additional sleep contributes O(||p||)
between bounded jump resolvents. At a minimal n no sequence with a
strictly smaller count vector can contribute a cancelling correction.

Thus for alpha_i>0, the asymptotic order of f_E(epsilon^alpha) is

    min_{n in H_E} alpha dot n,

with positive leading coefficient, and its Newton polyhedron is
conv(H_E)+R_{≥0}^V. This statement is only for one nonnegative event
probability. It does not transfer to D_v, a difference of two such
probabilities. Equations (V10)--(V11) show why that distinction matters.

### Aggregate guarantees are different from sure-filling rank

Let g_eta(C) be the smallest cardinality of a positive-probability terminal
IDLA support, with parking allowed only at C. For repeated copies of one
source x, put ell_C(x)=min_gamma |C intersect gamma|, minimizing over simple
x-to-sink paths. Then

    g_(k delta_x)(C)=min{k,ell_C(x)}.                          (V12)

If no particle is lost, k sites fill. If a particle is lost, its path
only encounters already filled C-sites; erasing loops shows that the
terminal support has size at least ell_C(x). Conversely all particles
can follow one minimizing simple path, filling its C-sites in order.
This positive-probability event attains the bound.

For example x branches equally to a1->a2->sink or b1->b2->sink, and
C={a1,a2,b1,b2}. Here D(x)={x}, so the full sure-filling matroid is the
same for one, two, or three inputs at x, and its rank on C is zero.
The minimum occupancies are nevertheless 1,2,2. With two inputs the exact
supports {a1,a2}, {a1,b1}, {b1,b2} have probabilities 1/4,1/2,1/4.

For general eta there is an exact recursion retaining the full support
S_eta(C) of the lower IDLA law. Let cut_x(U)=1 if every x-to-sink path
visits U, and zero otherwise. Then

    g_(eta+delta_x)(C)
       =min_(B in S_eta(C)) [|B|+cut_x(C without B)].          (V13)

Conditional on B, the added walk parks surely exactly when its vacant
capacity set is such a cutset; otherwise it can follow a finite avoiding
path and be lost. This needs more than the sure-filling matroid.

The added `matroid_neighborhood_audit` in `verify_sleep_valuations.py`
checks three fixed graphs: the three-site path, a cyclic four-site fork,
and the two-branch graph above. It verifies source union and the circuit
ideal update, circuit roots and source destruction, restriction/contraction,
and the exact branch laws. The fork's circuit trace has its sole leak 3/8.
The focused audit passed in the agent and main runs; the preceding proofs
establish the general statements.

### What this suggests, and the gaps

There are two separate combinatorial layers. Dominator-region capacities
give a laminar matroid and determine all sure-saturation faces exactly.
Its circuits force monomial support constraints for the signed response.
Marginal projection can then kill the whole first vector on a circuit
support, and extra successful sleeps can be needed before any response
is visible.

A useful next question is whether one can attach to each circuit the
minimal sleep events that prevent the released particles from returning
to v, and recover the remaining lower Newton faces from these data.
The cycle suggests this, but no general construction, completeness proof,
or sign theorem for those faces has been established here. Actual scalar
zero faces can be larger than the matroid faces, and cancellations can
occur again within a weighted face. Positive leading terms on individual
fixed rays still do not imply a uniform heterogeneous-cube theorem.
Only the particular three-site identity (V4) supplies that stronger
control in this note. General monotonicity remains open.


### Possible filling, exact supports, and reachability matchings

**Proved combinatorial characterizations, with one fixed exact audit.**
The sure-filling matroid uses dominators. Its possible-filling counterpart
uses ordinary reachability, including length-zero paths. Distinguish all
N source tokens of eta. A set C is filled with positive probability by
IDLA parking only on C iff C can be matched injectively to source tokens,
each token eligible for every vertex reachable from its source. These
sets form a transversal matroid, denoted M_eta^poss. Restriction of a
matching proves heredity; an alternating path between matchings proves
exchange. This does not describe exact terminal supports when parking
is allowed on a larger set.

The exact-support statement includes the missing constraint. For B subset C,
delete U=C\B from the directed graph entirely. A source in U then has
no reachability neighbors, even at time zero. The output B has positive
probability iff there is a perfect matching from all N source tokens to
one slot at each b in B and N-|B| distinct sink slots, using reachability
in G-U. The criterion is false if |B|>N. Equivalently there is one
matching between tokens and B that saturates B and every token unable
to reach the sink in G-U. No sink-accessibility assumption on G-U is made.

Necessity follows from the realized paths: a never-occupied U cannot be
visited, and distinct parked tokens occupy distinct B-slots. For
sufficiency, choose finite paths realizing a matching and process the
tokens in any fixed order. Follow the next assigned path to its first
vacant B-site or the sink. If it parks at c assigned to another token j,
give j the path consisting of its old path to c followed by the suffix
of the current path to its old endpoint. This endpoint may be a sink
slot. Remove the parked token and c. All remaining assignments are
again valid, with distinct vacant B-endpoints. A token reaching the
sink simply removes its sink slot. Concatenated paths may have loops;
they are finite allowed walks and can be loop-erased. The resulting
finite legal history has positive probability and ends exactly at B.
Putting B=C and assigning excess tokens to sink slots proves the
possible-filling criterion above.

This also gives a formula for the minimum possible occupancy g_eta(C).
For B subset C let T_(C\B) be the tokens unable to reach the sink after
deleting C\B. Then

    g_eta(C)=min{|B|: T_(C\B) can be matched into B
                         by reachability avoiding C\B}.    (V8)

To prove sufficiency for the contained-support event, assign the matched
tokens to B and all others to sink paths. If a token first hits an
unassigned vacant B-site, park it there and discard its old assignment;
otherwise use the same swap. The output is contained in B. Necessity
holds on every such output. Minimizing gives (V8), with no efficient
optimization algorithm asserted.

There is a useful comparison with the sure-filling matroid M_eta:

    rank(M_eta|C)<=g_eta(C),
    g_eta(C)=|C|-1 when C is a circuit of M_eta.              (V9)

For I subset C, use identical site-arrow stacks for the parking domains
I and C. Every legal C-toppling sequence is legal for I from the same
initial mass, since the I-capacities are smaller. Complete C first and
then continue to stabilize I. Abelianity gives u_C<=u_I, hence the
number of used sink arrows satisfies L_C<=L_I and the occupied counts
satisfy |occupied_C|>=|occupied_I|. Choose a sure-filling basis I of C.
For a circuit its rank is |C|-1, while dependence gives a deficient
output of positive probability. The two bounds give equality. This is
count order; containment of the two coupled outputs is not asserted.

For one source x, all possible supports over all particle counts form
the directed vertex-search antimatroid

    F_x={empty} union {B:x in B and every b in B is reachable
                                      from x inside G[B]}.

Its ground set is the vertices reachable from x. Union closure follows
by retaining witness paths; accessibility follows by removing a leaf
of an outward directed spanning tree. For k>=1, the exact outputs of
k delta_x, with parking on V, are precisely the nonempty B in F_x with
|B|<=k and either |B|=k or a path from x to the sink inside B plus the
sink. To realize B, first fill an outward spanning tree in parent-first
order and then send any surplus tokens to the sink inside B. Thus
union_(k>=0) supp(IDLA_(k delta_x))=F_x, but a fixed-k support family
need not itself be an antimatroid.

For full filling with exactly m=|V| tokens, this matching criterion
identifies the support of the polynomial Q in (A7a) explicitly:

    supp Q = supp product_(v in V)(sum_(a:a reaches v)b_a)
      = {eta>=0 integer: eta(V)=m,
                    eta(S)<=|Reach(S)| for every S subset V}. (V10)

These are Hall's inequalities on the source-token side. The coverage
function |Reach(S)| is submodular, so its base polytope describes these
integer inputs. Equivalently the Newton polytope is the Minkowski sum
of the coordinate simplices on the ancestor sets of each target.
For a direct symmetric exchange proof, take assignments realizing
eta,xi and draw a directed edge from the first assigned source to the
second for each target. Its divergence is eta-xi. A surplus source a
has a directed path to some deficit source b. Swap assignments on this
path to realize both eta-e_a+e_b and xi+e_a-e_b. Thus the M-convex
support conclusion also follows without any volume theorem. Equality
of supports is not equality of coefficients: the product counts
assignments, while Q weights IDLA outcomes.

`possible_support_audit()` in `verify_sleep_valuations.py` checks all
27 pairs B subset C on x->a, a->b or sink equally, b->x, with two tokens
at x. For C={a,b}, the outputs {a} and {a,b} have probabilities 1/2
each, explicitly requiring the sink-assignment swap. The audit checks
full filling, exact supports, contained supports, the minimum formula,
and sources in the forbidden set. Main reran this focused check; no
new parameter search was used.


## 33. Harmonic current tilts and completion of independent paths

**Proved identities, with an exact two-port certificate and independent
AI audit of (X1)--(X12), 2026-10-19.**
This uses the boundary-response object from Section 31. All exits are
absorbing and labeled. Splitting the original sink into several labels
changes the recorded currents but not the physical stabilization law.
No occupation comparison or literature-novelty claim is made.

### A harmonic change of routing leaves the sleep parameters unchanged

For positive exit weights z=(z_o:o in O), let

    h_i(z)=E_i[z_(exit)]
          =sum_{j in D} K(i,j)h_j(z)+sum_{o in O}K(i,o)z_o,
    chi_z(B)=product_{i in B} h_i(z),
    H_z=diag(chi_z(B): B subset D).

The expectation is for one ordinary killed random walk. Exit accessibility
gives a unique strictly positive h. Define a new jump kernel, with the
same support and the same sleep probabilities, by

    K^z(i,j)=K(i,j)h_j/h_i,    j in D,
    K^z(i,o)=K(i,o)z_o/h_i,    o in O.                         (X1)

Harmonicity says each row sums to one. For a stable B, one active
addition at i, terminal support C, and labeled exit counts k, the exact
change of measure is

    P^(K^z)_B(C,k | i)
       =P^K_B(C,k | i) z^k chi_z(C)/(h_i chi_z(B)).            (X2)

**Proof.** Extend the weight to an intermediate raw configuration by
product_x h_x^(number of particles at x), irrespective of activity,
and multiply by z to the accumulated exit counts. A sleep instruction
does not change this weight. Each jump multiplies it by h_j/h_i or
z_o/h_i. The ratios of the raw instruction probabilities under K^z
and K therefore telescope to the ratio in (X2). This also covers
ineffective sleeps and self-loops. Both kernels are sink-accessible,
so complete finite-input stabilization terminates almost surely.
Summing the nonnegative probabilities of complete histories proves
(X2), with no optional-stopping integrability gap.

Consequently the boundary matrices and binary all-active stabilization
matrix satisfy

    A_i^K(z) chi_z=h_i chi_z,
    H_z^(-1) A_i^K(z) H_z/h_i=A_i^(K^z)(1),
    P_p^K(z) chi_z=chi_z,
    H_z^(-1) P_p^K(z) H_z=P_p^(K^z)(1).                       (X3)

For P, the denominator is chi_z(B), since its input has one active
particle at each site of B. More generally an initial particle-count
vector eta has initial weight product_i h_i^eta(i). Formula (X2)
holds whether its particles are active or sleeping, whenever that
initial configuration is defined.

The same tilt extends the parameter transport of Section 29 to joint
exit counts. Write E_r^K(z)=P_r^K(z)^(-1). For r<=s, the physical operation
of independently waking current residents with probabilities 1-r_i/s_i
and stabilizing at s has generating matrix

    T^K(r,s;z)=E_r^K(z)P_s^K(z)
             =H_z T^(K^z)(r,s;1)H_z^(-1).                    (X3a)

To prove this, apply (X2) to the waking operation, then use the already
proved stochastic transport identity for K^z. The diagonal H_z does not
depend on r. Thus (X3a) retains nonnegative polynomial coefficients, not
just positive values at z>0. Its cocycle identity also records cumulative
labeled exits. At zero odds use the same limiting convention as Section 29.

For iid addition locations with fixed distribution rho, put
A_rho(z)=sum_i rho_i A_i(z). Its strictly positive eigenvector is chi_z
and its spectral radius is

    lambda_rho(z)=sum_i rho_i h_i(z).                         (X4)

Indeed conjugating by H_z and dividing by lambda_rho gives a stochastic
matrix, namely the driven K^z chain with new input distribution
rho_i^z=rho_i h_i/lambda_rho. This argument requires neither irreducibility
of the stable-state chain nor a simple dominant eigenvalue. The eigenvalue
in (X4) is independent of the sleep probabilities; the right eigenvector
is also independent of them. The left stationary eigenvector generally
depends on the sleep probabilities.

### Completed particle paths and current fluctuations

There is a pathwise explanation with a stronger uniform-in-time statement.
Start empty and label each added particle. Give it an independent complete
ordinary K-walk, stopping at its labeled exit. At a legal toppling choose
the moved particle by any rule using revealed information, then use its
next walk step. Unrevealed walk steps have the required conditional
K-law, so this realizes the usual ARW law. Independent sleep instructions
retain their usual probabilities. Sleeping interrupts paths, and waking
resumes them; labels are not copied or removed except at an exit.

After m inputs and complete stabilization, let H_m be the actual exit
count vector and Y_m the count of the preassigned endpoints of all m
particles. If R_m counts the endpoints of the particles still resident,
then simultaneously for all m,

    Y_m=H_m+R_m,    R_m>=0,    sum_o R_m(o)<=|D|.              (X5)

For a deterministic source word, the summands in Y_m are independent
with their ordinary source-dependent exit laws. For iid source law rho,
they are iid categorical with probabilities

    pi_o=sum_i rho_i P_i(ordinary walk exits at o).

Thus the interacting process differs from these independent exit counts
by at most the domain's storage capacity. Conditional on the revealed
stabilization history, the unused paths of the terminal residents are
independent fresh walks; their completion generating function is chi_z(C).
This also gives (X3) directly.

In particular, for iid driving,

    H_m/m -> pi almost surely,
    (H_m-m pi)/sqrt(m) => N(0, diag(pi)-pi pi^T).              (X6)

These follow from the corresponding elementary iid limit theorems and
the deterministic bound in (X5). The covariance may be singular, including
the compulsory zero fluctuation in the total-current direction. No
independence of the actual successive exit batches is asserted.

For any fixed positive z, set a=min(1,min_o z_o) and
b=max(1,max_o z_o). From a^|D|<=z^R_m<=b^|D|,

    b^(-|D|) E[z^Y_m] <= E[z^H_m] <= a^(-|D|) E[z^Y_m].

Consequently the current log moment-generating function per addition is

    lim_(m->infinity) m^(-1) log E[exp(theta dot H_m)]
       =log(sum_o pi_o exp(theta_o)),   theta in R^O.          (X7)

It is independent of p on every finite graph. With one unlabeled sink
these identities reduce largely to particle conservation; separate exit
labels retain the routing information. The terminal resident law and
finite-time currents still depend on p.

**Exact check.** `verify_exit_tilt.py` uses an asymmetric two-vertex graph,
two labeled exits, one self-loop, and p=(1/3,2/5). At exit weights
z_(-1)=2, z_(-2)=1/3, h=(13/9,7/6). Independently solving the raw instruction
chains under K and K^z verifies (X2) for all 64 nonzero joint entries of
the one-addition response, and 29 entries for all binary active inputs
plus eta=(2,1). A second exit-weight point checks the completion identity.
The driven eigenvalue for rho=(2/3,1/3) is 73/54. The certificate checks
the finite identities; the proofs above give (X5)--(X7) in general.

### The same transform counts internal jumps

Assign a positive weight w_(i,j) to each directed jump, including jumps
to labeled exits. Let W(i,j)=K(i,j)w_(i,j) for i,j in D, and
k_w(i)=sum_o K(i,o)w_(i,o). Assume the spectral radius of W is less
than one; this includes a neighborhood of all weights one. Then

    h=(I-W)^(-1)k_w
      =E_i[product of jump weights until ordinary-walk exit]

is finite and strictly positive. Replacing (X1) by

    K^w(i,j)=K(i,j)w_(i,j)h_j/h_i,
    K^w(i,o)=K(i,o)w_(i,o)/h_i                              (X8)

again gives a stochastic sink-accessible kernel. The same raw-history
likelihood proof gives (X2)--(X4), now weighting every jump. Unlike exit
counts, internal jump counts are unbounded; their generating entries
are convergent series on this domain, rather than polynomials. Their
finiteness follows directly by summing the transformed path probabilities.

Let J_m be the vector of actual directed jump counts after m iid inputs,
and let Z_1,...,Z_m be the full jump-count vectors of their independent
complete walks. Path completion gives

    sum_(a=1)^m Z_a=J_m+R_m^jump,    R_m^jump>=0.              (X9)

The correction is no longer deterministically bounded. Conditional on
the revealed history, however, it is the total future traffic of at most
|D| independent ordinary walks. On a fixed finite sink-accessible graph,
their lifetimes have a uniform exponential tail: a fixed positive chance
of exiting in a fixed number of steps gives a geometric block bound.
Thus R_m^jump has a uniform exponential tail, independent of m. The
strong law and multivariate central limit theorem for J_m are therefore
those of the iid vectors Z_a. For the strong law, apply the exponential
tail bound and Borel--Cantelli to R_m^jump/m; for the central limit theorem,
R_m^jump/sqrt(m) tends to zero in probability.

There is also an exact log moment-generating limit throughout the stated
domain. Let chi_w(C)=product_(i in C)h_i, and let chi_min,chi_max be its
minimum and maximum over the finite set of stable supports. Completing
the resident walks conditionally gives

    E[w^J_m chi_w(C_m)]=(sum_i rho_i h_i)^m.

Bounding chi_w between its two positive extrema proves

    lim_m m^(-1) log E[w^J_m]=log(sum_i rho_i h_i).            (X10)

Hence the limiting jump-current generating function depends on the ordinary
walk and input law, not on the sleep probabilities. These are counts per
added particle after complete stabilization, not activity per unit physical
time, raw failed-sleep counts, or avalanche-size independence.

This leads to the neighboring literature on
[geometric universality of currents, Chernyak--Chertkov--Sinitsyn, §§I,III.A](https://arxiv.org/pdf/1008.2749):
fixed routing probabilities can make long-run currents insensitive to
interactions that only change when particles move. Their setting uses
continuous-time open networks and conditions excluding accumulation under
the current tilt. We do not import that theorem into ARW. Here observations
are taken after full stabilization, and (X5), (X9), and the finite product
chi_w explicitly control the resident correction. The path geometry versus
holding-time distinction is the common mechanism; this is a connection to
an existing theory, not a new universality claim.

In the same two-vertex certificate, weight only the edge 0->1 by 3/2.
Then h=(7/5,6/5) and the driven eigenvalue is 4/3. A separate 12-state raw
Feynman--Kac absorption system checks all 30 nonzero terminal entries
against (X8). No jump-count cutoff is used.

### The full tilted spectrum comes from ordinary IDLA

Let D_i(w) be one IDLA addition, retaining the same jump weights, and
P_p(w) binary all-active ARW stabilization. Surplus IDLA and the Abelian
property preserve every directed jump count, so the ordinary intertwining
extends to

    P_p(w) A_i^p(w)=D_i(w)P_p(w).                             (X11)

There is no formal-inverse gap: (X8) gives
P_p(w)=H_w P_p^(K^w)(1)H_w^(-1), and the unweighted P on the right is
invertible by the complete-block theorem. Thus all driven matrices
A_rho^p(w)=sum_i rho_i A_i^p(w), at fixed K,rho,w, are similar to the
single matrix D_rho(w), for every heterogeneous sleep vector p.

IDLA can only enlarge its occupied support. In an ordering by support
size D_rho(w) is triangular. Its diagonal at B is the weighted probability
that an input walk is lost before finding a vacancy outside B:

    lambda_B(w)=sum_(i in B) rho_i
                      [(I-W_BB)^(-1)k_w|B]_i,    B subset D. (X12)

For B empty the value is zero; for B=D it is sum_i rho_i h_i.
These 2^|D| numbers, with multiplicity, are exactly the spectrum of
A_rho^p(w). In particular the entire spectrum is real and nonnegative
on the positive-weight convergence domain. Jordan blocks are permitted;
triangularity is not a diagonalizability claim. The simultaneous addition
algebra, not just its dominant eigenvalue, is independent of p up to
similarity. The conjugating matrices and resident eigenvectors depend on p.

At w=1 and when every site is accessible from the input distribution,
the largest eigenvalue below one is

    1-min_v P_rho(ordinary walk visits v before the sink).

Indeed every proper B omits some v, and lambda_B<=lambda_(D without {v});
the latter is the probability of avoiding v. This recovers the relaxation
scale already proved by
[Bristiel--Salez, Theorem 1](https://arxiv.org/pdf/2209.03274).
The underlying ARW/IDLA sampling relation is from
[Levine--Liang, §1.4](https://arxiv.org/html/2110.14008v3#S1.SS4).
Here (X11)--(X12) follow from the saved inverse and preservation of jump
counts; they are not presented as a literature-novelty claim.

The same two-site certificate verifies the full characteristic polynomials
by exact traces through degree four. The exit-weight spectrum is
{0,31/54,4/27,73/54}; the internal-edge-weight spectrum is
{0,1/3,1/6,4/3}, as predicted by (X12). This does not yield the occupation
sign: the spectral expansion of a target observable still has signed
coefficients, and the similarity is not a positive order isomorphism.


## 34. Addition algebras, sink-search antimatroids, and deformed cubes

**Derived arguments, 2026-10-19; independent AI audit passed for the
local algebras, polytope bridge, normalization, Jordan formulas, and
inverse-basis consequences (A10)--(A14).**
The construction below applies to every finite sink-accessible kernel,
without reversibility. The algebraic and polytope claims have proofs here;
the Lefschetz input is an identified external theorem. These are adjacent
structures, not an occupation-monotonicity proof. No novelty claim is made.

### The simultaneous addition algebra and its characters

Write k_i=K(i,s), M=I-K_V, and let D_i be ordinary IDLA addition. The
commuting D_i satisfy D_i^2=D_i(k_i I+sum_j K(i,j)D_j): when two particles
start at i, one surplus particle must jump. Define

    A_K=R[z_i:i in V]/<z_i(z_i-k_i-sum_j K(i,j)z_j)>.          (A1)

The squarefree monomials z_B form a basis, and the IDLA representation is
the regular representation. For spanning, apply the displayed relation
at a crowded site of any monomial. For each fixed initial degree there
are finitely many count configurations. The resulting surplus-IDLA
transient system has inverse I+Q+Q^2+...: the ordinary sink-hitting-time
potential bounds its expected number of jumps. Thus its finite linear
equations express every monomial as a combination of stable squarefree
ones. For independence, applying a proposed relation to the empty-state
row gives a relation between distinct rows delta_B, since
delta_empty product_(i in B)D_i=delta_B. This also identifies the regular
representation. By (X11), the ARW additions give the same algebra up to
similarity, at every allowed sleep vector.

For B subset V, let core(B) be the vertices in B having a path to the
sink whose internal vertices all lie in B. Call R feasible if core(R)=R.
These feasible sets form the vertex-search antimatroid of the reversed
routing graph, rooted at the sink. They contain the empty set and are
union-closed. Accessibility follows by choosing a sink-directed spanning
forest on R and removing a leaf. Equivalently,

    cl(C)=V\core(V\C)                                       (A2)

is a convex-geometry closure. For anti-exchange, if distinct a,b survive
deletion of C and deleting b cuts every a-to-sink path, deleting a cannot
also cut every b-to-sink path: take the suffix after b of a simple
a-to-sink path. In particular cl({d}) is exactly the dominator region
B_d of Section 32. Collective separators give the larger closures.
For the standard graph-search terminology see
[Nakamura--Okamoto, antimatroids from graph searches](https://doi.org/10.1016/S0166-218X(02)00581-4);
the arguments needed here are given above.

Every complex character of (A1) is real and corresponds to one feasible R:

    h^R_R=(I-K_RR)^(-1)k_R,    h^R_(V\R)=0.                 (A3)

Indeed, on the nonzero-coordinate support R, the character equations
reduce to the invertible linear system in (A3). Its inverse is
nonnegative, and every coordinate is positive exactly when R is feasible.
For this R put

    F_R={x outside R:k_x+sum_(y in R)K(x,y)>0},
    Y_R=V\(R union F_R),    m_R=|Y_R|.

The sets with core(B)=R are precisely R subset B subset R union Y_R.
One cannot include a frontier vertex; conversely a path from Y_R to
R or the sink must first pass through the excluded frontier. Therefore
the joint character (A3) occurs with multiplicity 2^m_R in the triangular
IDLA representation. This includes coincidences among scalar driven
eigenvalues without identifying their joint characters.

More precisely, the local algebra at this character is

    H_R=R[u_y:y in Y_R]/<u_y(u_y-sum_(z in Y_R)K(y,z)u_z)>.
                                                                  (A4)

To check this, in the local algebra the z_i with i in R are units, so
z_R=h^R_R+(I-K_RR)^(-1)K_R,Y_R z_Y_R after eliminating the frontier.
For f in F_R the other factor in its quadratic relation has a nonzero
constant term, hence z_f=0. A vertex of Y_R has no direct jump to R or
the sink, leaving exactly (A4). Conversely these substitutions satisfy
all equations. The algebra in (A4) is finite and supported only at zero,
so its positive-degree elements are nilpotent and the required units
are invertible. Hence A_K is the product of these local algebras.

For any transient substochastic kernel L on m vertices, write

    H_L=R[u_1,...,u_m]/<u_i((I-L)u)_i>.

It has squarefree basis and Hilbert series (1+t)^m. A direct proof uses
IDLA killed at its first loss: a crowded-site move keeps the degree,
but a jump out contributes zero. The same transient-system argument
gives spanning; the nonlosing addition matrices and their vacuum rows
give independence. In particular all degrees above m vanish. This proves
the finiteness used in the preceding local calculation. It also shows
that A_K is semisimple exactly when every k_i>0: in that case all Y_R
are empty; otherwise Y_empty is nonempty and its local factor contains
nonzero nilpotents.

### A polytope realizes each homogeneous factor

Put M=I-L, a nonsingular M-matrix. For b_i>0 consider

    P_b={x in R^m:x>=0, M^T x<=b}.                           (A5)

It is bounded since M^(-T)>=0 gives x<=M^(-T)b, and is full-dimensional.
Its paired lower and upper facets x_i=0 and (M^T x)_i=b_i do not meet:
at x_i=0 the latter left side is nonpositive. For every S subset [m],
there is a vertex x_(outside S)=0, x_S=(M_SS^T)^(-1)b_S>0; the remaining
upper inequalities are strict. Any vertex must choose one active facet
from each pair, since fewer than m active facets cannot determine a
vertex and both members of a pair cannot be active. The principal
submatrices are invertible. Thus P_b is a simple combinatorial cube,
and all b>0 have the same normal fan.

Let v_i,u_i be its lower and upper facet classes. Their normal vectors
are -e_i and column_i(M), respectively. The linear translation relations
give v_i=(Mu)_i, and the minimal empty-facet intersections are precisely
the pairs, giving u_i v_i=0. Its volume algebra is therefore H_L.
This presentation and the real-simple-polytope Lefschetz and
Hodge--Riemann theorems are stated in
[Bressler--Lunts, §4.1, Theorems 4.6, 4.8--4.9](https://arxiv.org/pdf/math/0302236),
reviewing Timorin's work. They apply to arbitrary real normals; no
rationality or approximation argument is needed. The ample class is
ell_b=sum_i b_i u_i. Consequently

    ell_b^(m-2d): (H_L)_d -> (H_L)_(m-d) is an isomorphism,
                                  0<=d<=floor(m/2).          (A6)

The corresponding primitive bilinear form, with sign (-1)^d, is positive
definite. This is a theorem application to (A5), not a new proof of the
polytope theorem. The original 1999 Timorin PDF was not obtained; the
quoted presentation and hypotheses were read in the accessible primary
Bressler--Lunts paper.

The volume normalization is useful probabilistically. Let Top pick the
coefficient of u_1...u_m in degree m. Then

    Top(ell_b^m)/m! = det(M) Vol(P_b).                        (A7)

For completeness, the volume polynomial F(b)=Vol(P_b) is annihilated by
partial_i (M partial)_i. Indeed allow lower supports c, so x>=-c;
translation gives F(b+M^T c), and the disjoint paired facets make the
mixed derivative in c_i,b_i zero. Since (H_L)_m is one-dimensional,
these relations determine its degree-m dual polynomial up to scale.
The full mixed derivative partial_1...partial_m F equals 1/det(M):
the corresponding finite difference is the volume of
b_i<(M^T x)_i<=b_i+epsilon_i for every i. Here x>=0 follows from
M^(-T)>=0 and b>0. Changing variables gives volume
product_i epsilon_i/det(M). The full mixed derivative of Top(ell_b^m)/m!
is one, proving (A7).

The vertices and support functions are linear in b on this fixed fan,
so P_(b+c)=P_b+P_c. This extends to b,c>=0 by continuity. Thus (A7) is
a mixed-volume polynomial in the bodies P_(e_i). If eta has total mass
m, then

    Top(product_i u_i^eta(i))
       =P_eta(ordinary IDLA fills every vertex without loss). (A8)

In particular Top(ell_b^m) is this full-filling probability when b is an
input probability vector and the m sources are sampled independently.
Write F_eta for the probability in (A8). The polynomial

    Q(b)=sum_(|eta|=m) binom(m;eta) F_eta b^eta                (A7a)

is Lorentzian: by (A7) it is a positive multiple of a convex-body volume
polynomial, to which
[Brändén--Huh, Theorem 4.1](https://web.math.princeton.edu/~huh/Lorentzian.pdf)
applies. Their Proposition 4.4, equivalently the mixed-volume
Alexandrov--Fenchel inequality here, gives for |eta|=m-2

    F_(eta+e_i+e_j)^2 >= F_(eta+2e_i) F_(eta+2e_j).          (A7b)

The same inequality holds with each chosen source replaced by any
nonnegative mixture of sources. Also Q^(1/m) is concave on the positive
orthant, and the positive coefficient support is M-convex. These are
claims about full filling with exactly m particles. Their coefficient
normalization is essential: (A7b) uses probabilities F_eta, not the raw
monomial coefficients of Q. None of these statements applies automatically
to a signed ARW occupation difference.

### Jordan blocks and a small exact example

For a driving law rho with rho_i>0 for every i, the local component of
sum_i rho_i z_i at R is lambda_R plus the positive linear form

    sum_(y in Y_R) [rho_y+
       sum_(i in R)rho_i((I-K_RR)^(-1)K_R,Y_R)_(i,y)] u_y.

By (A6) and the Hilbert series (1+t)^m_R, its Jordan blocks have sizes

    m_R-2d+1, with multiplicity
    binom(m_R,d)-binom(m_R,d-1),   0<=d<=floor(m_R/2).        (A9)

Here binom(m,-1)=0. This follows by decomposing into primitive Lefschetz
strings. The blocks for equal scalar eigenvalues from different R are
collected by direct sum. The same description holds for the driven ARW
chain by similarity. Without full-support driving the displayed linear
form may lie on the boundary of the ample cone, so (A9) is not asserted.

For 0->1->2 and 2->0,s with probabilities 1/2, the feasible R are
empty,{2},{1,2},V. Their characters are (0,0,0),(0,0,1/2),
(0,1/2,1/2),(1,1,1), and their local dimensions are 4,2,1,1.
At the zero character the relations are u_0^2=u_0 u_1, u_1^2=0.
For rho=(1/2,1/3,1/6), the blocks are [3,1] at 0, [2] at 1/12,
[1] at 1/4, and [1] at 1. `verify_addition_algebra.py` checks the
quadratic relations, regular basis, character fibers, and these block
sizes by exact ranks. This is one fixed example, not evidence replacing
the general algebra or polytope arguments.

For a biased nearest-neighbor interval, the homogeneous relations become
(q+1)u_i^2=q u_i u_(i-1)+u_i u_(i+1), with zero boundary variables.
This is exactly the q-Klyachko algebra studied by
[Nadeau--Tewari, §§1,3,5,9](https://www.math.utoronto.ca/vvtewari/qKlyachko_v1.pdf),
which already relates that algebra to IDLA and volume polynomials.
This identifies an existing mathematical neighborhood rather than a
claim that the IDLA/algebra connection is new.

For this same interval algebra, [Ruizhen Liu, Theorem 1.4 and §3](https://ruizhenliuu.github.io/qKlyachko.pdf)
constructs a projective toric fan from the q-deformed Cartan matrix and
obtains the Lefschetz and Hodge--Riemann package for positive rational q.
Its rays are the signed coordinate vectors and matrix columns, matching
(A5) up to a global sign and column normalization. The M-matrix argument
above supplies the present arbitrary-routing formulation; no claim about
its priority in the literature has been checked. The original
[Diaconis--Fulton paper, introductory algebra discussion](https://stacks.stanford.edu/file/xg052sf4924/xg052sf4924_EFS-NSF-374.pdf)
already constructs the commutative IDLA algebra on subset bases. Its
scanned introduction was obtained and read; this is not a claim to have
audited the whole paper.

### The inverse basis consists of opposite-facet combinations

The connection with the existing inverse is exact. In A_K set

    v_i=(Mz)_i-k_i,    w_i=z_i+r_i v_i.

Equation (A1) says z_i v_i=0. The complete-block formula says that the
stable-state basis vector at C is w_C=product_(i in C)w_i, expressed in
the all-active basis z_B by the row E_r(C,B). Its highest-degree part
lies in H_K, where the corresponding classes are u_i+r_i v_i with
v_i=(Mu)_i, exactly the two opposite facet classes of (A5).

In H_K, for C subset V, the vertex intersection number is

    Top(u_C v_(V\C))=det(M)/det(M_CC).                       (A10)

For a direct normalization check, permit all lower and upper supports
to vary. The full mixed derivative in the n facets meeting at that
vertex is the reciprocal absolute determinant of their normals, namely
1/det(M_CC). Equation (A7) changes integration to Top by det(M).
All determinants are positive. Hence

    Top(product_i (u_i+r_i v_i))
       =det(M) sum_(C subset V) r_(V\C)/det(M_CC)=Z_V.       (A11)

Thus the creator-sampler distribution nu_V(C) in Section 30 is exactly
the normalized list of terms in this opposite-facet intersection product.
The loop/forest determinant normalizers and the polytope intersection
normalizers describe the same numbers. The verifier checks all eight
vertex products on the single-leak three-cycle. It also checks nine
instances of (A7b) there and independently verifies (A7) by the shoelace
area of the asymmetric two-site polygon from Section 33.

These inverse classes are effective combinations of facets, but need not
be nef. On 0->1->sink, H_K has u_0^2=u_0 u_1 and u_1^2=0. Thus

    Top((u_0+r_0 v_0)^2)=1-r_0^2.                           (A12)

It is negative for p_0<1/2. This is a concrete obstruction to treating
all sleep-insertion factors as nef classes. The threshold resembles the
high-sleep arguments in Section 26, but no implication between those
arguments and (A12) has been proved.

There is also a flat algebraic interpolation:

    A_K(t)=R[t,z]/<z_i((Mz)_i-k_i t)>.                       (A13)

Weight every IDLA loss by t. The same finite absorption equations and
vacuum rows prove that the squarefree z_B form a free R[t]-basis.
More explicitly, z^eta=sum_B nu_eta(B)t^(|eta|-|B|)z_B, with nu_eta
the ordinary IDLA output law. Internal excursions preserve degree and
are t-independent; losses decrease particle count. Thus their finite
block resolvents introduce only constant real denominators, and the
displayed reduction is polynomial in t.
The specializations are H_K at t=0 and A_K at t=1. Equivalently start
with the fan face ring R[u,v]/<u_i v_i> and impose the linear equations
Mu-v=k t. This keeps the geometry of (A5) while deforming the ordinary
cohomology relations. No geometric identification of the target
occupation observable is inferred from this presentation.

One simple proposed bridge already fails. On A_K let Lambda extract
the coefficient of z_V; then Lambda(z^eta) is the ordinary IDLA
full-filling probability, with losses permitted if |eta|>|V|.
The pairing Lambda(ab) is nondegenerate: pass to highest degrees and
use Poincaré duality of H_K. Every target functional therefore has a
unique representer g_v with

    q_v(eta)=Lambda(g_v z^eta).

If g_v had nonnegative coefficients in the squarefree z_B basis,
IDLA addition monotonicity would prove the conjecture immediately.
But on 0->1->sink, for target 1, direct calculation gives

    g_1=p_1 z_0-p_0 p_1 z_1+p_0 p_1.                        (A14)

The coefficient of z_1 is strictly negative. This also rules out an
arbitrary positive mixture of ordinary IDLA completion inputs on this
same kernel: every positive polynomial reduces with nonnegative
coefficients to the squarefree basis. Formula (A14) is checked directly
on the four binary inputs, whose target probabilities are
0,(1-p_0)p_1,p_1,p_1; surplus IDLA then covers all inputs. This is a
failure of that particular representation, not a monotonicity failure.
The mixture here must be fixed across inputs. Countably many completion
inputs do not help: at eta=1_V every completion succeeds, so the total
mixture mass is at most one; Tonelli again reduces it to nonnegative
squarefree coefficients and contradicts the same unique representer.

### Capacities and full filling with surplus inputs

**Proved extension; independent AI audit passed.** Give vertex i capacity
d_i>=1 and put N=sum_i d_i. IDLA now parks whenever the current count
is below that capacity; kill the history at its first loss. Its graded
algebra and Hilbert series are

    H_(L,d)=R[u]/<u_i^d_i ((I-L)u)_i>,
    Hilb(H_(L,d),t)=product_i(1+t+...+t^d_i).                 (A15)

The stable monomials have 0<=eta_i<=d_i. Spanning follows from the
same finite killed-absorption equations as before, and applying to the
vacuum proves independence. The top monomial is u^d.

There is a direct polytope in dimension N. Retain x_i>=0 and introduce
y_(i,a)>=0 for 1<=a<=d_i-1, with inequalities

    (M^T x)_i+sum_a y_(i,a)<=b_i,    b_i>0, M=I-L.          (A16)

It is bounded and full-dimensional. Group i has d_i+1 facets: x_i=0,
the d_i-1 lower y-facets, and its upper inequality. All members of a
group cannot meet, since their equalities would force a nonpositive
left side to equal b_i>0. A vertex must omit exactly one facet from
every group. All such choices occur: let S be the groups omitting
x_i=0, set x_S=(M_SS^T)^(-1)b_S and x_(outside S)=0. In a group
omitting the upper facet, all its y's vanish. In a group omitting a
particular lower y-facet, only that y is nonzero, equal to
b_i-(M^T x)_i>0. Groups in S have all y's zero. The active normals
are independent by the same principal-submatrix calculation. Hence
the polytope is simple with the face structure of product_i Delta_(d_i),
and its normal fan is constant for b>0.

The upper facet and all lower y-facets in group i have the same class
u_i, by translation in y_(i,a). Translation in x gives the remaining
class v_i=(Mu)_i. The minimal nonfaces give v_i u_i^d_i=0, exactly
(A15). Thus this is a volume algebra, and positive sum_i b_i u_i is
a Lefschetz class. Writing Top_d(u^d)=1, the normalization remains

    Top_d((sum_i b_i u_i)^N)/N!
                       =det(M) Vol_N(P_b^d).                (A17)

For example, integrate the y-simplex fibers over x. Their volume is
product_i (b_i-(M^T x)_i)^(d_i-1)/(d_i-1)! on the base P_b.
Differentiating d_i-1 times in each b_i gives Vol(P_b); one more
derivative in every b_i gives 1/det(M), as in (A7). This fixes the
degree-N normalization. Equivalently choose the upper facet and all
lower y-facets in each group: they meet at x=M^(-T)b,y=0, with
normal determinant det(M), and their product class is u^d.

The following consequence concerns ordinary capacity-one IDLA on V,
where n=|V|, and any fixed number N>=n of inputs. For N>n, add a
buffer vertex o of capacity N-n, route every original sink jump into
o, and route a jump out of o to a new absorbing sink. The original
vertices retain capacity one. The enlarged kernel is transient, with
interior matrix [[K_V,k],[0,0]]. Give no initial input to o. A killed
history fills all capacities iff ordinary IDLA on V fills V: exactly
N-n particles then enter the buffer; if V remains deficient the
buffer overflows. This correspondence preserves all original routing
weights. For N=n no buffer is needed.

Therefore, with F_eta=P_eta(IDLA fills V) and |eta|=N,

    Q_N(b)=sum_(|eta|=N) binom(N;eta) F_eta b^eta             (A18)

is a convex-body volume polynomial up to a positive constant. It is
(A17) specialized to b_o=0, a nonnegative support limit. It follows
that Q_N is Lorentzian, Q_N^(1/N) is concave for b>0, and
F_(eta+e_i+e_j)^2>=F_(eta+2e_i)F_(eta+2e_j) for |eta|=N-2.
This concerns redistribution of a fixed total number of inputs.

It transfers exactly to full ARW occupation. Surplus IDLA gives

    P_eta(ARW occupies all V)=F_eta/Z_V.                     (A19)

Indeed a smaller IDLA support cannot produce n residents under ARW;
the full support does so with P_p(V,V)=1/E_p(V,V)=1/Z_V. Thus the
fixed-N input polynomial for full ARW occupation is also Lorentzian,
for every allowed heterogeneous sleep vector. As one consequence,
averaging an iid input distribution over automorphisms of K cannot
decrease full-filling probability: concavity of Q_N^(1/N) and symmetry
give this directly. If those automorphisms act transitively, uniform
input maximizes full filling. The sleep vector need not be symmetric
for this last ARW conclusion, since its factor Z_V is independent of
the input law. No corresponding assertion for a proper target set or
for fixed-input output marks follows from this argument.

At zero buffer intensity the class used for these probabilities need
not be Lefschetz. For one original site and N=3, the buffer algebra is
R[u,w]/<u(u-w),w^3>. Multiplication by u from degree one to degree two
kills u-w, although Top(u^3)=1 remains the correct filling probability.
The limiting volume and Lorentzian conclusions survive; a Lefschetz
isomorphism at that boundary is not part of the theorem.

There is also an exact-output version, still varying the inputs. Fix
nonempty B subset V, m=|B| and N>=m, and define Q_(N,B) as in (A18)
with coefficients P_eta(IDLA output exactly B). Inputs outside B make
these coefficients zero. Keep K_BB, kill every jump to V\B, and send
only original sink jumps k_B into a buffer of capacity h=N-m. The
nonlosing full-capacity event is precisely exact output B. Hence

    Q_(N,B)(b)/N! =det(I-K_BB)/h!
       integral_(x>=0,(I-K_BB)^T x<=b_B)(k_B dot x)^h dx.    (A20)

For h=0 use no buffer and read the integrand as one. This is again
a volume polynomial and therefore Lorentzian, with the zero-polynomial
convention. For N<m it is zero; for B empty only N=0 contributes.
For nonempty B and h>0 it is identically zero exactly when k_B=0.
Replacing k_B by the larger deficit 1-K_BB 1 would count forbidden
exits and invalidate the identity. This extension passed the same
independent capacity audit.

In contrast, fixed-input polynomials in output marks need not be
Lorentzian even for IDLA. On x->a or c equally, a->b->s and c->d->s,
three inputs at x give supports {x,a,b},{x,a,c},{x,c,d} with
probabilities 1/4,1/2,1/4. The homogeneous output polynomial is

    T(z)=(z_x/4)(z_a z_b+2 z_a z_c+z_c z_d).                 (A21)

Removing a from {x,a,b} cannot be repaired by adding c or d: b loses
its prerequisite. Thus this rooted-antimatroid level fails matroid
exchange. More directly, restricting the quadratic derivative T_x
to (z_a,z_b,z_c,z_d)=(u,u-2v,v,v) gives (u^2+v^2)/4, a
two-dimensional positive definite restriction. This violates the
Lorentzian Hessian condition. `output_mark_audit()` checks this single
five-site law and the restriction exactly; main reran it. The input
matching theorem fixes target slots, whereas varying the output changes
the forbidden vertices and the eligibility graph.

The focused capacity certificate in `verify_addition_algebra.py` uses
the asymmetric two-site kernel from Section 33 and four iid inputs with
law (2/3,1/3). Its full-filling probability is 631/648. This agrees both
with ordinary four-step IDLA and the 12-state capacity-(1,1,2) buffer
algebra. Independently integrating the two-dimensional base polygon
with its triangular buffer fibers gives four-dimensional volume
3155/46656; multiplication by 4! det(M)=72/5 gives the same probability.
The certificate also checks the capacity relations and commutativity.

### Source counts conditioned on the output

**Proved consequence; independent AI audit passed.** Let N>=1 and
take N iid input sources with a strictly positive distribution rho.
Condition on an exact IDLA output B of positive probability, and write
eta for the source-count vector. For Q=Q_(N,B) and its conditional
mean m_i=E[eta_i | output B],

    E[z^eta | output B]=Q(rho_1 z_1,...,rho_n z_n)/Q(rho),
    m_i=rho_i partial_i log Q(rho),  sum_i m_i=N.           (A22)

Concavity of Q^(1/N) bounds it by its tangent plane at rho. Hence
for c>=0,

    Q(c)/Q(rho)<=[sum_i (m_i/N)c_i/rho_i]^N.

In particular, for every real vector t,

    E[exp(t dot eta) | output B]
       <=[sum_i (m_i/N)exp(t_i)]^N.                       (A23)

The right side is the moment-generating function of a multinomial
vector with N draws and the same conditional mean. Differentiating
along any line at t=0 gives the matrix inequality

    Cov(eta | output B)<=diag(m)-m m^T/N                   (A24)

in positive-semidefinite order. Each coordinate has the corresponding
binomial Chernoff bounds, and arbitrary weighted source sums have the
multinomial Chernoff bound. This is a moment-generating-function order;
convex order and signs of individual off-diagonal covariances are not
asserted from (A24) alone.

For full ARW occupation, (A19) cancels the sleep-dependent constant in
Bayes' formula. Under **any** fixed prior on finite all-active inputs,
conditioning on full ARW occupation gives exactly the same posterior
input law as conditioning on full IDLA occupation, provided that event
has positive prior probability. This posterior is independent of the
sleep vector. Thus (A22)--(A24) also apply to full ARW occupation for
fixed-N iid input. No proper-target ARW analogue is claimed.

The remaining bridge to the target conjecture is specific: identify an
occupation-response functional or a cancellation identity controlled by
these positive/ample classes. The inverse P_p^(-1), the observable
P_p 1_(v occupied), and signed source differences do not yet have such
an interpretation. Merely knowing the addition algebra's spectrum,
Lefschetz property, or positive full-filling coefficients does not settle
their signs.


## 35. Poissonized IDLA as fluid growth with exponential capacities

**Proved representation, 2026-10-19; independent AI audit passed for
(F16)--(F21) and the conditional-current formula (F24).** This is an
IDLA representation. Its full-occupation probability
transfers to ARW by (A19); a corresponding representation of the whole ARW
occupation law is not asserted.

Let b_i>=0 be input intensities and let the initial particle counts be
independent Poisson(b_i). Equivalently, drive ordinary IDLA by independent
Poisson arrival processes of rates b_i and observe at time one. Set
M=I-K_V and k_i=K(i,s). Independently sample xi_i~Exp(1). Send deterministic
fluid into site i at rate b_i. A site stores fluid up to capacity xi_i;
once full it routes additional fluid according to the proportions K(i,j),
with sink mass lost. Let u(t) be total emitted fluid and mu(t) stored mass.
The stabilized fluid equations are

    mu=t b-M^T u,   u>=0, 0<=mu<=xi,
    u_i(xi_i-mu_i)=0,
    equivalently u=(t b-xi+K_V^T u)_+.                       (F16)

The solution is unique. For w=M^(-T)1>0, K_V^T w=w-1, so the positive-part
map is a contraction in the weighted sup norm with constant
max_i(1-1/w_i)<1. Its fixed point gives mu>=0 as well: at u_i=0,
mu_i=t b_i+sum_j K(j,i)u_j>=0; at u_i>0, mu_i=xi_i.

**Growth-process identity.** For each fixed b, the saturated set of this
fluid process has the same law, jointly in t>=0, as the occupied set of
Poisson-driven IDLA. To prove this, suppose the saturated set is B.
The derivative of emitted fluid there is

    dot u_B=(I-K_BB^T)^(-1)b_B,

and an unsaturated vertex y receives fluid at rate

    a_y(B)=b_y+sum_(x in B)b_x
                     [(I-K_BB)^(-1)K_(B,y)]_x.              (F17)

These rates are constant until the next saturation. Conditional on the
past saturation times and identities, the remaining capacity deficits
are independent Exp(1), by the exponential race and memorylessness.
Thus the next transition B->B+y has rate a_y(B). The same is true for
IDLA: an arrival at a vacant y parks immediately, while an arrival in B
walks until its first vacancy or the sink, giving exactly (F17). Lost
arrivals are self-transitions and do not alter the support generator.
Both finite pure-birth processes start empty. This proves equality of
their process laws. Zero rates cause no transition, and simultaneous
positive-rate saturations have probability zero. Saturation is defined
by mu_i=xi_i; at a fixed deterministic time its difference from u_i>0
is a null event.

Consequently full occupation has the particularly simple representation

    P(Poisson-IDLA fills V)
      =P(M^(-T)xi <= M^(-T)b).                              (F18)

If all sites are saturated, (F16) forces u=M^(-T)(b-xi)>=0. Conversely
this nonnegative u, with mu=xi, solves the same fixed-point equations,
so uniqueness gives full saturation. Threshold equalities are null.
The representation is a coupling of these marginal laws as b varies;
it is not asserted to coincide with the original joint arrow-stack and
Poisson-count coupling for every b simultaneously.

### Polyhedral chambers give every exact output probability

For B subset V write M_B=I-K_BB and define

    P_B(b_B)={u_B>=0:M_B^T u_B<=b_B}.

In the chamber where exactly B is saturated, u_(outside B)=0 and
xi_B=b_B-M_B^T u_B. Every y outside B has the inequality
xi_y>b_y+sum_(i in B)K(i,y)u_i. Integrating these independent exponential
tails and changing variables xi_B to u_B gives

    P(Poisson-IDLA output=B)
      =e^(-sum_i b_i) det(M_B)
                     integral_(P_B(b_B)) exp(k_B dot u_B) du_B. (F19)

Here k_B records only jumps to the original sink, not exits to V\B:
the latter occur in the tail factors and cancel from the exponent.
For B empty the integral and determinant are both one. The formula is
first derived at b>0; continuity gives all nonnegative intensities.
It is consistent with the fact that any input outside B prevents B
from being the exact IDLA output. The chambers partition the exponential
orthant, so summing (F19) over B gives one.

For B=V, (F19) also follows from the buffer-polytope construction.
With n=|V| and Q_N from (A18), integrate the buffer simplex fibers:

    Q_N(b)/N! =det(M) integral_(P_V(b))
                          (k dot u)^(N-n)/(N-n)! du, N>=n.  (F20)

Poisson summation is legitimate by nonnegativity and compactness of
P_V(b); summing the powers gives exp(k dot u). Changing variables
xi=b-M^T u gives (F18) directly. Thus the capacity-algebra and
fluid-generator derivations agree without identifying their sample paths.

For each fixed exact support B, (F19) is log-concave in the input vector b.
The integrand, including the indicator of u_B>=0 and M_B^T u_B<=b_B,
is a log-concave function jointly in (b,u_B). Its marginal is log-concave
by [Prékopa, Theorem 6](https://rutcor.rutgers.edu/Prekopa/pdf/SCIENT2.pdf).
If joint integrability is required, first restrict b to a bounded convex
box; u_B is then bounded too. Any two finite input vectors and their
segment fit inside such a box, giving the unrestricted inequality.
The factor exp(-sum b_i) is already included in this argument. In
particular full ARW occupation under independent Poisson all-active
inputs has probability (F18)/Z_V and is log-concave in b as well.
This does not give log-concavity for arbitrary sums of exact-support
events, or coefficientwise inequalities in the Poisson expansion.

`poisson_cell_audit()` in `verify_addition_algebra.py` checks all four
exact-support formulas on the same asymmetric two-site kernel. It
compares the finite birth-chain ODE solutions with exact exponential
integrals over the two singleton intervals and the full quadrilateral.
At b=(2/3,1/3), the full probability is
1+3 exp(-1)-2 exp(-2/3)-2 exp(-5/6). The check compares rational
coefficients of these exponentials, with no quadrature or rounding.

Surplus IDLA gives the exact remaining ARW mixture

    mu_b^ARW(C)=sum_B P(Poisson-IDLA output=B) P_p(B,C).       (F21)

The coefficients are positive, but sums of log-concave functions need
not be log-concave. More generally, the fluid set increases as b increases,
whereas no inclusion-monotonicity property of the continuation P_p is
available. This is the specific obstruction to importing the fluid
coupling as an ARW occupation proof.

### Proper-target input log-concavity fails even on a DAG

Take deterministic arrows 0->2, 1->2, 2->s and sleep probabilities
p_0,p_1,p_2 in (0,1). With two all-active input particles,

    q_2(2 delta_0)=q_2(2 delta_1)=p_2,
    q_2(delta_0+delta_1)=p_2(1-p_0 p_1).                     (F22)

A pair at one source must emit at least one particle, while separate
sources both retain their sole particle with probability p_0 p_1.
After any positive number of arrivals at 2, its last particle sleeps
with probability p_2. Thus the target's degree-two input polynomial,
restricted to sources 0 and 1, is

    p_2[b_0^2+2(1-p_0 p_1)b_0 b_1+b_1^2].

Its Hessian has two positive eigenvalues, or equivalently its normalized
mixed coefficient squared is smaller than the product of its pure
coefficients. It is not Lorentzian. This does not contradict the
fixed-output IDLA input statements or the full-occupation ARW statement.

The failure persists after Poissonization. With independent Poisson
inputs of means b_0,b_1 and no input at 2,

    q_2^Pois(b)=p_2[1-e^(-b_0-b_1)(1+p_0 b_0)(1+p_1 b_1)]. (F23)

At source i, no particle is emitted iff its input count is zero or is
one and that particle sleeps, of probability e^(-b_i)(1+p_i b_i).
Independence of the two source books proves (F23). If p_0=p_1=p,
then at (b_0,b_1)=(t+h,t-h), t>0, |h|<t,

    q_2^Pois(t+h,t-h)=q_2^Pois(t,t)+p_2 p^2 e^(-2t)h^2.

Equal symmetric endpoints therefore have larger values than their
midpoint, violating log-concavity. The target probability is nevertheless
coordinatewise increasing. This graph is already in the proved DAG
monotonicity class; the calculation rules out the proposed concavity
strengthening rather than providing a counterexample to monotonicity.

### Conditional exit counts form a Poisson mixture

Keep distinct absorbing exit labels o and write k_o(i)=K(i,o), so
k=sum_o k_o. The fluid construction also represents their IDLA counts.
At occupied support B, Poisson-driven IDLA has loss events of label o
at rate

    r_o(B)=sum_(i in B)b_i[(I-K_BB)^(-1)k_o|B]_i
          =(k_o|B) dot dot u_B.

These events do not change B. Construct the support chain from its growth
clocks and independently add unit-rate Poisson processes run at these
integrated intensities. The resulting joint generator for support and
loss counts is exactly that of Poisson-driven IDLA. Given the fluid
history through time one, its labeled losses L_o are therefore independent
Poisson variables with means k_o dot u(1). The same conditional statement
holds given only terminal u, since those means depend only on u.

In particular (F19) has the marked version

    E[1_(output=B) product_o z_o^L_o]
      =e^(-sum_i b_i) det(M_B)
          integral_(P_B(b_B)) exp(sum_o z_o k_o|B dot u) du. (F24)

Conditional on B with positive probability, u has density proportional
to exp(k_B dot u) on P_B(b_B). Given u, the labeled exit counts are
independent Poisson(k_o|B dot u). Thus the total loss count is a mixed
Poisson variable with mixing parameter k_B dot u; its conditional
variance is at least its conditional mean. No sign is asserted for
cross-covariances between different random label intensities.
Conditioning additionally on no exits cancels the exponential tilt:
u is then uniform on P_B(b_B), whenever that event has positive probability.

This gives a finite representation of Poisson-input ARW as well. Draw
the exponential capacities, solve (F16), take its saturated B, and
stabilize one active particle at each site of B using fresh ARW books.
The terminal law is exactly (F21). To retain labeled losses, independently
sample the Poisson variables in (F24) conditional on u and add the losses
of this final binary ARW stabilization. Surplus IDLA and fresh suffixes
justify the composition. This samples the Poisson-mixed law; it does
not supply a coupling conditional on a prescribed input-count vector.

The same two-site symbolic cell audit checks (F24) at sink weight z=2,
in addition to z=1. The full-support weighted value is
(3/4)exp(-1)+(1/4)exp(1)-(1/2)exp(-2/3)-(1/2)exp(-1/3).
The independently solved tilted birth ODE and the exponential polygon
integral agree coefficient by coefficient for all four supports.

### A Green-function polytope for the cover time

Let rho be an input probability distribution from which every vertex
is reachable, and let G=(I-K_V)^(-1). Put g_i=(G^T rho)_i>0 and let
a_i be column i of G divided by g_i. Along the rate-one drive rho,
the first full-occupation time has the representation

    T=max_i a_i dot xi,     xi_j independent Exp(1).         (F25)

This follows by applying (F18) to b=t rho for every t in the same fluid
construction. The vectors satisfy rho dot a_i=1. Thus T is the support
function of conv{a_i:i in V}, evaluated at an independent exponential
vector. Inclusion of these Green-column convex hulls is one sufficient
condition for pathwise ordering of the corresponding cover times under
this coupling. It is not an equivalence claim.

The ordinary hitting probabilities reappear directly. Write
h_j(i)=P_j(walk hits i before the sink), including h_i(i)=1. Then
G_(j,i)=h_j(i)G_(i,i), and H_i=P_rho(hit i)=g_i/G_(i,i). Hence
max_j (a_i)_j=1/H_i, and

    lim_(t->infinity) t^(-1) log P(T>t)=-min_i H_i.          (F26)

For the lower bound use the single summand xi_i/H_i. For the upper
bound, every exponential moment of T below min_i H_i is finite,
because exp(c max_i a_i dot xi)<=sum_i exp(c a_i dot xi).
This recovers the relaxation exponent already identified in (X12)
and the cited Bristiel--Salez result by a different representation.
No prefactor or polynomial correction is asserted in (F26).

With **full-support driving**, rho_i>0 for every i, the tail can be
sharpened. Put H=min_i H_i and I={i:H_i=H}. Then

    c_i=product_(j!=i)(1-h_j(i))^(-1),  C=sum_(i in I)c_i,
    P(T>t)=C exp(-Ht)+O(exp(-(H+epsilon)t))                 (F27)

for some epsilon>0 depending on K,rho. Every c_i is finite. Indeed if
j!=i and h_j(i)=1, hitting j implies hitting i, while a simple j-to-sink
path has an i-to-sink suffix avoiding j. Thus h_i(j)<1 and
H_i-H_j>=rho_i(1-h_i(j))>0. Such i cannot minimize H.

For critical i write H T_i=xi_i+S_i, with
S_i=sum_(j!=i)h_j(i)xi_j independent of xi_i. Its exponential moment
at one is c_i, and its moment at 1+delta is finite for some delta>0.
Conditioning on S_i yields
P(T_i>t)=c_i exp(-Ht)+O(exp(-(1+delta)Ht)). Noncritical forms have
an exponential moment at a parameter larger than H. For distinct
critical i,l, every coefficient of
H(T_i+T_l)/2=sum_j[(h_j(i)+h_j(l))/2]xi_j is strictly below one.
Their joint exceedance therefore also has a faster exponential rate.
The union and two-term Bonferroni bounds prove (F27).

There is a matching discrete statement. If tau is the number of iid
arrivals needed to fill V and n=|V|>=2, then, for some 0<beta<1-H,

    P(tau>N)=C(1-H)^N+O(beta^N).                           (F28)

The largest proper addition eigenvalue 1-H is semisimple, of
multiplicity |I|, for both IDLA and the similar ARW addition matrix.
To see this, (X12) gives the diagonal eigenvalue 1-H_A at support
V\A, where H_A=P_rho(hit A). It equals 1-H exactly for A={i}, i in I:
if A contains a critical i and another j, then
H_A>=H_i+rho_j(1-h_j(i))>H; if it contains none, H_A>H directly.
These critical coatoms communicate only with themselves and V.
Downward substitution in the triangular eigen-equations gives one
independent eigenvector per critical coatom and no Jordan coupling.
All remaining proper eigenvalues are smaller. The spectral projector
at 1-H has the same scalar coefficient in D^N and exp(t(D-I)); (F27)
identifies it as C. Lower Jordan factors are absorbed by beta. A
direct sink edge and n>=2 ensure 0<H<1. For n=1, tau=1 and T is
Exp(1), so (F28)'s ratio asymptotic is not asserted. This proof was
independently supplied and reviewed; it uses no new computation.

### Arrival-count moments and log-concavity

Keep rho strictly positive and normalized, n>=1, and let f_N=P(tau<=N).
For P=P_rho, take U uniform on P and write L=k dot U. Equation (F20)
gives

    a_h=f_(n+h)/binom(n+h,n)=f_n E[L^h],
    f_n=n! det(M) Vol(P)>0.                                (F29)

Here 0<L<1 almost surely and its support is [0,1]. In fact
u_*=M^(-T)rho is the coordinatewise maximum of P, k dot u_*=1,
and t u_* is interior for 0<t<1. Each level set of the nonconstant
linear form k dot u has zero volume. Therefore, with Delta the
forward difference, for all h,r>=0,

    (-1)^r Delta^r a_h=f_n E[L^h(1-L)^r]>0.                 (F30)

All shifted Hankel and localizing matrices
((-1)^s Delta^s a_(r+i+j))_(0<=i,j<=m) are positive definite:
their quadratic forms are f_n E[L^r(1-L)^s q(L)^2]>0 for nonzero q.
In particular a_h^2<a_(h-1)a_(h+1) for h>=1.

An additional geometric theorem gives the opposite curvature for
the unnormalized filling probabilities. The s=1 case of
[Fradelizi--Li--Madiman, Theorem 3.1](https://arxiv.org/pdf/1512.01490)
says that product_(j=1)^n(p+j) integral_P g(u)^p du is log-concave
for p>0 when g>=0 is concave on a convex body P. Apply it to
g=k dot u, extended by zero off P. This is an integrable 1-concave
function in the source's definition. Dominated convergence extends
the result to p=0. Consequently

    F(p)=f_n binom(n+p,n)E[L^p] is log-concave on [0,infinity),
    f_N^2>=f_(N-1)f_(N+1) for every N>=0.                  (F31)

The initial zeros N<n are included. Moment log-convexity alone would
not prove this; the normalized-moment theorem is essential. Combining
the two inequalities gives, for h>=1,

    1<=f_(n+h)^2/[f_(n+h-1)f_(n+h+1)]
      <(n+h)(h+1)/[h(n+h+1)].                              (F32)

The ratios f_(N+1)/f_N decrease to one, while a_(h+1)/a_h increase
strictly to one. Indeed f_N increases to one, since seeing every
vertex among the sources already fills V. Also a_h~n!h^(-n).
Strict log-concavity of f_N is not claimed: for one site f_N=1
for N>=1. These are arrival-count statements for a fixed source
distribution, not a target occupation order under prescribed inputs.

The same argument covers **every exact IDLA output B**, with the
input law rho_B left unnormalized. Write m=|B| and
P_B={u>=0:M_B^T u<=rho_B}, where M_B=I-K_BB. For nonempty B,

    f_(m+h)(B)=f_m(B)binom(m+h,m)E[L_B^h],
    f_m(B)=m!det(M_B)Vol(P_B)>0,
    L_B=k_B dot U,  U uniform on P_B.                      (F32a)

Here k_B contains only original-sink probabilities. Applying the same
normalized-moment theorem in dimension m proves log-concavity of
N->f_N(B). Its support is {m} if k_B=0 and all N>=m otherwise.
For B empty the support is {0}. These conclusions follow directly
from (F19); global sink accessibility makes M_B invertible but does
not imply accessibility of the original sink within B.

Put r_B=M_B^(-1)k_B and theta_B=rho_B dot r_B. Then
0<=r_B<=1 and 0<=L_B<=theta_B<=rho(B), with full support
[0,theta_B] when k_B is nonzero. Thus

    G_B(z)=sum_N f_N(B)z^N
          =f_m(B)z^m E[(1-z L_B)^(-m-1)], 0<z<1.          (F32b)

For a geometric prior P(N=j)=(1-z)z^j, conditioning on exact output B
gives an explicit surplus-count mixture. Tilt uniform U on P_B by
(1-z k_B dot U)^(-m-1), normalize, and conditional on U take
N-m to have negative-binomial probabilities
binom(m+h,m)(1-z k_B dot U)^(m+1)(z k_B dot U)^h.
The factors cancel to give z^(m+h)f_(m+h)(B)/G_B(z).
This specific mixture is log-concave, by the proved sequence theorem;
arbitrary negative-binomial mixtures need not have that property.
For proper B, f_N(B)<=rho(B)^N and tends to zero. When k_B is
nonzero its successive ratios decrease to theta_B<1, since moment
ratios tend to the upper endpoint theta_B. Exact-output sequences
are therefore generally unimodal rather than filling CDFs.

### Nearby mathematical formulations

Three primary-source comparisons clarify what has and has not been
imported. The deterministic fluid equations are exactly the oblique
Skorokhod problem S1--S4 in
[Debicki--Dieker--Rolski, Section 5](https://arxiv.org/pdf/math/0512119),
under W=xi-mu, initial state xi, driver -tb, reflection matrix M^T,
and regulator u. Their tree-specific explicit solution is not used.
The accessible general reflection equations alone do not identify the
exponential initial law with Poisson-IDLA. The original
Harrison--Reiman PDF was inaccessible and was not used.

On the deterministic path 1->2->...->n->s, the capacity region in
(F18) is exactly the
[Pitman--Stanley polytope, Section 1](https://math.mit.edu/~rstan/pubs/pubfiles/117.pdf):
xi>=0 and sum_(j<=i)xi_j<=sum_(j<=i)b_j. Our integral has exponential
weight; their principal object is ordinary volume.

For every accessible K, put G=M^(-T). In the notation of
[Xu, equation (3.4)](https://arxiv.org/pdf/0806.1127), which defines
exponential truncated powers as weighted linear pushforward densities,

    P(G xi<=G b)=E_(1,0)(G b | [G,I])
               =det(M) E_(0,1)(b | [M^T,I]).                (F33)

This follows by adjoining nonnegative slack variables and then setting
xi=b-M^T u. Both matrix cones are pointed: G>=0 with positive
diagonal, and M^(-1)1 is strictly positive on every column of [M^T,I].
Zero weights may require limits in Xu's separate generic-denominator
formula, which is not invoked here. These are exact substitutions
into existing definitions. No direct statement of the full stochastic
process identity was located in this bounded three-source reading;
no novelty claim follows from that absence.

### Association and convex avoidance events

**Proved; independent AI proof audit passed.** These statements concern
Poisson-input IDLA and its original exit counts, not ARW after the
binary-stabilization continuation. Put z=b-xi and, for B subset V,
G_B=(I-K_BB)^(-1). The emitted-fluid solution has the explicit
piecewise-linear form

    u_i(z)=max(0, max_(B contains i) [G_B^T z_B]_i).        (F34)

Indeed M^T u>=z and u>=0 imply
M_B^T u_B>=z_B+K_(V\B,B)^T u_(V\B)>=z_B, hence
u_i>=[G_B^T z_B]_i. At A={i:u_i>0}, complementarity gives
u_A=G_A^T z_A; this attains the maximum for i in A. For i outside A,
all displayed terms are nonpositive and the zero term attains it.
Thus u is increasing in every coordinate of z, and every u_i is
convex in z. All coefficients of the linear forms are nonnegative.

At any fixed b, the occupied indicator equals 1_(u_i>0) almost
surely. The possible exceptions are equality boundaries: conditional
on a candidate active set A not containing i, equality at i reads
xi_i=b_i+sum_(j in A)K(j,i)[G_A^T(b_A-xi_A)]_j. Each is a null
hyperplane under the continuous capacity law, and there are finitely
many A. This does not replace the saturated-site convention at random
jump times.

The occupied indicators are consequently coordinatewise decreasing
functions of independent xi_i. They are **positively associated**:
Cov(f,g)>=0 for bounded increasing functions of the occupied set.
The needed product-measure fact follows from the one-variable identity
2 Cov(f(X),g(X))=E[(f(X)-f(X'))(g(X)-g(X'))]>=0 and induction using
conditional covariance. Formula (F24) extends this result to the
joint fixed-time vector of occupied indicators and all original exit
counts. Given xi, realize each independent Pois(k_o dot u) count with
an independent inverse-CDF uniform. Every coordinate then increases
in the independent variables (-xi,U_1,...,U_l). In particular,

    Cov(L_o,L_l)=Cov(k_o dot u,k_l dot u)>=0, o!=l,
    Cov(1_(i occupied),L_o)>=0.                             (F35)

These are unconditional inequalities. Conditioning on a specified
occupied support is not asserted to preserve association. Fixed-input
IDLA, with deterministic total particle count, is also not covered.

Association does not give the FKG lattice inequality. On the existing
fork 0->2, 1->2, 2->s, take b_0=b_1=t>0 and b_2=0. The support {2}
has probability zero, whereas {0,2} and {1,2} each have probability
e^(-t)[1-e^(-t)(1+t)]>0. Their intersection and union therefore violate
the lattice inequality. By continuity the violation persists for
sufficiently small positive b_2.

For a target set R, write a_R(b)=P(no site of R is occupied). By (F34),
its capacity event is the convex set

    xi>=0, [G_B^T(b_B-xi_B)]_i<=0
    for every i in R and B containing i.                   (F36)

It is jointly convex in (b,xi). Integrating its indicator against
exp(-sum xi_i) and applying the same Prekopa theorem used for (F19)
proves that a_R is log-concave in b>=0. If needed, first restrict b
to an arbitrary bounded convex box, then exhaust the domain. This
avoids a global integrability requirement on the joint function.

For full-support driving rho, the first time tau_R when any site of R
is filled therefore has the representation

    tau_R=min_(i in R, B contains i)
          [G_B^T xi_B]_i/[G_B^T rho_B]_i.                 (F37)

All denominators are positive. This is the infimum of times at which
one of the corresponding emissions becomes positive; filling occurs
at that same boundary time. It is a concave, increasing, homogeneous
function of xi. Its survival probability a_R(t rho) is log-concave,
so its hazard rate is nondecreasing. The maximum-of-Green-columns
formula (F25) describes full cover; the minimum in (F37) describes
first arrival at a target set. These are different stopping times.

### Sure-hitting basins determine the first-hit tail polynomial

**Proved; independent AI proof audit passed.** Keep full-support rho
and a nonempty target R. Write h_j(R)=P_j(hit R before the sink),
H_R=sum_j rho_j h_j(R), and

    D_R={j:h_j(R)=1},  B_0=V\D_R,
    Y=D_R\R,  m=|Y|.                                      (F38)

The set D_R is exactly the closure cl(R) from Section 34: its
complement consists of vertices that can reach the sink while avoiding
R. The first-hit time from (F37) has the sharp form

    P(tau_R>t)~C_R t^m exp(-H_R t),  C_R>0.                (F39)

Here is a finite positive formula for the constant, requiring no
cancellation of spectral coefficients. For B not meeting R, let
a_y(B) be its Poisson-IDLA birth rate at y under drive rho, as in
(F17), and let h_B=sum_(y outside B)a_y(B)=P_rho(hit V\B).
For an ordered list pi, write pi_<j for its preceding sites. Put

    alpha_R=sum_(pi permutation of B_0)
       product_j a_(pi_j)(pi_<j)/(h_(pi_<j)-H_R),
    beta_R=sum_(sigma permutation of Y)
       product_j a_(sigma_j)(B_0 union sigma_<j),
    C_R=alpha_R beta_R/m!.                                 (F40)

Empty products and their permutation sums equal one. Every numerator
is positive because rho is strictly positive. Every denominator is
positive: for B subset V\R, h_B=H_R exactly when B_0 subset B.
Indeed, for A=V\B containing R, hitting A is equivalent to hitting R
if A subset D_R; otherwise a source j in A\D_R gives the strict
difference H_A-H_R>=rho_j(1-h_j(R))>0.

To prove (F39), kill the increasing support chain when it first hits R.
Its holding rate at B is h_B. The states of smallest holding rate H_R
form precisely the interval B_0 subset B subset V\R. Along an
increasing path, the maximum number of these critical states is m+1;
it is attained only if the path reaches B_0 before any Y-site or R.
The exponentially tilted entrance mass at B_0 is alpha_R:
integrating a birth from each earlier state contributes
a_y(B)/(h_B-H_R), and summing its orderings gives (F40).
Inside the critical interval every holding exponential has rate H_R.
Summing the m births that add all of Y gives the leading factor
beta_R t^m/m!. Paths entering the interval above B_0, or ending below
V\R, use fewer critical states and give lower polynomial degree;
paths never entering it decay faster. The finite acyclic path
expansion therefore proves the asymptotic and its positive constant.

Let tau_R^disc be the first iid addition that fills R at some site.
If H_R<1, put lambda=1-H_R. The same triangular expansion gives

    P(tau_R^disc>N)~C_R lambda^(-m) N^m lambda^N.           (F41)

The coefficient relation also follows from the common finite spectral
expansion: Poissonizing (N)_m lambda^N gives
(lambda t)^m exp(-H_R t). Lower degrees and lower eigenvalues cannot
change its leading coefficient. If H_R=1 then D_R=V and no particle
can leave the graph before visiting R; hence tau_R^disc<=m+1.
In that case the continuous survival is exactly exp(-t) times a
polynomial of degree m, still with leading coefficient (F40).

On the existing single-leak cycle with rho=(1/2,1/3,1/6), target
R={1} has D_R={0,1}, B_0={2}, H_R=11/12 and m=1. Formula (F40)
gives alpha_R=2, beta_R=7/12 and C_R=7/6. For R={2}, every site
hits R surely, m=2, and the survival is explicitly
exp(-t)(1+5t/6+7t^2/24). Thus the same sink-search closure that
organizes the addition algebra measures the tail's polynomial degree.

### First filling as optimal stopping on a reversed walk

There is a self-contained stopping interpretation of (F34), (F37).
Fix rho>0 and a known capacity vector xi>=0, and define

    h=(I-K^T)^(-1)rho, D=diag(h), Q=D^(-1)K^T D,
    kappa_i=rho_i/h_i, c_i=xi_i/h_i.                       (F42)

Since h=rho+K^T h, Q is substochastic with killing probability
kappa_i>0 in every state. This cemetery reverses source injection;
it is not the original K-sink. At live state i, choose either to
retire immediately for zero, or pay c_i and take one Q-step. If that
step is killed, receive reward t. Retirement occurs before the cost
and killing trial; each continued visit incurs its cost. Costs are
known to the controller; future transitions are not.

For a policy pi started at v, let p_v(pi) be its probability of
intrinsic killing before retirement and C_v(pi) its expected total
cost. Arbitrary randomized history-dependent policies are allowed.
The uniform positive lower bound on kappa makes their number of
continued visits geometrically bounded. The value solves

    y_i(t)=max(0,t kappa_i-c_i+sum_j Q(i,j)y_j(t)),
    y_v(t)=sup_pi[t p_v(pi)-C_v(pi)]=u_v(t rho,xi)/h_v.     (F43)

The Bellman map contracts in sup norm by at most 1-min_i kappa_i.
Its inequalities bound every finite-horizon payoff by y; geometric
tails let the horizon tend to infinity. Continuing exactly on
{i:y_i(t)>0} attains equality. This proves optimality without an
external stopping theorem. The value is continuous and 1-Lipschitz
in t, since each policy has slope p_v(pi) in [0,1].

For the policy that continues in B containing v and retires outside,
its Green matrix is D_B^(-1)G_B^T D_B. Therefore

    C_v(B)=[G_B^T xi_B]_v/h_v,
    p_v(B)=[G_B^T rho_B]_v/h_v>0,
    tau_v=min_(B contains v) C_v(B)/p_v(B)
         =inf_(pi:p_v(pi)>0) C_v(pi)/p_v(pi).              (F44)

The first minimum is (F37). At t=tau_v the value at v is zero, so
every policy has tau_v p_v(pi)-C_v(pi)<=0, while a minimizing domain
attains equality. Thus randomization or history dependence cannot
improve the ratio. This is a retirement threshold for first filling,
not a statement about original-sink hitting probabilities.

One canonical minimizing domain is

    B_v^*={v} union {i:tau_i<tau_v}.                        (F45)

At t=tau_v, continuation equality holds on the strict positive-value
set. It also holds at v by letting s decrease to t from s>t, where
y_v(s)>0. The value is zero outside B_v^*. Restricting (F43) to this
domain and inverting I-Q_(B_v^*,B_v^*) gives
C_v(B_v^*)=tau_v p_v(B_v^*). Ties cause no problem: continue at v
and retire at other equal-threshold vertices. Under exponential
capacities these thresholds form exactly the IDLA filling process.
Full cover is max_v tau_v, the different maximum in (F25). The
stopping proof was independently supplied and reviewed; no ARW
occupation-order statement or literature novelty claim is inferred.

### A deterministic semigroup with an exponential residual law

The fluid construction is more than a one-time representation. For
w,b>=0, let U_b(w) solve u=(b-w+K^T u)_+ and put

    Phi_b(w)=w-b+M^T U_b(w).

These are remaining capacities, with their zero coordinates already
filled. Complementarity gives 0<=Phi_b(w)<=w: if u_i>0 the residual
is zero, and if u_i=0 its value is w_i-b_i-(K^T u)_i<=w_i.
For all nonnegative inputs b,c,

    Phi_(b+c)=Phi_c composed with Phi_b,
    U_(b+c)(w)=U_b(w)+U_c(Phi_b(w)).                       (F46)

Indeed the summed emissions solve the combined linear equation.
If either summand at a site is positive, the final residual is zero,
because residuals only decrease. The combined complementarity
conditions and uniqueness prove (F46), including boundary inputs.
In particular the deterministic maps commute.

Let R_exp(B,dw) put zero residual on B and independent Exp(1)
residuals on V\B. Let D_i be the IDLA addition matrices. Then

    R_exp Phi_b=P_IDLA(b) R_exp,
    P_IDLA(b)=exp(sum_i b_i(D_i-I)).                       (F47)

On a phase with zero set B the emission rate is G_B^T b_B, giving
the birth rates (F17) at vacant sites. If their residuals are E_y
and rates a_y, the first event is min_y E_y/a_y. For event time t
and winning site j its joint density with surviving residuals v_y is

    a_j exp(-t sum_y a_y)dt product_(y!=j)exp(-v_y)dv_y.

This proves independent exponential residuals after the event,
independently of its time and label. The same holds conditional on
no event during a fixed interval. Zero rates cause no difficulty,
and initial zero sites are already routing sites. Induction through
at most |V| births proves (F47) and the residual law conditional
on the observed support history.

For any prescribed deterministic nondecreasing cadlag cumulative
input path b(t), starting at zero, (F46)--(F47) give

    W_t=Phi_(b(t))(W_0),  B_t={i:(W_t)_i=0},
    law(W_t | support history through t)=R_exp(B_t,.),
    P(B_t=C | support history through s)
       =P_IDLA(b(t)-b(s))(B_s,C), s<=t.                   (F48)

First apply the kernel identity on finite partitions; rational
observations and right continuity recover the full support history.
Thus this is IDLA driven by independent Poisson counts with cumulative
intensities b_i(t). A jump of b means a completely stabilized Poisson
batch, and may fill multiple sites at that time. No absolute
continuity of b is needed. This proof passed an independent audit.
It gives the process law along every prescribed path, not equality
with one original arrow-stack construction over all incomparable
input paths.

The ARW basis also fits this realization, with a precise sign
restriction. Define the mass-one signed face measures
Sigma_r=E_r R_exp. The addition intertwining gives

    Sigma_r Phi_b=A_r(b) Sigma_r,
    A_r(b)=exp(sum_i b_i(A_i^r-I)).                         (F49)

Faces with different zero sets are disjoint. Hence, for a stable-law
row mu, the measure mu Sigma_r is nonnegative exactly when
mu E_r>=0. This is exactly the cone of mixtures of binary all-active
stabilization laws: mu=(mu E_r)P_r, with the coefficients summing to
one. It contains every all-active input law by surplus IDLA and is
preserved by additions. Individual stable-state embeddings need not
be positive. On one site jumping directly to the sink,
Sigma_r({i})=(1+r)delta_0-r Exp(1). This describes the sign obstruction
in this particular exponential-face realization; it does not rule
out other positive representations.

There is a sharper limitation if the residual flow itself is kept.
For the direct-sink one-site model, write phi_b(w)=(w-b)_+ and
pi=(r,1)/(1+r). Every pair of mass-one finite signed measures with
H phi_b=exp(b(A-I))H has the form

    H_0=(1-a)delta_0+a Exp(1),
    H_1=(1+ra)delta_0-ra Exp(1), a real.                   (F50)

Indeed pi H is phi-invariant; its positive-axis tails must be
constant under translation and tend to zero, so pi H=delta_0.
The intertwining then gives H_j((b,infinity))=e^(-b)H_j((0,infinity)),
which determines the two displayed measures. Conversely these measures
intertwine because Exp(1) phi_b=e^(-b)Exp(1)+(1-e^(-b))delta_0.
For r>0 both rows are positive only when a=0, so both equal delta_0.
Even a positive embedding in an enlarged space projecting to phi must
put both state measures on the zero-residual fiber. If that fiber is
pointwise fixed, the two measures coincide. Thus a useful positive
extension must change the residual dynamics or carry evolving hidden
information on that fiber. This classification and (F49) passed an
independent proof audit; they do not exclude other representations.


## 36. Positive duality with complementary sleep parameters

**Proved, 2026-10-19; independent AI proof audit passed for (D1)--(D10).** This uses the
ordinary addition algebra (A1) and its coefficient functional Lambda.
It gives a nonnegative Markov duality and a stochastic intertwining with
a time-reversed chain. It does not identify a monotone dual representation
of the single-target indicator.

Write J_i=k_i+sum_j K(i,j)z_j and v_i=z_i-J_i, so z_i v_i=0.
For an odds vector r, let w_i^r=z_i+r_i v_i and w_C^r=product_(i in C)w_i^r.
These are the stable-state basis elements from the complete-block inverse.
Let

    F(C,D)=Lambda(z_C z_D)
          =P(IDLA from 1_C+1_D fills V),
    H_(r,s)(C,D)=Lambda(w_C^r w_D^s)
                =(E_r F E_s^T)(C,D).                        (D1)

The Frobenius pairing makes multiplication by every z_i self-adjoint.
In the corresponding stable-state bases this gives

    A_i^r H_(r,s)=H_(r,s)(A_i^s)^T,     every i.             (D2)

Both E's and F are invertible, so H is invertible as well. This is
duality for discrete additions; any common driving mixture or continuous
addition clocks obey the corresponding semigroup identity.

### Positivity from two-particle releases

Suppose r_i s_i<=1 for every i. Equivalently, the two sleep vectors
p_i=1/(1+r_i) and p'_i=1/(1+s_i) satisfy p_i+p'_i>=1. Then

    H_(r,s)(C,D)>=0 for every C,D.                           (D3)

For the proof put I=C intersection D and S=C symmetric-difference D.
The crucial local calculation is

    w_i^r w_i^s=z_i^2+r_i s_i v_i^2
       =(1-r_i s_i)z_i J_i+r_i s_i J_i^2.                  (D4)

Here z_i^2=z_i J_i and v_i^2=J_i^2-z_i J_i follow directly from
z_i(z_i-J_i)=0. Thus each overlap factor is a polynomial with
nonnegative coefficients in the z's, under the stated hypothesis.
The factors at nonoverlap sites combine into w_S^t, with t_i=r_i
on C\D and t_i=s_i on D\C; choose arbitrary positive t elsewhere.
For every active input eta,

    Lambda(w_S^t z^eta)
       =Z_V(t) P_t(full occupation from sleeping S plus active eta)>=0.
                                                                  (D5)

Indeed multiplication by z^eta is physical active addition in the
stable w^t basis, and Lambda(w_B^t)=Z_V(t)1_(B=V). All stabilizations
here terminate on the finite sink-accessible graph. Expanding the
overlap polynomials in (D4) and applying (D5) proves (D3).

There is a direct experiment behind the expansion. At each i in I,
with probability 1-r_i s_i retain one active particle at i and give
the other one independent K-jump; with the remaining probability give
both particles independent K-jumps. Sink destinations remove those
particles. Start with sleepers on S and stabilize the surviving active
particles at odds t. Its full-occupation probability times Z_V(t)
is exactly H_(r,s)(C,D). At complementary sleep probabilities,
s_i=1/r_i, every overlap releases both particles by a forced jump.
The use of sleeping inputs here only proves positivity of this auxiliary
kernel; no sleeping-input occupation monotonicity is assumed.

### Stationary normalization and reversal

Let pi_r=delta_V P_r, stationary for every A_i^r. Since the full IDLA
support stays full, F(V,D)=F(C,V)=1. Also E_r 1=1 and pi_r E_r=delta_V.
Consequently, for all finite positive r,s,

    pi_r H_(r,s)=1^T,       H_(r,s) pi_s^T=1.                (D6)

Under r_i s_i<=1 the matrix

    B_(r,s)=H_(r,s) diag(pi_s)                               (D7)

is therefore stochastic, and pi_r B_(r,s)=pi_s. Equivalently
pi_r(C)H_(r,s)(C,D)pi_s(D) is a coupling of the two stationary laws.
Each pi_s(D)>0: from one active particle at every site, route particles
outside D to the sink while retaining the distinguished particles of
D active, then use a sleep instruction at every site of D. A visiting
particle can be moved through a retained active site by a legal toppling.
This prescribes only finitely many positive-probability instructions.

Thus the time-reversed addition matrix

    A_i^(s,*)=diag(pi_s)^(-1)(A_i^s)^T diag(pi_s)

is well defined and stochastic. Equation (D2) becomes the positive
intertwining

    A_i^r B_(r,s)=B_(r,s) A_i^(s,*),       every i.           (D8)

In particular every allowed sleep vector has such a link to the reversal
of the chain at its complementary sleep probabilities. Self-duality gives
a nonnegative H_(r,r) when every p_i>=1/2. The sign condition is not
vacuous: on one site jumping directly to the sink,
H_(r,s)({i},{i})=1-r_i s_i, negative when p_i+p'_i<1.

Moreover B_(r,s)(empty,V)=1, since Lambda(w_D^s) vanishes unless
D=V, where it is Z_V(s), and pi_s(V)=1/Z_V(s). Therefore, for any
fixed input word, the forward r chain from empty, followed by B_(r,s),
has exactly the law of the reversed s chain from full under that word.
Also H_(r,s)(C,D)=0 whenever |C|+|D|<|V|, by degree. The stationary
coupling consequently satisfies |C|+|D|>=|V| almost surely.

The exact certificate uses the existing single-leak three-cycle with
p=(1/3,2/5,3/5) and p'=(2/3,3/5,2/5). `duality_audit()` in
`verify_addition_algebra.py` verifies all three addition dualities,
stationary normalizations, and stochastic reversal intertwinings.
Independently, for C={0,1}, D={1,2}, the overlap releases two particles
from 1, both jumping to 2. With sleepers {0,2} and stabilization sleep
vector (1/3,1/2,2/5), a 43-state raw absorption calculation gives full
probability 3/32. Its Z_V is 8, yielding the independently checked
pairing entry H(C,D)=3/4. No new graph or parameter sweep was used.

### The complementary link has a reachability support criterion

For s_i=1/r_i, let I=C intersection D, S=C symmetric-difference D,
and T=V\S. Make two tokens at every i in I, each required to take an
initial K-jump. Then

    H_(r,1/r)(C,D)>0 iff T matches injectively to these tokens,
    with eligibility given by a positive-length directed path
    from the token's original vertex to the matched target. (D9)

Paths cannot pass through the absorbing sink. A return to the original
vertex therefore requires a directed cycle, possibly a self-loop.
To justify (D9), first fix the forced destinations. For arbitrary
sleeping S and active eta, full ARW occupation is possible iff IDLA
parking on V\S can fill V\S from eta. Run the legal surplus-IDLA phase
first, treating S as already occupied. If it leaves fewer than |V|
particles, full occupation is impossible in the continuation. If it
fills V, sleeping all its remaining active singletons has positive
probability. Thus the positive-filling matching theorem of Section 32
applies. Allowing each positive-probability initial jump gives exactly
the positive-length reachability in (D9). Unmatched tokens may exit;
sink accessibility supplies the needed paths.

In particular on a DAG with no self-loops the complementary link is
the deterministic complement map:

    B_(r,1/r)(C,D)=1_(D=V\C).                              (D10)

If I is nonempty, choose an earliest vertex of I in a topological
ordering. No token from I has a positive-length path to it, so (D9)
fails. If I is empty, degree requires S=V. Stochasticity then fixes
the remaining entry to one. Thus off-complement randomness in this
link requires cycles; the general link is not asserted to impose
pointwise complementarity or an inclusion order.

### A sharp obstruction for scalar target payoffs

**Proved; independent AI proof audit passed for (D11)--(D12).**
Suppose the interior graph is strongly connected and n>=2. For every
pair of finite positive odds with r_i s_i<=1,

    B_(r,s)(C,D)>0 iff |C|+|D|>=n.                          (D11)

For complementary odds this follows from (D9): every vertex reaches
every vertex by a positive-length path, including itself. If I is
nonempty, all tokens are eligible for all targets and the matching
condition is just |T|<=2|I|, equivalent to the displayed inequality.
If I is empty, that inequality forces S=V. For general allowed r,s,
the event that every overlap takes the two-forced-jumps option in
(D4) has positive probability. Its support is unchanged by the
positive hybrid sleep probabilities. This gives positivity in the
same range; degree gives the converse.

For a nonempty proper R subset V, let f_R(C)=1_(R subset C). There is
no nonnegative scalar function g on dual supports satisfying

    B_(r,s)g=f_R.                                          (D12)

Choose v in R. The zero value at C=V\{v}, together with (D11), forces
g(D)=0 for every nonempty D. But C=R is proper, so B(R,empty)=0;
the equation at R could not give one. The obstruction also applies
if the proposed representation is required only after **all-active**
inputs: their binary stabilization rows form the invertible matrix
P_r, and equality on those rows implies (D12). For full occupation
the positive payoff does exist: g(D)=1_(D=empty)/B(V,empty).
Conversely, no nonzero nonnegative payoff represents an observable
which vanishes at full support, since B(V,D)>0 for every D.

This also rules out any inclusion-increasing payoff in vacancy
coordinates. Indeed L(C,E)=B(C,V\E) has L(empty,empty)=1; if Lh=f_R,
then h(empty)=0, so inclusion-increase would force h>=0 and contradict
(D12). Adding a constant to a signed payoff can make it positive but
introduces an additive subtraction in the observable and does not
repair a failed order. In the existing three-cycle certificate, for
the target f(C)=1_(0 in C), its unique transported
payoff g=L^(-1)f satisfies

    g({1})=g({2})=-29/598,
    g({1,2})=-65493339/194915225<g({1}).                    (D13)

Thus even inclusion-increase of the transported payoff fails; adding
a constant does not repair that. The exact certificate now checks
these values and the full support criterion on the same fixed cycle.

The duality retains source addition exactly and may still support a
reversed-history construction with more information than a terminal
scalar payoff. That application and the occupation conjecture remain
open; (D12) is an obstruction to this representation, not to the target
inequality.

### Compatibility with parameter transport

Put Pi_r=diag(pi_r). For a<=b, the Bayes reverse of the positive
sleep transport is

    T_hat_(b,a)=Pi_b^(-1) T_(a,b)^T Pi_a.                  (D14)

It is stochastic and sends pi_b to pi_a, since pi_a T_(a,b)=pi_b.
It cannot decrease particle count; the forward wake-and-stabilize
transport cannot increase it. For a<=b<=c the reversed kernels
satisfy T_hat_(c,b) T_hat_(b,a)=T_hat_(c,a).

If r<=t, s<=u and t_i u_i<=1, then all factors in

    B_(r,s)=T_(r,t) B_(t,u) T_hat_(u,s)                    (D15)

are stochastic. This is a direct cancellation of P_t E_t and
E_u^T P_u^T in E_r P_t E_t F E_u^T Pi_u
Pi_u^(-1) P_u^T E_s^T Pi_s. The identity itself holds for arbitrary
positive endpoints; the inequalities ensure probability kernels.
For C_r=B_(r,1/r) it gives, whenever r<=t,

    T_(r,t) C_t=C_r T_hat_(1/r,1/t),                       (D16)

both sides equal to B_(r,1/t). More generally every allowed B_(r,s)
factors through a complementary link C_t with r<=t<=1/s.

Symmetry of F also gives the Bayes reciprocity
Pi_s^(-1)B_(r,s)^T Pi_r=B_(s,r). Thus the stationary coupling at
(r,s) can be sampled from the one at (t,u) by conditionally independent
Bayes transport of its two coordinates. This is consistent along
ordered parameter chains; no simultaneous construction over every
incomparable parameter vector is asserted.

For the infinitesimal formula let G_x(s)=E_s partial_x P_s and
partial_x pi_s=pi_s G_x(s), as in Section 29. Then

    G_hat_x(s)=Pi_s^(-1)G_x(s)^T Pi_s
       -diag((pi_s G_x(s))(D)/pi_s(D)),
    partial_(r_x)B=-G_x(r)B,
    partial_(s_x)B=-B G_hat_x(s).                          (D17)

The reverse generator has nonnegative off-diagonals and zero row
sums. Its diagonal correction is needed because pi_s varies with s.
These algebraic formulas were independently derived and reviewed.

The links B themselves are not a cocycle. For nonempty V and positive
parameters for which all three links are stochastic,

    (B_(r,s)B_(s,t))(empty,empty)
       =Z_V(s)pi_t(empty)>0=B_(r,t)(empty,empty).           (D18)

The correct compatibility is (D15), involving forward transport on
one side and its Bayes reverse on the other. A Bayes step has the
finite rejection sampler C~pi_a, D'~T_(a,b)(C,.), accepted if D'=D;
the success probability is pi_b(D)>0. This proves existence of a
sampler but gives no useful efficiency estimate.

### Reachability order of the vacancy link

For complementary parameters write L_r(C,E)=B_(r,1/r)(C,V\E).
Say E<=_reach C if E can be matched injectively to distinct sites of C,
each matched source reaching its output by an interior directed path
of length at least zero. Then

    L_r(C,E)>0 implies E<=_reach C.                        (D19)

Equivalently, |E intersection U|<=|C intersection U| for every
ancestor-closed U. Such sets are ancestor-closed unions of strongly
connected components. To prove this put X=C\E,Y=E\C. Criterion (D9)
requires X union Y to match to two tokens at every x in X. A target
in U can only use a source in U, so |X intersection U|+
|Y intersection U|<=2|X intersection U|. For any A subset Y take
U=Anc(A). Then |A|<=|Y intersection U|<=|X intersection U|,
the last term being exactly the number of eligible removed sources
for A. Hall's theorem matches Y into X; match common sites to
themselves. The converse between reachability matching and the
ancestor-closed inequalities is the same Hall argument.

The communicating classes of L_r are exactly configurations with
the same particle count in every strongly connected component.
Every entry within such a class is positive. Counts on all
ancestor-closed component unions can only decrease along transitions;
equality on topological initial unions determines each individual
component count. This proves necessity for mutual communication.
For sufficiency, equal component counts give |X intersection Q|=
|Y intersection Q| in each component Q. If this is positive then Q
has at least two vertices and its 2|X intersection Q| tokens can
match all targets in that component by positive-length paths.
Combining these matchings proves (D9). Thus L_r is block triangular
in the induced partial order, with strictly positive diagonal blocks.

This is an order on a sampled output relative to its own input, not
monotonicity of the kernel as its input varies. It permits relocation:
on the single-leak three-cycle, L_r({2},{0})>0, using paths 2->0 and
2->0->1->2 for the two required targets. In particular the kernel
cannot be realized by deletion alone. Nor can ancestor-closed sets
be replaced by all dominator regions: B_0={0} on this graph, and
its count increases in this transition. The structural proof was
independently supplied and reviewed; it identifies support, not a
local transition sampler.

Iterating L_r defines a separate artificial chain; it is not the ARW
addition dynamics. Its absorbing states are exactly subsets of the
vertices lying on no interior directed cycle, counting self-loops.
For any occupied cyclic vertex i, (D9) gives
L_r(C,C\{i})>0 by one positive-length return. Repeatedly choosing such
deletions reaches an acyclic subset in at most |V| steps. Conversely,
if C is acyclic and a transition changes it, choose the earliest
component of X=C\E in the condensation order. Its required target
cannot be reached by any token from X, contradicting (D9). Thus it
is absorbing. Finiteness and a positive minimum over the chosen
deletion-path probabilities give a geometric absorption tail for
fixed K,r. The limiting acyclic support need not be C's original
acyclic part, since relocation may populate new sites.

If every vertex is cyclic, only empty is absorbing. With Q the
restriction of L_r to nonempty supports, its Green matrix
G=(I-Q)^(-1)=sum_(j>=0)Q^j is finite and nonnegative. Thus

    h_v(C)=E_C sum_(j<tau_empty)1_(v in C_j),
    h_v(empty)=0,   (I-L_r)h_v=1_(v occupied).              (D20)

This gives a positive history potential. It preserves order in the
reward and is excessive, L_r h_v<=h_v, but need not increase with
the starting support. On the maintained three-cycle, at target 0,

    h_0({0})=189979/74354,
    h_0({0,1})=40836185505115661/21542872302627536
              <h_0({0}).                                 (D21)

The same exact duality certificate checks the nonnegative Green
inverse and (D20), as well as this inequality. The general absorption
and potential arguments were independently supplied and reviewed.
On a reducible graph with nonempty absorbing supports, killing on
their union gives the analogous Green potential, but its right-hand
reward must vanish on that absorbing boundary.


## 37. Zero-loss ARW currents are the rejected loop soup

**Proved constructions, with independent AI audits of the current laws,
normal forms, coherent activation couplings, and block exchanges.**
The separate graphical suffix formula is marked locally pending audit.
This connects the earlier
periodic-loop construction of Sections 19--20 to actual ARW currents.
The first statements concern jump-count laws. The ordered-book
factorization and the subsequent coherent coupling across activation
sets have each passed separate proof audits.

For internal edge variables t=(t_ij), put W_ij=K(i,j)t_ij and
d_C(W)=det(I-W_CC), with d_empty=1. Extend the normalizer by

    Z(r;W)=d_V(W) sum_(C subset V) r_(V\C)/d_C(W).         (J1)

Here r_A=product_(i in A)r_i and p_i=1/(1+r_i). For a nonzero
internal edge-count vector nu, Sections 19--20 give coefficients

    ell_nu=[t^nu] log det(I-W)^(-1),
    c_nu=[t^nu] log g_p(t),
    a_nu=ell_nu-c_nu>=0,
    -log(p_V Z(r;W))=sum_(nu!=0)a_nu t^nu.                (J2)

The ordinary Markov-loop intensities are ell_nu; c_nu are the
intensities accepted by the periodically stabilizing-book test. Thus
a_nu are precisely the rejected intensities. They decrease as the
sleep probabilities increase and depend only on the vertices and
edges used by nu. Their sum is finite.

Start ARW with one active particle at every vertex, and condition on
final full occupation. This event allows no sink loss. If N is its
internal jump-count vector, then

    E[t^N | full]=Z(r;K)/Z(r;W)
       =exp(sum_(nu!=0)a_nu(t^nu-1)).                     (J3)

For the first identity, initially take strictly positive edge weights
near one with spectral radius(W)<1. Set h=(I-W)^(-1)k, where k is
the original sink vector. The weighted harmonic transform in Section
33 replaces W by diag(h)^(-1)W diag(h), with sink jumps k_i/h_i,
and preserves the sleep probabilities. The telescoping weight from
the full initial and final states is h_V/h_V=1. Principal minors
are invariant under the diagonal conjugation, so the weighted full
probability is 1/Z(r;W). Divide by its unweighted value 1/Z(r;K).
Continuity extends the generating-function identity to 0<=t_ij<=1.
Equation (J2) gives the second equality.

Consequently N is a finite compound-Poisson sum of rejected ordinary
Markov loops. Every atom is balanced and has connected Eulerian
support. This agrees with the independent normalization check

    P(N=0 | full)=p_V Z(r;K)
                  =exp(-sum_nu a_nu).                    (J4)

The no-jump event simply sleeps all n initial actives. In particular
the conditional current is infinitely divisible, and a common
ordinary loop soup with retention uniforms couples these currents
coordinatewise decreasingly as the sleep probabilities increase.
This comparison is conditional on full occupation.

There is a direct consequence for any finite all-active eta having
positive full-filling probability. Run surplus IDLA first, retaining
one active at each occupied site. Full ARW occupation is possible
only if this phase fills V; its fresh continuation is binary full
ARW and must suffer zero loss. Conditioning on that continuation's
full event multiplies every successful IDLA prefix by the same
factor 1/Z(r;K). Therefore the conditional current law splits into
independent parts:

    N_ARW | full =_law N_IDLA | full + N_rejected.          (J5)

All labeled sink exits belong to the IDLA part; the second part has
law (J3). This holds jointly with the input configuration under any
fixed all-active prior independent of the books. It strengthens the
posterior cancellation following (A24), but does not make the
IDLA-conditioned current itself infinitely divisible.

The creator's accepted-loop current has generating function g_p(t)/g_p(1).
Multiplying it by (J3) gives

    [g_p(t)/g_p(1)] [Z(r;K)/Z(r;W)]
       =det(I-K)/det(I-W).                                (J6)

Thus independent creator currents and conditional binary-full ARW
currents add to an ordinary Markov-loop current. This is an exact
decomposition in law of aggregate counts, not a concatenation claim
for arbitrary local word orderings.

### Full mixed activity gives increments of the same soup

Start with one particle at every site, sleeping on S and active on
T=V\S. Let r^S_i=r_i on S and zero on T; equivalently p^S_i=p_i
on S and one on T. Then

    P_r(full | sleeping S, active T)=Z(r^S;K)/Z(r;K).       (J7)

Indeed the full-occupation functional is
Lambda(w_S^r z_T)/Z(r;K). Expanding w_S and using the opposite-facet
pairing from (A10) gives
Lambda(w_S^r z_T)=d_V sum_(A subset S)r_A/d_(V\A)=Z(r^S;K).
Lower-degree terms cannot contribute to Lambda. The same weighted
harmonic transform applies because the initial and final particle
counts are again one per site. Thus the conditional internal current
has generating function

    E[t^N | full, S]=
       exp(sum_nu [a_nu(p)-a_nu(p^S)](t^nu-1)).            (J8)

All intensities are nonnegative by the monotonicity in (J2); the
auxiliary boundary probabilities p^S_i=1 are obtained by continuity.
Every atom intersects T, by locality. Activating more sites of this
full-support initial configuration increases the conditional current
under a common soup coupling: S shrinks, p^S increases and a(p^S)
decreases. At S=V the current is zero; at S=empty this is (J3).
The no-jump probability is p_T Z(r;K)/Z(r^S;K), as both formulas require.

These are positive loop descriptions of **zero-loss full-support
bridges**. They do not order arbitrary sleeping-input outcomes.

### A normal form for ordered zero-loss bridge records

**Proved; independent AI proof audit passed.** Use the buffer records
and least-action statement of Section 11. Let H be the monoid of
finite balanced internal bit-arrow books from Section 20, with
multiplicative weight wt. Let M contain empty and the books whose
first bit at every used vertex is one. Define R_* to consist of
books with no nonempty right divisor in M. Then there is a unique
factorization

    H=R_* M.                                               (J9)

For a direct proof, take h in H and append a terminal bit one at
each site after its finite word. Run binary full-input ARW on this
book. The whole vector of h's word lengths is a stabilizing upper
bound: it is balanced and the following bits are all one. Least
action therefore forces a finite stabilization using a prefix a of h.
All its arrows are internal, so it loses no particles and ends full;
its emitted prefix a is balanced. The remaining balanced suffix b
has first bits one at every used site, hence b is in M.

If a had a nonempty right divisor c in M, cutting just before c
would give an earlier full stable configuration. Thus a is in R_*.
Conversely, suppose h=a b with a in R_* and b in M. Cutting after
a is stable, so the actual prefix a' is at most a. Both are balanced.
If a'=a is false, write a=a'c; stability at a' makes the first bits
of the nonempty balanced c all one. This contradicts a in R_*.
Thus a is exactly the physical emitted prefix, proving uniqueness.

The weight series of H is f=1/d_V, and that of M is g=p_V Z/d_V.
Unique factorization and multiplicative weights imply that an
ordinary Boltzmann book h, of law wt(h)/f, splits into independent
factors with laws

    P(a)=p_V Z wt(a), a in R_*,
    P(b)=wt(b)/g, b in M.                                  (J10)

The first is exactly the conditional full ARW emission-book law in
this buffer realization.
Indeed each prescribed a in R_* occurs together with its terminal
sleep bits with iid-book probability wt(a)p_V; the full event has
probability 1/Z. Thus (J6) lifts to ordered books, in the order
**physical bridge prefix followed by creator prefix**. This does not
say that the physical prefix is obtained by merely deleting accepted
Lyndon factors from a chosen ordinary loop ordering.

An explicit sampler follows. Sample an ordinary cycle-popping balanced
arrow book, whose law is d_V times its arrow weight, and attach iid
Bernoulli(p_i) bits to all its records. Append terminal ones and
stabilize binary full input; retain only the emitted prefix a.
To complete a buffer book with the actual conditional law, discard the
remaining suffix, retain the terminal sleep bits, and independently
sample the unused terminal arrows and later iid records. The latter
were not constrained by the emitted-prefix event. This samples the
full zero-loss bridge without rejecting on the rare full-occupation
event. Its number of emitted arrows is bounded by the sampled
ordinary cycle-popping book length.

The language R_* is also closed under concatenation. For a,b in R_*,
replacing the terminal ones after a by the first bits of b can only
decrease sleep bits at the old terminal indices. The bit monotonicity
from Section 11 forces the new odometer to be at least a. By balance
of a, a remaining stabilizing vector is then a plus a stabilizing
vector for binary full input on b. Minimality for b forces it to be
at least a+b, which is already a stable upper bound. Thus the
concatenation ab is the actual emitted book. No trace-monoid presentation or unique
factorization into indecomposable R_* books is asserted here.

For mixed initial activity define H_S by requiring the first bit to
be one only at used sites in S. The same factorization restricts to

    H_S=(R_* intersection H_S) M.                          (J11)

The first factor inherits the prefix constraints; conversely an
M-suffix preserves them at every site not used by the prefix. Its
physical cylinder probability is p_T wt(a), because initial sleep
bits in S are fixed rather than sampled. Its normalized law is
therefore exactly the mixed zero-loss bridge law from (J7).
A finite sampler for H_S is the ordinary Boltzmann-book sampler
conditioned on its used first bits in S being one; its acceptance
probability is p_S Z(r^S;K)>0. Factor it as above. No efficiency
bound for this rejection step is claimed.

The two factors may not simply be reversed while preserving
the ordinary-book law. Already for one vertex with a positive self-loop,
the bit words 01 and
10 both map to 10 if their R_* and M factors are interchanged.
The current decomposition forgets that cut information. A joint
history exchange would have to retain or restore it with its exact
normalization.

### Coupling activation sets with common terminal records

**Proved; independent AI proof audit passed.** For S subset V, let
H_S have its preceding meaning and let N_S be the books with no
nonempty right divisor in H_S. Every h in H has a greatest right
H_S divisor rho_S(h), and

    H=N_S H_S uniquely.                                   (J12)

Use the trace-monoid realization of H from Section 20. Right divisors
of a finite heap are its upper order ideals. Their join is their union;
at each vertex its local word is the longer of the two suffix words.
Thus the join of two right H_S divisors remains in H_S: its first
used bit at a constrained site is the first bit of one of those
suffixes, hence one. There are finitely many divisors, proving
existence of the greatest one. If h=u b with b greatest, any nonempty
right H_S divisor c of u would make c b a larger one. Conversely,
in a proposed factorization with u in N_S, the greatest divisor b'
contains b, so b'=c b with c a right divisor of u. Since H_S is
closed under taking balanced left divisors, c is in H_S and must be
empty. This proves uniqueness.

Write R_S=R_* intersection H_S. If S_2 subset S_1, then

    R_(S_2)=(N_(S_1) intersection H_(S_2)) R_(S_1)
    with a unique factorization.                          (J13)

Indeed take the greatest H_(S_1) suffix b of a in R_(S_2).
Right factors of R_* remain in R_*, so b is in R_(S_1). The prefix
is in N_(S_1) and inherits H_(S_2)'s first-bit constraints. Conversely
both factors belong to R_* and H_(S_2), each a submonoid, so their
product is in R_(S_2). Uniqueness follows from (J12).

Multiplicative weights make the two factors independent under the
Boltzmann law on R_(S_2). In particular the suffix has precisely the
Boltzmann law on R_(S_1), which (J11) identifies with the conditional
mixed full-support bridge with sleeping set S_1. Therefore one sample
a from the all-active full-conditioned bridge yields, simultaneously
for every S,

    a_S=rho_S(a),
    law(a_S)=ordered local emission book conditioned full,
             with initial sleepers S and actives V\S.      (J14)

For S_2 subset S_1, the local words of a_(S_1) are suffixes of those
of a_(S_2). These projections are coherent:
rho_(S_1)(rho_(S_2)(a))=rho_(S_1)(a), because every H_(S_1) divisor
is an H_(S_2) divisor. Thus activating more sites prepends balanced
records, leaving the terminal portions intact. Append the same
independent unused terminal arrows and later iid records to every
book, with terminal sleep bits one, to obtain a coupling with common
future tails as well. The claim concerns ordered local records, not
a common global chronological toppling order.

For one new activation x, put S_1=S_2 union {x}. Every nonempty
connected component of the extra prefix in (J13) must contain x.
Otherwise that component obeys H_(S_1)'s first-bit constraints and,
being disjoint from the other components, is a commuting right
H_(S_1) divisor, forbidden by N_(S_1). Hence a nonempty extra prefix
is connected and contains x. This is a statement about its record
graph, not a claim that its sampled law is independent of K or p.

The focused book certificate uses the existing two-site kernel with
a self-loop at 0, and the single balanced table with one edge each
0->0, 0->1, 1->0. It checks all 16 decorated local books, their
unique physical-prefix/creator-suffix factorizations, and 64 greatest
suffix projections with nested coherence. On the same table, starting
asleep at 0 and active at 1, the local orders (0,1) and (1,0) at
vertex 0 have conditional probabilities 3/5 and 2/5 when p=(1/3,2/5).
Uniformly shuffling the current table would therefore lose the bridge
law. This finite check supplements the general proof.

Every marginal in (J14) is conditioned on final full occupation.
The construction does not compare the probabilities of distinct
terminal supports or establish the unconditioned target conjecture.

### Exact exchanges of activation blocks

The partition functions in the preceding factorization are explicit:

    z_S=sum_(a in R_S)wt(a)=Z(r^S)/(p_(V\S)Z(r)),
    kappa_(S_2,S_1)
       =Z(r^(S_2))/(p_(S_1\S_2)Z(r^(S_1))), S_2 subset S_1. (J15)

The second is the extra-prefix normalizer in (J13), and is at least
one since the empty book contributes one. To check the first directly,
summing bits on a balanced arrow book of support A in H_S gives
p_(A intersection S). Expanding
product_(i in S)[p_i+(1-p_i)1_(i outside A)] and using total arrow-book
weight 1/d_(V\B) when B is avoided gives
sum_(h in H_S)wt(h)=p_S Z(r^S)/d_V. Divide by g=p_V Z(r)/d_V.
Product weights and (J13) now give independent successive prefixes
along every fixed activation chain. Independence of increments for
incomparable sets is not asserted.

For x in S define U_x(S)=N_S intersection H_(S\{x}) and
kappa_x(S)=Z(r^(S\{x}))/(p_x Z(r^S)). Activating x and then y gives
two independent blocks u in U_y(S\{x}) and v in U_x(S), in book order
u v. Multiplication bijects this product space with
N_S intersection H_(S\{x,y}). Refactor the same book w=u v by

    v'=rho_(S\{y})(w),  w=u'v',
    u' in U_x(S\{y}), v' in U_y(S).                       (J16)

This is an invertible exchange between the two activation orders.
It preserves every ordered paired record and its total weight. The
independent product laws match exactly because

    kappa_y(S\{x})kappa_x(S)
       =kappa_x(S\{y})kappa_y(S)
       =Z(r^(S\{x,y}))/(p_x p_y Z(r^S)).

For the multiplication bijection, extract the greatest H_(S\{x})
suffix of w. It inherits N_S by right-factor closure; the remaining
prefix belongs to N_(S\{x}) and inherits H_(S\{x,y}). Conversely
both products lie in the required intersection, since H_S and N_S
are submonoids. For completeness, N_S closure follows directly from
heaps: if a,b have no right H_S divisor and d is such a divisor of
ab, intersect d's upper ideal with the lower ideal for a. That
intersection is a right divisor of a and a left divisor of d, hence
lies in H_S and must be empty. Then d lies entirely in b, also
impossible. No trace presentation of N_S is used.

An old suffix in H_S can be kept unchanged in (J16); normal-form
uniqueness gives the same refactorization across its boundary. Iterating
over distinct activation labels gives a unique factorization for each
tagged activation order. Therefore the adjacent exchange maps obey

    sigma_j^2=id,  sigma_i sigma_j=sigma_j sigma_i (|i-j|>1),
    sigma_j sigma_(j+1) sigma_j
       =sigma_(j+1) sigma_j sigma_(j+1).                  (J17)

These are maps between activation-order-tagged factor spaces, with
the sleeping base updated at every swap. Both sides preserve the
same complete book and end in its unique factorization for the new
order. The product normalizer is Z(r^(S\D))/(p_D Z(r^S)), independent
of the order listing the activated set D. This proof and all the
normalizations passed an independent audit. The maps exchange block
ownership; they do not identify block order with a physical chronology.

For one activation x, the extra prefix is empty iff the more-active
book's first x-bit is one, with its appended terminal one used when
x is otherwise unused. Its probability is 1/kappa_x(S). A nonempty
prefix is connected, contains x, and has first x-bit zero.

### A finite reachability graph computes all constrained suffixes

**Proved; independent AI proof audit passed.** Represent a balanced
book h as a heap of decorated simple-cycle occurrences. Form a graph
on those occurrences and one cemetery vertex. At each site i, direct
an unconditional edge from each occurrence using i to its next such
occurrence. For i in the sleeping constraint set S, also add an edge
from every occurrence whose i-bit is zero to its preceding occurrence
at i, or to the cemetery if there is no preceding occurrence. Then

    rho_S(h) consists exactly of the occurrences
    with no directed path to the cemetery.                (J18)

An upper heap ideal is a right divisor precisely when it is closed
under the unconditional edges. Its first selected i-bit is one
precisely when each selected zero-bit occurrence at i forces its
immediate predecessor to be selected; a first zero occurrence is
forbidden. Thus an admissible H_S suffix is an outward-closed set
avoiding the cemetery. The largest such set is exactly the set of
vertices not reaching the cemetery, proving (J18).

Activating i removes the i-labeled backward edges, so retained
occurrences can only increase. For the all-active bridge a in R_*,
all occurrences reach the cemetery when S=V. Applied to an ordinary
book h at S=V, the doomed occurrences form its physical bridge prefix
and the surviving ones its creator suffix. This graph encodes local
record order; it is not the original routing graph. The existing
16-book certificate compares this reachability construction with
exhaustive balanced-cut suffixes in 128 cases, including the ordinary
books and their physical prefixes.

Keep occurrence identities and edge labels, even for repeated cycles
or parallel constraints. Reconstruct the book in its original heap
order: the constraint graph may have directed cycles. If applying the
rule to an extracted suffix, rebuild its local chains; a zero bit that
becomes first needs a cemetery edge. An induced subgraph alone can miss
this constraint.

### Sleep coupling, bounded currents, and the first bits

Several further consequences passed a separate proof audit. Sample
one ordinary cycle-popped arrow book, attach independent uniforms to
its records, and at parameter p set each bit to 1_(U<=p_i). Append
terminal ones and use the sampler from (J10). The upper book lengths
stabilize every threshold assignment. Sleep-bit monotonicity gives,
simultaneously for all positive p<=p', nested emission vectors and
used arrow prefixes. The bits themselves change; nested paired
records, raw instruction histories, or regenerated full books are not
asserted in this sleep-parameter comparison.

For every mixed full-support bridge, (J8) bounds its loop intensities
by those of an ordinary loop soup. With G=(I-K)^(-1),

    E[sum_j N_ij | full,S]<=G_ii-1,
    E[sum_(i,j) N_ij | full,S]<=tr(G)-|V|.                 (J19)

Indeed an ordinary soup has expected edge count K(i,j)G(j,i);
summing gives (KG)(i,i)=G(i,i)-1. This is uniform in the sleep
parameters and initial activity set, for one particle per site.
A conditioned surplus-IDLA phase contributes its own additional
current as in (J5).

As every p_i tends to zero, at arbitrary relative rates, the mixed
bridge with active set T converges in total variation of its current
vector to ordinary loops confined to T:

    E[t^N | full,S] -> d_T(K)/d_T(W).                      (J20)

Use g_p(W)=sum_C p_C(1-p)_(V\C)/d_C(W). Then g_p tends to one
and g_(p^S) tends to 1/d_T. Coefficient continuity and domination by
the summable ordinary loop intensities give convergence of the total
intensities, hence a Poisson coupling in total variation. The
all-active case gives the whole ordinary soup and makes (J19) sharp.
This is a limit of positive-sleep conditional laws, not conditioning
on a full event at zero sleep.

Let B be the sites whose first bit is one in the all-active
full-conditioned buffer book. If L is an ordinary loop range and P
an independent Bernoulli(p_i) field, then

    B =_law (V\L) union (L intersection P),
    P(S subset B | full)=p_S Z(r^S;K).                     (J21)

The normal form preserves the first-bit field of its original book
when unused sites have their appended one: if the physical prefix
does not use a site, the creator suffix's first bit there is one.
Before factorization, bits on the loop range are independent with
probabilities p_i. This proves the law and its normalization. The
complement field L intersection Bernoulli(1-p_i) is an increasing
function of independent Poisson loop counts and independent bits;
the product covariance inequality proves association of it and of B.
No FKG lattice condition is asserted. These are first-record bits
conditioned on full occupation, not unconditioned occupation marks.

### The physical bridge monoid is not a trace monoid

The preceding normal forms do not give a trace/Lyndon factorization
of R_* itself. There is an explicit positive-weight obstruction on
the wired path 2--0--1--3: vertices 0,1 choose their two neighbors
equally, and each endpoint chooses its neighbor or the sink equally.
For records (bit,destination), form three two-cycle books:

| Site | m | a | b |
|---|---|---|---|
| 0 | (1,1) | (0,2) | empty |
| 1 | (1,0) | empty | (0,3) |
| 2 | empty | (0,0) | empty |
| 3 | empty | empty | (0,1) |

Here m is in M, while a,b,ma,mb are distinct atoms of R_*. The
balanced left divisors of ma are exactly empty,m,ma; those of mb
are empty,m,mb. Their nonempty right divisors are respectively
a,ma and b,mb, none in M. The simple cycles a,b have no proper
balanced factors. Because a and b have disjoint supports,

    (ma)b=(mb)a=mab.                                      (J22)

The two atomic factorizations have different atom multisets, so no
choice of commutation relations can connect them. This rules out
even an abstract trace-monoid isomorphism for R_*. The complete
balanced-prefix list of mab is empty,m,ma,mb,mab: balance on each
edge of the tree forces paired orientations, and local record order
forces either outer cycle to follow the central one. Thus no factor
has been missed. Its weight is 2^(-6)p_0p_1(1-p_0)(1-p_1)
(1-p_2)(1-p_3)>0 for every interior sleep vector. This hand proof
received an independent audit. It leaves (J9)--(J18) intact, but
prevents treating the physical bridge atoms as the old trace alphabet.

### Activation prefixes are conditional auxiliary avalanches

**Proved; independent AI proof audit passed.** In (J13) put
A=S_1\S_2 and T_1=V\S_1. Set auxiliary sleep probabilities
q_i=p_i on S_1 and q_i=1 on T_1. Force every emitted bit at T_1
to one, writing F(u) for the resulting book. Balanced right cuts
depend only on arrows, and for each such suffix c,
c belongs to H_(S_1) iff F(c) belongs to M. Therefore F maps
N_(S_1) intersection H_(S_2) onto exactly the emission-book support
of the auxiliary process started active on A and asleep elsewhere,
conditioned to end full. Its fibers are all bit assignments at T_1.

For every book a in that auxiliary support, its mixed-input cylinder
has probability p_A wt_q(a). The auxiliary full probability is
Z(r^(S_2))/Z(r^(S_1)). Hence its conditional probability is
wt_q(a)/kappa_(S_2,S_1). Replacing its emitted T_1 bits by independent
Bernoulli(p_i) bits gives exactly

    P(U=u)=wt_p(u)/kappa_(S_2,S_1),
    U in N_(S_1) intersection H_(S_2).                    (J23)

The weight identity is wt_p(u)=wt_q(F(u)) times the independent
decoration probability; summing the decorations gives one. No unused
initial-sleeper or terminal-arrow factor is missing. Probability-one
sites use the deterministic buffer rule, without an infinite-rate
raw-stack convention. At each such site every arriving particle
displaces one resident, so its emission equation is m_i=I_i(m).
For A={x}, the extra prefix is thus a conditional one-active avalanche
among sleepers, with previously active sites acting as occupied unit
capacities. Its arrows and S_1 bits are independent of the original
sleep parameters on T_1; those parameters enter only in the decorations.
This identifies finite buffer records, not original raw histories.

### Open input paths and the full terminal boundary

**Proved; independent AI proof audit passed.** Let eta be all active,
N=|eta|>=|V|, and F_eta>0 its ordinary IDLA full-filling probability.
Let J_eta contain the finite paired books with
eta+I(ell)-ell=1_V. They have exactly N-|V| exit arrows. Let Q_eta
be the local arrow prefixes of complete IDLA stabilizations ending
full, with arbitrary bit decorations. Complete stabilization includes
all surplus exits, even after the last vacancy fills. Then

    J_eta=Q_eta H=R_eta M,
    R_eta=Q_eta R_*, all uniquely,                         (J24)

where R_eta forbids a nonempty balanced right M divisor. These are
right-module factorizations; J_eta need not be a monoid.

Append terminal ones to j in J_eta. Both IDLA and ARW stabilize below
its lengths. Each uses at most N-|V| exit arrows, while stable capacity
forces at least that many, so both end full. Extract the complete IDLA
prefix q. Its remainder h is balanced and internal; conversely qh
replays that same IDLA stabilization. Every balanced right suffix of
j lies after q: removing it leaves another full IDLA stable cut,
which is at least the least IDLA cut q. Thus the balanced right
divisors of j are precisely those of h, and (J9) proves (J24).
The R_eta factor is the actual ARW prefix by the same least-action
argument as (J9).

For the index check, write shifted bits and arrows at i as their
original entries q_i+j. Since eta+I(q)=q+1, the ARW fixed-point
equation at m=q+r becomes
r_i=1+I_i^q(r)-B_i(q_i+1+I_i^q(r)). The least shifted solution is
therefore the binary full-input odometer. This is an equality of local
emission records, without asserting a chronological IDLA phase in ARW.

Complete IDLA prefix cylinders have total weight F_eta. Consequently

    W(J_eta)=F_eta/d_V, W(R_eta)=F_eta/(p_V Z),
    P_eta(full)=F_eta/Z,
    law(a | full)=(p_V Z/F_eta)wt(a), a in R_eta.           (J25)

Under a=qb the two factors are independent, with laws wt(q)/F_eta
and p_V Z wt(b). This is the ordered form of (J5), retaining every
paired record and labeled exit. Sample IDLA conditioned to fill,
decorate its emitted arrows with independent p bits, and concatenate
an independent binary bridge from (J10). Keep terminal ones and sample
unused arrows and later tails freshly. If eta>=1_V, F_eta=1: the
IDLA prefix is sampled by a fixed order of independent surplus walks
to the sink. For arbitrary eta, no efficient conditional-IDLA sampler
is asserted. The ordered surplus identity also refines the existing
CB Section 5 reduction; no novelty claim is made.

### Formal stable cuts with arbitrary terminal supports

**Proved; independent AI proof audit passed.** For a finite all-active
eta and B subset V, call j a formal stable-cut book if
eta+I(ell)-ell=1_B, the last bit at each used i outside B is zero,
and terminal ones are appended at ell_i+1 for i in B. Its weight is
p_B wt(j). Let Q_eta(B) be the sum of these weights. The appended
bits make ell a stable upper bound, but the actual stabilization below
ell may end at a different support C.

Define a remainder book b from C to B by
1_C+I(|b|)-|b|=1_B, first bit one at every used site in C, and last
bit zero at every used site outside B. Empty local words have neither
test. Define the nonnegative matrix

    R_p(C,B)=(p_B/p_C) sum_(b:C to B)wt_p(b).

Legal stabilization of j gives a unique physical emitted prefix a
ending at C and such a remainder b. At a used site of C, b starts
with the terminal bit already read by the physical stabilization.
Conversely concatenate a physical eta-to-C prefix with any such b:
the legal stabilization replays a, and the combined book has the
specified formal B cut. The weights factor as
p_B wt(ab)=[p_C wt(a)][(p_B/p_C)wt(b)]. Hence

    Q_eta(B)=sum_C P_eta(C) R_p(C,B).                      (J26)

Here P_eta(C) is the unconditioned ARW terminal probability. The
identity also holds with multiplicative marks on all arrows, replacing
P_eta(C) by its joint current transform. R_p(C,B)=0 if |B|>|C|;
equal cardinality does not force C=B or inclusion.

All sums are finite. Choose h>0 and a<1 with Kh<=a h, and put
b=max_i k_i/h_i. In a book from eta to B the exit count is fixed
at k=|eta|-|B|. Scaling an internal edge i->j by h_j/h_i and a
sink edge by 1/h_i bounds the total arrow weight at length L by
h^(eta-1_B) binom(L+|V|-1,|V|-1) binom(L,k) a^(L-k)b^k.
This is summable; bit restrictions and fixed boundary factors preserve
finiteness. The same argument applies to R.

On the diagonal the remainder is balanced and internal. Summing its
bits gives a factor alpha_i per used site, where alpha_i=p_i on C
and alpha_i=1-p_i outside C. Therefore

    R_p(C,C)=g_alpha(K),
    R_p(V,V)=p_V Z(r;K)/d_V(K).                           (J27)

The full column of (J26) has only C=V and reduces to (J25). Partial
columns mix different terminal supports through positive open remainders.

A two-site example makes this mixing explicit. Let 0 jump to 1 and
let 1 jump to 0 or the sink with probability one half. Start one active
at 0 and prescribe the paired words
j_0=((0,1),(0,1)), j_1=((1,0),(0,0)), with a terminal one at 0.
The formal length vector (2,2) ends at B={0}, but the particle jumps
once from 0 to 1 and sleeps: C={1}, emissions (1,0). The remaining
book is b_0=((0,1)), b_1=((1,0),(0,0)), carrying net flow C-B.
The focused certificate checks all eight bit assignments on this
arrow table compatible with the formal B cut, their decompositions,
and exact weights. Six end at {0}, two at {1}; these are counts, not
probabilities. It also checks the positive factorization on one site
with independently marked internal and sink arrows. This is a failure
of identifying a formal partial cut with the actual outcome, not a
counterexample to monotonicity.

There is an exact reciprocity for the remainder matrix. Reverse
each local word and complement its bits, keeping every arrow direction
and each bit-arrow pair. This carries a C-to-B remainder at p to a
B^c-to-C^c remainder at 1-p, with unchanged weight. For the product
Bernoulli mass mu_p(C)=p_C(1-p)_(V\C), the prefactors give

    mu_p(C) R_p(C,B)
       =mu_p(B) R_(1-p)(B^c,C^c).                        (J28)

The boundary tests exchange exactly: first ones on C become last
zeros outside C^c, and last zeros outside B become first ones on B^c.
This reversal concerns finite remainder words, not physical time reversal.
The positive factorization, finiteness, diagonal formula, and reciprocity
passed an independent proof audit, including the unused-site boundaries.

### The remainder matrix is the complementary-sleep duality

**Proved; independent AI proof audit passed.** Write q=1-p,
r=q/p, s=p/q, d=det(I-K), and mu_p(B)=p_B q_(V\B). The pairing
H_(r,s) from Section 36 has the exact record interpretation

    R_p(C,B)=[mu_p(B)/d] H_(r,s)(C,B^c).                 (J29)

Here is a direct constant-term proof, including its convergence
convention. Put J_i(z)=k_i+sum_j K(i,j)z_j and a_i(z)=J_i(z)/z_i.
Choose h>0 with Kh<=theta h, theta<1, and take the torus
|z_i|=t h_i with t sufficiently large that |a_i(z)|<1. Use the
absolutely convergent geometric expansions on this torus. Summing
the unrestricted emitted bits, except the last zeros outside B, gives

    Q_eta(B)=p_B CT z^(eta-1_B)
       product_(i in B)(1-a_i)^(-1)
       product_(i outside B)(1-p_i a_i)/(1-a_i)
       =mu_p(B) Phi(z^eta w_(B^c)^s),
    Phi(f)=CT f(z)/product_i[z_i-J_i(z)].

This is a specified large-torus constant term, not an expansion at
the origin. For every polynomial g, Phi(z_i(z_i-J_i)g)=0:
after cancellation, every remaining denominator expansion has only
nonnegative powers of z_i, and the numerator contains an extra z_i.
Thus Phi descends to the addition algebra. For squarefree z_C with
|C|<|V|, every expansion term has total degree
|C|-|V|-number_of_sink_arrows<0, so Phi(z_C)=0. For C=V only
balanced internal arrow books contribute, with total weight 1/d.
Therefore Phi=Lambda/d on the squarefree basis, without assuming
semisimplicity. On binary input rows this says
Q(A,B)=[mu_p(B)/d](F E_s^T)(A,B^c). Multiplying Q=P_r R_p by
E_r proves (J29).

Let pi_s=delta_V P_s be the stationary law at sleep probabilities q.
The complemented stochastic link L_r(C,B)=B_(r,s)(C,B^c) from
Section 36 is consequently

    L_r=R_p diag(f),
    f_B=d pi_s(B^c)/mu_p(B)>0,   R_p f=1.                 (J30)

Equivalently, for fixed C, assign a remainder b from C to B the
probability f_B(p_B/p_C)wt_p(b). Its B marginal is exactly L_r.
This realizes the old duality by finite positive records; it does not
change the proper-target payoff obstruction in Section 36.
The algebraic identification and stochastic normalization are stated
here for the original unmarked kernel. The independent book identity
(J26) and reciprocity (J28) retain arbitrary arrow marks.

The focused certificate sums all remainder books on the same two-site
single-leak cycle exactly. At site 1, a length m word with k forced
sink exits has binom(m,k) orders; after the first three lengths its
bit factor is constant, so its tail is a geometric derivative.
All 16 resulting entries agree with (J29), all 16 reciprocity entries
agree with (J28), and every normalized row in (J30) sums to one.

### Single activation currents and fixed-book Boolean responses

Two scoped limitations clarify the earlier activation coupling. For
x in S let R=S\{x} and write p^S for p on S and one elsewhere.
The extra block from (J15) has directed-current transform
kappa_x(S;W)/kappa_x(S;K), with

    kappa_x(S;W)=g_(p^R)(W)/g_(p^S)(W),
    lambda_(nu,x)(S)=a_nu(p^S)-a_nu(p^R)>=0.              (J31)

Thus its currents are compound Poisson, with atoms visiting x. For
explicit means, put F_S=sum_(D subset S)r_D/d_(V\D), and give D
probability r_D/[d_(V\D)F_S]. If Gamma^D is the killed Green matrix
on V\D, extended by zero, then
E J_ij=K_ij(E_R Gamma^D_ji-E_S Gamma^D_ji). This follows by
differentiating log kappa; it asserts no order of the two D laws.

The law of this same labeled block need not increase or decrease
when other sites were activated earlier. On x->y->x with return
probability c, leakage 1-c at y, and both p=1/2, let z mark one
complete circuit. The two unnormalized block transforms are

    kappa_x({x,y};z)=(1-cz/2)/(1-3cz/4),
    kappa_x({x};z)=1/(1-cz/2).                            (J32)

Their circuit-m intensities are c^m[(3/4)^m-(1/2)^m]/m and
c^m(1/2)^m/m. The first difference is positive at m=1 and negative
at m=2. At c=3/4, their probabilities of at least one circuit are
3/10 and 3/8, but of at least two are 27/160 and 9/64. Hence the
two block laws are not stochastically comparable, consistently with
the positive ordering of total bridges along a fixed activation chain.

For a fixed decorated heap h, (J18) gives another exact description.
For an occurrence e, let P_e be the inclusion-minimal sets of labels
appearing on directed e-to-cemetery paths. Activating A retains e
iff A meets every member of P_e. Its minimal activation sets are
therefore the minimal transversals of this finite path-label family.
These are grouped edge cuts: all constraints with one label disappear
together. The resulting Boolean response need not be submodular or
supermodular. For a realizable example, take a cycle 0->1->2->0
with bits (1,0,0), followed by a zero-bit self-loop at 0. For the
self-loop occurrence, the response is

    f(A)=x_0 OR (x_1 AND x_2), x_i=1_(i in A).            (J33)

Use K(0,0)=K(0,1)=1/2, K(1,2)=1, K(2,0)=K(2,s)=1/2.
The book has positive weight and belongs to R_*. Its cemetery paths
have minimal label sets {0,1},{0,2}. The choices A={1},{2} violate
submodularity, and A={0},{1,2} violate supermodularity. Minimal
activators {0} and {1,2} are not matroid bases; even adjoining the
empty set does not make the successful-set family accessible. For a
fixed A, however, all admissible suffixes form a distributive lattice:
contract the constraint graph's strongly connected components, discard
those reaching the cemetery, and take upper ideals of the remaining
poset. This fixed-graph statement does not impose an antimatroid on
activation labels. The formulas and these hand examples were checked
in the main thread after separate agent derivations; no sweep was used.

### A constructive bridge between remainders and the pairing

**Proved; independent AI proof audit passed.** There is also a direct
finite-word proof of (J29), yielding a conditional remainder sampler.
Fix C,B and put I=C\B, H=B\C, U=C intersection B,
O=V\(C union B), S=U union O. Set alpha_i=1-p_i on O and alpha_i=p_i
elsewhere. At every i in I a remainder has at least two records:
its first bit is one and its last bit zero. Remove those two records,
retaining their arrows as an ordered first/last pair. Let eta count
the internal destinations of all these forced arrows; retain sink
labels separately. At every i in O reverse the local word and
complement its bits. The resulting h satisfies

    eta+1_S+I(h)-|h|=1_V,
    first used bits on S are one.                        (J34)

Thus h is a formal full stable-cut book for auxiliary input sleeping
S plus active eta at sleep vector alpha. Conversely insert the two
removed records at the ends of each I-word and undo the O-word
operation. This restores every remainder, including empty intervening
words. Arrow directions never change, and the first/last labels give
no factor 2!. If d denotes the forced-arrow tuple and w(d) its product
K-weight, the pointwise identity is

    (p_B/p_C)wt_p(b)
       =p_H(1-p)_I w(d) wt_alpha(h)
       =r_I w(d)[p_(I union H) wt_alpha(h)].              (J35)

The bracket is the mixed-input formal-full weight: initial sleep bits
on S are fixed, while the terminal factor outside S is p_(I union H).
This remains correct when additional actives immediately wake sleepers.

For each forced tuple, the full normal form is h=a m, where a is the
auxiliary physical prefix ending full and m belongs to M_alpha. The
stable upper bound and particle count force the physical prefix to end
full; its unused suffix is balanced with first bits one. Conversely
every such pair concatenates. Writing
F(d)=P_alpha(full | sleeping S plus active eta(d)), summing (J35) gives

    R_p(C,B)=r_I g_alpha E_d F(d)
       =[mu_p(B) Z_alpha/d_V] E_d F(d).                  (J36)

The expectation uses independent K-arrows at the two forced positions
per i in I. Impossible forced tuples simply have F(d)=0. In the
pairing H_(r,1/r)(C,B^c), the overlap I contributes exactly J_i^2,
and nonoverlap sites U,O are sleepers with odds r,1/r respectively.
Section 36 therefore gives H=Z_alpha E_d F(d), proving (J29) again
by this explicit record transformation.

If R_p(C,B)>0, sample the forced arrows and auxiliary ARW run jointly
conditioned on full occupation. Independently sample a creator word
m with law wt_alpha(m)/g_alpha, concatenate its finite records after
the auxiliary emission book a, and invert the transformation. This
has exactly the conditional remainder law. The forced-arrow marginal
is w(d)F(d)/E F; conditioning only the run after sampling an unbiased
fixed tuple would be incorrect. No efficiency bound is claimed.
Only emitted records and m enter the inverse map. Auxiliary terminal
arrows, later tails, and the artificial appended terminal ones are
not inserted as extra records.

### Finite sleep mixtures, parameter corners, and uniform tails

**Proved; independent AI derivation and main-thread proof check.** Keep
I,H,U,O from (J34). For A subset U union O, let c_A(K) be the total
weight of ordered arrow books with flow C to B, at least two emissions
at every i in I, and used sites among U union O exactly A. There is
no emission requirement on H: those sites may receive without emitting.
Summing the remainder bits gives exactly

    R_p(C,B)=p_H(1-p)_I
       sum_(A subset U union O)c_A(K) p_(A intersection U)
                                      (1-p)_(A intersection O). (J37)

At an I-site the first one and last zero are distinct and, after the
boundary ratio, contribute 1-p_i. At an H-site only its terminal
factor p_i remains. A used U-site contributes p_i, a used O-site
contributes 1-p_i, and an unused U- or O-site contributes one. This
proves the formula, including all zero-length cases. Each c_A is finite
by the geometric bound after (J26).

Writing alpha_i=p_i on B and 1-p_i outside B, this is the polynomial
alpha_(C triangle B) sum_A c_A alpha_A, with nonnegative coefficients
and degree at most one in each coordinate. Hence R_p(C,B) increases
in p_i for i in B and decreases outside B. These are unnormalized
weights, not monotonicity claims for their conditional laws.
The polynomial extends to every closed-cube corner p=1_T. There
a nonzero entry requires B triangle T subset C triangle T. Thus
R_(1_T) is triangular when states are ordered by distance to T, with

    R_(1_T)(C,C)=1/d_(V\(C triangle T)),
    det R_(1_T)=product_(D subset V)d_D^(-1).              (J38)

This follows from (J27): only balanced arrows confined to the sites
where C agrees with T remain. It concerns the polynomial extension,
without conditioning a zero-sleep ARW event.

For R_p(C,B)>0, its conditional ordered-arrow law is the finite mixture
of the c_A-normalized arrow laws, weighted proportionally to
c_A p_(A intersection U)(1-p)_(A intersection O). Each component is
independent of p; the common p_H(1-p)_I cancels. In particular the
arrow law is independent of the parameters on C triangle B. Given
the arrows, endpoint bits are fixed and all other bits are independent
with their original probabilities, so the full paired law need not
be independent of p.

The geometric bound supplies a common radius z_0>1 for all arrow-length
series c_A(z). There are only finitely many positive components, so
for total record length L,

    sup_(0<p<1) E[z_0^L | C-to-B remainder]
       <=max_(A:c_A(1)>0) c_A(z_0)/c_A(1)<infinity.        (J39)

Thus rare remainder endpoints do not cause unbounded conditional
record lengths as sleep varies on a fixed graph. The arrow laws lie
in a finite-dimensional compact simplex in total variation. A parameter
face has a unique conditional limit if the reduced mixture denominator
is nonzero; otherwise relative vanishing rates can matter. For example,
take C={0}, B=empty, arrows 0->1,2,s with probabilities 1/3 each,
and 1,2->0 deterministically. The empty neutral signature is impossible,
but both singleton signatures have positive weight. Sending
(1-p_1,1-p_2) to (epsilon,epsilon^2) or to its reversal gives different,
disjoint-support arrow limits. Uniform exponential tails still hold.

When B=C^c, there are no neutral sites U,O, so the entire conditional
arrow law is independent of p. There is an IDLA sampler for this case.
Jointly sample two forced outgoing arrows at every i in C and complete
IDLA from their internal destinations, conditioned on filling V. Keep
the forced arrows and IDLA emission prefix. Concatenate an independent
ordinary balanced arrow book, decorate these remaining records with
iid p bits, and restore the first-one/last-zero records at C. Indeed
(J34) now has S empty, and its formal full book factors as Q_eta H
by (J24). Its arrow mass is F_eta/d_V, so conditioning the forced
tuple and IDLA jointly gives exactly (J37)'s sole component. This
corollary uses the audited factorizations; it does not claim an
efficient conditional-filling algorithm.

### Stable cuts form lattices, but open remainders lose distributivity

**Proved; independent AI derivation and exact finite check.** Remainder
books form a category: concatenate a C-to-D book and a D-to-B book
sitewise to obtain a C-to-B book. If the first word is unused at
i in C, balance forces i into D, so a used second word starts with
one; the analogous argument preserves last zeros outside B. Boundary
weight ratios telescope, and empty books are identities.

Fix a C-to-B book h of lengths ell, append terminal ones on B, and
put beta_i(0)=0. Its formal stable cuts are precisely the fixed points
in 0<=m<=ell of

    F_i(m)=1_C(i)+I_i(m)-beta_i(1_C(i)+I_i(m)).           (J40)

Every index is specified by h and its terminal bits. Since
n-beta_i(n) is nondecreasing, F is monotone, with F(0)=0 and
F(ell)=ell. A fixed point has intermediate support
D=1_C+I(m)-m and corresponds exactly to a factorization into C-to-D
and D-to-B remainders. For two fixed points, iterate F upward from
their coordinate maximum to obtain their least fixed upper bound;
iterate downward from their coordinate minimum for the greatest fixed
lower bound. Finite termination proves a complete lattice of cuts.

This lattice need not be distributive. Take sites
(v_1,v_2,v_3,e_12,e_23,e_31), initial and final support
C=B={v_1,v_2,v_3}, and the following twelve paired records:

| Site | First record | Second record |
|---|---|---|
| v_1 | (1,e_12) | (0,e_31) |
| v_2 | (1,e_23) | (0,e_12) |
| v_3 | (1,e_31) | (0,e_23) |
| e_12 | (1,v_1) | (0,v_2) |
| e_23 | (1,v_2) | (0,v_3) |
| e_31 | (1,v_3) | (0,v_1) |

Give each site probability 1/3 to each of its displayed neighbors
and 1/3 to the sink, with arbitrary interior sleep probabilities.
The book has positive weight, every site sends and receives two arrows,
and its endpoint conditions hold. At a fixed cut every local length
is zero or two. Writing x_i=m_(v_i)/2 and y_ij=m_(e_ij)/2, (J40) becomes

    x_i=y_ij OR y_ik,   y_ij=x_i AND x_j.                 (J41)

Thus the selected v-indices are empty, one of the three pairs, or all
three. In the displayed site order the five cuts are zero, ell=(2,...,2),
and
(2,2,0,2,0,0), (0,2,2,0,2,0), (2,0,2,0,0,2).
Any two distinct middle cuts have meet zero and join ell: one iteration
from their coordinate minimum or maximum gives those bounds. This is
the five-element nondistributive diamond M_3. The verifier independently
checks every one of the 729 possible cuts, finding exactly these five.
No minimality claim is made for the example.

For C=B=V, in contrast, every cut is balanced: no sink arrows are
available and cardinality forces its intermediate support to be V.
Coordinate maxima and minima of balanced prefixes remain balanced by
monotonicity of incoming counts and equality of their total masses;
their next bits remain one. Hence the full-support cut lattice is
distributive. In (J41) the whole book is balanced, but the three middle
cuts transfer occupation to different supports and are not balanced
heap divisors. No net loss is needed for the distinction.

The same example also rules out a global least-common-extension rule
for remainder prefixes. Keep its first two middle prefixes fixed and
interchange the two arrow destinations at e_31 in the full book, leaving
bits and all other records unchanged. This site is unused in both
prefixes. Both full books are valid common extensions, with exactly
the same fixed-cut lattice, so each is a minimal common extension.
They have different first e_31 arrows and are not equal; a common
prefix dividing both cannot be their required full join. Thus no
least common right multiple exists in the whole remainder category,
although cuts inside either fixed book form a lattice. This last
consequence follows directly from the displayed words and fixed cuts.

More generally, take any finite directed multigraph with two incoming
and two outgoing internal arc occurrences at each site. Give each site
records (1,d_1),(0,d_2), and let C be the OR sites, with all others AND.
Every fixed coordinate is zero or two, and the Boolean equation at
that site is the OR or AND of its two incoming source values. This
is an order isomorphism between stable cuts and network fixed points.
Parallel arcs and self-loops are allowed, counted with multiplicity;
the two local arrow orders do not affect this fixed-point correspondence.
Setting K(i,j) to one third of the arc multiplicity and giving each
site sink probability 1/3 realizes every such book with positive weight.
No assertion about arbitrary fan-in/fan-out networks is made.

Already four sites give the nonmodular pentagon N_5. Take C={a,b,c},
local arrows (a,d) at a, (b,a) at b, (c,d) at c, and (b,c) at d,
always with bits (1,0). The equations are

    x_a=x_a OR x_b, x_b=x_b OR x_d,
    x_c=x_c OR x_d, x_d=x_a AND x_c.                    (J42)

Their solutions are 0000,1000,1100,0010,1111. Indeed x_c=0 forces
x_d=0 and x_b<=x_a, while x_c=1 forces x_a=x_b=x_d. The first two
nonzero solutions form a chain and 0010 is incomparable to both.
This is N_5, with one vacancy at both endpoints and no emitted sink
arrows. A separate exhaustive check of its 81 candidate cuts verifies
the five solutions. Thus even this near-full class need not have a
modular cut lattice. The network correspondence and both examples
were independently derived and reviewed.

### Association on diagonal remainders and its boundary failure

For C=B, (J27) identifies the joint arrow-range/current law as ordinary
Markov-loop soup tilted by alpha_L, where L is its range and
alpha_i=p_i on C, 1-p_i outside. All alpha_i lie in (0,1). On the
independent Poisson loop-count lattice, the indicator that a given site
is visited is submodular: its value at a join is the maximum, and at
a meet is at most the minimum. Therefore log(alpha_L) is supermodular.
The product-Poisson mass is log-modular, so the tilted counts satisfy
the FKG lattice inequality. The finite-lattice association step uses
[Fortuin--Kasteleyn--Ginibre, Proposition 1, p. 91](https://math.bme.hu/~balint/oktatas/perkolacio/percolation_papers/fortuin_kasteleyn_ginibre.pdf#page=3),
whose statement was read in the primary PDF; no stronger version is
imported.

For completeness truncate both the number of loop types and every
count. Apply that proposition on the finite product of chains, then
remove the count bounds by dominated convergence. The total loop
intensity is -log d_V<infinity, so only finitely many loop instances
are present almost surely; removing the type truncation then converges
pointwise. The tilt lies between alpha_V and one and has positive
normalizer d_V g_alpha. Thus under every diagonal remainder law,

    all range indicators and directed edge counts
    are jointly positively associated.                  (J43)

This persists under further vertex suppression fields in [0,1],
including hard absence conditions by limits. It also persists after
fixing finitely many latent loop counts: sites already covered by
those loops have effective remaining tilt one. Unbounded current
observables follow by truncation and finite transient-walk moments.
The proof and its limits were checked after an independent derivation.

Association fails for general endpoints, already on (J39)'s fork
with C={0}, B=empty. The neutral used set A is never empty, but
both singleton possibilities have positive weight. With q_i=1-p_i
and D=5(q_1+q_2)+22q_1q_2, exact arrow-word sums give

    P(A={1})=5q_1/D, P(A={2})=5q_2/D,
    P(A={1,2})=22q_1q_2/D,
    Cov(1_(1 used),1_(2 used))=-25q_1q_2/D^2<0.          (J44)

For one branch the arrow sum is sum_(m>=1)(m+1)/3^(m+1)=5/12;
the total is 8/3, leaving 11/6 for both branches. Their bit factors
are respectively q_0q_1, q_0q_2, q_0q_1q_2, proving the formulas.
The same fork shows the limits of (J43): diagonal loop ranges can
be empty, {0,1}, {0,2}, or {0,1,2}, but not {0}. Their range law
therefore fails the lattice inequality, despite association; conditioning
site 0 to be visited produces strictly negative correlation between
the other two range indicators. Neither arbitrary presence conditioning
nor FKG for the range image is asserted. These are record-range
statements, not assertions about ARW terminal occupations.

### An absolutely convergent inverse from stable-cut topology

**Proved; independent AI proof audit passed.** Let L(h) be the stable-cut
lattice of a typed remainder h:C->B, and let mu_h be its Mobius number
from zero to its full length vector; the identity book has mu_h=1.
Write W_p(h)=(p_B/p_C)wt_p(h), and mark every emitted arrow by t.
Then

    [R_p(t)^(-1)](C,B)
       =sum_(h:C->B) mu_h W_p(h) t^|h|.                 (J45)

The sum is absolutely convergent per ordered paired book for
|t|<1/rho(K), with infinite radius when rho(K)=0. In particular it
converges at the original weights t=1. The following elementary bound
supplies the convergence, without a spectral assumption on R-I.

First, cuts ending at a fixed support D form a distributive sublattice
of L(h), even when the whole lattice is nonmodular. For two such cuts
m,n, let u=max(m,n), v=min(m,n) coordinatewise. Monotonicity gives
I(u)>=u+1_D-1_C and I(v)<=v+1_D-1_C. Both original cuts use exactly
k=|C|-|D| sink arrows. Summing these inequalities gives at most k sink
arrows for u and at least k for v; prefix inclusion gives the reverse
bounds. Thus equality holds in every coordinate. The relevant next-one
or last-zero bit at each new cut comes from one of the old cuts at
that site, so u,v are stable with endpoint D. The fiber need not be
convex: its operations are coordinatewise, not a constraint on cuts
of other supports lying between its elements.

There is at most one lattice atom with any endpoint D different from
C: two distinct atoms would have a nonzero stable coordinate minimum
with endpoint D. Atoms ending at C have pairwise disjoint positive
coordinate supports, since their coordinate minima must be zero.
There are at most |V| of them. Therefore every book, independently of
its length, has at most

    A_n=2^n-1+n atoms, n=|V|,
    |mu_h|<=M_n:=2^(2^n-1+n).                            (J46)

For the second assertion use the elementary crosscut identity
mu(0,x)=sum_(S subset atoms: join S=x)(-1)^|S|. Summing its right
side over all y<=x gives (1-1)^(number of atoms below x), which is
one at zero and zero otherwise, exactly the Mobius recurrence.
This proves the uniform bound; no sharpness is claimed.

Particle counts sharpen this estimate. If |C|=c and |B|=b, every
intermediate endpoint has size between b and c. The same argument
therefore bounds the atom count by n-1+sum_(j=b)^c binom(n,j).
In particular a one-particle zero-exit book has at most 2n-1 atoms,
without any restriction on its word lengths or bits.

To prove (J45) formally, factor h into k nonempty remainders. Such
factorizations are exactly strict chains of k intervals from zero to
the full cut, with all intermediate supports determined by the cuts.
The boundary factors telescope, so their alternating contribution in
the formal inverse is W_p(h)mu_h. There is no global chronological
interleaving count. The cut lattice of a prefix is its lower interval,
and that of a suffix is the corresponding upper interval.

For any radius z<1/rho(K), choose a>rho(K) with az<1 and v>0
satisfying Kv<=a v. The fixed-exit geometric estimate after (J26)
then proves convergence of R_p(z). By (J46),

    sum_(h:C->B)|mu_h|W_p(h)|t|^|h|
       <=M_n R_p(|t|)(C,B).                              (J47)

This holds locally uniformly in the stated disk. The convolution
with R is absolutely bounded by M_n R_p(|t|)^2. Regrouping by the
concatenated book is therefore legitimate, and each nonempty book
cancels by sum_(m in L(h))mu(0,m)=0. The suffix recurrence gives
the opposite product. Both products equal the identity, proving the
actual inverse throughout the disk. Absolute convergence of the
unreduced factor-count Neumann series at t=1 is not asserted.

For full-support remainders, distributivity gives mu_h in {0,1,-1}:
the crosscut sum vanishes unless the lattice is Boolean, in which
case it is (-1)^(number of atoms). In the M_3 book, however,
mu_h=-1+3=2. Its individual contribution is
2*3^(-12)*product_i p_i(1-p_i)*t^12. This is one ordered-book term,
not the whole coefficient of t^12. The exact cut certificate checks
Mobius values 2 for M_3 and 1 for N_5 and verifies the endpoint-fiber
closure on both examples.

Equivalently, for a nonidentity book mu_h is the reduced Euler
characteristic of the order complex of L(h) with its endpoints
removed: a k-piece factorization has k-1 proper cuts and contributes
the required alternating simplex sign. The one-piece chain supplies
the reduced -1 term. Thus (J26) can be inverted through an explicit,
absolutely convergent topological cancellation over finite records.
The sign needed for occupation monotonicity remains unproved.

The formal part is categorical Mobius inversion: local concatenation
defines the arrows, and positive total length bounds every strict
factorization. This is a Mobius category in the sense discussed in
[Leinster, *Notions of Mobius inversion*, Section 1 and Theorem 6.1](https://arxiv.org/pdf/1201.0413).
Those sections were read; no novelty is claimed for the formal
inversion principle. The length-independent bound and convergence
at stochastic weights are supplied by the argument above.

The inverse remains regular at every sleep-parameter face. Indeed,
(J45)--(J47) give entrywise

    |R_p^(-1)(C,B)|<=M_n R_p(C,B).                       (J48)

The polynomial extension (J37) is uniformly bounded on the closed
cube. Hence its interior inverses are uniformly bounded, which forces
every boundary matrix to remain nonsingular by passage to the limit
in ||v||<=L||R_p v||. The inverse extends continuously and obeys
(J48) on the whole cube. Its determinant is positive everywhere,
since at p=1 it is diagonal with entries 1/d_C.

There is also a direct uniform series argument. For any fixed arrow
book, cancel each initial forced-one factor against its denominator
in p_C. The remaining paired-book weight is a product of p_i,1-p_i
and arrow weights, without denominators. Summing all free-bit assignments
gives at most the arrow weight, uniformly on the closed cube. Apply
the bound M_n before this bit sum, then the polynomial-times-geometric
length bound after (J26). Thus the absolute tails of (J45), grouped
by total length, tend uniformly to zero over the entire closed cube.
This avoids paying an exponential factor for the free decorations.

For an explicit bound, put N=2^n. The first-cut interpretation (J71)
below gives R_p(C,B)<=R_p(B,B)<=1/d, with d=det(I-K), and the
maximum 1/d is attained at sleep corners. In the maximum-row-sum
norm, ||R_p||<=N/d and ||R_p^(-1)||<=M_n N/d. These
graph-dependent bounds are coarse and concern R, not the binary
stabilization inverse P_r^(-1), which can diverge at zero sleep.
The boundary and convergence consequences passed a separate audit,
using (J45)--(J46) as the already proved inputs.

### A one-particle Boolean subclass has only unit Mobius weights

**Proved; independent AI derivation and main-thread proof check.** In
the two-in/two-out Boolean construction after (J41), take a single OR
site o, so C=B={o}. All other gates are AND. First suppose the internal
arc multigraph is strongly connected. Every nonempty fixed selected
set contains o, since a nonempty predecessor-closed set of AND sites
cannot exclude a vertex in a strongly connected graph. A proper fixed
set has exactly one incoming boundary arc, entering o; Eulerian balance
then gives exactly one outgoing boundary arc. Its intermediate particle
is at the head of that outgoing arc.

Let u_j->o, j=1,2, be the two incoming arc occurrences. Define A_j={o}
if u_j=o; otherwise let A_j consist of o and the ancestors of u_j in
G-o, including u_j itself. Each A_j is fixed: its selected nonroot
vertices contain all their predecessors, its chosen root input is
present, and its at most one outgoing boundary arc cannot activate an
outside AND gate. Every nonempty fixed set contains one A_j. Two
proper A_j are incomparable and join to V; if only one is proper,
it is the unique atom. The crosscut formula therefore gives

    mu_h=1_(A_1!=V)+1_(A_2!=V)-1 in {-1,0,1}.           (J49)

When neither incoming arc is a self-loop this becomes
mu_h=1-1_(u_1 reaches u_2 in G-o)-1_(u_2 reaches u_1 in G-o),
with length-zero reachability included. Thus two ancestor searches
determine the inverse coefficient in this subclass.

A finite Eulerian digraph has no arcs between distinct strongly
connected components: a source component has no entering arcs and
therefore no leaving arcs, and removal gives induction. Components
away from o have only AND gates and each contribute a two-element
factor to the fixed-cut lattice. Their Mobius values multiply by -1.
Consequently |mu_h|<=1 holds without strong connectivity as well.
This covers the specified two-record Boolean books only. Longer words
allow partial emissions, and the argument does not establish the
same bound for arbitrary one-particle remainders; (J64) below gives
a counterexample with longer words.

### Computing an endpoint fiber by IDLA and heap implications

**Proved; independent AI proof audit passed.** A fixed C-to-B book h
admits the following representation of every stable-cut fiber ending
at D. First run complete ordinary IDLA from 1_C+1_(V\D) on h's arrows,
ignoring its bits. Any D-flow cut would be a full stable upper bound
for this input, and its exit budget forces the legal IDLA run to
end full below that cut. Thus failure to finish full within the finite
book rules out the fiber. Otherwise its emitted vector q_D is the
least D-flow cut, with

    1_C+I(q_D)-q_D=1_D.

It need not satisfy the sleep-bit conditions. Every other D-flow cut
is q_D+a for a balanced internal prefix of the remaining tail. Even
though that tail need not be balanced, popping available head cycles
finds its unique greatest balanced prefix g. Indeed balanced prefixes
are closed under coordinate maximum by the incoming-count/exit argument
above; if any prefix extended beyond the popped g, its residual would
expose another head cycle. Exhausted heads and sink arrows act as roots
and are not popped. Possible a are exactly the lower ideals of g's
decorated simple-cycle heap.

At site i list its heap occurrences e_1,...,e_m in local order. A
selected prefix length k is allowed precisely when

    beta_i(q_D(i)+k+1_D(i))=1_D(i).                      (J50)

Use the original bits, including beta_i(0)=0 and the actual terminal
ones on B. A boundary bit can lie in q_D or beyond g; it is not reset.
Every required index is available: q_D+g is a D-flow cut below h,
so exhausting a site occupied in D forces it to be occupied in B.

Build a directed implication graph on the cycle occurrences, with
a forced source and a forbidden cemetery. Add each later-to-previous
edge e_(j+1)->e_j. For every forbidden k in (J50), add source->e_1
if k=0, e_k->e_(k+1) if 0<k<m, and e_m->cemetery if k=m>0.
If m=0 and its sole cut is forbidden, add source->cemetery. The
feasible selections are exactly outward-closed sets containing the
source and avoiding the cemetery. They exist iff the source does not
reach the cemetery; their minimum is the source-reachable occurrences,
and their maximum is the occurrences not reaching the cemetery.
Keep the original heap order when reconstructing the words.

This also finds all global lattice atoms. For D!=C retain the minimum
of each nonempty fiber. For D=C, q_C=0 and the empty cut is allowed;
its minimal nonempty selections are sink SCCs of the safe implication
graph. Remove candidates strictly dominating another candidate. The
survivors are exactly the global atoms. Compute their subset joins
by the upward iteration in (J40), then apply the crosscut formula.
There are at most 2^n fiber queries and 2^(A_n) atom subsets; for
fixed n this is polynomial in explicit book length, with no practical
or optimal dependence on n asserted.

The focused verifier compares this representation with direct cut
enumeration for all 64 endpoints of the M_3 book and all 16 endpoints
of the N_5 book. All 80 fiber queries agree. It enumerates closed
selections only inside these tiny extracted heaps to check the
representation; the general procedure uses reachability and SCCs.

The ordered surplus decomposition also extends to every terminal support,
without a full-filling assumption. For arbitrary all-active eta, let
Qcal_(eta,S) denote complete IDLA emission prefixes ending at S, with
unrestricted emitted-bit decorations. If Jcal_(eta,B) denotes formal
stable books and Acal_(eta,C) actual ARW emission books, then uniquely

    Jcal_(eta,B)=disjoint_union_S Qcal_(eta,S) Jcal_(1_S,B),
    Acal_(eta,C)=disjoint_union_S Qcal_(eta,S) Acal_(1_S,C),
    Jcal_(eta,B)=disjoint_union_(S,C)
       Qcal_(eta,S) Acal_(1_S,C) Rcal(C,B).               (J51)

Here Rcal is the remainder language, and products concatenate local
paired records. A formal ARW stable cut is an IDLA stabilizing upper
bound, so its complete IDLA prefix q lies below it. Crucially q_i=0
outside its final support S: IDLA never empties a visited site.
Subtracting balances gives the formal binary-S suffix. If that suffix
is unused at a site outside B, its balance forces the site outside S,
so q_i=0 too; no last-zero constraint leaks into the arbitrary q bits.

Every actual ARW cut is also an IDLA stable upper bound. Shifting by
q makes its buffer equation exactly that of all-active binary S.
For zero shifted arrivals, i is outside S and q_i=0, preserving the
beta_i(0)=0 convention. Least action in both directions therefore
identifies the original emission vector as q plus the shifted one.
This proves the same ordered triple whichever normalization is done
first. There is no chronological IDLA phase asserted inside ARW.

Weights multiply without an extra p_S factor at the IDLA interface.
If nu_eta is the ordinary IDLA output law, then Q_eta=nu_eta Q and
P_eta=nu_eta P; the marked versions retain the emitted currents of
each factor. This is a paired-record refinement of the existing CB
Section 5 reduction, with no novelty claim. The complete-support
extension, empty-site indexing, and weight accounting passed an
independent proof audit.


### Determinant geometry on the closed sleep cube

**Proved; independent AI proof audit passed.** Fix K, put n=|V|>=1,
N=2^n, m=N/2, d=det(I-K), and let J be the complement permutation
on subsets. Define the polynomial matrix and its determinant by

    U_p=diag(p_C)E_(q/p),  A(p)=det U_p.

The C row represents product_(i in C)(z_i-q_i J_i(z)) in the
squarefree addition-algebra basis, so this definition extends to
the whole closed cube and U_1=I. Taking determinants in (J29) gives

    det R_p=c_K A(p)A(1-p),
    c_K=det(F)det(J)/d^N>0,
    A(p)>0 for every p in [0,1]^V.                       (J52)

Indeed each coordinate occurs m times in product_C p_C. The same
count for the column factors mu_p(B)=p_B q_(B^c) cancels all
denominators. At p=1/2, (J48) implies c_K>0. The polynomial identity
and nonsingularity of R on the closed cube then show that A has no
zeros there; connectedness and A(1)=1 fix its sign.

Fix i and interior values of the other sleep coordinates. At p_i=1,
let W_i leave configurations without i unchanged; from a support
containing i, force its resident to make one K-jump and then stabilize
against the other sleeping residents. Let Q_i be W_i's block whose
rows and columns both contain i. This is a substochastic matrix of
size m. With s=r_i and E(s) the corresponding inverse basis, the
parameter-transport derivative gives

    E'(0)E(0)^(-1)=I-W_i,
    E(s)=[I+s(I-W_i)]E(0).

The second identity is exact because every row is affine in s.
Ordering supports without i first, the relevant block of the
bracket is (1+s)I-sQ_i; all other diagonal entries are one.
Since p_i=1/(1+s), row normalization therefore gives

    A(p)=A(p_i=1) det(I-(1-p_i)Q_i).                    (J53)

Here Q_i depends only on the other probabilities. A prescribed
simple path from i to the sink, routing one moving particle and
restoring residents along the way, leaves i vacant with probability
at least its positive arrow weight times the product of
min(p_j,1-p_j) over its intermediate vertices. This bound works for
every initial support containing i. Thus rho(Q_i)<1 and

    log det(I-tQ_i)=-sum_(k>=1) t^k tr(Q_i^k)/k,
    d^2/dt^2 log det(I-tQ_i)
       =-sum_(k>=2)(k-1)t^(k-2)tr(Q_i^k)<=0, 0<=t<=1.

Nonnegativity of Q_i suffices; symmetry is not assumed. Equations
(J52)--(J53) prove that log det R_p is separately concave. Positivity
and continuity extend the concavity inequalities to all boundary
slices. If i lies on an internal directed cycle, starting from full
support and routing the forced particle around that cycle shows
Q_i(V,V)>0. Hence the k=2 term makes the concavity strict when the
other probabilities are interior. If i lies on no internal cycle,
the forced activity cannot return to i, so Q_i=0 and the determinant
is independent of p_i, including on boundary slices by continuity.

All corners have the value in (J38). Successive one-coordinate
concavity inequalities therefore give

    det R_p >= product_(D subset V) det(I-K_DD)^(-1)
             for every p in [0,1]^V.                    (J54)

This is a global lower bound from separate concavity. Joint
log-concavity and maximality at p=1/2 are not asserted.


### A fixed-book target sign already fails with one particle

**Proved; independent AI audit and exact certificate passed.** Let
O_v(C,C)=1_(v in C). The conjugated target projector G_v=R^(-1)O_v R
has an absolutely convergent book expansion with coefficient

    gamma_v(h)=sum_(m in L(h)) mu(0,m) 1_(v in D(m)),
    G_v(C,B)=sum_(h:C->B) W_p(h) gamma_v(h).             (J55)

This follows by concatenating the prefix in R^(-1) and suffix in R;
absolute convergence is bounded by M_n R(|t|)^2 with an arrow mark t.
A tempting stronger sign rule is gamma_v(h)<=0 when v is initially
absent and gamma_v(h)>=0 when it is initially present. It fails.

Take four sites 0,1,2,3, C={0}, B={3}, and paired rows

    h_0=((1,1),(0,2),(0,3)),
    h_1=((1,0),(0,1)), h_2=((1,0),(0,2)), h_3=empty.

Append the terminal one at site 3. The four stable cuts are exactly

    emitted vector   endpoint   mu(0,m)
    (0,0,0,0)        {0}         1
    (2,2,0,0)        {2}        -1
    (2,0,2,0)        {1}        -1
    (3,2,2,0)        {3}         1.

Indeed m_1,m_2 are either zero or two. Root arrival counts one,
two,three yield root emissions zero,two,three. With exactly one
nonroot site selected its self-loop and the relevant root arrow
supply its two arrivals; with both selected the root emits all
three arrows. These are the displayed Boolean diamond and exhaust
the possibilities. Thus (gamma_0,gamma_1,gamma_2,gamma_3)=(1,-1,-1,1),
and the initially absent target 3 has a positive coefficient.

For a concrete accessible kernel, let site 0 choose 1,2,3,sink with
probability 1/4 each, sites i=1,2 choose 0,i,sink with probability
1/3 each, and site 3 jump to the sink. The book's positive weight is
p_3 q_0^2 p_1 q_1 p_2 q_2/5184. At all p_i=1/2, however, the full
sum satisfies G_3({0},{3})=-5625/238336. This exact value follows
from (J29) and the singleton block of R; the focused verifier checks
both all 36 candidate cuts and this rational matrix entry.

Consequently this obstruction concerns cancellation within one fixed
book, even with one particle. It does not disprove a sign after
averaging books, and G_v's sign by itself has not been identified
with the original add-one inequality.

In fact the offending term has an equal-weight partner. Exchange the
two zero records at site 0, giving h_0=((1,1),(0,3),(0,2)). The cuts
are now the chain (0,0,0,0), (2,2,0,0), (3,2,2,0), with endpoints
{0},{3},{3} and Mobius values 1,-1,0. Its target coefficient is -1.
These are the only two admissible local orderings of these paired
record multisets, so the grouped target-3 coefficient is zero.
The focused certificate checks the 36 cuts of each book. This
suggests grouping local permutations; no general grouped sign is
claimed.

For comparison, [Fomin's loop-erased path switching, Theorem 6.1
and its two-path proof](https://arxiv.org/pdf/math/0004083), cancels
signed walk families by exchanging tails at a canonically recoverable
intersection with a loop-erased path. Sections 5--6 and Theorem 7.1
were read. Adapting that choice to paired records would additionally
have to preserve their stable-boundary bits and recover the selected
cut. This is a possible mechanism, not an application of that theorem.


### One-particle joins and a scoped unit-Mobius theorem

**Proved; independent AI derivation checked in the main thread.** Let
h:{r}->{r} be a zero-exit remainder, with arbitrary words, and write
D(m)={v(m)} for the singleton endpoint of a stable cut. Call an atom
balanced if its endpoint is r, and transporting otherwise.

For two cuts a,b, their coordinate maximum z has nonnegative mass
at each site: use the old cut supplying z_i and monotonicity of
incoming counts. Conservation makes that mass a singleton. Every
empty site's last bit is zero, inherited from a cut supplying that
coordinate. The only possible instability is therefore at the one
occupied site. Follow zero-bit arrows from there until the next bit
is one. The full book bounds this walk, and the resulting cut is
a join b by least action. Thus

    join = coordinate maximum followed by one zero-bit walk.    (J56)

In particular this correction never emits a new one-bit record.

If there is at most one transporting atom, then

    mu_h is either 0 or (-1)^(number of atoms).                 (J57)

Here is the atom-sum proof. Balanced atoms have pairwise disjoint
positive coordinate supports, and their joins are coordinate unions.
For any cut j ending at v!=r and a balanced atom b, if j_v>=b_v,
the coordinate minimum min(j,b) has mass at most zero off r and
at most one at r. At v, use the coordinate supplied by b; at other
nonroot sites either old cut supplies a zero-mass bound. Conservation
forces equality, and the inherited bits make this a stable r-cut.
It is zero or b, since b is an atom. Consequently a balanced atom
not below j can intersect j only if its support contains v. At most
one balanced atom has that property.

If every atom is balanced, only their full collection can join to
the top, proving (J57). Otherwise let a be the unique transporting
atom. It intersects at most one balanced-atom support by the preceding
observation and atom minimality. Thus a cannot lie below the union
of balanced atoms: it would lie below one of them. Every top-joining
atom subset must include a.

Start j=a and declare a required. If an unrequired balanced atom is
already below j, toggling it pairs all top-joining subsets with
opposite signs, so mu_h=0. Otherwise, when v(j)!=r, every unrequired
atom except possibly the unique one containing v(j) is disjoint from
j, avoids v(j), and can be added without changing its endpoint.
A top-joining subset must include that exceptional atom, since the
top ends at r. If it does not exist, there is no such subset. If it
exists, declare it required and replace j by their join. When j
ends at r, all remaining atoms are disjoint from j and from each
other. Every one is necessary to reach the top. The finite procedure
therefore either pairs the alternating sum, finds no top-joining
subset, or leaves just the full atom set. This proves (J57).

A concrete sufficient condition allows arbitrary words away from r.
Suppose r has exactly two records with bits (1,0), and its first
arrow follows a simple return path to r using first arrows whose
bits are zero at all nonroot sites; a root self-loop is allowed.
Make r absorbing, release both root arrows, and stabilize the other
sites. Their full book is a formal stable upper bound ending empty,
so this has a least emitted vector q within the given records.
Following the protected first return before releasing the second
particle proves that at least one particle reaches r. Hence the
total absorbed is one or two. Restoring root emission count two
makes (2,q) a stable original cut, ending respectively away from r
or at r. Every root-positive cut has root coordinate two and
dominates this one. A root-zero cut necessarily ends at r. There
is therefore at most one transporting atom, and (J57) applies.

The protected return is an explicit hypothesis. Without it the
auxiliary two-particle run can leave two sleepers and return neither
particle, so (2,q) need not be an original stable cut. The general
one-particle unit-Mobius bound is false by (J64). If the root is unused,
has only its forced first one-bit record, or has two one-bit records,
every cut ends at r and the distributive-fiber argument already
gives the unit bound.


### Reversal of the inverse and the exact target commutator

**Proved; independent AI proof audit passed.** The local reversal and
bit complementation h->h^dagger from (J28) reverse concatenation:
(ab)^dagger=b^dagger a^dagger. A cut m becomes ell-m, with endpoint
D(m)^c. Thus L(h^dagger) is the order dual of L(h), its full Mobius
number is unchanged, and (J28) holds for the inverse as well:

    mu_p(C) R_p^(-1)(C,B)
       =mu_p(B) R_q^(-1)(B^c,C^c).                      (J58)

No arrows change direction; this is formal record reversal, not
physical reversal of a legal history. Put M_p=diag(mu_p) and
define D_p(X)=M_p^(-1) J X^T J M_p. This map reverses products.
Since D_p(O_v)=I-O_v, the target operator in (J55) satisfies

    G_v(p)=I-D_p(checkG_v(q)),
    checkG_v(q)=R_q O_v R_q^(-1).                       (J59)

The opposite sandwich is essential. Equivalently, if gamma_v^+(h)
uses upper-interval weights mu(m,ell) in place of mu(0,m), then
gamma_v(h^dagger)=1_(h is an identity)-gamma_v^+(h).
These identities hold also with common arrow mark t throughout the
absolute-convergence disk; transpose is not complex conjugation.

Here is the precise interface with the original inequality. Let
A_x^p be physical addition at x, f=R_p^(-1)1 the positive column
in (J30), and Q_eta=P_eta R_p the formal stable-book row. Define
Xi_x=R_p^(-1) A_x^p R_p. Addition duality gives
Xi_x=D_p(A_x^q)>=0 and Xi_x f=f. Therefore

    q_v(eta+delta_x)-q_v(eta)
       =Q_eta (Xi_x-I)G_v f
       =Q_eta [Xi_x,G_v]f.                              (J60)

The second equality uses Xi_x f=f. Also G_v^2=G_v,
Q_eta G_v=P_eta O_v R_p>=0, and Q_eta(I-G_v)>=0.
These are positivity statements on represented physical laws;
they do not make G_v entrywise positive or settle the commutator.
Replacing f by the all-ones column is generally incorrect.

By (J51), Q_eta=nu_eta Q with Q=P R_p, so (J60) also equals
nu_eta(D_x-I)P 1_(v occupied), where D_x is ordinary IDLA addition.
This retains the all-active scope. Normalizing R to the stochastic
link L=R diag(f) recovers the scalar payoff L^(-1)1_(v occupied),
whose proper-target positivity and inclusion-order obstructions
were already established in Section 36. Reversal does not remove
them.


### Grouping a whole local inventory still need not fix the target sign

**Proved; independent AI audit and exact certificate passed.** The
two-book cancellation after (J55) does not extend to arbitrary local
inventories. Start with the six-site M_3 book (J41), add a fresh
target 6, and append a zero-bit arrow 0->6 at site 0. Explicitly,

    C={0,1,2}, B={1,2,6},
    arrow rows=((3,5,6),(4,3),(5,4),(0,1),(1,2),(2,0),()),
    bit rows=((1,0,0),(1,0),(1,0),(1,0),(1,0),(1,0),()).

The displayed book still has an M_3 cut lattice: its three proper
atoms are unchanged, and its top now emits three arrows at 0 and
ends at B. Target 6 occurs only at the top, so gamma_6=2.

Now fix only each site's arrow multiset and its numbers of zero and
one bits, allowing their pairing to change. Boundary constraints
force exactly the displayed bit words. There are 3!*2^5=192 local
arrow orders, all valid remainders with the same weight W. At the
five other used sites a stable emission coordinate is zero or two,
so both arrows are used together and their order is irrelevant.
Only the order of the three arrows at 0 matters.

If 0->6 is last, the cut lattice and target coefficient remain M_3
and 2. If it is first or second, root emission two sends one arrow
to 6 and one to an internal edge-site. The fixed cuts are a Boolean
diamond with a new top above it. Target 6 occurs at one atom, their
join, and the new top, whose Mobius values are -1,1,0. Hence its
coefficient is zero. This also follows directly from the other
sites' OR/AND equations: choosing the first internal edge gives
two root-positive proper cuts, while the old atom avoiding root 0
remains the other atom.

There are 64 orders with 0->6 last and 128 others. Consequently

    sum_(h in this inventory) gamma_6(h)=128>0,
    sum_(h in this inventory) W_p(h)gamma_6(h)=128W>0.   (J61)

Any accessible K positive on the displayed arrows and any interior
sleep vector makes W positive. The focused verifier exhausts this
one 192-book inventory and checks the coefficient in every order.
Even if paired records are kept intact, the only two admissible
orders have coefficients 2 and 0, so that narrower grouping fails
as well. This rules out the proposed initially-vacant target sign
at either of these finite groupings. It does not determine the sign
of the full sum G_v, or of the add-one commutator (J60).

The grouping has an exact commutative-algebra formulation. Give each
paired type (i,b,j) a variable x_(i,b,j), and let Rscript(C,B) sum
the corresponding inventory monomials over remainder books. Then
Rscript(0)=I and

    [x^N]Rscript^(-1)(C,B)=sum_(h with inventory N)mu_h,
    [x^N](Rscript^(-1)O_v Rscript)(C,B)
       =sum_(h with inventory N)gamma_v(h).             (J61a)

Distinct words, not labeled permutations of repeated records, are
counted. No extra factorial belongs in these formulas. The physical
matrix is diag(p_C)^(-1) Rscript(x) diag(p_C) after setting
x_(i,b,j)=p_i^b q_i^(1-b)K(i,j). To forget bit-arrow pairing first,
substitute x_(i,b,j)=u_(i,b)w_(i,j); (J61) is a positive coefficient
even after this coarser substitution.

There is a compact constant-term expression. Put z_sink=1,
a_i=sum_j x_(i,1,j)z_j/z_i, b_i=sum_j x_(i,0,j)z_j/z_i. Then

    Rscript(C,B)=CT_z z^(1_C-1_B)
       product_(i in C\B)a_i b_i
       product_(i in C intersect B)(1-b_i)
       product_(i outside C union B)(1-a_i)
       / product_i(1-a_i-b_i).                         (J61b)

Expand in record variables before taking the constant term. The
four local word series are respectively ab/(1-a-b), 1/(1-a-b),
(1-b)/(1-a-b), and (1-a)/(1-a-b) on C\B, B\C, their intersection,
and their joint complement. Flow rules out an unused C\B site.
This proves (J61b) directly and retains ordering multiplicities.

The four-site pair after (J55) also separates the inverse and target
coefficients: its full Mobius numbers sum to 1+0=1 in the off-diagonal
singleton inverse entry, while its target coefficients sum to zero.
Thus a coefficientwise nonpositive singleton inverse with independent
paired-record markers is already impossible. Numerical inverse signs
after summing all inventories are a different question.


### An explicit complex neighborhood without zeros

**Proved; independent AI proof audit passed.** The finite-book bound
also controls joint complex changes of the sleep parameters. Fix
the original real nonnegative transient K. For p in C^V put
q_i=1-p_i, c_i=|p_i|+|q_i|, and C_c=diag(c_i). Then

    rho(C_c K)<1  implies  det R_p != 0.                (J62)

Here R_p denotes its polynomial extension (J37), not a stochastic
law at complex p. The inverse series (J45) converges absolutely in
this region.

For a fixed arrow book, cancel the first forced-one factor against
the denominator at every used site of C. At unused sites in C,
flow balance forces membership in B, so the two boundary factors
cancel directly. Thus each individual normalized book weight is
a polynomial in p_i and q_i, even at p_i=0. Summing its absolute
value over all admissible bit assignments is at most

    product_(i in B\C)c_i * product_i c_i^(ell_i)
       * (the original arrow weight).

For example, a used site in C\B has distinct forced first-one
and last-zero records: after cancellation its bit sum is
|q_i| c_i^(ell_i-2)<=c_i^ell_i. A site in B\C contributes its
terminal |p_i| times the free-bit sum. The other cases have no
larger factor. We used c_i>=1 throughout.

The last display is bounded by product_i c_i times the arrow-book
weight with internal matrix C_c K and sink weights c_i k_i. Its
total is finite under (J62), by the fixed-exit geometric estimate:
the number of sink arrows is |C|-|B|<=n. Multiplying by the uniform
bound M_n gives absolute convergence of the inverse series.
Multiplicativity of the normalized weights holds polynomially,
since it holds where every p_i is nonzero. The two convolution
cancellations from (J45) therefore remain valid, proving (J62).

In particular, if max_i dist(p_i,[0,1])<=epsilon, then
c_i<=1+2epsilon. Thus a uniform zero-free neighborhood of the
whole real cube is given by

    (1+2epsilon)rho(K)<1.                              (J63)

For an acyclic K, C_c K remains nilpotent for all finite p, so this
argument imposes no complex-sleep restriction. These are sufficient
zero-free regions, not claims of multivariate real stability.


### One particle can already have an M_3 cut lattice

**Proved; independent AI construction, main proof check, and exact
certificate.** The general extension of (J49) or (J57) fails even
when the root still has only two records. On five sites take
C=B={0} and

    h_0=((1,1),(0,2)),
    h_1=((1,3),(1,4),(0,0),(0,0)),
    h_2=((1,1),(0,2)),
    h_3=((1,1),(0,3)), h_4=((1,1),(0,4)).

The incoming and outgoing counts are both (2,4,2,2,2), and every
nonroot last bit is zero. Append the root's terminal one. The
root and leaf emission coordinates are zero or two. Write a cut
as (2R,H,2b_2,2b_3,2b_4). The hub's arrival-to-emission function
is phi(0),...,phi(4)=(0,0,1,3,4). The cut equations are exactly

    R=1_(H>=3), b_2<=R, b_3<=1_(H>=1), b_4<=1_(H>=2),
    H=phi(R+b_2+b_3+b_4).

If R=0, the only possibility is H=0 and all leaves zero: H=1
would require two hub arrivals, but only leaf 3 could return.
If R=1, either H=3 with two selected leaves, or H=4 with all
three. Thus the complete cut list is

    cut                    endpoint    mu(0,m)
    (0,0,0,0,0)            {0}           1
    (2,3,0,2,2)            {2}          -1
    (2,3,2,0,2)            {3}          -1
    (2,3,2,2,0)            {4}          -1
    (2,4,2,2,2)            {0}           2.             (J64)

The lattice is M_3 and all three atoms transport the unique particle.
They use the same first return arrow to the root, so assigning atoms
to the root's two incoming arrow occurrences cannot bound their
number by two. Any accessible kernel positive on the displayed
arrows realizes this positive-weight book. The verifier exhausts
all 405 candidate cuts; no minimality is claimed.

This is precisely the two-sleeper case excluded by the protected
return condition after (J57): release the two root arrows with the
root absorbing, and both particles immediately sleep, at the hub
and leaf 2. The formal alternatives return through any two of the
three leaves. The hub's second one-bit prevents a smaller balanced
vacuum cut, so there is no balanced atom to peel off. The Boolean
two-record subclass and the one-transporting-atom theorem remain
valid.

The construction extends to every M_k, and the transporting-atom
bound n-1 is sharp. For k>=2 use vertices 0,1,...,k, root records
((1,1),(0,2)), and leaf records ((1,1),(0,i)) at each i=2,...,k.
At the hub 1 use, in order,

    (1,3),...,(1,k),(1,0),(1,1),(0,1),(0,0).

The initial range is empty for k=2. Its bit word is 1^k00, and
its two self-loops lie between its two root-return arrows. This
balanced book has lengths (2,k+2,2,...,2), with C=B={0}.
Write a cut (2R,H,2b_2,...,2b_k). The root is selected iff
H>=k-1, leaf 2 requires R=1, and leaf i>=3 requires H>=i-2.
The hub input is

    N=R+sum_i b_i+1_(H>=k)+1_(H>=k+1).

Its emission function is zero at N=0, N-1 for 1<=N<=k, and N
for N=k+1,k+2. If R=0, at most H leaves can return, ruling out
positive H. If R=1, the only possibilities are: H=k-1 with all
leaves selected; H=k+1 with exactly one leaf unselected; and H=k+2
with all selected. These give zero, k incomparable atoms, and top.
The first atom ends at the hub and the others at their missing
leaves. Thus

    L(h)=M_k, number of transporting atoms=k=n-1,
    mu_h=k-1=n-2, total record count=3n-1.             (J64a)

Distinct transporting atoms cannot have the same nonroot endpoint
by the fiber-minimum lemma, proving sharpness of their number.
No universal sharp bound on Mobius values or record minimality is
claimed. For k=3 this gives the four-site book
h_1=((1,3),(1,0),(1,1),(0,1),(0,0)), with the root and leaf rows
as above. Its atoms are (2,2,2,2), (2,4,0,2), and (2,4,2,0),
and its top is (2,5,2,2). The symbolic family received an independent
AI derivation and main proof check.

These lattice examples do not depend on permitting self-loops.
Replace every occurrence of an arc u->v by u->w, with w new,
and give w the same number N of records, all (0,v). Keep w out
of both endpoints. If c(t) counts replaced occurrences among the
first t records at u, every new stable cut satisfies

    m_w=c(m_u), endpoint at w vacant.                  (J64b)

Substitution recovers exactly the old cut equations. Projection
and m->(m,c(m_u)) are inverse order isomorphisms, preserving all
interval Mobius numbers and the atom types. This applies also
to u=v. Redirect K(u,v) to w and set K(w,v)=1; accessibility of
the sink is preserved. Any interior p_w gives these N zero bits
positive probability. One new site per loop vertex therefore
removes recorded self-loops without changing the lattice. This
is a statement about these specified books, not equality of the
two unrestricted ARW laws. The subdivision proof was independently
derived and checked.


### Component factors and exact directed-cycle determinants

**Proved; independent AI derivation, main proof check, and focused
exact certificates.** Let S_a be the strongly connected components,
n_a=|S_a|, n=|V|. In each local kernel K_a=K_(S_a,S_a), treat all
outgoing transitions as exits. Then

    A_K(p)=product_a A_(K_a)(p_(S_a))^(2^(n-n_a)),
    det R_K(p)=product_a det R_(K_a)(p_(S_a))^(2^(n-n_a)). (J65)

Order components topologically and filter supports by their vectors
of component counts. Counts in every initial segment cannot increase,
so R is block triangular. Preserving the whole count vector forces
zero cross-component transfers and zero losses, successively from
the first component onwards. Its diagonal blocks are therefore tensor
products of the local fixed-count blocks. Taking determinants and
summing spectator dimensions gives the exponent 2^(n-n_a). The
same argument for physical P, followed by inversion and the row
normalization defining U, proves the A formula. Polynomial continuity
includes the faces. A loopless singleton contributes one, although
it doubles the exponent of every other component factor. A singleton
with self-loop a contributes A=1-aq and det R=(1-ap)(1-aq)/(1-a)^2.

Now take a directed simple cycle of n>=2 sites, with forward weights
a_i=K(i,i+1)>0, all other transitions exiting. Put kappa=product_i a_i<1
and d=1-kappa. In the homogeneous algebra x_i^2=a_i x_i x_(i+1),
let B_k(p) be the degree-k matrix taking x_C to
product_(i in C)(x_i-q_i a_i x_(i+1)), and Delta_k=det B_k. Then

    A_K(p)=product_(k=0)^n Delta_k(p),
    Delta_0=1, Delta_1=1-kappa product_i q_i,
    Delta_n=d+kappa product_i p_i.                     (J66)

Every intermediate B_k also has an explicit rule. Split C into maximal
occupied cyclic runs. In a run [a,...,b], either do nothing with
weight one, or choose j in the run, remove j, insert the following
vacancy b+1, and give this choice weight

    -q_j product_(t=j)^b a_t product_(t=j+1)^b p_t.

Multiply the choices from distinct runs. Indeed, in the expansion,
let j be the first selected forward term in a run. Reduction follows
the same arrows j,...,b regardless of subsequent selected terms;
summing their signs and q weights gives the displayed product of p's.
Different runs end at different original vacancies and do not interact.
For a full C, every nonempty selected set instead makes one complete
turn, giving the last formula in (J66).

For equal-size supports C,B define the nonnegative integer height
vector h(C,B) by

    h_i-h_(i-1)=1_C(i)-1_B(i), min_i h_i=0,
    H_k(C,B)=product_i a_i^h_i, |C|=|B|=k.

Then the count-preserving remainder block is

    R_K[k,k]=I+[kappa/d] H_k diag_(|B|=k)(p_B q_(B^c)).  (J67)

Equal counts forbid exits. All arrow-length vectors are h+t*1. If
C!=B, a rise from height zero to one occurs at a site of C\B;
its forced first one and last zero exclude t=0. Every t>=1 is
admissible, and summing its free bits gives exactly p_B q_(B^c).
The arrow weights form a geometric series. For C=B, t=0 supplies
the identity and the same tail gives the diagonal. This proves
(J67), including k=0,n, without ARW stabilization calculations.

Writing b_k=binom(n,k), the resulting determinant identities are

    det H_k=d^(b_k-1), Delta_k(0)=d for 1<=k<=n,
    det R_K[k,k]=Delta_k(p)Delta_(n-k)(q)/d,
    det R_K(p)=A_K(p)A_K(q)/d^(n+1), A_K(0)=d^n.        (J68)

To verify the nontrivial height identity, the homogeneous pairing
F(C,B^c) equals H_k(C,B): it is the probability of the unique least
full-filling IDLA prefix h. The degree blocks of (J29) therefore give
det R[k,k]=d^(-b_k) det H_k Delta_k(p)Delta_(n-k)(q), with all
sleep boundary factors canceling. For 0<k<n, set p=1 and use
(J67) to get det H_k Delta_(n-k)(0)=d^b_k. Gauge the forward weights
to (1,...,1,kappa), which conjugates both matrices. Their determinants
are polynomials in kappa. At kappa=1, H_k is all ones, so its
determinant is divisible by d^(b_k-1). Also B_(n-k)(0) has zero row
sums there: only the last site in each run can move, and summing
run choices gives product_runs(1-1). Thus its determinant has a
factor d. The product identity forces both remaining factors to be
constant; at kappa=0 the acyclic triangular matrices have determinant
one. This proves (J68), with extreme sectors immediate. The common
corner determinant is 1/d. Redistributing forward weights with the
same kappa only conjugates these blocks.

For example, on a three-cycle put P=product_i p_i, Q=product_i q_i,
L=Q+p_1q_2q_3+q_1p_2q_3+q_1q_2p_3. Then
A=(1-kappa Q)(1-kappa L-kappa^2 P Q)(d+kappa P).
Intermediate sectors need not have a single linear factor in kappa.

Self-loops can also be removed exactly, including their normalization.
If sigma_i=K(i,i), let ell_i=1-sigma_i and let barK divide all other
routing and exit weights at i by ell_i. Put b_i(p)=1-sigma_i q_i,
barp_i=p_i/b_i(p), and hatp_i=q_i/b_i(q). In the identical addition
algebra, z_i-q_i J_i^K=b_i(p)[z_i-(1-barp_i)J_i^barK], so

    A_K(p)=product_i b_i(p)^(2^(n-1)) A_barK(barp).

For a cycle after removal, bar kappa=product_i a_i/ell_i, and

    det R_K(p)=product_i [b_i(p)b_i(q)/ell_i^2]^(2^(n-1))
       * A_barK(barp)A_barK(hatp)/(1-bar kappa)^(n+1).   (J69)

The effective parameters barp and hatp are generally not complementary.
The constant in (J52) scales by product_i ell_i^(-2^n), since F
is unchanged while det(I-K) acquires product_i ell_i. This proves
(J69). The focused certificate checks all five count sectors of one
four-cycle, their height determinants and zero-sleep factors, and
the SCC factorization on two cycles joined by a directed edge.


### A positive first-cut kernel and single-state elimination

**Proved; independent AI proof audit passed.** Fix a support D. In
this subsection F denotes a first-cut kernel, not the Frobenius matrix
of Section 36. Let F(C,D) be the total normalized
weight of C-to-D remainders with no proper D-cut, and set F(D,D)=1
for the empty book. Any book possessing a D-cut has a unique least
one, since its D-fiber is closed under coordinate minimum. Splitting
there gives a first-D prefix followed by an arbitrary D-to-B remainder.

Conversely, concatenate such a prefix a and an arbitrary suffix b.
Any D-cut m of ab has min(m,ell(a)) as a D-cut inside a. Minimality
forces this minimum to be ell(a), so ell(a)<=m. Thus concatenation
introduces no earlier incomparable D-cut, and the factorization is
bijective. Multiplicativity of weights now gives

    F(C,D)=R(C,D)/R(D,D),
    R(C,B)-R(C,D)R(D,B)/R(D,D)
       =sum_(h:C->B with no D-cut) W_p(h)>=0.           (J70)

If C=D or B=D, the right side is zero. Otherwise this is exactly
the Schur complement eliminating index D. In particular its diagonal
entries are at least one, from the empty book, and every two-index
principal minor satisfies
R(C,C)R(D,D)-R(C,D)R(D,C)>=max(R(C,C),R(D,D))>0.
The factorization also proves positivity before commutative grouping;
absolute convergence justifies the displayed numerical sums.

The normalized remainder law (J30), assigning mass W_p(h)f_B to
(B,h:C->B), makes F an actual cut-visitation kernel:

    F(C,D)=P_C(the sampled book admits a stable D-cut),
    0<=F(C,D)<=1-f_C<1 for C!=D, F(C,C)=1.             (J71)

Indeed the total suffix mass is (Rf)_D=1. Conditional on visiting D,
the first-D prefix has law W_p(a)/F(C,D), independently of the
suffix with its full normalized remainder law starting from D.
The strict bound follows by summing the avoidance term in (J70)
against f and retaining the identity C-to-C book of mass f_C>0.

Consequently every column is maximized at its diagonal:
R(C,B)<R(B,B) for C!=B in the interior. Diagonal books are balanced
arrow books with additional bit restrictions, so R(B,B)<=1/d with
d=det(I-K). Both bounds extend to the closed cube by continuity.
At the corner p=1_B, that diagonal has all its effective range
weights equal to one, giving R(B,B)=1/d. This proves the explicit
uniform maximum used after (J48).

There is consequently a directed distance on this support space:
put delta(C,D)=-log F(C,D), with -log 0=infinity. Equation (J70)
gives F(C,D)F(D,B)<=F(C,B), so delta satisfies the triangle inequality
and is positive off the diagonal. On every class of mutually positive
entries its symmetrization is the finite metric

    [delta(C,D)+delta(D,C)]/2
       =(1/2)log[R(C,C)R(D,D)/(R(C,D)R(D,C))].

This is the transition inequality and logarithmic-distance construction
in [Chebotarev, *The graph bottleneck identity*, Sections 2--3,
Theorem 1 and its footnote](https://arxiv.org/pdf/1003.3904), which
were read. Here equality in (J70) means every positive-weight C-to-B
book admits a D-cut. An equivalence with bottlenecks in an ordinary
directed graph has not been established.

The conditional splitting above is not a chronological Markov-chain
representation. In the four-site diamond (J55), split at the cut
ending at {2}. Both factors avoid {1}, but their concatenation has
the incomparable cut ending at {1}. Thus books avoiding one support
are not closed under concatenation. Iterating (J70) inside an avoided
language would require another argument. This is also why visiting
a union of endpoint states need not have a uniquely earliest cut.


### Singleton inverse geometry and a cycle target-sign regime

**Proved; independent AI derivations and main proof checks.** The
singleton remainder block has a useful exact reduction. Put
M=I-K, d=det M, Z=M^(-1), g_i=Z(i,i), and let D_p,D_q be the
site diagonals of p,q. Write P_(p,k), E_(p,k) for the particle-count
k diagonal blocks of binary stabilization and its inverse, with
sleep vector p. Reindex the complementary (n-1)-block by holes:
barP_q(i,j)=P_(q,n-1)(V\{i},V\{j}), and similarly barE_q.
If D_mu(i,i)=p_i product_(l!=i)q_l, then

    E_(p,1)=D_p^(-1)(I-D_q K),
    Hhit(i,j)=Z(i,j)/g_j,
    R_1=d^(-1) E_(p,1) Hhit barE_q^T D_mu,
    R_1^(-1)=d D_mu^(-1) barP_q^T diag(g) M P_(p,1).    (J72)

Indeed Hhit is the homogeneous pairing of x_i and x_(V\{j}):
park the latter particles and the surplus walks until hitting j
or exiting. Only degrees 1 and n-1 can contribute in (J29), proving
the displayed factorization. All count-block determinants of E
are positive by invertibility and continuity from p=1; det Hhit
=1/(d product_i g_i)>0. Thus det R_1>0.

For n<=3, R_1^(-1) is a nonsingular M-matrix. For n=3 and distinct
i,j,k, the cofactor formula and (J70) give

    R_1^(-1)(i,j)
       =-R_1(k,k) U_{ {k} }({i},{j})/det R_1<=0,

where U_D is the no-D-cut weight in (J70). Diagonal cofactors are
positive by the same result; n<=2 is immediate. Its inverse R_1
is nonnegative, completing the assertion. In the three-site case
an off-diagonal entry is strictly negative whenever R_1(i,j)>0:
take a positive-weight arrow book i->j, set its first source bit
to one and every other emitted bit to zero. This keeps its required
outer boundary and forbids a stable {k}-cut.

For larger n, (J72) isolates the unresolved inequality. Put
A=I-D_q K. Since M A^(-1)=I-D_p K A^(-1), the desired off-diagonal
sign is precisely, for i!=j,

    g_j barP_q(j,i)
       <=sum_z g_z barP_q(z,i) p_z [K A^(-1)](z,j).

The nonnegative left multiplier barP_q^T in (J72) need not preserve
M-matrix signs by general matrix algebra. The nonclosure after
(J71) also prevents an automatic iterated-Schur proof. The existing
four-site example (J55) passes the numerical inverse sign, but no
general n>=4 conclusion is claimed.

On a directed cycle, the height formula gives a different sufficient
condition in every fixed-count sector. For 0<k<n, use (J67), write
H=H_k, M_mu=diag(mu_B), tau=kappa/(1-kappa), and
s_k=sum_(|B|=k)mu_B<=1. For any C,D,B of size k,
h(C,D)+h(D,B)=h(C,B)+m*1 with an integer m>=0. Therefore

    H(C,D)H(D,B)=kappa^m H(C,B)<=H(C,B).

If A=tau H M_mu, this gives A^2<=tau s_k A. When tau s_k<1,
the row-sum bound ||A||<=tau s_k permits the paired Neumann series
I-(I+A)^(-1)=sum_(j>=0)A^(2j)(A-A^2)>=0. Continuity, applied first
to lambda A for lambda<1, includes tau s_k=1. Thus

    tau s_k<=1 implies R[k,k]^(-1)(C,B)<=0 for C!=B,
    (2O_v-I)(R[k,k]^(-1)O_v R[k,k]-O_v)>=0 entrywise.   (J73)

For an initially absent target the latter follows by summing inverse
entries over supports containing it; for an initially present target
subtract the complementary sum from I. In particular (J73) holds
in all sectors and at all sleep vectors when kappa<=1/2. The extreme
sectors have G_v=O_v. At sleep corners the nonzero correction to I
in (J67) has rank one, and its explicit inverse gives the same signs
for every kappa<1. These are sufficient regimes for fixed-count
blocks, not the lower-count entries of the full G_v or add-one
monotonicity.

The cycle's inverse height matrix also has an exact exclusion form.
For 0<k<n let B_0=B_k(0). At a support C choose any subset T of
the pairwise disjoint edges i->i+1 with i occupied and i+1 vacant,
and move their particles simultaneously. The coefficient of the
resulting support is (-1)^|T| product_(i in T)a_i. This is B_0 by
the zero-sleep run formula. The same edges move the complemented
output to the complemented input, so the pairing identity implies

    H^(-1)=B_0/d, Q=B_0+kappa M_mu,
    R^(-1)=I-kappa Q^(-1)M_mu,
    G_v-O_v=[kappa/d]Q^(-1)[B_0,O_v]H M_mu.             (J74)

This signed parallel-move operator is not stochastic: simultaneous
two-particle moves have positive coefficients, whereas one-particle
moves have negative coefficients. Its sparse commutator alone does
not have the target sign. The later chain argument (J107) removes
(J73)'s return-mass restriction for loop-free directed cycles,
including all particle-loss sectors.


### Prefix probabilities, stopping cuts, and overlapping visits

**Proved.** The normalized law in (J30) has a stronger prefix
property than the first-support statement. For a fixed remainder
a:C->D, let [a] be the event that the sampled book has a as a
local paired-record prefix, with its full length vector a stable
cut. Splitting at that vector gives a unique arbitrary suffix.
Consequently

    P_C([a])=W_p(a),
    law(suffix | [a])=P_D.                            (J75)

Indeed its mass is W_p(a) sum_(B,b:D->B)W_p(b)f_B=W_p(a).
The conditioning statement includes the terminal support. Thus
these cylinder probabilities multiply under concatenation. More
generally, any family of prefixes whose cylinders are disjoint
can be selected with their weights W_p(a); conditional on the
selected prefix, the suffix again has its normalized law.

This resembles the cylinder chain rule for Markov measures in
[Abbes, *Markovian Dynamics of Concurrent Systems*, introduction,
equation (1)](https://arxiv.org/pdf/1505.05536), which was read.
That paper works with trace-monoid actions. Our category need not
have unique minimal common extensions, so its trace-monoid
realization and states-and-cliques results have not been applied.
Here the books are finite and their terminal masses f are essential.

For prescribed supports D_1,...,D_m, let E(D_1,...,D_m) be the
event that the book admits a weakly increasing chain of stable
cuts with those successive endpoints. Choose the least D_1-cut,
then the least D_2-cut above it, and so on. If any such chain
exists this greedy procedure succeeds: by induction its current
cut is below the corresponding cut in every competing chain.
Repeated first-cut splitting therefore gives

    P_C(E(D_1,...,D_m))
       =F(C,D_1) product_(j=1)^(m-1)F(D_j,D_(j+1)).    (J76)

Consecutive equal supports allow the same cut twice, consistently
with F(D,D)=1. This formula does not count strict return chains.

Simultaneous visits also have a canonical stopping cut. For a
nonempty finite collection T of supports, suppose all are visited, let
d_D be the least D-cut for each D in T, and take their join j in
the stable-cut lattice. Define J_T(C,A) as the sum of weights of
prefixes C->A for which j is the full cut. Such a prefix remains
the selected prefix under every extension: each d_D remains the
least D-cut by (J70), and their join remains j because the cut
interval below j is unchanged. Hence these prefix cylinders are
disjoint, and

    sum_A J_T(C,A)=P_C(every D in T is visited)<=1,
    sum_(h:C->B visiting all T)W_p(h)
       =sum_A J_T(C,A)R(A,B).                         (J77)

The join is not necessarily the coordinate maximum and need not
end at a support in T. These are stopping-cut decompositions;
they do not put all stable cuts in one chronological order.

**Two-support elimination; proved, independent AI proof audit
passed.** Fix distinct D,E, let T={D,E}, and index the following
matrices by O=2^V\T. Let S be the Schur complement eliminating T
from R, and let U(C,B) sum books avoiding both D and E. Define
X(C,A) to sum prefixes whose least D- and E-cuts have full join,
with A outside T. This is the restriction of J_T in (J77), so

    sum_A X(C,A)<=1-f_C,
    S=U-XS=(I+X)^(-1)U
      =sum_(m>=0)(-1)^m X^m U.                       (J78)

All these numerical series converge absolutely. The excluded
identity book has mass f_C>0, giving ||X||_infinity<1 on the finite
index set. The suffix normalization here sums over ALL terminal
supports, including D and E; X itself is restricted to O.

Here is a proof that also identifies every coefficient. Write
a=F(D,E), b=F(E,D), so ab<1. The inverse of the normalized T
block [[1,a],[b,1]] expands into alternating first-cut chains.
For a fixed book h:C->B, let L_D(h),L_E(h) be the lengths of the
greedy alternating chains, beginning respectively at the least
D- or E-cut; a missing starting cut gives length zero. Every
switch strictly increases the cut. Least-cut splitting shows
that the Schur coefficient of this book is

    s(h)=1-1_(L_D(h) odd)-1_(L_E(h) odd).              (J79)

The matrix expansion is absolutely convergent since ab<1.
The greedy chains are longest for their respective starting
types; deleting the first cut proves |L_D-L_E|<=1. Thus s is
zero for unequal lengths and (-1)^m when both lengths equal m.

If both least cuts d,e exist, let j=d vee e. If j has endpoint
D, it is the least D-cut above e. The E-cuts above d are exactly
those above j, so the E-starting chain is one longer; the case
with endpoint E is symmetric. These books have coefficient zero.
If j has endpoint A outside T, split h=xy there. Every E-cut
above d lies above e and hence above j, and conversely; the same
holds with D,E exchanged. The upper cut interval is the suffix
lattice, giving

    L_D(h)=1+L_E(y), L_E(h)=1+L_D(y), s(h)=-s(y).

The prefix x belongs to X and is uniquely recoverable from h.
Its first cuts and join survive every extension by (J77).
Books visiting neither support contribute U; books visiting only
one contribute zero. This proves S=U-XS without multiplicities.
Iteration gives the stronger disjoint-class identity

    (X^m U)(C,B)
       =sum_(h:C->B with L_D(h)=L_E(h)=m)W_p(h),
    sum_(m>=0)(X^m U)(C,B)<=R(C,B).                   (J80)

The four-site diamond (J55), with D={1}, E={2}, has both lengths
one and contributes -W_p(h) to S-U. In the one-particle M_3 book
(J64), take two atom endpoints for D,E. The whole book z is an
X(C,C) prefix; z^m has coefficient (-1)^m by canonical splitting.
Thus either sign occurs before averaging books. No sign of the
numerical S or S-U, and no general inverse-M-matrix theorem,
follows from (J78).

Dagger supplies a compatible last-cut decomposition. Write Y_p
for the weights of suffixes whose greatest D- and E-cuts have
meet zero, with both outer endpoints in O. For q=1-p eliminate
the complemented pair on the q side, and write Phi_p for weighted
complement-transpose, as in (J58), restricted to these index sets.
Order reversal gives

    U_p=Phi_p(U_q), S_p=Phi_p(S_q), Y_p=Phi_p(X_q),
    (I+X_p)S_p=S_p(I+Y_p)=U_p,
    X_p U_p=U_p Y_p, X_p^m U_p=U_p Y_p^m.             (J81)

The first line follows from avoidance, Schur complementation,
and exchanging least joins with greatest meets. Applying Phi to
the q-side identity proves the right factorization; subtract S
and multiply by I+X to obtain the positive intertwining. Each
layer in (J80) therefore has two canonical factorizations: m
first-crossing prefixes followed by an avoiding suffix, or an
avoiding prefix followed by m last-crossing suffixes. There is
no record reordering. The weighted column bound

    sum_C mu_p(C)Y_p(C,B)
       <=[1-f_q(B^c)]mu_p(B)

gives absolute convergence of the right-hand alternating series.
The matrix Y need not be row-substochastic or equal X. Since R
and its eliminated principal block are invertible, S is invertible;
so is U=(I+X)S, and Y=U^(-1)XU. In particular complementary
first-crossing kernels are isospectral and have the nonnegative
duality matrix U_p diag(mu_p)^(-1)J, where J complements supports.
This does not assert positive definiteness or a target sign.


### The normalized remainder is a greatest cut in iid stacks

**Proved; independent AI proof audit passed.** Use infinite iid
paired stacks, with independent Bernoulli(p_i) bits and K-routing
arrows, and condition their first bits at C to be one. Denote this
product law by Q_C^stack, to distinguish it from the formal-cut
matrix Q. For a remainder a:C->D of length vector m, the event E_a
fixes its emitted records and requires beta_i(m_i+1)=1 on D.
It implies every conditioned initial one: at an unused site in C,
balance forces membership in D. Thus

    Q_C^stack(E_a)=p_D wt_p(a)/p_C=W_p(a).             (J82)

After conditioning on E_a and shifting past only the emitted
records, the tails have exactly law Q_D^stack. The next arrows
are still fresh, even where their paired bits are terminal ones.
The events E_a are precisely the finite stable cuts of the field.

These cuts have a greatest element almost surely. First, on a
probability-one arrow event, every edge of positive probability
occurs infinitely often in its source stack. Every legal sequence
with finitely many particles then has finitely many emissions:
otherwise an infinitely emitting site would send infinitely many
particles along every edge of a path to the sink, contradicting
the finite exit budget. This holds after every finite shift and
for every bit field. For two stable cuts m,n, their coordinate
maximum u satisfies F_C(u)>=u for the buffer map (J40). Its
residual input eta_u=1_C+I(u)-u is nonnegative. Legally stabilizing
this residual in the shifted stacks produces the least fixed
point above u, by least action, hence their join. If eta_u(i)=0,
F_C(u)>=u ensures the old bit at u_i is zero, so the shifted
zero-input bit convention agrees. This verifies the shift even
at empty sites and proves directedness without assuming an upper
cut in advance.

Second, Tonelli and (J82) give expected cut count sum_D R(C,D),
which is finite. Hence there are finitely many cuts almost surely;
they contain zero and are directed, so have a unique greatest cut.
Let H_C be its emitted book, and put g_D=Q_D^stack(no nonzero cut).
Given E_a, a is greatest precisely when its shifted tails have
no nonzero cut: directedness rules out an incomparable exception.
Consequently P(H_C=a)=W_p(a)g_D. Summation gives Rg=1, whence
g=f by invertibility. We have proved

    law(H_C)=P_C from (J30),
    f_C=Q_C^stack(zero is the only stable cut),
    E_(Q_C^stack)(number of cuts)=sum_D R(C,D).        (J83)

This supplies one iid field behind all the prefix probabilities
and stopping-cut rules. Its selected cut is the greatest; ordinary
legal stabilization selects the least, which here is zero because
the initial bits on C were conditioned to one.

For a general finite all-active eta, use the unconditioned iid field.
The same argument applies: the least legal cut exists, any two cuts
have a finite join, and their expected number is sum_B Q_eta(B)<infinity.
Writing H_eta for the greatest emitted book, a formal eta-to-B book
a has probability

    P(H_eta=a)=p_B wt_p(a)f_B,
    law(endpoint H_eta)=Q_eta diag(f)=P_eta L.         (J83a)

After any specified stable prefix its tails have Q_B^stack law,
so the same no-further-cut factor f_B applies. The statement holds
simultaneously for all finite integer inputs on a common probability-one
event. It does not equate the greatest endpoint with the physical one.

Intersections of prefix events now have a finite product formula.
For a:C->D_a and b:C->D_b, reject inconsistent emitted records
and any terminal-one requirement meeting an emitted zero. Otherwise
let u be their merged raw word, of coordinatewise maximum length,
and let S contain sites with a remaining terminal-one requirement
at position |u_i|+1. All other terminal constraints are already
checked inside u. The raw word need not be a stable cut. If
MCE(a,b) denotes the globally minimal common right extensions, then

    sum_(h in MCE(a,b))W_p(h)
       =Q_C^stack(E_a intersect E_b)
       =wt_p(u)p_S/p_C.                               (J84)

Conflicting constraints give zero. Every common extension contains
the join of a,b in its cut lattice; this is globally minimal.
Two distinct minimal extensions cannot share a common extension,
since the join inside it would lie below both. Their cylinders
are therefore disjoint and cover the intersection. Alternatively,
given the merged constraints, legally stabilize eta_u in the tails,
whose only forced bits are their first ones on S. This generates
the minimal completion. Here S only contains sites with eta_u>=1.
Fresh arrows and the termination argument above ensure completion.
Thus (J84) sums possibly infinitely many completions exactly.
Compatible cylinder events are positively correlated: their
intersection probability is the product of their probabilities
divided by the probability of their shared constraints, under
the already conditioned product law.

The possibility of infinitely many completions is real. In (J41)
take the two atoms using respectively (v_1,v_2,e_12) and
(v_2,v_3,e_23). Keep all other rows fixed, and replace the e_31
row by N copies of (1,e_31), then (1,v_3),(0,v_1), for N>=0.
These are distinct common completions h_N:C->C. A common divisor
must use both records at the other five sites. If it uses j<=N
records at e_31, that site's residual count is two; if j=N+1,
the residual is one but the next bit is zero. Only the full cut
is stable. Hence every h_N is minimal. They have positive weight
in one fixed accessible kernel with a positive e_31 self-loop.
This rules out a finite completion list even for a fixed pair of
finite prefixes. Formula (J84) does not require such a list.


### An arrow envelope and adjoint input maps

**Proved; independent AI proof audit passed.** One finite random
arrow book bounds every stable cut for every binary input and
every bit field. First drain one particle per site completely,
without sleeping, using the iid arrows. Let d be its emission
vector. In the unused stacks pop all cycles of current head arrows
until the heads form an acyclic graph, and let c be the popped
vector. Both operations terminate on the common arrow event used
above. Put Delta=d+c. Then

    Delta=1_V+I(Delta),
    m<=1_V+I(m) implies m<=Delta.                     (J85)

Here is the maximality argument. Balanced arrow prefixes are
closed under coordinate maximum: for their maximum u, I(u)>=u,
and summing forces equality and zero exits. Every nonempty balanced
prefix exposes a head cycle on its used vertices. Thus the completed
cycle-pop vector c is the greatest balanced prefix of the shifted
stacks. Least action bounds the drain d by any full-drain solution;
the remaining suffix is balanced, proving Delta is the greatest
solution of its displayed equation. More strongly, for any m as
in (J85), u=max(m,Delta) satisfies u<=1_V+I(u). It has at least
|V| sink arrows, since Delta does, so summing forces equality.
Then u-Delta is balanced in the acyclic-head tails, hence zero.
This proves (J85) without assuming fresh tails at a selected cut.

For C subset V, F_C(Delta)<=Delta and every queried bit index
1_C(i)+I_i(m), m<=Delta, is at most Delta_i. Descending iteration
of F_C from Delta therefore computes the greatest stable cut,
using only its first Delta_i bits, with at most |Delta| strict
decreases. Ascending iteration from zero gives the physical least
cut. This works simultaneously for all bits, including deterministic
ones, and requires no initial-one conditioning for the algorithm.

Writing d_K=det(I-K), the consumed envelope h has the exact law

    P(envelope=h)=d_K wt_ar(h),
    sum_(h: |h|=1_V+I_h(|h|))wt_ar(h)=1/d_K.          (J86)

For a fixed full-drain word h, it is the envelope iff its first
unused arrows are acyclic. Their probability is d_K by directed
cycle inclusion-exclusion. This proves (J86), including books
that are not themselves legal drains. The next-arrow forest is
independent of h, with its ordinary forest law; only the arrows
after that forest are fresh iid tails. Equivalently h factors
uniquely into its physical complete drain and a balanced suffix,
independently distributed under (J86). In particular the envelope
law is invariant under reversing every local arrow word.

**Input maps; proved under (J83a), independent derivation and main
proof check.** On one iid field, let L_eta be the stable-cut lattice
for finite all-active eta, and take eta<=eta'. Since F_eta<=F_eta',
define

    alpha(m)=least eta'-cut above m=lim F_eta'^k(m),
    beta(n)=greatest eta-cut below n=lim F_eta^k(n).    (J87)

The first iteration increases and is bounded by legal stabilization
of the residual D_eta(m)+eta'-eta in the shifted stacks. The second
decreases from the finite vector n. The shifted first ones on
D_eta(m) represent existing sleepers; they are retained, not redrawn.
These extremal definitions give the adjunction

    alpha(m)<=n iff m<=n iff m<=beta(n).               (J88)

Thus alpha preserves joins, beta preserves meets, beta alpha is
a closure, and alpha beta is an interior operator. They carry
the lower least cut to the upper least cut and the upper greatest
cut to the lower greatest cut, respectively. For eta<=eta'<=eta'',
the maps compose exactly along the chain. Indeed, between any
eta-cut m and eta''-cut n with m<=n, upward F_eta' iteration
stays between them and supplies the needed intermediate cut.
This proves alpha composition; uniqueness of right adjoints gives
beta composition. In particular the formal backward deletions commute.

These maps can lose information in both directions. On one site,
prescribe arrows (self,sink,self,sink,sink) and the first four bits
all one. For inputs 1 and 2 the complete cut sets are {0,1} and
{2,3}, respectively. Then

    alpha(0)=alpha(1)=2, beta(2)=beta(3)=1.            (J89)

Both maps fail injectivity and surjectivity. Beta of the larger
physical cut is a nonphysical cut for the smaller input. Therefore
beta is not an inverse physical addition, nor does the adjunction
imply endpoint occupation order. A formal choice of cut may also
bias its tails; the deterministic addition interpretation of alpha
does not provide fresh-stack conditional probabilities by itself.

The focused finite-field certificate uses a two-site full-drain
book and its local reversal, with acyclic next heads. All 128 bit/input
fields have the predicted least and greatest cuts; independent buffer
execution agrees with the least cut, and all 400 adjunction cut-pair
checks pass. This is a finite check supplementary to the arguments.

The arrow envelope cannot be held fixed in a target comparison.
The old excursion obstruction has a direct realization here. Take
V={0,1,2} and the arrow words

    h_0=(2,s,s,s), h_1=(0), h_2=(0,0).

They have unit divergence at every vertex. Under
K(0,2)=K(0,s)=1/2, K(1,0)=K(2,0)=1, the next root arrow being
s makes the next heads acyclic. This is exactly an envelope event
of probability 1/32, independent of all seven sleep bits.
For C={0}, C'={0,1}, direct buffer stabilization gives

    P(0 occupied | C,h)=p_0[1+(1-p_0)(1-p_2)],
    P(0 occupied | C',h)=p_0[1+p_1(1-p_0)(1-p_2)].     (J89a)

If the extra sleeps at 1, the lower run remains. Otherwise it
reaches 0; the final root bit is its second or third iid bit,
according to whether the excursion sleeps at 2, so its mean is
just p_0. The lower particle has an additional return-and-sleep
opportunity. Thus the conditional increment is strictly negative,
-(1-p_1)p_0(1-p_0)(1-p_2). All 128 bit fields at p=1/2 give
80 lower occupations and 72 upper occupations, a deficit 1/16.
This is a checked version of the SX fixed-path obstruction in the
new envelope space, not a counterexample to the original conjecture.
Any sufficient comparison here must average different arrow books.


### Monotone cut couplings and stationary densities

**Proved; independent AI proof audit passed.** Fix the source C,
use common arrows and uniforms, and set beta_i^p(j)=1_(U_i(j)<=p_i),
except that the first bits on C are pinned to one for every p.
The arrow envelope bounds all cuts simultaneously for every bit
field, so there is no exceptional-set issue across parameters.
If p<=p', then F_p>=F_p' and descending iteration from the common
envelope gives

    m_C^+(p)>=m_C^+(p').                              (J90)

Every directed-arrow current, including each labeled exit count,
therefore decreases; the greatest endpoint's total particle count
increases. If only p_v increases, every old stable cut occupied
at v remains stable: its queried v-bit was already one. If the
old greatest cut is occupied there, (J90) forces the new greatest
cut to be exactly the old one. Its endpoint and selected arrows
remain unchanged; bits inside that prefix may still change.
Thus L_p(C,B) increases in p_v for each exact B containing v.
The unconditioned greatest cut for a fixed all-active eta has the
same property, for its law P_eta L. These do not assert occupied-set
inclusion under general parameter changes.

For comparison, the physical least cut decreases when p increases.
If only p_v changes and the new least endpoint is vacant at v,
its queried bit is zero and was also zero before the change. It
is therefore an old fixed point, forcing equality of the two least
cuts. This proves physical own-p_v occupation monotonicity, and
P_eta(B) decreases in p_v for each exact B avoiding v. No particle
is added in this parameter comparison.

Conditioning on the arrow envelope makes these finite Bernoulli
functions. For greatest occupation at v, or an exact greatest
endpoint B containing v, each unpinned bit of row v has nonnegative
influence. Differentiating gives the convergent pivotal formula

    partial_(p_v) E X
       =sum_(j unpinned) P(bit (v,j) is pivotal)>=0.

For the conditioned law omit j=1 when v is in C. Finite expected
envelope length justifies differentiation: its physical-drain part
has finite Green-function mean, and its balanced part has finite
mean by the exponential tail in (J26). For general finite eta the
same drain/cycle construction with a dominating input gives a finite
mean bound. The focused certificate checks 256 individual-bit
comparisons, including both freezing implications.

The no-remainder probabilities have two orders:

    p<=p' implies f_C(p)<=f_C(p'),
    C subset C' implies f_C'(p)<=f_C(p).              (J91)

The first follows from f_C=P(m_C^+=0) and (J90). For the second,
pin precisely the first bits on each respective source, using common
underlying uniforms. Add i to C and put n=I_i(m). At n=0 both
map coordinates are zero. At n>=1 the old and new coordinates
are n-b_i(n) and n+1-b_i(n+1), whose difference is nonnegative.
The other coordinates agree. Thus F_C' >= F_C, and their greatest
cuts are ordered in the same direction; their zero-cut events are
nested in reverse. These are conditioned remainder fields, not
the two physical all-active laws.

The first inequality extends to the closed cube with first bits
pinned, including faces where literal conditioning has probability
zero. At p=1, R is diagonal with R(C,C)=1/d_C, so f_C(1)=d_C.
There is also a positive lower bound uniform in C and p. Choose
an acyclic routing forest toward the sink, with edge probabilities
alpha_i>0, and let t_i count vertices in its subtree at i. Force
the first t_i+1 arrows at i to follow that forest edge. This event
has probability gamma=product_i alpha_i^(t_i+1)>0; its envelope
has lengths t. Every remainder cut is zero by induction from
forest leaves, since used arrows follow the forest and initial
first bits are pinned. Consequently

    0<gamma<=f_C(0)<=f_C(p)<=d_C.                     (J91a)

In particular the crossing contraction in (J78) can be bounded
uniformly over the closed sleep cube for a fixed K.

Here is a stationary-law consequence. Write pi_t for the full-input
stationary law at sleep vector t, and mu_t(A)=t_A(1-t)_(A^c).
The identity (J30) and (J91) give

    pi_t(A)/mu_t(A)=f_(A^c)(1-t)/d,
    A subset B implies pi_t(A)/mu_t(A)<=pi_t(B)/mu_t(B).
                                                               (J92)

Thus its density relative to product Bernoulli(t) is increasing.
Equivalently, every full one-site conditional obeys

    pi_t(v occupied | occupied set outside v is A)>=t_v,
       A subset V\{v}.

This follows by comparing the two density values at A and A union
{v}; averaging gives the bound under any conditioning on other
sites. Sequential revelation with common uniforms therefore couples
pi_t above product Bernoulli(t). The density also decreases in
every sleep coordinate t_i, by the first inequality of (J91).
These statements do not assert the two-pi FKG lattice inequality
or positive association. An exact heterogeneous three-cycle check
passes all twelve density cover comparisons. No new-literature
claim or add-one occupation conclusion is made.

The changing reference product law matters: pi_t itself need not
increase in an off-target sleep parameter. For K(0,1)=1 and
K(1,0)=K(1,s)=1/2, set p_0=a, p_1=t. Direct two-site inverse
algebra gives

    pi_(a,t)(0 occupied)
       =a[2+(1+a)t+(1-a)t^2]/[(1+a+(1-a)t)(1+at)],
    derivative in t
       =a(1-a)^2(t^2+2t-1)/[(1+a+(1-a)t)^2(1+at)^2]. (J92a)

For every interior a its minimum is at t=sqrt(2)-1. At a=2/3,
increasing t from 1/10 to 1/2 lowers this probability from 217/272
to 35/44, a change -7/2992. To verify the formula, put
Z=(1+a+(1-a)t)(1+at); the full-input probabilities at
empty,{0},{1},{0,1} are respectively
(1-a)(1-t)/Z, a(1-t)[2+(1-a)t]/Z, 2(1-a)t/Z,
and 2at/(1+at). They solve pi E=delta_V for the two-site
complete-block inverse. This symbolic calculation and the three
rational evaluations were checked separately. It concerns parameter
variation, not the add-one conjecture.

The conditional lower bound in (J92) has an exact strictness criterion:
for every interior t and every outside configuration A,

    pi_t(i occupied | outside i equals A)>t_i
       iff i lies on a positive-probability directed cycle.     (J92b)

Here a self-loop counts. First, the full stationary law factors over
SCCs: process components topologically, and use full-input saturation
to see that the next component's law is unchanged by its random
incoming active batch. Determinants and Bernoulli weights factor,
so f_C also factors over SCCs. An acyclic singleton has local f=1
for both sources, giving equality and independence at that coordinate.

For strictness, fix i outside C on a simple directed cycle Gamma.
Choose an acyclic positive routing forest F to the sink. Let the
first-arrow graph H follow Gamma on its vertices and F elsewhere;
its only cycle is Gamma. If a counts its incoming arrows and T is
the internal routing matrix of F, put t=(I-T^T)^(-1)a and delta=1+t.
Prescribe first arrows H, then arrows F through delta+1. The
length-delta word has unit divergence and forest next heads, so
bounds all cuts. Make all relevant bits at C union {i} one, and
first bits on Gamma outside that set zero. This has positive
probability, including the old unpinned first bit at i.

For the old source C, any stable endpoint contains C; conservation
forces it to equal C, with no exits. No arrow can arrive at i,
whose prescribed one bit would add an occupant. Thus any nonzero
cut would be balanced and expose a head cycle avoiding i, impossible.
For the new source C union {i}, one emission around Gamma is a
valid nonzero balanced cut. The nested zero events therefore give
f_C>f_(C union {i}); substitution in the two-state conditional
ratio proves (J92b). This finite-event proof received an independent
derivation and main check. Strictness can disappear on sleep faces.


### Envelope inventories remove the last-exit restriction

**Reduction proved; independent AI proof audit passed. The proposed
coefficient signs remain open.** The envelope construction works
with any finite integer drainage input rho>=0: drain rho, then pop
all remaining head cycles. Its greatest flow solution satisfies
Delta=rho+I(Delta), and the same maximum-and-exit argument proves
m<=rho+I(m) implies m<=Delta. Its ordered arrow book still has
law d_K wt_ar(h), on all books with divergence rho. For eta<=rho,
all queried bit indices are at most Delta_i. A zero-length site
has rho_i=eta_i=0 and queries only the convention beta_i(0)=0.

Fix C'=C union {x}, x outside C, and rho>=1_(C'). For an edge
inventory J with divergence rho, write m_i=sum_j J(i,j), including
sink destinations. Every local ordering of these arrow multisets
is allowed. Give its m_i bits exactly k_i ones, and sum the
physical target increment over this pool, calling it c_env^rho(J,k).
The pool size and each atom's probability are

    product_i [m_i!/product_j J(i,j)!] binom(m_i,k_i),
    d_K K^J product_i p_i^k_i(1-p_i)^(m_i-k_i).        (J93)

The next-arrow forest has been integrated out; it restricts none
of these orders. Both physical stabilizations stay within the same
book, so their target difference is determined by this finite
pool. Nonnegative c_env for every such inventory would imply the
binary conjecture, hence the full all-active conjecture by the saved
IDLA reduction. The conditioned-first-one remainder law is not
used for these physical observables.

There is a precise relation to the earlier one-family pool in
[DI Sections 1--2 and 6](drain-inventory.md). Use absolute edge
variables Y, put W=(Y_ij) on internal edges, b_i=Y_is,
h=(I-W)^(-1)b, H_r=product_i h_i^r_i, and d(Y)=det(I-W).
Let C_env^rho(Y,p) and C_1(Y,p) be the unnormalized signed series
over envelope and DI atoms, respectively. Then

    C_env^rho(Y,p)=[H_(rho-1_(C'))(Y)/d(Y)] C_1(Y,p).
                                                               (J94)

To prove it, first use positive weights with rho(W)<1 and h>0.
The harmonic tilt K^Y_ij=Y_ij h_j/h_i, K^Y_is=Y_is/h_i is stochastic,
with determinant d(Y). Every divergence-rho book has weight
Y^J/H_rho under this tilt. Its normalized envelope expectation
is the ordinary physical increment DeltaQ^Y. Hence
C_env^rho=H_rho DeltaQ^Y/d(Y), whereas DI gives
C_1=H_(1_(C')) DeltaQ^Y. This proves (J94). The series converge
absolutely near zero, giving the formal identity without inverting
the zero-constant-term h_i formally.

The multiplier P_r=H_r/d has nonnegative integer coefficients
a_r(L), counting spectator walk families and balanced arrow books.
Padding adds m_i(L) unused coins, so the exact coefficient relation is

    c_env^rho(J,k)=sum_(J0+L=J) a_r(L)
       *sum_k0 c_1(J0,k0) product_i binom(m_i(L),k_i-k0_i),
       r=rho-1_(C').                                  (J94a)

This follows by inserting (p_i+(1-p_i))^m_i(L) and using the fixed
degree Bernstein basis. No additional edge-allocation factorial
occurs. The relation concerns signed marginal sums; the path-coin
and site-coin couplings need not agree atom by atom. In particular

    c_1>=0 implies c_env^(1_(C'))>=0 implies c_env^(1_V)>=0.

The first padding consists only of balanced loops; the second
also adds spectator walks. These are weaker sufficient conditions
in this implication sense; no strictness or converse is asserted.

The inventory containing (J89a) has a complete small cancellation.
Move its sole root arrow to 2 through positions t=1,2,3,4 among
the three sink arrows. Put X_1,...,X_4 for root bits, Y for the
site-1 bit, and Z,Z_2 for site-2 bits. If Y=1 all target-0 increments
are zero. For Y=0 and Z=0, the increments at t=1,2 are respectively

    X_3-(X_1 OR X_2), (X_2 OR X_3)-X_1.

Pair t=1 with t=2 by (X_1,X_2,X_3)->(X_3,X_1,X_2), using the
inverse rotation in the reverse direction. In all other Y=0
cases the increment is X_2-X_1; keep t and swap X_1,X_2. Fix
Y=1 states. This is an involution preserving every site's coin
count and the arrow inventory, and reversing the target-0 sign.
Thus that coefficient is zero in every pool. The other two target
increments are pointwise nonnegative: target 1 gains Y, while
target 2 gains (1-Y)X_1 Z at t=1, (1-Y)(1-X_2)Z at t=2, and
zero otherwise. All 512 decorated books and 30 coin-count sectors
passed an exact check, supplementing the explicit involution.
This proves one inventory, not general (J93) positivity.

A second fixed inventory also has nonnegative coefficients. For
C={0}, C'={0,1}, take row multisets (s,s,1), (s,2), (0,0), with
lengths (3,2,2). All six arrow orders and 128 bit fields are included.
If P_v(x,y,z)=sum_k c_env(J,k)x^k_0 y^k_1 z^k_2, exact enumeration gives

    P_0=x(1+x)y(1+y)(1+z),
    P_1=(1+x)^2 y(1+z)[(1+x)+2(2+3x)(1+y)(1+z)],
    P_2=(1+x)^2(1+y)z(1+z)(3+3x+y).                 (J94b)

Every coefficient is nonnegative. At equal sleep one half the
signed totals are (8,336,112) among 768 decorated books. Individual
losses still occur, including 24 at target 0 and 32 at target 1.
The focused verifier checks all 36 sector populations and these
three factorizations exactly. This is evidence for a second
inventory, not a general sign theorem.

There is also a finite random-set interpretation of (J92). Fix
underlying arrows and bits at 1-t, and for each C pin its first
bits to one as in (J91). Let I be the family of C whose greatest
cut is zero. It is downward closed, and may be empty. Its
complemented family U={A:A^c in I} is upward closed. The envelope
decides the whole family using finitely many records. Independently
sample A from product Bernoulli(t). Then

    P(A in U)=d_K,
    law(A | A in U)=pi_t.                             (J95)

Indeed the joint acceptance mass at A is mu_t(A)f_(A^c)(1-t),
which equals d_K pi_t(A). Conditional on the field, acceptance
requires A to meet every minimal forbidden source set of I.
Thus the stationary law is a mixture of product laws conditioned
on these upward events, with the fields tilted by their acceptance
probabilities. These families need not be matroids, as shown below;
no FKG structure is asserted. Equation (J95) is an exact rejection representation,
not an efficiency claim; d_K can be small.

The source ideal factors pathwise over SCCs, I=product_S I_S,
when outgoing arrows are killed in each local problem. A global
nonzero cut has an earliest nonzero component, giving a local cut
with no external arrivals. Conversely embed a local cut with zero
emissions elsewhere; it is postfixed for the global buffer map.
Upward iteration below the common envelope completes it, changing
only downstream components. Thus local and global no-cut events
factor before averaging, with no fresh-tail assumption at the chosen cut.

The finite envelope permits an exact comparison with burning tests,
developed in (J99)--(J101) below. It requires a cyclic processor and
careful treatment of its boundary state; an empty reset state is
not automatically a recurrent local phase.

The no-cut source ideal need not satisfy matroid exchange. There
is an explicit positive-probability three-site field:

| Site | Consumed arrows | Bits | Next arrow |
|---|---|---|---|
| 0 | (1,s,1) | 101 | s |
| 1 | (2,s,2,s,2) | 10101 | s |
| 2 | (0,1,0,1) | 1010 | 0 |

Its lengths Delta=(3,5,4) have unit divergence and acyclic next
heads, so every binary-input cut is bounded there. All first bits
are already one. Write C's indicators c_i, and any cut as
m=(2a,2b,2c); alternating bits force even coordinates. Its endpoint is

    (c_0+c-2a, c_1+a+c-2b, c_2+b-2c),
    0<=a<=1, 0<=b,c<=2.

For a nonzero cut, c=0 forces a=b=0, while c=2 makes the last
endpoint coordinate negative. Thus c=1. The last coordinate
requires b>=1, and the second rules out b=2. Hence b=1, c_2=1,
a=c_0, and c_0+c_1>=1. Conversely these conditions give a cut.
Therefore

    I={empty,{0},{1},{2},{0,1}}.                       (J96)

The independent sets {2} and {0,1} violate augmentation; the
minimal forbidden sets {0,2},{1,2} violate circuit elimination.
Choose each row's two displayed destination types with probability
one half. The twelve arrows and three next heads have probability
2^(-15); the twelve specified bits have positive probability for
all interior sleep vectors. Later records cannot alter the ideal.
The full hand classification and all 144 even-cut/input candidates
were checked independently; the focused verifier retains the latter.

There is an exact energy interpretation in a nearby deterministic
model. At every site repeat a fixed ordered pair of outgoing arrows,
and repeat bits 1,0. Assume the resulting directed graph reaches the
sink. If A_ij counts its two internal outgoing occurrences, put
M=2I-A^T. Then M is a nonsingular M-matrix, with M^(-1)>=0.
For any integer input eta>=0, stable cuts have m=2z and endpoint
eta-Mz in {0,1}^V. In fact the greatest cut is

    m^+=2z*, z*=max {z>=0 integer : Mz<=eta}.          (J97)

The feasible set is finite because z<=M^(-1)eta, and closed under
coordinate maximum because M has nonpositive off-diagonal entries.
Its greatest vector therefore exists. If B=eta-Mz* has B_i>=M_ii,
adding e_i to z* remains feasible, impossible. Thus B_i<M_ii<=2,
so 2z* is a stable cut. No further nonnegative integer script s
can leave B-Ms nonnegative. This is z-superstability, which tests
arbitrary firing multiplicities, not only subsets.

It is the unique such representative in the whole nonnegative
integer class eta+M Z^V. Indeed, if g=eta-Mw>=0 for an arbitrary
integer w, the same coordinate-maximum argument makes max(w,z*)
feasible and nonnegative. Hence w<=z*, so g=B+Ms for s>=0.
Unless s=0, this g admits the nonzero script taking it to B.
Moreover

    ||M^(-1)g||_2^2-||M^(-1)B||_2^2
       =||s||_2^2+2 s^T M^(-1)B>0 when s!=0.          (J98)

Thus the greatest-cut endpoint uniquely minimizes this energy over
the NONNEGATIVE class. This agrees with the general theorem in
[Guzman--Klivans, *Chip-firing and energy minimization on M-matrices*,
Definition 4.3 and Theorem 4.6, with energy (3.1)](https://arxiv.org/pdf/1403.1635),
whose definitions, conventions, and proof were read. The direct
bridge above passed an independent AI proof audit. Self-loops give
M_ii=1 and force the greatest endpoint vacant there; the net matrix
diagonal, rather than the two emitted records, determines that bound.

For binary C all first bits already equal one, so its no-cut
condition is exactly z-superstability of 1_C. The source ideal in
(J96) is the example with arrow pairs (1,s), (2,s), (0,1).
Ordinary chi-superstability or G-parking-function identifications
require additional hypotheses in the directed setting; none is
substituted here. The general periodic field is a deterministic
comparison model. The finite cylinder (J96) establishes that its
nonmatroid behavior also occurs with positive iid probability.

For a general unit-divergence arrow inventory J, write m_i>=1 for
the row lengths, D_m=diag(m), and use J also for its internal part.
Thus m-J^T 1=1. Periodize any ordering of its arrows and any bit
word beta_i(1),...,beta_i(m_i), retaining all m_i labelled phases.
Here beta_i(0)=beta_i(m_i), in contrast to the empty reset convention.
At local phase a modulo m_i, the resident is beta_i(a) and the
emission head is e=a-beta_i(a). One arrival advances the phase
and emits 1+beta_i(a)-beta_i(a+1) consecutive arrows. Over m_i
arrivals every arrow is emitted once. This is a finite unary
abelian network with total kernel D_m Z^V and

    P=J^T D_m^(-1), M=(I-P)D_m=D_m-J^T,
    M1=1, Pm=m-1, rho(P)<1.                         (J99)

The last inequality follows in the maximum norm weighted by m.
Thus the network halts on every input. The burning input is 1_V:
a phase is globally recurrent iff adding 1_V returns it to itself.
There are exactly det M recurrent phase vectors, and the critical
group is Z^V/M Z^V. These applications use
[Bond--Levine, *Abelian networks III*, Theorems 2.9, 3.3, 3.11,
3.18, and 5.5](https://arxiv.org/pdf/1409.0170); the hypotheses,
kernel conventions, counting theorem, and burning proof were read.
Local irreducibility holds because each labelled phase map is a cycle.

To represent the original reset processor, use phase zero only if
beta_i(m_i)=0. Otherwise add a transient state with empty resident
and head zero, whose first arrival enters phase 1 and emits
1-beta_i(1) arrows. For every positive r, reset input r emits
r-beta_i(r) arrows beginning with arrow 1. The transient state is
not cyclic phase zero and cannot be globally recurrent. Binary
all-active inputs still reproduce the original envelope calculation:
every queried positive bit index is at most m_i. Iid full-input
saturation is not being transferred to these periodic networks.

A cyclic state with resident support C and all emission heads zero
exists exactly when the first bit is one on C and the last bit is
zero outside C. Its phase vector is 1_C modulo m. Adding 1_V there
equals the reset calculation with formal input 1_V+1_C; the full
emission vector m is a stable upper cut with endpoint C. For a
cut u<=m with endpoint B, reverse the local arrow words and
complement-reverse their bits. No arrow direction is reversed.
Then n=m-u is a cut from source C^c with endpoint B^c, since

    Ibar(n)=m-1-I(u),
    1_(C^c)+Ibar(n)-n=1-B.                           (J100)

The bit indices are complementary: the reverse index is
m_i+1-(u_i+B_i). An index zero on the reverse side corresponds
only to u_i=m_i,B_i=1,i in C; its forward bit m_i+1 repeats the
required first one. All other indices lie in the original words.
This proves the cut anti-isomorphism including the boundaries.
The burning operation returns its phase iff its least emission
cut is m: a return requires each positive arrival count to be a
multiple of m_i, and the upper bound forces exactly m_i. Hence
recurrence is equivalent to absence of a nonzero reverse cut in
[0,m]. Appending acyclic unused heads makes the reverse word an
envelope, so this is also the full no-cut condition. A periodic
continuation of that reverse word need not have this envelope bound.

This comparison applies to every pinned-source test in (J95).
If S is the pinned source, indices queried outside S are at most
m_i-1. Set those unused last bits to zero without changing any
cut. Reverse-complement now produces a boundary-compatible cyclic
state with C=S^c; its recurrence is exactly the original no-cut
event. This is a deterministic equivalence, but the padding changes
coin counts and does not preserve an arbitrary fixed-(J,k) pool.

There is nevertheless an exact count in the boundary-compatible
ensemble. Fix bit counts k_i and let W contain all local arrow
orders and bit orders. Put

    N=|W|=product_i [m_i!/product_dest J_(i,dest)!]
                        * binom(m_i,k_i).

Then

    #(w,C compatible)=N,
    #(w,C compatible and recurrent)
        =N det(M)/product_i m_i
        =N det(I-D_m^(-1)J).                         (J101)

To prove this, take w uniform in W and a uniform labelled phase a.
Set C_i=beta_i(a_i), e_i=a_i-C_i modulo m_i, and rotate both local
words by e_i. This gives a compatible pair with new phase C_i.
Every compatible pair has exactly product_i m_i preimages: choose
the offsets e_i, undo the rotations, and set a_i=e_i+C_i. This
remains true for words of smaller period because the phases are
labelled. Rotation is a network isomorphism, so it preserves
recurrence. Each fixed w has det M recurrent phases among
product_i m_i choices, proving both counts. The empirical routing
matrix D_m^(-1)J need not equal the original kernel K. The marginal
on words is weighted by the number of compatible supports; it is
not uniform. The construction and boundary checks passed an
independent AI audit. For the inventory in (J94b), the focused
verifier checks 9,216 word/phase burning tests (ten recurrent
phases out of twelve for each word), all 768 compatible pairs
(640 recurrent), and 27,648 individual cut reversals. These
results give a finite burning model
and a normalizer, not the two-input target sign.

There is an exact addition action on this finite ensemble. Starting
from a recurrent compatible pair (w,C), stabilize an input b, let
ell be its emission vector, and let D be its final resident support.
Rotate both local words by ell_i, so their new head is zero. Write
this operation as A_b. Since the final phase is ell_i+1_D(i), the
rotated word is compatible with D. It preserves both inventories.
On recurrent pairs,

    A_b A_c=A_(b+c), A_1=id,
    A_(e_x)^(-1)=A_(1-e_x).                          (J102)

For precision, use labelled pairs (w,a) and quotient by rotations
(w,a)->(R_e w,a-e), where e lies in product_i Z/m_i Z. This action
is free, even for words of smaller period, and each orbit has one
head-zero representative. Stabilization commutes with these network
isomorphisms. Its addition maps commute, preserve recurrence, and
the burning input is the identity there; hence they descend to
(J102). Uniform recurrent compatible pairs are therefore invariant
under every single-site addition. The retained-word process is
stationary; its resident projection need not be Markov. Cyclic word
orders are conserved, so uniqueness on the whole ensemble is not
asserted.

For fixed J, now give a compatible pair weight
N_arrow(J)^(-1) product_i p_i^k_i q_i^(m_i-k_i), where q=1-p.
This is a probability law with resident marginal mu_p: first-one
compatibility costs p_i and last-zero compatibility costs q_i.
Its recurrence probability is d_J=det(I-D_m^(-1)J), independently
of p, by (J101). Conditional on recurrence, the counts k_i remain
independent Binomial(m_i,p_i); given k, the pair is uniform among
the recurrent compatible pairs. Let pi_(J,p) be its resident law.

The original stationary law has the mixture representation

    pi_p(B)=sum_J D_K(J) pi_(J,p)(B),
    D_K(J)=N_arrow(J) K^J d_J.                       (J103)

Here D_K is exactly the arrow-inventory law of an ordinary no-sleep
drain from 1_V. Indeed a unit-divergence ordered book is such a
drain iff its last arrows form a sink forest. Necessity follows
by ordering last departure times. For sufficiency, an unconsumed
suffix after drainage would be balanced without exits; its last
arrows would have a directed cycle. Under uniform local orders,
the last destinations are independent with kernel D_m^(-1)J, and
their forest probability is d_J. Thus D_K sums to one.

To identify the resident law, first use the envelope inventory law
d_K N_arrow(J)K^J. Recurrence conditioning tilts it by d_J/d_K,
giving exactly D_K. Given resident B, reverse the arrows and
complement-reverse the bits. The reverse source is B^c, with first
ones there and last zeros outside. Those last zeros are irrelevant
to its cut test and can be regenerated as iid Bernoulli(q) bits.
The recurrence event is precisely no-cut, whose joint mass at B
is mu_p(B)f_(B^c)(q)=d_K pi_p(B), proving (J103). This does not
identify the cyclic word with the actual drain ordering; only
their inventory laws agree. It also does not replace retained
periodic words by fresh iid instructions. The addition action,
mixture weights, and normalizations passed an independent AI audit.

Stationary occupation also has an exact pivotal-current formula.
Keep the no-cut parameter p, and let C~mu_p independently of the
envelope and its free bits. For i fixed, let D~mu_(p,-i) be a
subset of V\{i}. Differentiating sum_C mu_p(C)f_C(p)=d_K gives

    J_i(p):=E_C partial_i f_C(p)
        =E_D [f_D(p)-f_(D+i)(p)] in [0,1],
    pi_(1-p)(i occupied)-(1-p_i)
        =[p_i(1-p_i)/d_K] J_i(p).                  (J104)

The currents on both sides have explicit nonnegative events.
For a fixed envelope of lengths m, use m_j-1 independent free
coins Z_(j,r) at each site and set beta_j(r+c_j)=Z_(j,r), with
c_j=1_C(j), forcing the first bit on C. The common no-cut map is

    F_(C,Z)(u)_j=0                        if I_j(u)=0,
                =I_j(u)+c_j-Z_(j,I_j(u))  otherwise. (J104a)

Thus adding source i increases this map exactly by e_i when its
incoming count is positive. Raising a free coin decreases the map.
If Phi_C is the no-cut indicator, partial_i f_C is the expected
number of free row-i coins pivotal from Phi_C=0 to 1. Their
absolute indices are 1+c_i,...,m_i-1+c_i; neither a forced first
bit nor an unused last bit contributes. The right side of (J104)
is the probability that Phi_D=1 and Phi_(D+i)=0 in this shifted
coupling. The envelope has an exponential length tail, and the
conditional derivative is bounded by m_i-1; differentiation and
the pivotal sums are justified, including one-sided limits at the
closed-cube faces. The probability formula for pi follows directly
by writing d_K=(1-p_i)E_D f_D+p_i E_D f_(D+i). This proof received
an independent AI audit. No bijection between marked coin-pivotal
and source-pivotal configurations has yet been proved.

The shifted coupling also yields the exact cross-face identity
f_D(p_i=0)=f_(D+i)(p_i=1), with other parameters unchanged:
the row-i map is I_i(u) in both cases. In particular its comparison
survives arbitrary changes to that entire row's unforced bits
when the source is inserted: n-beta_old(n)<=n+1-beta_new(n+1)
for n>0, and both maps vanish when n=0.

This entire-row comparison also works without replacement. Fix
(J,k), let theta_i=k_i/m_i, and let nu_(J,k) be the resident law
of uniform recurrent compatible pairs. The unconditioned compatible
pair has resident prior mu_theta, since its two local counts are
binom(m_i-1,k_i-1) and binom(m_i-1,k_i). After reversal, source
S=C^c has m_i-k_i one bits, with first one on S or last zero off S.
Let g_S be its no-cut probability under these conditional uniform
rows and uniform arrow orders. For a feasible source insertion,
share ALL arrows and the other bit rows; couple the changed bit
row arbitrarily. The last inequality above orders the two cut
maps for every such pair of rows, so g_(S+i)<=g_S. Consequently

    nu_(J,k)(C)/mu_theta(C)=g_(C^c)/d_J
        is increasing in C,
    nu_(J,k)(i occupied | exact other residents)>=k_i/m_i. (J105)

The density is defined on the product-prior face: k_i=0 or m_i
forces the corresponding resident bit. The conditional bound is
asserted only for positive-probability backgrounds. Zeros in g
cause no difficulty: for h_1>=h_0>=0 the conditional probability
theta_i h_1/[(1-theta_i)h_0+theta_i h_1] is at least theta_i.
Averaging gives the same bound under partial conditioning, hence
stochastic domination of the product prior by sequential uniforms.
This proof passed an independent AI audit. It compares source tests
within one fixed-inventory ensemble, not different k ensembles,
and asserts neither FKG nor physical add-one occupation order.

For the distinct question of stationary log-submodularity, exact
finite evidence is positive in two prescribed cases. With
Delta_ij(A)=pi(A)pi(A+i+j)-pi(A+i)pi(A+j), all six elementary
minors are negative on 0->1->2->0/s, with the last choice uniform
and sleeps (1/3,2/5,3/5). On the four-site kernel of (J55), with
all sleeps one half, twelve minors are negative and the twelve
involving the independent sink-pointing vertex are zero. The
focused verifier checks all thirty minors rationally. Since the
product Bernoulli factors are modular, the equivalent f inequality
is f_C f_(C+i+j)<=f_(C+i)f_(C+j), at complementary sleep parameters.
These are finite checks, not a general lattice inequality or a
consequence of the increasing density in (J92) or (J105).

The recurrent ensemble has a concrete forest correspondence. Fix
one periodic arrow/bit word with unit divergence. For a cyclic
phase a, let c=beta(a) be its residents and e=a-c its emission
heads. Release c particles into the pure rotor system at heads e,
and let D_w(a) be the final heads. Here a rotor head denotes the
LAST-used occurrence; it advances before the next departure.
For every input b>=0,

    D_w(T_b^buffer a)=T_b^rotor D_w(a).              (J106)

Indeed interpret pending buffer letters x and stored residents c
as x+c rotor particles. A buffer step emits 1+c_i-c'_i<=1+c_i
arrows while at least 1+c_i rotor particles are available there.
Those departures are legal, including the zero- and two-emission
cases and self-loops. The remaining rotor particles are the final
residents. Draining them completes the same abelian routing. The
emission vectors add exactly, not only modulo the periods.

Both cyclic networks have critical group Z^V/M Z^V and burning
input 1_V. Thus (J106) maps recurrent states equivariantly between
two free transitive group actions and is a bijection there. The
recurrent rotor heads h are precisely those whose last occurrences
j_i(h_i) form a sink forest. To see this directly, drain 1_V from h:
its odometer is at most m and each coordinate is positive. Return
requires exactly m. A completed drain orders last-departure times
strictly along internal last arrows, excluding cycles. Conversely,
an unconsumed suffix before m would be balanced without exits,
forcing a cycle among its last arrows. Parallel occurrences remain
distinct states; selected self-loops cannot lie in the forest.

There is an explicit inverse, including the reset qualification:

    D_w^(-1)(h)=stabilize 1_V in EMPTY reset buffers
                    with emission heads at the rotor forest h. (J106a)

Write I for the incoming-count map in the words shifted to h.
The pure rotor drain of 1_V uses exactly m. If w is any stable
buffer cut for reset input 2*1_V, its binary endpoint forces
w>=1_V+I(w). Rotor least action therefore gives w>=m. Periodicity
and unit divergence then give

    Fix(F_(2*1_V))=m+Fix(F_(1_V)).

No reset-index ambiguity occurs: the shifted bit indices on the
right are 1+I_i(w-m)>=1. If u is the physical least cut for 1_V,
the second full addition therefore emits exactly m and returns
the cyclic phase. The first output is recurrent. Completing its
residents by pure rotor drainage returns h, proving (J106a).
The forest hypothesis is essential to this proof.

Quotienting by simultaneous local word rotations gives a bijection
at each fixed (J,k):

    complete drain book d
      -> reset stabilization of 1_V, emissions u, endpoint C
      -> recurrent compatible pair (R_u d,C).       (J106b)

The inverse releases C into pure rotors from head zero, with
emissions r, then rotates by r. The resulting last arrows form a
forest, hence the new book is a complete drain. In both directions
u+r=m, so the rotations cancel and all arrow and bit counts are
preserved. Words of smaller period cause no multiplicity issue.
This correspondence and its inverse passed an independent AI proof
audit. The focused verifier checks 27,648 addition intertwinings,
7,680 forest/reset inverse pairs, and the complete bijection between
640 decorated drain books and recurrent compatible pairs in (J94b).

This also identifies the conditional laws in the original iid
construction. Define J by the no-sleep drain from 1_V on its shared
arrows. Given J, that drain book is uniform among its last-forest
orders, while its bits remain iid; given also k, their orders are
uniform. Physical ARW from 1_V uses only this finite book. Thus
its conditional resident law is exactly nu_(J,k), and (J105) is an
occupation bound conditional on the actual auxiliary drain inventory
and coin counts. Averaging gives (J103) directly. Replacing complete
drain words by all envelope words would change this law.

The proposed stronger physical sign after averaging only local
paired-word rotations is false by (J110) below. Two earlier bounded checks
were positive.
The 768 words in (J94b) split into 96 rotation orbits: 64 of size six
and 32 of size twelve. Every orbit has a nonnegative signed increment
at every target for {0}->{0,1}. Separately the single word in (J96)
has a 30-word orbit, with row periods (3,5,2). Its eight input totals
are exactly (20c_0,18c_1+9(1-c_1)c_2,(15+9c_1)c_2), so all 36
cover-target sums are nonnegative. The focused verifier retains
these physical least-cut checks. Neither recurrence nor the forest
bijection supplies the general rotation-orbit sign.

For a loop-free directed simple cycle the entire stable-cut order
is simpler than the general lattice. Let a_i>0 be its forward
probability, ell_i=1-a_i its leak, and kappa=product_i a_i<1.
For fixed source C and fixed stacks, the cut equation at i is
m_i=c_i+S_(i-1)(m_(i-1))-beta_i(c_i+S_(i-1)(m_(i-1))), where
S_i counts forward arrows. Its right side is nondecreasing in the
previous coordinate. Ordering any two cuts at one coordinate
therefore propagates around the cycle: ALL its cuts form a chain.
This is not a comparison between different sources.

Let A(C,B) sum the normalized weights of nonempty remainders having
no proper nonzero cut. Every nonempty remainder splits uniquely
at its first nonzero cut into such an atom and a remainder. Conversely
concatenation cannot introduce a cut incomparable to that atom.
The weights multiply, and distinct atom cylinders are disjoint.
Hence, on the full support space,

    R=I+AR, R^(-1)=I-A,
    A>=0, A1=1-f, rho(A)<1.                         (J107)

Thus R^(-1) is a nonsingular M-matrix, strictly row diagonally
dominant with excess f_C in row C. For O_v=diag(1_(v in C)),

    G_v=R^(-1)O_v R=O_v+(O_v A-A O_v)R,
    (2O_v-I)(G_v-O_v)>=0 entrywise.                (J107a)

A target-vacant row is -sum_(D containing v) A(C,D)R(D,.).
A target-occupied row is its identity row plus the analogous sum
over D avoiding v. These formulas include all particle-loss blocks
and impose no restriction beyond kappa<1. They extend to the closed
sleep cube using the proved continuous inverse and positive lower
bound on f. This sharpens (J73); it does not assert the same statement
for original uncollapsed self-loop books. The chain and factorization
passed an independent AI audit; the operator consequences received
a separate proof check.

The cycle also has a routing-only formula for ALL remainder sectors:

    R_p=I+Hcal diag(mu_p).                           (J108)

Every nonempty remainder has positive forward current at every
site. Indeed zero incoming current forces a zero local word: its
only possible positive length, a 1-to-0 row of length one, would
require its sole bit to be both one and zero. Propagate around
the cycle. For positive currents, each row's summed normalized bit
weight is p_i when its terminal site is occupied and q_i otherwise.
In a 1-to-0 row, the positive incoming current guarantees distinct
first and last bits. This proves that Hcal is independent of p.

For k=|C|-|B|>=0, sum over exit vectors e>=0 with |e|=k. Let h be
the integer height vector with min h=0 and
h_i-h_(i-1)=1_C(i)-1_B(i)-e_i. Then

    Hcal(C,B)=sum_e product_i a_i^h_i ell_i^e_i
               * sum_(s>=1) kappa^s
                    product_i binom(s+h_i+e_i,e_i). (J108a)

For k<0 the entry is zero. The positive forward currents are exactly
h+s*1, s>=1; the binomial factors count all local arrow orders,
without a chronological interleaving factor. Their product is a
degree-k polynomial in s with nonnegative coefficients. Each s^r
can be replaced by (kappa*d/dkappa)^r[kappa/(1-kappa)], keeping
the prefactor outside that derivative. Thus this is a finite rational
formula with denominator dividing (1-kappa)^(k+1). The common bit
factor and routing enumeration passed an independent audit and
agreed in all 64 entries with the complete-block/IDLA-pairing formula
on the prescribed heterogeneous single-leak three-cycle.

These stronger operator signs do not give an entrywise positive
add-one commutator. Put g=G_v f=(I-A)1_v, Xi=R^(-1)A_x^p R,
Xihat=diag(f)^(-1)Xi diag(f), and u=g/f coordinatewise. Then
Xihat is stochastic, u_C<=0 when v is absent, and u_C>=1 when
it is present. For the represented probability row nu_eta=Q_eta diag(f),

    q_v(eta+e_x)-q_v(eta)
        =nu_eta (Xihat-I)u.                         (J109)

The required subharmonic average on these particular rows remains
open. In fact h=[Xi,G_v]f has both signs whenever a particle from x
can reach v. Full-input saturation gives the strictly positive
invariant row Q_full=pi_p R and Q_full h=0. For empty input, Q_0
has only the empty-support entry R(empty,empty), so
h(empty)=q_v(e_x)/R(empty,empty)>0. Some other entry must be
negative. This obstruction uses no cycle-specific hypothesis; it
rules out entrywise commutator positivity on every accessible graph
with x reaching v and interior sleeps. It does not rule out the
physical averages (J109).

### Rotation averages fail; exchanging cyclic orders can repair them

**Exact obstruction and proved repair family, 2026-10-20.** Local paired
rotations preserve too much information for a general occupation comparison.
For two vertices, with sink denoted s, take

    arrows_0=(0,s,s), bits_0=010;
    arrows_1=(0),     bits_1=0.

The lengths are (3,1), with unit divergence at both vertices. Compare
reset inputs {0} and {0,1}, targeting 0. The three rotations at 0 give
the following least emission vectors and target indicators:

    shift   lower emissions   upper emissions   lower/upper target
      0         (1,0)             (3,1)                1/0
      1         (0,0)             (2,1)                1/0
      2         (1,0)             (2,1)                0/1.

Thus the orbit-averaged increment is -1/3. This conditioning is not
the iid ARW law. If the arrow order and the position of the single
one bit are instead chosen independently, all nine words have the
same inventory and coin counts. Four occupy 0 from each nonempty
input; the full-order increment is zero. The nine words split into
three rotation orbits with signed sums -1,+1,0. Even permuting the
same three paired records joins the negative and positive orbits
and cancels their deficits.

There is a loop-free family on a genuine directed cycle of every
size N>=2. Put L=N+1. Row 0 has one forward arrow to 1 and N sink
arrows; row i, 1<=i<N, has i+1 copies of i->i+1 modulo N. The
lengths have unit divergence. All non-target bits are zero. Put
the sole one at position d+1 in row 0, whose forward arrow is at
position one, with 1<=d<=N. All initial particles reach 0, and its
one forward arrow permits at most one further circuit.

For r=|C| in {1,...,N}, let t be the forward-arrow position and
beta the target bit row. The physical target indicator is

    B_r=beta_r                         if t>r-beta_r;
        beta_(r+1)                     if t<=r-beta_r.

Indeed r-beta_r is the candidate departure count without a circuit;
it either stops before the forward arrow or produces exactly one
additional arrival. Under the L paired rotations, the one bit is
cyclically d positions after that arrow. A one at r survives exactly
when r<=d; a one at r+1 is reached exactly when r>=d. Consequently

    q_rotation(0)=0,
    q_rotation(r)=[1+1{r=d}]/L, 1<=r<=N.             (J110)

Every cover from size d to d+1, d<N, has increment -1/L, including
full-minus-one when d=N-1. Thus removing self-loops or restricting to
cycles does not rescue rotation averaging.

The entire family's full-order inventory sign nevertheless holds.
Allow any fixed target coin count k, 0<=k<=L, and average arrow and
bit orders independently. First average the displayed B_r over t:

    E_t B_r=beta_r+(r-beta_r)(beta_(r+1)-beta_r)/L.

Exchangeability of adjacent bits gives, for every nonempty input,

    q_all(r)=k/L+k(L-k)/[L^2(L-1)].                  (J110a)

This is independent of r. If the bit b attached to the forward arrow
is also fixed, averaging all permutations of the paired records gives
instead k/L+1{b=0}k/[L(L-1)], again independent of r: terms with t<r
cancel by exchangeability, and only t=r with b=0 contributes a
correction. These are repairs for this family, not a general paired-
permutation or inventory theorem. The all-orders envelope coefficient
and the original occupation conjecture remain open.

The focused verifier checks the two-vertex table by legal queue
execution and least-cut iteration. It also checks (J110) on four
vertices for d=1,2,3 and (J110a) for every target coin count there,
using every binary input. Separate agent calculations additionally
enumerated all bounded cuts of these witnesses. The formulas above
give the arbitrary-N proof.

### Recurrent completion and an exact full-minus-one criterion

**Proved identities; the comparison remains open.** Put R=V\{x}. In a
complete decorated drain book d of unit divergence, let v(d) be the
pure-rotor odometer for input R, and rotate both local rows by v(d).
Call the result Phi_x(d). This permutes the drain books at fixed arrow
inventory and coin counts: rotor addition of R permutes recurrent
heads, and its inverse is addition of x. More explicitly, the remaining
suffix has divergence e_x and last-arrow forest, so it is the complete
x-to-sink walk with odometer m-v(d). The two rotations cancel. This
also works for words with smaller periods.

Let L(d) and U(d) be the physical reset states from R and V, retaining
their heads and any unvisited empty-reset states. Then

    residents U(Phi_x(d)) = residents T_V L(d),
    T_x L(d)=U(d),  T_x T_V L(d)=U(d).                (J111)

For the first equality, every stable cut w for input R+V satisfies
w>=R+I(w), because its endpoint is binary. Rotor least action gives
w>=v(d). Subtracting v gives a bijection with the stable cuts for
input V in the shifted rows. All shifted queried bit indices are
1+I_shift>=1, so there is no reset-index exception. Least cuts
correspond. The shifted book remains a drain book, hence its full
output is recurrent by (J106a). Commuting additions and T_V U=U give
the second line. Thus one full burn makes these particular lower
outputs recurrent; no statement about arbitrary transient states is
needed.

Uniform averaging over all drain books at fixed (J,k) consequently gives

    E[U_v-L_v]=E[(T_V L)_v-L_v].                     (J111a)

This is the actual conditional law of the iid full-drain construction
in (J106b). The recurrent completion and lower state coalesce after
addition at x, but their resident sets need not be ordered.

Even averaging over an entire recurrent rotor class is insufficient.
On four vertices with sink s, take paired rows

    arrows=( (s,1,3), (2,0,s), (3,1,s), (0,2,s) ),
    bits  =(   111,     010,     100,     111   ).

Among the 3^4 paired rotations, exactly 45 have last arrows forming a
sink forest. For target 1, full input occupies 20 of these books;
input V\{0} occupies 21. Thus the conditional increment is -1/45.
All vertices are in one loop-free SCC. Averaging larger pools, always
restricted to last-arrow forests, repairs this example:

    pool                              books   upper/lower target totals
    fixed paired cyclic orders           45            20/21
    all paired-record permutations      720           332/326
    independent arrow/bit orders        6480          3004/2512.

The latter two increments are 1/120 and 41/540. The focused verifier
checks these totals by legal buffer execution and least-cut iteration.
This rules out a proof confined to each recurrent class of paired cyclic
orders; it does not refute (J111a) after all orders or its iid mixture.

There is a separate exact product-measure formulation. Write calR for
the remainder matrix, f=calR^(-1)1 and gamma_v=calR^(-1)1_v. Put

    h_x=0, h_i=P_i(sink before x) for i!=x,
    sigma_x=K(x,s)+sum_i K(x,i)h_i>0,
    Z_x=product_(i!=x)(p_i+q_i h_i),
    theta_i=p_i/(p_i+q_i h_i), i!=x.

Let D be a product Bernoulli(theta) subset of R and d=det(I-K). Then

    q_v(V)-q_v(R)
      =(p_x Z_x/d) E[gamma_v(D+x)-sigma_x gamma_v(D)],
    E[f(D+x)]=sigma_x E[f(D)].                       (J112)

The averaging law is explicit; h_i=0 simply forces i into D.
To prove the identities, use the squarefree addition algebra, with
Lambda extracting its full-support coefficient. IDLA with only x
vacant gives, for every polynomial a,

    z_R a(z)=a(h)z_R+[a(1)-a(h)]z_V.

The constant-term formula (J29) yields
(Q_V-Q_R)(B)=mu_p(B)w_(B^c)^(p/q)(h)/d. The factors w_i(h)
equal h_i for i!=x, whereas w_x(h)=-(p_x/q_x)sigma_x.
For B=D+x and B=D the two coefficients are therefore
(p_x/d)mu_(p,R)(D)h_(R\D) and -sigma_x times that number.
Normalize their common product weight by Z_x, then multiply by
calR^(-1)1_v or f. This proves both formulas without semisimplicity.

Their sign still needs averaging over the whole product law. On
0->1, 1->0/s equally, with sleeps one half, gamma_1(empty)=0 but
gamma_1({0})=-8/49. Even averaging just the target coordinate
need not repair every remaining background. For the bidirected
four-cycle K(i,i-1)=K(i,i+1)=K(i,s)=1/3, p_i=1/2, x=0,v=2,
fix D\{2}={1,3}. Here sigma_x=5/7 and theta_2=7/12; the
target-averaged bracket is
-87744420438699/6964520196491312. The actual full-minus-one
increment is nevertheless 6235498643/991955589872>0.
These are obstructions to smaller averages, not to (J112)'s sign.
An agent's exact algebraic checks verified (J112), including its
normalization, on two four-site kernels, one with mandatory portals;
the displayed derivation is the general proof. The maintained verifier
now checks 192 formal coefficients, 48 target identities and 12 balances,
including the uniform-sleep obstruction just displayed.

### One-coordinate release interpolation

**Exact identities; two sufficient sign conditions remain unproved.**
Fix a vertex i and all other sleep probabilities. Let P_0 be binary
stabilization with p_i=1, and let A_j,J_j be its addition and forced-jump
operators. The kernel W_i from (J53) fixes supports without i; from B
containing i it removes that resident, forces its K-jump, and stabilizes:
delta_B W_i=delta_(B\i)J_i. Put s=(1-p_i)/p_i and z=s/(1+s).
The p_i=1 boundary uses the continuous extension: failed sleeps when
two or more particles meet are skipped before their compulsory jump.
Inverting (J53)'s affine identity gives exactly

    P_s=P_0[I+s(I-W_i)]^(-1)
       =(1-z) sum_(k>=0) z^k P_0 W_i^k.             (J113)

The series converges in matrix norm for z<1. Thus the required comparison
is a geometrically weighted comparison of release rounds, with all
other parameters unchanged.

The absorbing projection is explicit. Let V_i remove i from a sleeping
support. Then

    W_i A_i=A_i,   lim_(k->infinity)W_i^k=A_i V_i.   (J113a)

For i in B, the first identity follows from
delta_(B\i)J_i A_i=delta_(B\i)A_i J_i=delta_B J_i=delta_B A_i;
on i-vacant rows it is immediate. The last equality uses p_i=1 and
the quadratic addition relation. The block Q_i on i-occupied states
is transient by (J53). Hence a W_i-harmonic kernel is determined by
its values on i-vacant states. A_i V_i is harmonic and equals the
identity there, proving the limit.

When i belongs to both initial inputs, the p_i=1 endpoint reduces to
ARW on V\i with the ordinary trace kernel: retain the initial resident
at i and route every later arrival through it. At p_i=0, include also
the ordinary walk of that initial particle until its first visit to
V\i or the sink. Thus both endpoints of a common-root comparison
reduce to a smaller graph. Positivity between the endpoints does
not follow just from this reduction.

Fixed release times already fail. On the two-site kernel
K(0,1)=K(1,0)=K(0,s)=K(1,s)=1/2, with p=(1,1/2), compare
{1}->{0,1}, target 1, and release i=0. In support order
empty,{1},{0},{0,1}, the exact kernel is

    W_0=[[1,0,0,0], [0,1,0,0],
         [5/8,1/4,1/8,0], [1/16,5/8,19/112,1/7]].

After one release the target probabilities are 9/16 and 107/196,
giving increment -13/784. In fact, writing delta for the difference
of the two P_0 rows,

    delta W_0^k 1_1=(3/2)(1/8)^k-(10/7)(1/7)^k.

Every k>=1 has a negative increment, but the cumulative sums stay
above their limit 1/21. The actual geometrically averaged increment is
2(s+2)/[(6s+7)(7s+8)]>0. This demonstrates necessary cancellation,
not an ARW counterexample.

A possible sufficient statement is

    delta sum_(k=0)^M W_i^k 1_v>=0 for every M>=0,   (J113b)

where delta=(delta_(C+x)-delta_C)P_0 and i is a common initial root.
Abel summation would give the main inequality:
(1-z)sum z^k delta W_i^k 1_v
=(1-z)^2 sum z^M delta sum_(k=0)^M W_i^k 1_v.
The terms grow at most linearly, so this rearrangement is justified.
The fixed-time obstruction does not refute (J113b).

There is a distinct finite coefficient condition. Let L_i=I-W_i,
D_i(s)=det(I+sL_i), and H_k=[s^k]adj(I+sL_i). The row Laplacian's
forest formula gives H_k>=0 and H_k1=d_k1, with
d_k=[s^k]D_i(s)>0 for 0<=k<=m=2^(|V|-1). Its forests live on
sleeping-support states, with off-diagonal weights W_i(B,B'), not on
the original graph. They have k arcs and send each state to its root.
The coefficients follow the finite recursion
H_0=I, d_k=tr(L_i H_(k-1))/k, H_k=d_k I-L_i H_(k-1).
See [Chebotarev--Agaev, §§2,4–5, Theorem 3' and Proposition 4](https://arxiv.org/pdf/math/0508178);
these definitions, statements and proofs were read.

Consequently

    D_i(s)[q_v(C+x)-q_v(C)]
       =sum_(k=0)^m s^k delta H_k 1_v.              (J113c)

Nonnegativity of every coefficient is another sufficient condition,
not a consequence of H_k>=0 because delta is signed. The last kernel
H_m/d_m is (J113a)'s absorbing projection. These forest-degree
conditions and the cumulative-time condition (J113b) are different.
Neither is asserted in general; the release identity itself imposes
no such signs.

### Integrated release: a two-site theorem and a numerator reduction

**Proved for all two-site kernels and horizons.** First remove self-loops.
Let i=0 be the common root, K(0,1)=a, K(1,0)=b, with the remaining
probabilities going to the sink, ab<1, and p_1=p in (0,1). The only
nontrivial common-root cover is {0}->{0,1}. Define

    rho=ab(1-p), tau=p/(1-rho), sigma=ab tau.

With root 0 permanently occupied, the other initial particle occupies
1 with probability tau. Thus the difference of the two P_0 laws is
tau(delta_{0,1}-delta_{0}). In order empty,{0},{1},{0,1}, the release
row from {0} is (1-rho-ap,rho,ap,0). From {0,1} it is

    (a(1-b)^2(1-p), c, e, sigma),
    c=ab[1-tau+(1-b)(1-p)], e=1-a+a(1-b)p.

For the latter row, a jump to 1 wakes its resident; its first departure
either fills 0 or exits, leaving a one-particle continuation. These two
possibilities give the displayed probabilities. Solving the two
transient-state reward recurrences and summing release prefixes yields

    sum_(M>=0) A_M(1) z^M
      =tau(1-az)(1-abz)/[(1-z)^2(1-rho z)(1-sigma z)],
    sum_(M>=0) A_M(0) z^M
      =tau ab[1-b(1-p)]z/[(1-z)(1-rho z)(1-sigma z)],
    A_M(v)=delta sum_(k=0)^M W_0^k 1_v.              (J114)

Both series have nonnegative coefficients. For the first, use
(1-cz)/(1-z)=1+(1-c)sum_(j>=1)z^j for c=a,ab; the other factors
are positive geometric series. The second is immediate. This proves
(J113b) for every loop-free two-site kernel, with either common root
by relabelling. Individual rounds need not be positive.

Self-loops can be removed for this question on any graph. At a non-root
j the usual elimination replaces p_j by p_j/[1-(1-p_j)K(j,j)] and
normalizes the non-self routing. At the common root, alpha=K(i,i)<1
gives P_0^old=P_0^new and W_i^old=alpha I+(1-alpha)W_i^new.
Indeed a forced self-jump restores the resident immediately; other
self-jumps are removable waiting steps. Therefore

    A_M^old(v)=sum_(j=0)^M b_(M,j) A_j^new(v),
    b_(M,j)=binom(M+1,j+1)(1-alpha)^j alpha^(M-j)>=0. (J114a)

Expand the lazy powers, then sum by parts. The coefficient of
(W_i^new)^j before summation by parts is
sum_(k=j)^M binom(k,j)alpha^(k-j)(1-alpha)^j; its difference from
the next coefficient is b_(M,j). This proves preservation of the
prefix condition and completes the two-site theorem with self-loops.

For a general graph, let T,A be the i-occupied and i-vacant states,
write W_i=[[I,0],[R,Q]], and restrict delta to T. Since i is a common
root, it has no A component. The exact prefix generating function is

    sum_(M>=0) A_M(v)z^M
      =delta_T adj(I-zQ)[(1-z)1_(v,T)+z R1_(v,A)]
          /[(1-z)^2 det(I-zQ)].                     (J114b)

Absorbed states keep contributing their target reward at later rounds;
the term R1_(v,A) must be retained. The inverse determinant is always
coefficientwise nonnegative:

    1/det(I-zQ)=exp(sum_(k>=1) z^k tr(Q^k)/k).

This is a formal identity near zero, and Q>=0 makes every coefficient
on its right nonnegative. No real-eigenvalue or reversibility hypothesis
is needed. Positivity of the numerator after the two partial sums
represented by (1-z)^(-2) would suffice, but the final agent certificate
below disproves this stronger condition even after restricting to reachable
states. Determinant compensation cannot generally be discarded.

There is an exact all-horizon check beyond two sites. On
0->1->2, 2->0/s equally, with p=(1,1/2,1/2), common root 0,
C={0,2}, x=1, v=2, the transient Q in order
{0},{0,2},{0,1},{0,1,2} is

    [[1/8,0,0,0],
     [15/112,1/7,1/28,0],
     [15/112,1/7,9/56,0],
     [139/1120,13/70,117/560,1/5]].

Put D(z)=(1-z/8)(1-17z/56+z^2/56)(1-z/5)=det(I-zQ).
The prefix series is exactly
(16-15z+z^2)/[280(1-z)D(z)]. Its numerator divided by 1-z
has coefficients 16,1,2,2,..., proving every prefix positive.
In contrast its individual round increments begin
2/35,-173/9800. The agent also checked the release construction
against independent raw absorption with failed sleeps removed.
The maintained sleep-transport verifier checks the two-site formulas at
five parameter choices, 250 exact prefixes, and the full polynomial
identity for the displayed three-site generating function. These results
support the integrated inequality; they do not establish it on general
graphs.

### Forest-degree comparisons and two commuting boundary families

**Proved reductions and two-site coefficient theorem; general sign open.**
In the block notation of (J114b), put M=I-Q, m=2^(|V|-1), and
F_k=[s^k]adj(I_m+sM), with F_(-1)=F_m=0. Block inversion gives the
full-state coefficients from (J113c):

    H_k=[[d_k I_m,0],[F_(k-1)B,F_k]],
    F_0=I_m, F_k=d_k I_m-MF_(k-1).                  (J115)

The off-block term cannot be dropped. Since B1=M1, its row sums
together with F_k are d_k. Let T=M^(-1)B be the absorption kernel
on outside supports. For a target v!=i, with outside indicator h,
the occupied block's payoff is

    F_k h+F_(k-1)B h=d_k T h+F_k(I-T)h.             (J115a)

Here T is an actual addition kernel on the graph traced through i:
with i occupied at p_i=1, add at i and then remove its permanent
resident. Its added particle follows an ordinary walk from i to its
first outside vertex or the sink. This also follows from (J113a).
Neither M and T nor their conjugates are assumed to commute.

For two loop-free vertices i,j write K(i,j)=a, K(j,i)=b, p_j=t,
q=1-t, alpha=t/(1-abq), e=1-abq and f=1-ab alpha. Both e,f
are positive. The only common-root cover has delta equal to
alpha(delta_{i,j}-delta_i), and D_i(s)=(1+es)(1+fs).
Direct substitution of the release rows in (J114) gives

    k       delta H_k 1_i                  delta H_k 1_j
    0              0                              alpha
    1       alpha ab(1-bq)                  alpha(2-a-ab)
    2              0                       alpha(1-a)ef.     (J115b)

Every coefficient is nonnegative. The top coefficient also follows
from the one-site trace: at zero sleep at i, the lower initial
particle reaches j with probability a, whereas the upper already
has a particle there. Self-loop elimination changes the other sleep
parameter as in (J114a); a root self-loop of weight ell multiplies
degree k by (1-ell)^k. Thus the coefficient criterion is proved for
every accessible two-site kernel, not only for the evaluated resolvent.

More generally, if i lies on no internal cycle after removing its
self-loop, then Q=0 and W_i is a projection. Hence

    H_k=binom(m-1,k)I+binom(m-1,k-1)W_i.             (J115c)

Each comparison is a positive combination of the two endpoint
comparisons on the smaller trace graph. This proves the criterion
on DAGs and gives a conditional reduction whenever that smaller
graph's all-input theorem is known. It is not asserted for cyclic i.

There are further positive identities which retain the boundary
information. Remove the root self-loop and write physical addition
at x!=i, at p_i=1, in blocks:

    A_x=[[U_x,V_x],[0,tilde A_x]],  A_i=[[0,I],[0,T]].

U_x is addition killed on any visit to i: once visited, i stays
occupied. The kernel tilde A_x is physical addition on the trace
graph. Commutation with A_i and the first forced jump give

    tilde A_x=U_x+V_x T,   tilde A_x T=T tilde A_x,
    B=K(i,s)I+sum_(j!=i)K(i,j)U_j,
    Q=sum_(j!=i)K(i,j)V_j,  B=(I-Q)T.               (J115d)

All U_x commute with B, since the killed additions commute. All
tilde A_x commute with T. The two families generally do not commute
with each other. These are concrete positive boundary decompositions
for a possible forest/path exchange, not a proved coefficient sign.

Conjugating M into the trace input basis does not by itself supply an
M-matrix argument. On the bidirected four-cycle with equal probabilities
1/3 to each neighbour and the sink, root i=0, and other odds
(r_1,r_2,r_3)=(10,1/3,3), let P_c be binary stabilization on the
trace graph. Then
(P_c M P_c^(-1))({3},{1})=56259620/21199632657>0.
This rules out that particular off-diagonal sign premise, while
leaving (J113c)'s physical comparisons open.

The maintained [forest verifier](verify_release_forest.py) implements the
block recurrence, endpoint/trace identities and two-site formulas, including
self-loops. Its focused suite passed in the main thread: 1,740 degree
coefficients on two prescribed kernels and an independent 127-transient-state
raw calculation checking all 256 entries of P_0 and 128 nontrivial entries
of W on the bidirected four-cycle. The agent's nine-kernel calculation gave
21,876 nonnegative coefficients (21,102 positive, 774 zero); `--full`
reproduces that prescribed suite. These are finite certificates.

### Release obstruction, paired records, and verification status

**Verification status, 2026-11-08.** The numerator obstruction (J116)
was independently reproduced in the main thread by the maintained
`verify_sleep_transport.py::numerator_shortcut_audit`: seven reachable
states and the exact values of A_0, A_1 and the negative coefficient below.
The correlated-routing arguments remain saved agent reports awaiting
independent audit. Neither settles the general occupation conjecture.

For (J114b), restrict Q to supports reachable from either physical input
law, set h=1_(v,T), beta=R1_(v,A), and write

    D(z)=det(I-zQ),
    N(z)=delta_T adj(I-zQ)[(1-z)h+z beta].

The sufficient condition `[z^m]N(z)/(1-z)^2>=0` fails. The agent's
exact rational certificate uses rows with the sink last:

    K=[[0,1/2,1/2,0,0],
       [0,0,0,1,0],
       [1/2,0,0,0,1/2],
       [1/2,0,1/2,0,0]],
    p=(1,1/5,4/5,1/2), i=0, C={0,3}, x=1, v=2.

The reachable supports are
{0},{0,3},{0,2},{0,2,3},{0,1},{0,1,3},{0,1,2}, with respective
signed initial masses
(-17323/970256,-737713/1455384,-17323/121282,3648/8663,
-41269/2910768,1178/8663,1080/8663).
Here tr(Q)=14775529/9702560, A_0=48869/121282 and
A_1=9220242269537/17651188228800>0, yet

    [z]N(z)/(1-z)^2=A_1-tr(Q)A_0
                   =-13280929/145538400<0.          (J116)

The agent also verified gcd(N,D)=1 and checked P_0,W against a
99-transient-state raw absorption solver. All prefixes through M=24
were nonnegative. Thus this refutes only the stronger numerator
condition. Exact numerator certificates did establish all-horizon
positivity for all 48 common-root-0 cover/target comparisons on
0->1->2->3, 3->0/s equally, p=(1,1/2,1/2,1/2).
The original temporary reproduction files are not available in the current
session. The kernel and obstruction above, and the maintained focused
verifier, are the durable evidence. The gcd, 99-state raw audit, and
48-comparison all-horizon claim remain the original agent's reported checks;
they were not rerun for the results compilation.

**Paired-pool bridge; signs still open.** Fix a unit-divergence multiset
of paired records (arrow,bit) at each vertex. Average physical reset
endpoints over all distinct permutations of those pairs. The proposed
cover sign, and its full-minus-one version restricted to last-arrow
forests, remain open; paired rotation obstructions do not refute them.
For iid pairs put w_i(0,j)=q_i K_i^0(j), w_i(1,j)=p_i K_i^1(j),
q_i=1-p_i and Kbar_i=q_i K_i^0+p_i K_i^1, with accessible sink.
An empty arrival sleeps with probability p_i, reserving its K_i^1
arrow, or emits a K_i^0 arrow. At an occupied arrival, emit the old
reserved K_i^1 arrow and independently process the new bit in the
same way. The pending arrow is independent of the exposed history.
Fixed-stack commutation and a dominating ordinary drain for Kbar
give almost-sure halting.

Let A_i denote physical addition and
J_i^b=K_i^b(s)I+sum_j K_i^b(j)A_j. Processing two arrivals gives
A_i^2=A_i(p_i J_i^1+q_i J_i^0). On an i-vacant support B, its
first-arrival law is p_i delta_(B+i)+q_i delta_B J_i^0.
In the corresponding quotient algebra, b_i=(z_i-q_i J_i^0(z))/p_i;
the products b_B form the stable-support basis, and their coefficient
matrix E in the squarefree z-basis satisfies EP=I. Ordinary IDLA
with kernel Kbar supplies multiplication in that z-basis. Do not
interpret b_i as a stochastic operation on already occupied supports.

Applying the arrow-envelope theorem to Kbar and then conditioning
the bits on their arrows gives decorated-word probabilities

    det(I-Kbar_internal) product_(i,t) w_i(bit_i(t),arrow_i(t)).

For an ordinary full drain, replace the determinant by the indicator
that the last arrows form a sink forest. Products are constant on
each paired-permutation pool, giving the uniform conditional law.
Consequently universal paired-pool positivity would imply binary
occupation comparisons for every such correlated iid model.
A negative correlated-model comparison would force a negative pool.
The converse implication is not justified. The marginal inventory
law uses Kbar; conditional bits are arrow-dependent, not iid p_i.

The agent supplied two scoped proofs. On a DAG, target v's arrival
count N_v is independent of its own row permutation. If that row
has k_v one bits among m_v records, its mean occupation is
(k_v/m_v)Pr(N_v>0). Least-cut monotonicity increases N_v when input
is added, proving every paired-pool sign there.
For a loop-free two-site correlated model, let
p=p_0,t=p_1,a=K_0^0(1),b=K_1^0(0),c=K_0^1(1),d=K_1^1(0),
u=(1-t)b, Delta=1-(1-p)au and H=1-puc-(1-p)a(td+u).
Writing Q_0(C) for occupation of 0, direct first-step equations give

    Q_0({0})=p/Delta, Q_0({1})=up/Delta,
    Q_0(V)-Q_0({0})=p(1-p)atd(1-u)/(Delta H)>=0.

Indeed, with L=p/Delta and U=Q_0(V), the full-input equation is
U=p(1-u)+pu[cU+(1-c)L]+(1-p)[1-a(td+u)]uL
+(1-p)a(td+u)U. Subtract HL and use L Delta=p.
Here H>=pt>0; symmetry gives every binary cover. This is an iid
statement, not finite two-site inventory positivity.
The agent's separate nine-state raw calculation at
(p,t,a,b,c,d)=(2/5,3/7,2/3,3/5,1/4,4/5) gave
Q_0(V)=9574/18271 and gain 1104/18271.
Its bounded exact checks found no negative sign in 81 paired pools
on the bidirected four-cycle, nor in the 81 corresponding iid models
with p_i=1/3, deterministic K_i^1 and K_i^0 uniform on the other
two left/right/sink destinations. These are agent-reported finite
checks. The displayed argument is preserved here; the additional finite checks retain their original reported status.

**Cycle arm and continuation.** The cycle first-chance resolvent was
already in imported N3 evidence; it was not discovered anew. Its
killed return kernel is exactly the occupied block Q_0 of W_0.
For inputs containing 0, its entrance law is the occupied-output
block of P_0. This joins that arm to (J113). On every directed three-cycle
with arbitrary forward leaks, p_0=1 and interior off-root sleeps, Q_0 now
has an exact nonnegative Bernstein certificate for preservation of
nonnegative increasing functions. `verify_cycle_return.py` independently
reconstructs all six polynomials and 2,232 coefficients; its separate
28-state raw check also passes. The proof and domain are collected in
[the reader](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-monotonicity.pdf), the three-cycle killed-return theorem.
General cycle lengths and entrance-law domination remain separate open
conditions. Physical trace
addition itself need not preserve inclusion, so that shortcut is
unavailable. The old temporary cycle report is not a maintained source.
The same collection gives the positive boundary block interpolation for
physical additions, assembled from the saved release identities; the
independent raw audit checks 256 addition and 128 binary-stabilization
entries on two 44-state chains. The remaining occupation forcing is signed.
On resumption, prioritize a release forest/path exchange retaining
both boundary families in (J115d), or cycle-specific cancellation
retaining determinant compensation; paired pools offer a distinct
direction through the correlated-routing bridge. First audit the
remaining paired-routing arguments above. Retain the earlier exact-pair
remainders in (127f); their omission remains disproved.


## Exact audit

Run `python3 -B verify_response_renewal.py`.
It uses exact rational raw-ARW absorption, with no history truncation,
on a two-site path, a three-site path, a heterogeneous directed
three-site graph, and CB's four-site graph. It checks (1)–(2) for every terminal payoff basis vector,
(5) entrywise, (6) on every formal response basis vector, and (4) with
strictly positive defects. A further exact audit checks (7) for four
all-active inputs, including a multiplicity and a target-free input,
and verifies the conditional-order counterexample above. These checks
audit the formulas and their link
to the saved signed M; the preceding arguments prove the general theorem.
