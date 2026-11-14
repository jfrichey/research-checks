# Sink reinjection and stationary rank order

Selected working proof from the ARW negative-association record (20 October 2026), retaining equation labels (50)–(58). The saturated-stationarity input is the credited Levine–Liang result. The harmonic-eigenfunction argument used below is in the [selected negative-association technical record](../arw-negative-association/report.md), under “Harmonic polynomial identities.” These are ordinary working proofs, not Lean formalizations. The full derivative and finite transform are tested by `../arw-negative-association/weak_interaction.py --reinjection-flow`; these finite checks do not establish configuration-wise order or negative association.

**Global reinjection: exact flow and all-rate rank monotonicity
(complete working proof, 2026-10-20).** Fix any finite killed kernel
`K`, finite positive site sleep rates, and a probability distribution
`rho` on its non-sink vertices. Write `k_u=K(u,sink)`, and, for `0<=z<1`,
replace each would-be absorption by a fresh jump with distribution `rho`
with probability `z`. Thus the internal kernel and its killing vector are

`K_z=K+z k rho^T`, `k_z=(1-z)k`.

Every `K_z` is still killed: successive reinjections give a geometric
number of finite original killed walks. Let `nu_z` be its stationary
stable law, let `T_z=sum_v rho_v T_(v,z)` be one addition with drive
`rho`, and let `N` denote multiplication of a measure by its particle
count: `(N eta)(S)=|S| eta(S)`. No reachability hypothesis on `rho`
is needed, since the saturated stationary law is invariant under every
addition. The full vector identity is

`d nu_z/dz=(N nu_z-T_z N nu_z)/(1-z)`.                    (50)

To prove it, let `F_(u,z)` be the saturated stabilization law marked by
the number of absorbing departures from `u`, divided by `k_z(u)`;
omit vertices with zero killing. The marked identity (45), with a
returned particle sent instead to vertex `v`, gives response
`(T_(v,z)-I)F_(u,z)` for a unit increase in `K_z(u,v)` replacing killing.
This follows by delaying the marked particle until the unmarked
stabilization is complete; the unused instruction suffixes are fresh.
Finite absorbing-chain differentiation gives the identity even for
trajectories with several possible marks, by summing over the marks.

There are initially `n=|V|` particles. Hence conservation gives the
exact marked-law identity

`sum_u k_z(u) F_(u,z)=(n I-N)nu_z`.

The perturbation speed in row `u` is
`k_u rho_v=k_z(u)rho_v/(1-z)`. Sum the responses and use
`T_z nu_z=nu_z` to obtain (50).

Let `(X,Y)` be a stationary pair consisting of an initial sleeper
configuration with law `nu_z` and the result of one independent
`rho`-addition. Put `A=|X|`, `B=|Y|`; these have the same marginal
distribution. For any rank threshold `1<=r<=n`, let `c=r-1/2`.
Equation (50) and stationarity imply

`d Pr_(nu_z)(N>=r)/dz`
` =E[(A-c)(1_(A>=r)-1_(B>=r))]/(1-z)>0`.                 (51)

Indeed, the derivative first equals
`(E[B 1_(B>=r)]-E[A 1_(B>=r)])/(1-z)`.
Replace the first term by `E[A 1_(A>=r)]` and subtract the zero term
`c E[1_(A>=r)-1_(B>=r)]`. The displayed integrand is zero without a
threshold crossing and strictly positive on either direction of crossing.

Strictness needs no irreducibility of the driven chain. The saturated
law has full support: route the unwanted initial particles to the sink,
keep one active at each desired site, then put these remaining particles
to sleep, a finite legal prescription of positive probability. Choose
`v` with `rho_v>0` and an `(r-1)`-site set not containing `v`.
From this positive-mass state, immediate sleep of the added particle at
`v` produces an upward crossing with probability at least
`rho_v lambda_v/(1+lambda_v)>0`.

Consequently the stationary rank law is strictly stochastically
increasing in the reinjection probability throughout `0<=z<1`.
All nonconstant increasing functions of rank increase strictly, by
finite threshold decomposition. In particular,

`d E_(nu_z) N/dz=E[(B-A)^2]/(2(1-z))>0`.                (52)

This follows from (50) and `E[A^2]=E[B^2]`. Integrating (51) compares
any two reinjection probabilities below one, with the same original
kernel, rates, and destination distribution.

If killing is possible at only one port `u`, choose `rho=delta_u`.
Then this perturbation only replaces killing at `u` by its self-loop.
The general response `U` in (42) is therefore strictly positive on all
nonconstant increasing functions of total particle count, at every
positive collection of sleep rates. The conversion is
`d nu_z/dz=k_u U` at the current kernel. This supplies an all-rate rank
case of the marked-departure order, not its full configuration order.
Neither ULC, coordinatewise stochastic monotonicity, nor NA follows from
rank monotonicity alone.

**Finite reinjection transform (complete working proof).** The same family
admits an exact transformation using only the original addition matrices.
Write `T_0=sum_v rho_v T_(v,0)`, and define a matrix `M_z` by its columns:

`B_z=(I-zT_0)/(1-z)`,
`M_z e_S=B_z^(n-|S|) e_S`.

Then

`M_z nu_z=nu_0`,
`M_z T_(v,z)=(1-z)T_(v,0)(I-zT_0)^(-1) M_z`.             (56)

Its inverse is a stochastic matrix with columns

`C_z e_S=[(1-z)I+zT_z]^(n-|S|) e_S`, `C_z=M_z^-1`.       (57)

In particular, `nu_z=C_z nu_0`. Formula (57) has a direct meaning:
given `S`, take a binomial number `Bin(n-|S|,z)` of independent
rho-additions, stabilizing with the reinjected kernel.

First suppose every original vertex has positive direct killing.
Use the joint harmonic eigenbasis from the earlier harmonic-polynomial
section. For `J subset V`, its active weights `a^J` equal zero outside
`J` and, within `J`, the probability of reaching the sink before leaving
`J`; put `a_rho^J=sum_v rho_v a_v^J`. If `b^J` are the associated
stable-site weights, the eigenfunction is
`h_J(S)=product_(i in S)b_i^J`.

Restarting an originally killed walk at rho multiplies every active
weight, and every first-step harmonic weight, by the same scalar
`r_J=(1-z)/(1-z a_rho^J)`. Indeed the return-to-rho escape probability
satisfies its one-variable renewal equation; a restart outside `J`
counts as failure. The stable weights also scale by `r_J`, so
`h_(J,z)(S)=r_J^|S| h_J(S)`. Since `h_J T_0=a_rho^J h_J`,
the column definition gives

`h_J M_z=r_J^-n h_(J,z)`.

The harmonic expectations of both saturated stationary laws vanish
for every proper `J` and equal one for `J=V`. This proves the first
identity in (56). The eigenvalue of `T_(v,z)` is `r_J a_v^J`;
the same row identity proves the second. The joint eigenfunctions form
a basis under positive direct killing, as proved earlier.

Averaging the intertwining over `v` shows
`M_z[(1-z)I+zT_z]=B_z^-1 M_z`. Apply its `(n-|S|)`th power to `e_S`
and use the column definition of `M_z`. This proves `M_z C_z=I`.
Each column in (57) is a probability measure, so `C_z` is stochastic.
For a general killed kernel, approximate by `(1-epsilon)K`, which has
positive direct killing everywhere. Both families of addition matrices
are continuous at the original killed kernels by their finite absorption
resolvents. Taking limits preserves all displayed identities, including
the stochastic inverse; no diagonalizability assumption is needed in
the limiting case.

For `0<=a<b<1`, the kernel at `b` is obtained from the kernel at `a`
by reinjection probability `(b-a)/(1-a)`. For `delta=(b-a)/(1-a)`, let
`B_(a,b)=(I-delta T_a)/(1-delta)`. Averaged intertwining gives the
explicit matrix identity

`M_(0,a) B_(a,b)=B_b B_a^-1 M_(0,a)`.

The two matrices B_a,B_b commute, being rational functions of T_0.
Applying the `(n-|S|)`th power to each column proves
`M_(0,a) M_(a,b)=M_(0,b)`, hence
`C_(a,b) C_(0,a)=C_(0,b)`.
This argument respects the column-dependent powers; it does not commute
M with the count operator. Rebasing gives the full cocycle
`C_(b,c)C_(a,b)=C_(a,c)`. These are stochastic transition
matrices for a Markov evolution in the reinjection parameter. Its
generator, acting on columns, is

`Q_z=(T_z-I)(nI-N)/(1-z)`.                                (58)

Thus a configuration receives a rho-addition, using the current kernel,
at rate `(n-|S|)/(1-z)`. On the stationary family, `Q_z nu_z` equals
(50). This exact continuation does not give an NA-preservation theorem.
The special input and the dependence of `T_z` on `z` remain essential.

In fact, the vector field in (50) does not preserve the entire ULC cone.
At `K=0`, uniform addition refreshes a uniformly chosen coordinate to
Bernoulli `p=lambda/(1+lambda)`, so its rank kernel `T_0` has transition
probabilities `p(n-j)/n` from `j` to `j+1` and `(1-p)j/n` to `j-1`.
For a binomial input `r_j=binom(n,j)t^j(1-t)^(n-j)`, the relative
vector field `((I-T_0)N r)_j/r_j` is quadratic in `j` with quadratic
coefficient `(t-p)/(n t(1-t))`. Writing `w_j=r_j/binom(n,j)`,
every interior normalized LC gap, initially zero, therefore has derivative

`2(p-t)w_j^2/(n t(1-t))`.

For `n=2,p=1/2,t=3/4`, the rank input is `(1/16,3/8,9/16)`, the
vector field is `(-3/32,-3/8,15/32)`, and the LC-gap derivative is
`-3/64`. This is an exact cone-preservation obstruction, checked by
`complete_graph.py --renewal`. It does not concern the actual stationary
trajectory, whose initial density is `t=p`; a continuation proof must
use additional stationary identities.

The exact `reinjection_flow` certificate compares (50) with the
nilpotent derivative of the IDLA law and independent raw addition
kernels. It uses a directed three-site kernel with self-loops and
unequal rates at `z=0,2/5`, and the endpoint three-path at common rate
one and `z=0,1/3`, with reinjection at its sole killed port. Each case
checks both raw toppling orders, every addition-stationarity equation,
the complete signed derivative vector, all rank tails in (51), and (52).
It also verifies the finite transform (56), its stochastic inverse (57),
and the intertwining for every individual addition, using independent
raw addition matrices at both the original and reinjected kernels.

