## Lattice Supersymmetry and Order-Disorder Coexistence in the Tricritical Ising Model

In general, this paper was the first real introduction into the language of one dimensional physics.What it sought to show was that the interacting one dimensional ising model with an interaction term of $\gamma_i \gamma_{i+1} \gamma_{i+3} \gamma_{i+4}$ shows some interesting physics, in particular:
1. Self-duality requires no parameter tuning
2. Coexistence of phases with ordered and disordered pairing - in other words, coexistence of topologically ordered ground states with topologically disordered ground states.
3. Tricritical point with couplings of order 1 (difference between parameter strength is same order of magnitude) vs 250x times larger int/non-int of other models
4. Is written in the same form as that of a super-symmetric field theory, which opens the question for different models that include super-symmetry built into the theory, not just as an emergent property.
Using the Jordan Wigner transformations, we can show that the Majorana interaction term of $\gamma_i \gamma_{i+1} \gamma_{i+3} \gamma_{i+4}$ is mapped onto a three spin interaction term of $$\gamma_i \gamma_{i+1} \gamma_{i+3} \gamma_{i+4} \rightarrow \sigma^x_{i}\sigma^z_{i+1}\sigma^z_{i+2} + \sigma^z_i \sigma^z_{i+1}\sigma^x_{i+2},$$
The full Hamiltonian of the system studied is reduced down to:
$$H = \lambda_1 \text{pair interaction} (H_1) + \lambda_3 \text{three-spin interaction} (H_3) + \text{Ground state energy} E_0$$
with parameters $\lambda_1, \lambda_3$ . This can be rewritten as a sum of squares of 2 operators, which are analogues to supersymmetry generators. A SUSY generator acts as follows:
![[Pasted image 20250717164816.png]]
The simulations of DMRG seeked to reproduce some well known CFT constants for both the 2D Ising model as well as the tricritical generalized Ising model. It was able to do so successfully and showed that the nontrivial critical point in $H_3$ is in the tricritical Ising universality class
- This means that the long range interaction of the point in question is explainable by the tricritical Ising CFT
## Emergent Supersymmetry from Strongly Interacting Majorana Zero Modes

Paper studied the Hamiltonian of $$H = it \sum_j \gamma_i \gamma_{i+1} + g \sum_j \gamma_i \gamma_{i+1} \gamma_{i+2} \gamma_{i+3} $$ which corresponds to nearest neighbour interactions as well as a a 4 block of interactions. It is shown that for $g > 0$, this hamiltonian contains supersymmetric phases. 

Important things discovered in this paper is the general methodology of DMRGs with CFTs. 

## Critical properties of the Majorana chain with competing interactions

This paper went on to study the Hamiltonian of $$\mathcal{H} = it \sum_a \gamma_a \gamma_{a+1} - g \sum_a \gamma_a \gamma_a+1 \gamma_a+2 \gamma_a+3 -f \sum_a \gamma_a \gamma_{a+1} \gamma_{a+3} \gamma_{a+4},$$ where the g interaction was studied in the previous paper. The paper was interested in the interaction between f and g. To study this, t was set to 1. The phase diagram of 

## Conformal data and renormalization group flow in critical quantum spin chains using periodic uniform matrix product states

In general, this paper studied the 