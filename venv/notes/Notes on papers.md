
## Introductory papers

### Lattice Supersymmetry and Order-Disorder Coexistence in the Tricritical Ising Model

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

Between thsi paper and the next paper, the most important things to learn are the techniques employed. 
### Emergent Supersymmetry from Strongly Interacting Majorana Zero Modes

Paper studied the Hamiltonian of $$H = it \sum_j \gamma_i \gamma_{i+1} + g \sum_j \gamma_i \gamma_{i+1} \gamma_{i+2} \gamma_{i+3} $$ which corresponds to nearest neighbour interactions as well as a a 4 block of interactions. It is shown that for $g > 0$, this hamiltonian contains supersymmetric phases. 

Important things discovered in this paper is the general methodology of DMRGs with CFTs. We want to use DMRG to find the ground states of the hamiltonian, and with that calculate specific values to then verify against theory. Some of the values that are predicted by CFT that we wish to verify numerically include:
- Central charge c
- Various energy difference ratios, in the form of $$\frac{E_1 - E_0}{E'_1 - E'_0}$$
The identification of these CFTs not only help verify that our method worked, but also helps to determine the phase of what we are working under. For example, the phase diagram of this study looks like:
![[Pasted image 20250718155150.png]]
As we can see, we are primarily concerned with identifying the phases based on the central charge. THe way that they went about doing this in the paper is by splitting up the work into 2 different boundary conditions - antiperiodic and periodic. Then they split up the hamiltonian into an odd and even section. They looked at the ground state and first excited state of these antiperiodic/periodic boundary conditioned and even/odd hamiltonians. 

With these different energies, they calculated using CFT what they expect, and then numerically found the same numbers at hte critical point - this proved that at the critical point, the system they were studying was at a TCI phase!
- How did they find the critical point? They simply simulated for different values of $t_c$ and checked tot see which one got closest to theoretical predictions.
### Critical properties of the Majorana chain with competing interactions

This paper went on to study the Hamiltonian of $$\mathcal{H} = it \sum_a \gamma_a \gamma_{a+1} - g \sum_a \gamma_a \gamma_{a+1} \gamma_{a+2} \gamma_{a+3} -f \sum_a \gamma_a \gamma_{a+1} \gamma_{a+3} \gamma_{a+4},$$ where the g interaction was studied in the previous paper. The paper was interested in the interaction between f and g. To study this, t was set to 1. The phase diagram of this Hamiltonian presented *nine* different phases (very high compared to what we are expecting). In general, the paper focuses on the following phases, which I might expect to come out of my own Majorana simulations:
- Ising-1/Ising-2: critical phases that behave in the Ising universality class (described by the CFT that describes the Ising model)
- Floating phase + Ising: This kind of phase is special. We see both incommensurate LL behaviour as well as Ising criticality. In other words, we have 

### Conformal data and renormalization group flow in critical quantum spin chains using periodic uniform matrix product states

In general, this paper looked less on new physics and more on discovering a new method. In particular, it saw that a Block-state ansatz of puMPS is a good option for 1D quantum spin chain analysis - this is in effect what I will be doing in my research so I will need to learn quite a bit from the master :))

It backs up most of it's claims with a description of the accuracy and it's reproduction of certain CFT variables, such as central charge and energy gap ratios. The paper is quite short but it features quite a hefty additional supplemental material for those looking to use their methods. They used these methods to study the "O'Brien-Fendley model" which is the physical model that studies the 4 particles over 5 sites interaction of the first paper in this review.

To extract the most information, I go straight to the supplemental material. It references a paper by Ganahl et al. called (Continuous MPS for Quantum Fields: an Energy Minimization Algoirthm). We work with a gradient descent method, and we are trying to minimuze the energy function:
$$ E = \frac{\bra{\Psi(A_L)}H\ket{\Psi(A_L)}}{\bra{\Psi(A_C}\ket{\Psi(A_C))}}$$
Minimizing this function is difficult since the left canonical tensor of our MPS is not necessarily linear - we use instead an auxiliary tensor as a replacement and subsitute it back in afterwards. To minimize, we try to find the gradient of the central tensor $A_C$ and change $A_L$ in that same direction. By then representing this as a tensor network:
![[Pasted image 20250718171830.png]]
We can contract our tensor network and get a local gradient. This is repeated until we minimize the norm (in theory, the norm of the gradient should vanish but in our case due to small errors it simply gets very very small)

We can extract things such as scaling dimensions and central charge - the conformal data that we desire, even though the infinite MPS of a critical system has an artificially finite correlation length that scales with bond dimension. 


Some other techniques that might be useful - preconditioning, look it up as i go along

We want to find this conformal data because conformal data usually tells us about the phasee that the system is in for some point in parameter space. When there is a jump, we know we are in a critical region! 
### Phase diagram of the interacting Majorana chain model

As before, we study a particular hamiltonian here:
$$ \mathcal{H} = \sum_j \left[ it\gamma_i \gamma_{i+1} + g \gamma_i \gamma_{i+1} \gamma_{i+2} \gamma_{i+3} \right] $$
This is the four majorana over four sites interaction - the relatively basic one that started it all. Compared to the other paper discussed, this one exhibits SUSY only on the scale of 250. 
### Topological and quantum critical properties of the interacting Majorana chain model

They looked at a similar hamiltonian of 

* Ask Natalia about one of the papers - why is there some kind of difference? Oh never mind she already brings it up in the paper hahahaha
Double checking 


## Second round of papers 