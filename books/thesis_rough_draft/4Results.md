Results

Guiding questions:

- What was found? How did you find it?

General Results:

1. Phase Diagram
2. List of phases and phase transitions (written below)
3. 

**Ising transition**

The Ising transition for low $k$ remain

---

**G-1 to Floating-1transition:**

The first transition that was studied in depth was the AFM first gapped phase transition and the Floating transigion for h < 1. 

-------------

**Incommensurate Luttinger Liquid:**

The incommensurate Luttinger Liquid phase was identified using a multitude of methods. Firstly, and the most straightforward, we use a well known formula from CFT to find the fit of the central charge from the entanglement entropy of the system. The Calabrese Cardy formula \cite{calabrese cardy} predicts the following relation between the entanglement entropy of a 1+1D system in a critical phase as 

\begin{equation}
S_A(x) = \frac{c}{6} \log{(\frac{2L}{\pi a} \sin{\frac{\pi x}{L})} }+ \tilde{c}_1'.

\label{eq: calabrese-cardy}

\end{equation}

We can define a conformal length $ l_C = \log{(\frac{2L}{\pi a} \sin{\frac{\pi x}{L})} }$ that turns equation\ref{eq:calabrese-cardy} into a linear fit. Examples of the fit are given in Figure \ref{fig:fit_c_ll}. We expect the following values for the central charge based on the CFT that describes it:

\begin{table}[H]
    \centering
    \begin{tabular}{|c|c|}
    \hline
    CFT & $c$ \\
    \hline
    TFIM & 0 \\
    \hline
    TFIM, critical & 1/2 \\
    \hline
    Luttinger Liquid & 1 \\
    \hline
    Luttinger Liquid + TFIM, critical & 3/2 \\
    \hline
    \end{tabular}
    \caption{Caption}
    \label{tab:placeholder}
\end{table}

The incommensurability of the luttinger liquid phase is identified via identification of the Friedel oscillations in the local magnetization $\sigma^X_i$. riedel oscillations are the effect of boundary perturbations to the order parameter in the bulk, leading to a modulation of the local magnetization. In our example,the perturbations are simply a finite size effect. Using boundary CFT (bCFT) \cite{natalia topological paper} , we expe24:18ct the modulation to appear as

\begin{equation}

$$\langle \sigma_i \rangle \propto \frac{\cos{(qi + \phi)}}{\left[(N/\pi) \cos{(\pi i/N)}\right]^K}$$,

\label{eq:friedel_oscillations}

\end{equation}

where $K$ is the Luttinger liquid parameter and $q$ is the real-space wavevector. For commensurate systems, $q$ is a multiple of $\pi$; for incommensurate $q \neq $n\pi$. 

----------------

**Self Duality: Transition at h=1**

The existence of self duality in our Hamiltonian provides an exciting opportunity to study the
transition along the self-dual line. We note that we identify two primary transitions - the first
being the expected Ising transition with a central charge of c = 1/2 along the gapped phases of
low k. This transition occurs approximately between 0 <  k < 0.15 when h = 1. Secondly, the transition
between two Luttinger Liquid phases across k > 0.15. 

The transition between the critical .

Since the disappearance of the Ising criticality for the floating phases is quite different from what was expected with the $f$-term, we also studied the non-symmetric case where $k_L \neq k_R$. This was done in the following way, by defining the ratio $r_k = \frac{k_L}{k_R}$ and requiring $k_L + k_R = 2k$, we obtain:

\begin{equation}

k_R = \frac{2k}{r_k+1}, k_L = \frac{2k * r_k}{1+r_k}

\end{equation}

The idea is that maybe the lack of a transition at h = 1 is an exception in the symmetric case of $r_k = 1$. Otherwise, we consider

---

**ED results**

Although the primary tool for investigating low energy physics is DMRG, ED still has various uses for the investigation of the model. 

The main reason i wanted to use ED is to fill in the knowledge gaps that DMRG left behind. Namely, the information of states above the ground state. ED makes it easy to examine the energy spectrum as well as isolate specific terms in relatively quick simulations for small system sizes. The energy spectrum was used to determine the natural degeneracy of the ground state, and what polarized boundary conditions would be needed.



Additionally, the $k$-term was studied in isolation to determine when and where it would be used.

---

**Competing Interactions** 

As an additional extension to the primary investigation of this project, that being the phase diagram for the $h-k$ model, we studied the $k-f$ model
