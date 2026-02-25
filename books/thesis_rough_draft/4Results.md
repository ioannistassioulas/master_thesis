Results

Guiding questions:

- What was found? How did you find it?

General Results:

1. Phase Diagram
2. List of phases and phase transitions

First order Transition:



Incommensurate Luttinger Liquid:

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

The incommensurability of the luttinger liquid phase is identified via identification of the Friedel oscillations in the local magnetization $\sigma^X_i$. Named after Jacques Friedel who predicted them in 1952 \cite{} Friedel oscillations are the effect of boundary perturbations to the order parameter in the bulk, leading to a modulation of the local magnetization.Using boundary CFT (bCFT) \cite{natalia topological paper} , we expect the modulation to appear as

\begin{equation}

$$\langle \sigma_i \rangle \propto \frac{\cos{(qi + \phi)}}{\left[(N/\pi) \cos{(\pi i/N)}\right]^K}$$,

\label{eq:friedel_oscillations}

\end{equation}

where $K$ is the Luttinger liquid parameter and $q$ is the real-space wavevector. For commensurate systems, $q$ is a multiple of $\pi$; for incommensurate $q \neq $n\pi$. 



