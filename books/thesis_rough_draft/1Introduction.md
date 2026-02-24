This is just a nice rough draft so that I can write about my thesis without worrying about damaging my LaTeX document.

Introduction

Guiding questions:

- Hook
- Background and context of research (why is this research underway? What world context do you see this research embarking on?)
- Briefly discuss the topic of the research and what was done before. 
- Discuss the problem that you wish to remedy in this landscape and what you want to add
- Define research question and exact goals of thesis, and the approach that will be taken in understanding this

The study of one dimensional systems is one that has garnered a lot of attention in the past years among the theoretical condensed matter physics community. Due to the reduced dimensionality, interactions between particles takes on a noticeably different character. One of the basic theoretical models of interacting fermions, Fermi liquid theory, famously breaks down in 1-D and requires a radically different theory to be described, Luttinger liquid theory. Similarly, the method of bosonization is an interesting theoretical field theory that treats systems of interacting fermions as bosonic, and only describes one dimensional systems well. More specifically to this thesis, the discovery of Majorana zero modes in a simplified toy model of a 1-D nanowire, a phenomenon that opens the door for fault tolerant topological quantum computing architectures. These are all examples of the interesting physics associated with one-dimensional systems. 



The primary focus of this thesis will be in studying an interacting Majorana chain model, that seeks to extend the toy model introduced by Kitaev. The context in which this model will be studied is primarily in quantum phase transitions and the various quantum phase. Research of this variety has seen significant interest in the past (source all the papers of similar work), although with one glaring interaction not previously studied in any depth, a five-site, quartic Majorana interaction with an off center gap,

\begin{equation}

\gamma_i\gamma_{i+1}\gamma_{i+2}\gamma_{i+4} + \gamma_i\gamma_{i+2}\gamma_{i+3}\gamma_{i+4},

\end{equation}

for some interaction coefficient $k$. This $k$-term represents an adaption of the more commonly studied $f$-term, where the particle gap finds itself exactly in the middle. A diagram of all the different interactions can be found in Figure \ref{fig:interacting-terms}.



The exact goal of the thesis will be to study the $k$-term as in Figure \ref{fig:interacting-terms}. This is a gap in the literature of the interacting Majorana chain model and calls for further investigation, especially considering the wealth of critical phenomena discovered in similar investigations.



This thesis will be organized as follows. After this brief introduction to the topic, a more complete literature review will proceed covering the theoretical groundwork as well as the general state of research. The methodology section will explain in detail in the numerical methods applied, including an explanation and justification of all the simulation parameters for both the DMRG and ED computations. The results section will begin with a sketch of the phase diagram, followed by more detailed subsection describing each quantum phase and quantum phase transition. All novel discoveries, critical analysis and implications are expanded upon in the discussion section. Finally, the takeaways, limitations, and extensions of this investigation will be concluded and presented in the conclusion. All source code and more developed derivations are kept to the appendix.

