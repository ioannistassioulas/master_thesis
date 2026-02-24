Methodology



Guiding questions:

- What sort of methods are being employed? theoretical, numerical, etc

This thesis will investigate the research question using primarily numerical methods. 

- Why were these methods chosen? what are their advantages?
- Give rough breakdown on the exact specifics of the DMRG algorithm i.e. what sort of L values did you take? did you increase chi with each sweep etc. and why was this done. Discuss boundary conditions and which ones were chosen, and why (reference lit review explaining why this is done)

DMRG was used as the primary numerical method in this thesis. The advantages and limitations were already covered in depth in the literature review, and will not be repeated here. Two different implementations were used throughout the investigation of the thesis. The first one involves a source code written by Pietro Richelli, available on his github page \cite{Pietro's github page}. The second DMRG simulation used was built on the TeNPy library \cite{Tenpy}, a Python library for general tensor network methods. In both cases, the Hamiltonian of my Model had to be derived and constructed. The source code for both methods can be found in the appendix. The primary difference between both implementation of DMRG lies in how the tensor network contractions are saved. In Pietro's code, the MPS and the contractions were all written to disk and fetched as they are required in each step of the sweeps, while TeNPy keeps all the tensors in memory for speed. This decision between storing the MPS and tensor contractions on disk storage vs volatile memory depends on the nature of the system. Disk storage allows for larger systems, both in terms of number of sites as well as bond dimension, to be studied without the fear of out-of-memory errors. The price to pay for this is the computational slowdown from constantly reading and writing to storage. Keeping the MPS and the contractions in memory, on the other hand, removes this bottleneck and speeds up the simulation time, albeit only applicable for mid-memory size systems. Memory requirements of the DMRG algorithm are as follows:

\begin{table}[H]

MPS : L\chi^2d : Full dimensions of MPS tensor network

Left/Right contractions: l L \chi^2 \chi_{MPO} : Virtual bonds of <psi|H|psi>

\end{table}

Therefore, general simulations with TeNPy were run with a maximum system size of $L\approx 300$ and a maximum MPS bond dimension of $\chi < 300$, in order to keep memory usage below the limit of 12GB, which was the maximum amount of memory available for each DMRG run in the DelftBlue cluster. Only when the precision required to understand the physics of the system demanded a larger system size or bond dimension was the disk-storage method chosen.

Along

- Give an explanation of why ED was chosen, and for what value it served over DMRG

The second numerical method employed in the investigation of the $k$-term interacting Majorana chain model was exact diagonalization. As the name suggests, exact diagonalization solves the 
