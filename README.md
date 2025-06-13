# Universal-Description

**"A Universal Description of Stochastic Oscillators "** [[PDF]](https://www.pnas.org/doi/reader/10.1073/pnas.2303222120)

at [[PNAS]](https://www.pnas.org/doi/abs/10.1073/pnas.2303222120)
by A. Pérez-Cervera, B.Gutkin, P.J.Thomas and B. Lindner

## Theoretical background

Next, we briefly summarize few details of the theoretical background. For the complete information, please, check our paper and its [[Supplemental Information]](https://www.pnas.org/doi/suppl/10.1073/pnas.2303222120/suppl_file/pnas.2303222120.sapp.pdf)

We assume our stochastic oscillator can be described by a Langevin equation

$$\frac{d\textbf{x}}{dt}=\textbf{f}(\textbf{x}) + \textbf{g}(\textbf{x})\xi(t), \qquad \textbf{x} \in \mathbb{R},$$

where $\textbf{f}$ is a is an n-dimensional vector, $\textbf{g}$ is an $n \times n$ matrix, and $\xi$
is $n$-dimensional white noise with uncorrelated components, so it satisfies, $\left\langle \xi_i(t)\xi_j(t') \right\rangle  = \delta(t-t')\delta_{i,j}$.

The process described by the previous Langevin equation is a $n$-dimensional Markov process, that is uniquely determined by the transition probability density $P(\textbf{x},t | \textbf{y},s)$ (for $t>s$). This central statistics satisfies both the forward Kolmogorov (or ``Fokker-Planck'') equation 

$$\frac{\partial P}{\partial t}=\mathcal{L}[P]=-\nabla_\textbf{x}\cdot\left( \textbf{f}(\textbf{x}) P \right)+ \sum_{i,j}\frac{\partial^2}{\partial x_i x_j}\left(D_{ij}(\textbf{x})P\right), $$

$where ~ D=\frac12 gg^\intercal$. 

Here the functional $\mathcal{L}^\dagger$ acts with respect to the $\textbf{x}$ coordinates. $P(\textbf{x},t | \textbf{y},s)$  also obeys the backward Kolmogorov equation 

```math
-\frac{\partial P}{\partial s} = \mathcal{L}^\dagger [P]=\textbf{f}(\textbf{y})\nabla_\textbf{y}\left(P \right)+\sum_{i,j}D_{ij}(\textbf{y})\frac{\partial^2 P}{\partial y_i y_j} 
```
where the operator $\mathcal{L}^\dagger$ acts with respect to the $\textbf{y}$ coordinates.

Here we assume a discrete set of eigenvalues $\lambda=\mu_\lambda+i\omega_\lambda$ with corresponding forward ($\mathcal{L}$) and backwards ($\mathcal{L}^\dagger$) eigenfunctions
```math
\mathcal{L}[P_\lambda]=\lambda P_\lambda,\qquad\mathcal{L}^\dagger [Q^*_{\lambda}]=\lambda Q^*_{\lambda}.
```



## How to use the code?

The attached code is ready to generate the spectral decomposition of the "noisy Stuart-Landau" model in Eq. (11) in the main text with parameters as in Fig. 1 in the main text. 

Next, we explain how to modify it in order to obtain the discretisation of the Kolmogorov operators associated to a different planar SDE:

Whereas the library "kolmogorovTools.py" contains the functions computing the terms composing $\mathcal{L}^\dagger$, the script "diagonalizeKolmogorov.py" builds the operator for the model of interest and provides the set of eigenfunctions and eigenvalues for $\mathcal{L}$ and $\mathcal{L}^\dagger$. That is why, we just need to modify "diagonalizeKolmogorov.py".
- Make sure both functions are at the same directory  
- Open "diagonalizeKolmogorov.py", edit the functions $f_x (x,y)$ and $f_y (x,y)$ and its associated parameters. Then, choose the domain you aim to discretise and how many points you want in each direction (params N,M)
- Then, execute it and enjoy your eigenfunctions :)
- As tip: chose a domain large enough so the probability of trajectory reaching the boundaries is low.
- The matrices in the code are written following sparse SciPy libraries. This way one can consider big number of points to discretise without running out of memory.
- How many points to use when discrtising? Well, try and error will teach you to find a good compromise (whereas few points may generate weird results, choosing too many points may cause the code to be too slow)
- We also encourage the reader to read Section 5 in the SI, where complementary details about the procedure can be found.

## What produces?

By running the code you will generate the files containing the eigenvectors and eigenvalues of the constructed matricial representation of the forward and backward operators. The eigenvectors you obtain correspond to the eigenfunctions evaluated on the chosen grid. As the eigenvalues and eigenfunctions are computed via the SciPy routine "eigs", we refer to [[scipy.eigs documentation]](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html) for information about how they are stored.

