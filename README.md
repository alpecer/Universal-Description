# Universal-Description

**"A Universal Description of Stochastic Oscillators "** [[PDF]](https://www.pnas.org/doi/reader/10.1073/pnas.2303222120)

at [[PNAS]](https://www.pnas.org/doi/abs/10.1073/pnas.2303222120)
by A. Pérez-Cervera, B.Gutkin, P.J.Thomas and B. Lindner

## Theoretical background

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
\mathcal{L}[P_\lambda]=\lambda P_\lambda,\qquad\mathcal{L}^\dagger [Q^*_{\lambda}]=\lambda Q^*_{\lambda},
```
where the smallest eigenvalue is $\lambda_0$, corresponds to the stationary state, which we denote $P_0(\textbf{x})$ and assume to be unique.




Our code is efficient as it relies on
- a diagonalisation of the linear system yielding the coeficcients $K_{\alpha, m-\alpha}(\theta)$ by means of the Floquet Normal Form
- automatic differentiation techniques

## How to use it?

Our code considers a 3D neuron model (see Eq. 59) but it is easy to adapt for any system of interest. To that end, open 'obtainKs.py' and
- Define its vectror field equations at "myVectorField" 
- Define the jacobian matrix at "myJacobianMatrix"
- introduce the period, $T$
- introduce an initial condition on the limit cycle, $x_0$
- write the automatic differentiation expression for your system (see Appendix B for more details)

There are extra parameters $b_1, b_2$ controlling the norm of $K_{10}$ and $K_{01}$ (see Remark 3.2 in the manuscript for more details)

## What produces?

By running the code you will generate three files "kx.dat", "ky.dat" and "kz.dat". Each file will have as many rows as monomials are generated and as many columns as points used to discretise $\theta$. Monomials are stored following the sum
$$\sum_{m=0}^{L} \sum_{\alpha=0}^{m} K_{\alpha, m-\alpha}(\theta)$$

Note the code prints the maxError of each monomial (the error is computed via equation 48 in the manuscript)
