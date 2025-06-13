# Universal-Description

**"A Universal Description of Stochastic Oscillators "** [[PDF]](https://www.pnas.org/doi/reader/10.1073/pnas.2303222120)

at [[PNAS]](https://www.pnas.org/doi/abs/10.1073/pnas.2303222120)
by A. Pérez-Cervera, B.Gutkin, P.J.Thomas and B. Lindner

## Theoretical background

We assume our stochastic oscillator can be described by a Langevin equation

$$\frac{d\textbf{x}}{dt}=\textbf{f}(\textbf{x}) + \textbf{g}(\textbf{x})\xi(t), \qquad \textbf{x} \in \mathbb{R},$$

where $\textbf{f}$ is a is an n-dimensional vector, $\textbf{g}$ is an $n \times n$ matrix, and $\xi$
is $n$-dimensional white noise with uncorrelated components, so it satisfies, $\left\langle \xi_i(t)\xi_j(t') \right\rangle  = \delta(t-t')\delta_{i,j}$.

The process described by the previous Langevin equation is a $n$-dimensional Markov process, that is uniquely determined by the transition probability density $P(\textbf{x},t | \textbf{x}_0,s)$ (for $t>s$). This central statistics satisfies both the forward Kolmogorov (or ``Fokker-Planck'') equation 

$$\frac{\partial P}{\partial t}=\mathcal{L}[P]=-\nabla_\textbf{x}\cdot\left( \textbf{f}(\textbf{x}) P \right)+ \sum_{i,j}\frac{\partial^2}{\partial x_i x_j}\left(D_{ij}(\textbf{x})P\right), $$

$where ~ D=\frac12 gg^\intercal$. 

Here the functional $\mathcal{L}^\dagger$ acts with respect to the $\textbf{x}$ coordinates. $P(\textbf{x},t | \textbf{x}_0,s)$  also obeys the backward Kolmogorov equation 

$$-\frac{\partial P}{\partial s} = \mathcal{L}^\dagger [P]=\textbf{f}(\textbf{x})\nabla_\textbf{x_0} + D_{ij}$$

where the operator $\mathcal{L}^\dagger$ acts with respect to the $\textbf{x}_0$ coordinates.

```
with $T$ the period of the limit cycle and $\lambda_i$ its correspondent Floquet exponents.

To find such a parameterisation $x=K(\theta, \sigma)$ one has to solve the following invariance equation
```math
\frac{1}{T}\frac{\partial}{\partial\theta}K(\theta, \sigma) + \sum_{i=1}^{d-1} \lambda_i \sigma_i \frac{\partial}{\partial \sigma_i}K(\theta, \sigma) = X(K(\theta, \sigma))
```

The code we next present, solves the above mentioned invariance equation for the case $d=3$ and provides an expression for the parameterisation $x=K(\theta, \sigma)$ in the following form
$$K(\theta, \sigma) = \sum_{m=0}^{\infty} \sum_{\alpha=0}^{m} K_{\alpha, m-\alpha}(\theta) \sigma^{\alpha}_1 \sigma^{m-\alpha}_2$$

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
