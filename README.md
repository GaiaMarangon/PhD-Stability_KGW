# Read me

This folder contains a `MATLAB` implementation of the radial linearized perturbation problem for the Klein-Gordon-Wave system. It includes the following `.m` files:
1. `main_spectralSolver`, solving the spectral problem 
2. `main_spectralElaborate_nFixed` performing some analysis at fixed $n$ (excitation index of the eigenstate serving as a background for the perturbations)
2. `main_spectralElaborate_nVary` performing some analysis across different values of $n$ (excitation index of the eigenstate serving as a background for the perturbations)

## Problem
The Klein-Gordon-Wave problem reads:
```math
\begin{align}
\begin{cases}
     \square \, u  = (\mu^2 + 2\phi) \,u \\
   \square \, \phi  =  u ^2 
\end{cases}\,.
\end{align}
```

Consider small variations of the $n$-th stationary solution $(\mu^2_n,u_n,\phi_n)$:
```math
\begin{align*}
\begin{cases}
    u(t,x) \equiv u_n(x) + \delta u(t,x) \\
    \phi(t,x) \equiv \phi_n(x) + \delta \phi(t,x) \\
\end{cases}
\end{align*}
```

Consider perturbations in the harmonic form:
```math
\begin{align}
    \delta u (t,x)\equiv U(x) e^{i \omega t}\,, \qquad
    \delta \phi (t,x)\equiv \Phi(x) e^{i \omega t}\,,
\end{align}
```

The linearized perturbation problem reads:
```math
\begin{align}
    \begin{bmatrix}
        -\triangle+\mu^2_n +2\phi_n & 2u_n\\
        2 u_n & -\triangle \\
    \end{bmatrix}
    \begin{pmatrix}
        U\\
        \Phi\\
    \end{pmatrix}
    = \omega^2
    \begin{pmatrix}
        U\\
        \Phi\\
    \end{pmatrix}
    \,.
\end{align}
```

In the radial case, consider the convenient change of variables:
```math
\begin{align*}
    V(r) \equiv r U(r), \qquad W(r) \equiv r \Phi(r)\,.
\end{align*}
```

The radial part of the linear perturbation problem reads:
```math
\begin{align}
    \begin{bmatrix}
        -\partial^2_r+\mu^2_n +2\phi_n & 2u_n\\
        2 u_n & -\partial^2_r \\
    \end{bmatrix}
    \begin{pmatrix}
        V\\
        W\\
    \end{pmatrix}
    = \omega^2
    \begin{pmatrix}
        V\\
        W\\
    \end{pmatrix}
    \,.
\end{align}
```

We solve it with boundary conditions appropriate for unstable modes:
```math
\begin{align}
    \begin{cases}
        V(0) = W(0) = 0 \\
        V(\infty) = W(\infty) = 0
    \end{cases}
\end{align}
```

##  The `main_spectralSolver` file
This file solves the spectral problem (4) with BCs (5) for perturbations around the $n$-th eigenstate (called background). 

The problem is solved for different refinements of the radial domain ($N$ is the number of discretization points in the domain), and the results are saved to `.mat` files.  
Each result file refers to a single background excitation index $n$ and includes the results for different numbers $N$ of discretization points.
When the code is executed, for each $N$ it first checks the result file, and the problem is solved only for new combination $(n,N)$.

The solving procedure consist in defining a uniform discretization of the domain, computing the spectral matrix accordingly and applying the boundary conditions (Dirichlet BCS are handled by removing first and last rows/columns of each block). Then, the discretized problem is solved with the `eig` built-in function, and results are sorted by growing eigenvalues. Spurious modes are optionally filtered (excluding largest eigenvalues) and the eigenvectors are normalized ($L^2$-normalization), taking care of the BCs (inserting missing first/last zero components). Finally, some analysis is performed on stability results: eigenvalues are classified as stable (positive), unstable (negative) or zero (vanishing, within a given tolerance), and the index and total number of each class of eigenvalues (stable, unstable, zero) is saved. 



## The `main_spectralElaborate_nFixed` file

This file performs some analysis of the results for fixed background excitation index $n$. 

The spectrum of the radial linear perturbation problem includes $n+1$ negative eigenvalues $\{\omega_k^2\}_{k=0}^n$ for each $n$-th background stationary solution, conventionally sorted in descending order.
The analysis includes the following plots:
- plot of the eigenvalues $\omega_k^2$ against their index $k$;
- plot of the spectrum for different values of $N$; 
- plot of the relative differences of the eigenvalues across different values of $N$: $\Delta_k (N_j) = \frac{\omega_k^2(N_j)-\omega_k^2(N_{j-1})}{\omega_k^2(N_0)}$;
- plot of the eigenmodes $V(r)$. They are optionally compared to the corresponding eigenfunctions $f_n(r)$ of the original stationary Klein-Gordon-Wave problem having the same number of nodes.
- plot of the eigenmodes $W(r)$. 

Additionally, the stability results are printed to file.

## The `main_spectralElaborate_nVary` file

This file performs some analysis of the results for fixed background excitation index $n$.

It includes the following plots:
- plot of the eigenvalues $\omega_k^2$ against their index $k$, across different values of $n$.