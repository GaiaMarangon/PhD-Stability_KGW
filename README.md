# Read me

This folder contains an implementation of the radial linearized perturbation problem for the Klein-Gordon-Wave system. It includes the following files:
1. `main_spectralSolver` 
2. `main_spectralElaborate_nFixed`
2. `main_spectralElaborate_nVary`

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
