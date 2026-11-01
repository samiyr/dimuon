# Decay grid generator

This tool can be used to generate the decay grids necessary for running the dimuon code. The decay grids are tabulations of the integral
```math
\begin{equation}
\frac{\pi \tau}{m_h}\int \mathrm{d}\rho \, \rho E_h^2 \int \mathrm{d}(\cos\theta) \, d(w)\Theta(\rho E_h-E^{\text{min}}),
\end{equation}
```
where
- $\tau$ is the lifetime of the resonance
- $E_h$ is the energy of the hadron
- $\Theta$ is the Heaviside theta function
- $E^{\text{min}}$ is the energy cut
- $w=\rho E_h (E_h-\sqrt{E_h^2-m_h^2})/m_h^2$

and the decay function is parametrized as
```math
\begin{equation}
d(w)=Nw^\alpha(1-\gamma w)^\beta \Theta(0\leq w \leq 1/\gamma),
\end{equation}
```
where $N$, $\alpha$, $\beta$, and $\gamma$ are fitted parameters. More information can be found in e.g. [[2405.12677]](https://doi.org/10.1007/JHEP09(2024)043).

The code calculates eq. (1) for all values of $E_h$ in a grid. The grid uses Chebyshev nodes and is done in log-space (i.e. the eventual polynomial interpolation used in `libdimuon` is based on $\log E_h$). The parameter values are read from a file, which must contain at least one row of values in the format
```
N,alpha,beta,gamma
```

This input file can of course contain multiple rows, each line following the same format.

The calculations are parallelized using `OpenMP`.

The executable takes the following options:
- `minE`: $E^{\text{min}}$ (required)
- `input`: input file, i.e. the file containing parameter values (required)
- `output`: output directory (default: current working directory)
- `bin_min`: lower $E_h$ grid edge, in $\mathrm{GeV}$ (default: `minE`)
- `bin_max`: upper $E_h$ grid edge, in $\mathrm{GeV}$ (default: 300)
- `bin_count`: number of Chebyshev nodes (default: 30; values above 40 can cause numerical problems in the GSL implementation of polynomial interpolation used in `libdimuon`)

Note that the output files follow the naming scheme defined in `dm::decay_grid::filename`.