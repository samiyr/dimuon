# Numerical evaluation

## Structure of the cross section

Let's start by writing down the cross section we want to calculate (while dropping irrelevant details):
```math
\begin{equation}
\sigma_{\mu\mu}=\sum_h \overbrace{\int\mathrm{d}z \, \sigma_{h\mu}(\mathcal K, z) B_h(\mathcal K, z)}^{\equiv \sigma_{\mu\mu}^h},
\end{equation}
```
where 
- $\sigma_{\mu\mu}$ is the dimuon cross section
- $\sigma_{h\mu}$ is the hadron production cross section 
- the sum over $h$ runs over all intermediate hadrons
- $\mathcal K$ represents a kinematical point $(x, y, E)$
- $B_h(\mathcal K, z)$ is the energy-dependent branching ratio

The hadron-production cross section can be written as a linear combination of structure functions as
```math
\begin{equation}
\sigma_{h\mu}(\mathcal K, z)=\sum_{k=1}^3 F_k^h(\mathcal K, z)\omega_k(\mathcal K)
\end{equation}
```
where
- $F_k^h(\mathcal K, z)=F_k^h(x, z, Q^2)$ is a structure function
- $\omega_k(\mathcal K)$ is some simple function of the kinematical variables

The structure functions themselves are written as convolutions of coefficient functions:
```math
\begin{equation}
F_k^h(\mathcal K, z)=\sum_{\alpha\beta} \int \frac{\mathrm{d}\xi}{\xi} \, \frac{\mathrm{d}\zeta}{\zeta} \, C_k^{\beta\alpha}(\xi, \zeta)f_\alpha(\chi/\xi)D^h_\beta(z/\zeta)
\end{equation}
```
where
- $C_k^{\beta\alpha}$ is a coefficient function of the channel $\alpha\to \beta$
	- expanded in a perturbation series: $C_k^{\beta\alpha}=\sum_p (\frac{\alpha_s}{2\pi})^p C_k^{(p),\beta\alpha}$
	- each coefficient $C_k^{\beta\alpha}$ has a corresponding charge conjugate $C_k^{\bar\alpha \bar\beta}$, which is given by $C_k^{\beta\alpha}$ up to a possible change of sign [[2506.19926](https://inspirehep.net/literature/2938806): §3.3.1, §3.3.2]
		- `sidis::partonic_channel::conjugate_type`
		- `coupling_sign` in `sidis::evaluator::channel_factor`
	- $C_k^{\beta\alpha}$ is flavor-independent, i.e. it only matters whether $\alpha$ and $\beta$ are quarks, antiquarks, or gluons
		- this is true only for the formula of the coefficient, but the actual numerical output depends on final-state masses via e.g. the slow-rescaling variable
		- at NNLO this gets a bit more complicated, since a $q\to q$ coefficient is actually one of multiple coefficients (which at lower orders don't exist) depending on the flavors, but each of these coefficients is then independent of the flavors
	- the coefficient of a physical channel $\alpha\to\beta$ is built up of multiple contributions
		- the analytical expressions for the coefficients match the ones given in [[2506.19926](https://inspirehep.net/literature/2938806): §3.3.1, §3.3.2], wherein multiple terms $C_k^{\beta\alpha,s}$ combine to give a total contribution $C_k^{\beta\alpha}=\sum_s |V_{q_s q'_s}|^2 C_k^{\beta\alpha,s}$, where $V$ is the CKM matrix, and $q_s$ and $q'_s$ are the quarks connected to the $W$ boson in the subcontribution $s$ (`evaluator::channel_factor`)
		- we only care about the physical channels, so no need to track each subcontribution separately $\to$ `sidis::physical_channel`
- $f_\alpha$ is the PDF of parton $\alpha$ (factorization scale left implicit)
- $D_\beta^h$ is the FF for the transition $\beta\to h$ (fragmentation scale left implicit)
- $\chi$ is the slow-rescaling variable (note that the final state and therefore $\chi$ depends on the subcontribution $s$ of a channel, which is hidden away in the index $s$ in eq. (8))

For each $\alpha\to\beta$, the structure function contribution reads
```math
\begin{equation}
F_k^{h,\beta\alpha}(\mathcal K, z)=\int \frac{\mathrm{d}\xi}{\xi} \, \frac{\mathrm{d}\zeta}{\zeta} \, C_k^{\beta\alpha}(\xi, \zeta)f_\alpha(\chi/\xi)D^h_\beta(z/\zeta)
\end{equation}
```
so that $F_k^h=\sum_{\alpha\beta} F_k^{h,\beta\alpha}$. The code doesn't actually compute $F_k^{h,\beta\alpha}$, but rather
```math
\begin{equation}
\chi z F_k^{h,\beta\alpha}(\mathcal K, z)=\int \mathrm{d}\xi \, \mathrm{d}\zeta \, C_k^{\beta\alpha} \frac{\chi}{\xi} f_\alpha(\chi/\xi) \frac{z}{\zeta} D_\beta^h(z/\zeta).
\end{equation}
```
Thus, "$f_\alpha$" and "$D_\beta^h$" in this document (and in the code) actually refer to the "$xf(x)$" and "$zD(z)$" values. The extra factors of $\chi$ and $z$ are divided out at the end for $F_1$ and $F_3$. For $F_2$, only $z$ is divided out, since $F_2$ naturally appears with the $\chi$ factor. **This normalization convention, which differs from the usual one for $F_k$, is used throughout this document and the code.**

## The structure of coefficient functions

The coefficient function $C_k^{\beta\alpha}$ is built up of products of
- regular functions in $\xi$ and $\zeta$
- the delta functions $\delta(1-\xi)$ and $\delta(1-\zeta)$
- the plus distributions $\left[\frac{\log^n(1-\xi)}{1-\xi}\right]_+$ and $\left[\frac{\log^m(1-\zeta)}{1-\zeta}\right]_+$ for some values $n$ and $m$

The regular functions don't change the arguments of the PDF/FF, but the distributions do, yielding terms proportional to $f_\alpha(\chi)$ and $D_\beta^h(z)$. In addition, the quark-initiated SACOT coefficients are written in terms of $\zeta'$, which depends on $\xi$, $\zeta$, $Q^2$, and the mass of the fragmenting parton (and thus on $\beta$ of the channel $\alpha\to \beta$), and use $D_\beta^h(z/\zeta')$ instead of $D_\beta^h(z/\zeta)$.

Thus, after convolving the distributions, $C_k^{\beta\alpha,s}$ is bilinear in the PDF/FF in the sense that
```math
\begin{equation}
f_\alpha(\chi/\xi) D_\beta^h(z/\zeta) C_k^{\beta\alpha,s}(\xi, \zeta)=\sum_{\mathfrak a=0}^1 \sum_{\mathfrak b=0}^2 f_{\alpha, \mathfrak a} D_{\beta, \mathfrak b}^h C_{k, \mathfrak a\mathfrak b}^{\beta\alpha,s}(\xi, \zeta, \chi, z)
\end{equation}
```
where
```math
\begin{equation}
\vec f_\alpha=\begin{bmatrix} f_\alpha(\chi) \\ f_\alpha(\chi/\xi) \end{bmatrix} \quad\text{and}\quad \vec D_\beta^h=\begin{bmatrix} D_\beta^h(z) \\ D_\beta^h(z/\zeta) \\ D_\beta^h(z/\zeta') \end{bmatrix}.
\end{equation}
```
$C_{k,\mathfrak a\mathfrak b}^{\beta\alpha,s}$ picks up dependence on $\chi$ and $z$ from e.g. the remainder terms of the plus distributions.

Putting everything together,
```math
\begin{equation}
\sigma_{\mu\mu}=\int \mathrm{d}z \, \mathrm{d}\xi \, \mathrm{d}\zeta \sum_k \sum_{\alpha\beta} \sum_{\mathfrak a \mathfrak b} \sum_p \left(\frac{\alpha_s}{2\pi}\right)^p \sum_s f_{\alpha, \mathfrak a}^s \left(\sum_h B_h(\mathcal K, z) D_{\beta, \mathfrak b}^h\right) C_{k, \mathfrak a\mathfrak b}^{(p),\beta\alpha,s}(\xi, \zeta, \chi, z) \omega_k(\mathcal K).
\end{equation}
```
As a reminder, the indices are
- $h$: hadron
- $k$: structure function ($k=1,2,3$ for $F_k$)
- $\alpha\beta$: the physical channel (this is regarded as a single index)
- $\mathfrak{a}$: component of $\vec f_\alpha$, which contains the two possible PDF values at a given integration point
- $\mathfrak{b}$: component of $\vec D_\beta^h$, which contains the three possible FF values at a given integration point
- $p$: perturbative order
- $s$: subcontribution, from which the (NNLO) coefficients for a given physical channel are built

## Evaluating the coefficients

Numerically, one cannot use the coefficient functions directly, but instead work with the convolved version where the distributions have been convolved against the PDF/FF under the integral sign, e.g.
```math
\begin{equation}
\delta(1-\xi)\frac{1}{[1-\zeta]_+}f_\alpha(\chi/\xi) D_\beta^h (z/\zeta)\to \frac{f_\alpha(\chi)}{1-\chi}\frac{D_\beta^h(z/\zeta)-D_\beta^h(z)}{1-\zeta}+\frac{f_\alpha(\chi)}{1-\chi}\frac{D_\beta^h(z)}{1-z}\log(1-z)
\end{equation}
```

Additionally, the coefficients are written in such a way that the decomposition in terms of $\mathfrak a$ and $\mathfrak b$ (i.e. which terms are proportional to which one of the $2\times 3=6$ combinations of PDF/FF arguments) is not directly available. These could in principle be determined analytically, but one can also extract $C_{k, \mathfrak a'\mathfrak b'}^{(p),\beta\alpha,s}$ numerically from $C_k^{(p),\beta\alpha,s}$ by setting $f_{\alpha,\mathfrak a}=\delta_{\mathfrak a, \mathfrak a'}$ and $D^h_{\beta,\mathfrak b}=\delta_{\mathfrak b,\mathfrak b'}$. 

## Vectorization and masking

The calculation of the coefficients can be vectorized, i.e. instead of computing each possible product of the PDF and FF values, these can be computed at once as components of a single vector. 

First, the 2D index structure $\mathfrak a, \mathfrak b$ has to be unraveled into a 1D index structure $[\mathfrak c]$, where
- $[\mathfrak c]=3\mathfrak a+\mathfrak b$
- $\mathfrak a=\lfloor [\mathfrak c]/3\rfloor$
- $\mathfrak b=[\mathfrak c] \ \mathrm{mod} \ 3$

What we want then is essentially the Kronecker product $\otimes$ (`vectorization::kronecker_product`), since
```math
\begin{equation}
\vec f_\alpha\otimes \vec D_\beta^h=
\begin{bmatrix}
f_\alpha(\chi) D_\beta^h(z) \\
f_\alpha(\chi) D_\beta^h(z/\zeta) \\
f_\alpha(\chi) D_\beta^h(z/\zeta') \\
f_\alpha(\chi/\xi) D_\beta^h(z) \\
f_\alpha(\chi/\xi) D_\beta^h(z/\zeta) \\
f_\alpha(\chi/\xi) D_\beta^h(z/\zeta')
\end{bmatrix}
\end{equation}
```
This gives us all $2\times 3=6$ combinations that we need for the coefficients. The coefficient can now be regarded as a 6D vector $\mathbf C_k^{(p),\beta\alpha,s}$ with entries given by $C_{k,[\mathfrak c]}^{(p),\beta\alpha,s}$ (or $2\times 3$ matrix with entries $C_{k,\mathfrak a \mathfrak b}^{(p),\beta\alpha,s}$), so that 
```math
\begin{equation}
f_\alpha(\chi/\xi) D_\beta^h(z/\zeta) C_k^{(p),\beta\alpha,s}(\xi, \zeta)=\mathbf 1_6\cdot (\vec f_\alpha \otimes \vec D_\beta^h)*\mathbf C_{k}^{(p),\beta\alpha,s},
\end{equation}
```
where $*$ is the elementwise (Hadamard) product. Note that $\mathbf 1_n \cdot \mathbf v = \mathbf 1_n \cdot (v_1, \ldots, v_n)=\sum_{i=1}^n v_i$, where $\mathbf 1_n=(1, 1, \ldots, 1)$ is the identity vector with respect to the Hadamard product.

For the moment, let's drop the extra indices from eq. (11) and rewrite it in a more familiar way as
```math
\begin{equation}
\mathbf 1_6\cdot (\vec f \otimes \vec D)*\mathbf C=(\vec f \otimes \vec D)\cdot \mathbf C.
\end{equation}
```
Using the [mixed-product property](https://en.wikipedia.org/wiki/Kronecker_product#Relations_to_other_matrix_operations) $(A\otimes B)(C\otimes D)=(AC)\otimes(BD)$, the above can be further written as
```math
\begin{equation}
\begin{aligned}
(\vec f \otimes \vec D)\cdot \mathbf C&=(\vec f \otimes \vec D)^T\mathbf C \\
&=(\vec f^T \otimes \vec D^T)\mathbf C \\
&=(I_1 \otimes \vec D^T)(\vec f^T \otimes I_3)\mathbf C \\
&=\vec D^T(\vec f^T\otimes I_3)\mathbf C \\
&=(\vec f^T\otimes I_3)\mathbf C \cdot \vec D,
\end{aligned}
\end{equation}
```
where $I_n$ is the $n\times n$ identity matrix. The advantage of eq. (13) over eq. (11) is that the dependence on $\vec f_\alpha$ and $\vec D_\beta^h$ is separated: the PDF dependence (i.e. $(\vec f^T\otimes I_3)\mathbf C$) can be calculated separately from $\vec D$. This is what the code actually does. First, it calculates $\mathbf C$. Then, for every PDF row, the code computes $(\vec f^T\otimes I_3)\mathbf C$. Finally, for every FF row, it computes the full dot product. PDF and FF rows will be explained later as they depend on the run mode.

For a practical implementation, we need the components of the vector $(\vec f^T\otimes I_3)\mathbf C$. Consider the general case $(\vec a^T\otimes I_d)\mathbf b$, where $\vec a$ and $\mathbf b$ are $d_1$- and $d_2$-dimensional column vectors, and $d=d_2/d_1$ (with the assumption that $d_1|d_2$). Writing for convenience $i=dk+l$ and $j=dm+n$, the components are given by
```math
\begin{equation}
\begin{aligned}
\left[(\vec a^T\otimes I_d)\mathbf b\right]_i &=(\vec a^T\otimes I_d)_{ij} b_j \\
&=\left[\vec a^T\right]_{km} \delta_{ln}b_{dm+n} \\
&=a_m b_{dm+i},
\end{aligned}
\end{equation}
```
where $k=0$ (and thus $i=l$) since $\vec a^T$ is a row vector. The index structure used above can be found on e.g. [Wikipedia](https://en.wikipedia.org/wiki/Kronecker_product#Definition).

The actual coefficients, however, are not written in vector form. Rather, they look like
```math
\begin{equation}
f_\alpha(\chi)D_\beta^h(z)(\cdots)+f_\alpha(\chi)D_\beta^h(z/\zeta)(\cdots)+\cdots.
\end{equation}
```
 If one replaces the originally-scalar PDF/FF values (which are $f_\alpha(\chi)$, $f_\alpha(\chi/\xi)$, $D_\beta^h(z)$, $D_\beta^h(z/\zeta)$, $D_\beta^h(z/\zeta')$) in the actual coefficients with vectors, one gets the Hadamard product, not the outer product. The next step then is to figure out the proper masking vectors to extract the elements of the coefficient $\mathbf C_k^{(p),\beta\alpha,s}$.

 *Side note: the reason why the Kronecker product cannot easily be used is that it would require the coefficient implementations to distinguish between $\text{scalar}\times \text{vector}$ and $\text{vector}\times \text{vector}$ products.*

The Kronecker and Hadamard products are connected via another [mixed-product property](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)#The_mixed-product_property):
```math
\begin{equation}
(\vec A\otimes \vec B)*(\vec C\otimes \vec D)=(\vec A*\vec C)\otimes (\vec B*\vec D).
\end{equation}
```
From this it follows that
```math
\begin{equation}
\vec A\otimes \vec B=(\vec A\otimes \vec 1_{n_B})*(\vec 1_{n_A}\otimes \vec B),
\end{equation}
```
where the right-hand side is the form we want: it's a Hadamard product of two vectors (i.e. the elementwise product of the vectorized PDF/FF values in the coefficients). Here $n_A$ and $n_B$ are the dimensions of $\vec A$ and $\vec B$.

Therefore the masking vectors should be constructed as follows: starting from the definitions
```math
\begin{equation}
\vec v_0=\begin{bmatrix} 1 \\ 0 \end{bmatrix} \quad\quad \vec v_1=\begin{bmatrix} 0 \\ 1 \end{bmatrix}
\end{equation}
```
and
```math
\begin{equation}
\vec u_0=\begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} \quad\quad \vec u_1=\begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} \quad\quad \vec u_2=\begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix},
\end{equation}
```
the masking vectors are given by $\vec v_{\mathfrak a}\otimes \vec 1_3$ and $\vec 1_2 \otimes \vec u_{\mathfrak b}$ (`sidis::vectorization::make_masking_vector`). For example, if one sets
```math
\begin{equation}
\mathbf f_\alpha(\chi)=\vec v_0\otimes \vec 1_3=\begin{bmatrix} 1 \\ 1 \\ 1 \\ 0 \\ 0 \\ 0 \end{bmatrix} \quad\text{and}\quad \mathbf D_\beta^h(z)=\vec 1_2 \otimes \vec u_0=\begin{bmatrix} 1 \\ 0 \\ 0 \\ 1 \\ 0 \\ 0 \end{bmatrix},
\end{equation}
```
then
```math
\begin{equation}
\mathbf f_\alpha(\chi) * \mathbf D_\beta^h(z)=\begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{bmatrix}.
\end{equation}
```
Looking back at $\vec f_\alpha \otimes \vec D_\beta^h$, the above vector corresponds exactly to $f_\alpha(\chi)D_\beta^h(z)$, as we wanted. Put another way: by replacing the scalar values $f_\alpha(\chi)$ and $D_\beta^h(z)$ in the coefficients with the above vectors, we end up with the vector
```math
\begin{equation}
\begin{bmatrix}
C^{(p),\beta\alpha,s}_{k,0} \\ 0 \\ 0 \\ 0 \\ 0 \\ 0
\end{bmatrix}.
\end{equation}
```
Summing over all the combinations results in the full vector $\mathbf C_k^{(p),\beta\alpha,s}$. Thus these give precisely the masking vectors we need, but which can be obtained from the Hadamard product (and not the Kronecker product) of the vectorized PDF and FF values.

It's important to understand the difference between e.g. $\vec f_\alpha$ and $\mathbf f_\alpha(\chi)$. $\vec f_\alpha$ is the 2D vector of the two possible values $f_\alpha(\chi)$ and $f_\alpha(\chi/\xi)$ at any given integration and kinematical point. These vectors make the rows of the matrix returned by `parton_values::pdf_rows`. $\mathbf f_\alpha(\chi)$ is the 6D vector that replaces the scalar value $f_\alpha(\chi)$ in the actual coefficient files (i.e. `p.pdf`, where `p` is an instance of `sidis::evaluation_parameters`). All in all, the correspondence goes

- $\mathbf f_{\alpha}(\chi)$: `p.pdf`
- $\mathbf f_{\alpha}(\chi/\xi)$: `p.pdf_hat`
- $\mathbf D_\beta^h(z)$: `p.ff`
- $\mathbf D_\beta^h(z/\zeta)$: `p.ff_hat`
- $\mathbf D_\beta^h(z/\zeta')$: `p.ff_hat_p`

The complications (vectorization and masking) described in this section exist for performance reasons – without them, the code simply runs slower (either due to repeated work or the lack of SIMD vectorization described below).

### SIMD

Vectorization exploits the architecture of modern CPUs. Instead of working with a single number at a time, SIMD keeps multiple distinct numeric values in a single register. These SIMD registers are a realization of vectors: a SIMD register is essentially an $n$-dimensional vector. In SIMD parlance, a vector component is called a lane, and the size of the vector (in bits) is called the width. In addition to SIMD registers, there are SIMD instructions that manipulate these vectors directly by e.g. performing basic arithmetic. 

The width of SIMD registers depends on the hardware. The hardware used in the development of this code supports the AVX-512 instruction set. As the name suggests, the SIMD registers are 512 bits in width. Such a register can fit e.g. 8 double-precision floating-point numbers (each of which takes 64 bits). We only need 6 doubles, so they fit within a single register. Thus, all 6 PDF/FF combinations can be manipulated all at once with AVX-512, instead of one-by-one. Note that the older AVX2 instruction set uses 256-bit registers (and even older instruction sets use even smaller registers), which can only hold 4 doubles – not enough for all 6 values. Therefore the performance of this code can depend on the available vector instruction set, but only slightly (the code is not the most optimal in terms of vectorization to keep the code more readable).

The basic vector template type is `vectorization::vec`, defined as the typealias
```cpp
template <std::size_t lanes>
using vec = std::experimental::fixed_size_simd<double, lanes>;
```
Thus the vectors are just instances of `std::experimental::fixed_size_simd`, which is based on the `C++ parallelism TS v2`, hence the `experimental` namespace. `C++26` will include a proper SIMD library `<simd>`, but this isn't supported by compilers as of now. In the future, `std::experimental::fixed_size_simd` can be replaced by `std::simd::vec`.

## Computing integrals (shared code path)

The rough execution flow is the following:
```
 sidis::instance                                                                         
┌───────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                                                                                   │
│         (entry point) ─────────► lepton_pair_cross_section_xy                                     │
│                                               │                                                   │
│                                               ▼                                                   │
│         (entry point) ─────────► lepton_pair_cross_section_xQ2                                    │
│                                               │                                                   │
│                                               ▼                                                   │
│                                        integrate_common                                           │
│                                               │                                                   │
└───────────────────────────────────────────────┼───────────────────────────────────────────────────┘
                                                │                                              
 sidis::evaluator                               │                                              
┌───────────────────────────────────────────────┼───────────────────────────────────────────────────┐
│                                               ▼                                                   │
│                                      evaluate_observable                                          │
│                                               │                                                   │
│                                               ▼                                                   │
│                                     evaluate_coefficients                                         │
│                                               │                                                   │
│                                               ▼                                                   │
│                                      evaluate_integrand                                           │
│                                               │                                                   │
└───────────────────────────────────────────────┼───────────────────────────────────────────────────┘
                                                │                                              
 sidis::coefficients                            │                                              
┌───────────────────────────────────────────────┼───────────────────────────────────────────────────┐
│                                               ▼                                                   │
│    f{k}_{hw,gc}::{lo,nlo,nnlo}::{quark_to_quark,quark_to_gluon,gluon_to_quark,gluon_to_gluon}     │
│                                                                                                   │
└───────────────────────────────────────────────────────────────────────────────────────────────────┘
```
The entry point to the calculation is `sidis::instance`. The two entry-point functions, `lepton_pair_cross_section_{xQ2,xy}`, correspond to $\mathrm{d}\sigma_{\mu\mu}/\mathrm{d}x \, \mathrm{d}Q^2$ and $\mathrm{d}\sigma_{\mu\mu}/\mathrm{d}x \, \mathrm{d}y$, respectively. They differ by just a Jacobian: the code always computes $\mathrm{d}\sigma_{\mu\mu}/\mathrm{d}x \, \mathrm{d}Q^2$ and then multiplies that by a Jacobian if necessary.

`lepton_pair_cross_section_xQ2` then calls into `integrate_common`. This function is responsible for setting up the integrand and then calling the integration function (provided by `lepton_pair_cross_section_xQ2`) to calculate the necessary integrals up to some given perturbative order. It also enforces a possible energy cut on the primary muon (not the decay muon), which generates a phase space cut for $y$. 

### Setting up the integrand

Let us remind ourselves that we are trying to calculate
```math
\begin{equation}
\sigma=\int \mathrm{d}z \, \mathrm{d}\xi \, \mathrm{d}\zeta \sum_k \sum_{\alpha\beta} \sum_{\mathfrak a \mathfrak b} \sum_p \left(\frac{\alpha_s}{2\pi}\right)^p \sum_s f_{\alpha, \mathfrak a}^s \left(\sum_h B_h(\mathcal K, z) D_{\beta, \mathfrak b}^h\right) C_{k, \mathfrak a\mathfrak b}^{(p),\beta\alpha,s}(\xi, \zeta, \chi, z) \omega_k(\mathcal K),
\end{equation}
```
which can be written as
```math
\begin{equation}
\sigma=\sum_{\alpha\beta}\sum_p \left(\frac{\alpha_s}{2\pi}\right)^p \sigma_{\alpha\beta,p}
\end{equation}
```
where
```math
\begin{equation}
\sigma_{\alpha\beta,p} = \int \mathrm{d}z \, \mathrm{d}\xi \, \mathrm{d}\zeta \sum_k \sum_{\mathfrak a\mathfrak b} \sum_s f_{\alpha, \mathfrak a}^s \left(\sum_h B_h(\mathcal K, z) D_{\beta, \mathfrak b}^h\right) C_{k, \mathfrak a\mathfrak b}^{(p),\beta\alpha,s}(\xi, \zeta, \chi, z) \omega_k(\mathcal K).
\end{equation}
```

`integrate_common` essentially calculates $\sigma_{\alpha\beta,p}$. It returns a `sidis::result` object, which contains
- a 3D tensor `channel_decomposition`
	- the first two dimensions are related to kernel mode and are described in more detail below – in direct mode, these two dimensions are trivial (i.e. the tensor has dimensions $1\times 1 \times n$)
	- the third dimension indexes the physical channel $\alpha\to\beta$
		- the channels are stored in `physical_channels`
- the sum $\sum_{\alpha\beta} \sigma_{\alpha\beta,p}$ as `value`
- the index $p$ for the perturbative order is handled by storing everything as a `perturbative_quantity`, which stores the LO, NLO, and NNLO values separately

`integrate_common` constructs the integrand, which in turn calls `sidis::evaluator::evaluate_observable` with the actual coefficient functions. These coefficient functions depend on the supplied parameters; the selection is handled by `with_coefficient_functions`. The integrand, which has the signature

```cpp
void(const std::span<const double> input, const std::span<double> components)
```
is fed to `integrator_mc`, where
- `input` is the integration point (if the integration point is $(\xi, \zeta, z)$, then `input[0]` is $\xi$, `input[1]` is $\zeta$, and `input[2]` is $z$)
- `components` is an accumulator that tracks individual components
	- the integrator can integrate multiple components at once, using the same grid adaptation
	- this way all physical channels can be computed at once, with nearly no CPU overhead
		- memory usage is increased, as all components must be kept in memory at the same time

### Filling components

When the integrand function is called, the integrand is expected to fill `components` appropriately. When `integrate_common` constructs the integrand function, the filling of `components` is handled by a lambda:
```cpp
[](
	const distribution_index pdf_row, 
	const distribution_index ff_row, 
	const channel_index physical_channel_idx, 
	const double value
) {
	components[layout(pdf_row, ff_row, physical_channel_idx)] += value;
}
```
This lambda (which is called further down the line in `evaluator::evaluate_observable`) is given the appropriate index `(pdf_row, ff_row, physical_channel_idx)`. The layout is determined by
```cpp
const shims::layout_right::mapping layout{
	shims::dextents<std::size_t, 3>{pdf_row_count, ff_row_count, physical_channels.size()}
};
```
The `shims::layout_right::mapping` determines how the multi-dimensional tensor data is stored as one-dimensional data.

*Side note: the `shims` namespace is used for some standard library features that are not available on `GCC 15.2`. They are defined in `utility/shims.hpp` and work using feature macros. If the compiler has support for the library feature, then that implementation is used. If not, a reference implementation fallback is used instead.*

### Channels

Before moving on, we have to consider channels. A partonic channel (`sidis::partonic_channel`) is the fundamental channel type. It describes a partonic channel through the following properties:
- `pto`: perturbative order (i.e. the perturbative order this channel is associated with), used for deciding the correct coefficient function
- `type`: channel type
- `initial_line`: initial-state parton (connected to the PDF)
- `fragmentation_line`: fragmenting parton (connected to the FF)
- `nonfragmentation_lines`: list of final-state partons that don't fragment (e.g. gluon-initiated channels where there are two or more final-state quarks)
- `mass_scheme`: mass scheme
- `true_fragmentation_lines`: if the fragmenting parton isn't a charm quark, then it will first have to split into a charm-anticharm pair, which changes the number of final-state charm quarks and thus the final-state mass 

At LO, there are $3\times 2\times 2=12$ partonic channels:
- 3 possible initial-state partons ($d$, $s$, $b$ in neutrino scattering)
- 2 possible final-state partons ($u$, $c$ in neutrino scattering)
- the conjugates of each channel (e.g. the conjugate of $s\to c$ is $\bar c\to\bar s$; see also `sidis::partonic_channel::conjugate_type`)

At NLO, there are $(3\times 2+3\times 2+3\times 2)\times 2=36$ channels, including the equivalent channels as at LO, gluon fragmentation channels, and gluon-initiated channels.

At NNLO, there are over 300 channels. At this point the classification gets significantly more complicated because of the different types of diagrams with different allowed quark flavors. 

The different channel types, from which all these channels are built by iterating over all allowed flavors, are defined in [[2506.19926](https://inspirehep.net/literature/2938806), §3.3.1, §3.3.2] and implemented in `sidis::partonic_channel::channel_type`. The construction of the channels is handled in `sidis::channel_factory`.

However, there are two kinds of redundancies in these channels. 

For one, we don't actually care about the various subcontributions that make up the contribution $\alpha\to\beta$. That is, we don't need to know exactly how the individual coefficients contribute to e.g. the $s\to c$ channel, we just want to know how much $s\to c$ contributes to the full cross section. A channel that only cares about the initial- and final-state partons ($\alpha\to\beta$) is called a physical channel within the codebase (`sidis::physical_channel`). At LO, there isn't really a difference between partonic and physical channels, but at NLO and NNLO one physical channel can be thought to be made up of multiple partonic channels. Thus, when computing the channel decomposition (i.e. integrating components), it's done over physical channels and not partonic channels. This reduces the amount of components that need to be held in memory. The physical channels are computed from the list of partonic channels supplied to `sidis::instance` using the constructor `physical_channel(const partonic_channel &)` in `sidis::physical_channel`.

For another, some coefficient functions are the same for certain kinds of channels. The meaning of "certain kinds of channels" is implemented in `sidis::coefficient_channel`, but it boils down to what properties of the channel the coefficients depend on. The coefficients are then independent of whatever properties remain, and all channels sharing those remaining properties can use the same coefficient values. Grouping the channels this way leads to `partonic_channel_group<sidis::partonic_channel>`, where the actual grouping is done by `sidis::evaluator::group_by_coefficients`. The remaining properties are actually (almost) orthogonal to physical channels: the coefficients depend on everything else except the initial- and final-state flavors, since those just affect the PDF/FF values. The "(almost)" is that the coefficients differ by at most a sign between a channel type and its conjugate. Thus, conjugate channels are included in a group.

### Evaluating the integrand

`sidis::instance` is responsible for setting up the integrand, while `sidis::evaluator` is responsible for actually evaluating it.

The `sidis::evaluator` constructor takes in a list of partonic channels and a list of physical channels, and uses these to construct the channel groups (as described in *Channels*).

The integrand calls `evaluate_observable`, which takes in 
- the integration point
- parameters
- $3\times 4=12$ functions
	- these are the structure functions $F_1$, $F_2$, and $F_3$, with 4 different combinations depending on the channel type: $q\to q$, $q\to g$, $g\to q$, and $g\to g$
- an `observable` callable, which is responsible for combining the $F_1$, $F_2$, and $F_3$ values into an observable
	- the observable for the cross section is `cross_section_observable`, which essentially implements the sum over $k$ (remember: $k$ indexes the structure functions as $F_k$)
	- this observable callable provides a customization point for calculating e.g. only a specific structure function
- a `fill_component` callable, which handles the filling of the integration component vector as described in *Filling components*

The first thing `evaluate_observable` does is check for threshold points:
- $\xi=1$
- $\zeta=1$
- $\xi=\zeta$

The coefficients can evaluate to `NaN` at these points. For `gc` coefficients, the bound $\zeta\geq z$ is enforced already here. For `hw` coefficients, the check is in the coefficient function implementations, since the SACOT coefficients legitimately go from $\zeta=0$ to $\zeta=1$. For NNLO coefficients, the regions
- $\xi=1$
- $\zeta=1$
- $\zeta=\xi$
- $\zeta=1-\xi$
- $(1-\xi)(1-\zeta)=0$

are excluded within a finite tolerance (i.e. the excluded regions have a non-zero measure) to keep the coefficients finite and to avoid any numerically problematic areas. The coefficients diverge numerically in these regions, but the divergences are either spurious ones that cancel or genuine integrable divergences (like $\xi=1$). In principle these cuts introduce a small systematic error to the calculation, but for the genuine divergences there isn't much that can be done about it, save for some serious analytical work or arbitrary-precision arithmetic.

Next, the PDF/FF values are computed by constructing a `sidis::parton_values` object. The constructor of this object does the actual computation of the PDF/FF values.

Then, an evaluation cache `sidis::evaluation_cache_gc` is constructed. It contains values shared by all (or most) coefficient functions. This is mostly relevant at NNLO, where ~ 70 transcendental function values can be computed just once and then shared across all NNLO coefficients.

After the caches have been created, the structure function weights $\omega_k$ are projected out. This is done by setting $F_k=1$ and the others to zero for each $k$. Because the observable is assumed to be linear in the structure functions (i.e. $\sigma=\sum_k F_k \omega_k$), then setting $F_k=\delta_{kl}$ in the sum gives $\sum_k \delta_{kl} \omega_k = \omega_l$.

Then the real work begins. First, each partonic channel group is iterated over:
```cpp
for (const auto &[coefficient_channel, group] : coefficient_channels) { ... }
```
All channels within `group` share the same coefficient functions (which depend only on the properties of `coefficient_channel`). These are then evaluated only once per channel group with `evaluate_coefficients`. Evaluating these yields a vector (`sidis::vectorization::distribution_product`) for each $k$.

The responsibility of `evaluate_coefficients` is to determine which of the four ($qq$, $qg$, $gq$, $gg$) coefficient types should be used based on the current coefficient channel. 

The actual call to the integrand function selected by `evaluate_coefficients` is done by `evaluate_integrand`. This function enforces the cut $\xi\geq\chi$.

After the coefficients have been evaluated and control returns to `evaluate_observable`, the PDF/FF values are applied for every channel within the group by looping over them:
```cpp
for (const partonic_channel_group<sidis::partonic_channel>::indexed_channel &indexed_channel : group) { ... }
```
In this loop, the $k$–indexed contributions are combined with $\omega_k$ using the weights obtained from the `observable` callable into a single vector. Then, all PDF/FF values are applied by computing $\mathbf 1_6\cdot (\vec f_\alpha \otimes \vec D_\beta^h)*\mathbf C_{k}^{(p),\beta\alpha,s}$ (see eq. (11)). This computation is done for each PDF/FF row, and the results are filled into the output tensor by finally calling the `fill_component` lambda. Note that the actual computation uses eq. (13). Again, in direct mode there is only one row of PDF and FF values, but in kernel mode there are more.

The logic of how `evaluate_observable` works is a slight simplification of how the code actually works. The call to `evaluate_coefficients` is actually deferred until a non-empty channel (i.e. one where at least one PDF/FF value is non-zero) within `group` is found.

#### Integration methods

In direct mode, adaptive (quasi) Monte-Carlo integration methods are used. There are two options: Vegas and Suave, both implemented by the [Cuba](https://feynarts.de/cuba/) package (default is Suave). Vegas has constant memory usage, whereas the memory used by Suave is proportional to the number of evaluations (or subdivisions, since Suave keeps track of those, which is why its memory consumption is not constant).

In kernel mode, a simple non-adaptive quasi-MC method is used instead. Since all kernels are computed in one go, the integrands are all different and there is no point in trying to adapt the integration grid for such a diverse set of integrals.

### Coefficient implementations

There are two kinds of coefficient implementations, `hw` ("handwritten") and `gc` ("generated code"). As the names suggest, the `hw` coefficients are the original hand-implemented coefficients and `gc` coefficients have been generated by a Mathematica program. 

Available configurations:

| 	   			| `hw` 				| `gc` 				|
| ------------- | ----------------- | ----------------- |
| LO (ZM/SR)	| yes *(default)* 	| yes				|
| LO (SACOT)	| yes *(default)*	| yes				|
| NLO (ZM/SR)	| yes *(default)*	| yes				|
| NLO (SACOT)	| yes *(default)*	| no				|
| NNLO (ZM/SR)	| no				| yes *(default)*	|
| NNLO (SACOT)	| no				| no				|

The namespace structure is the same for all kinds of coefficient functions:
```cpp
sidis::coefficients::f{k}_{hw,gc}::{lo,nlo,nnlo}::{quark_to_quark,quark_to_gluon,gluon_to_quark,gluon_to_gluon}
```
Implementation details for `hw` coefficients are hidden in the `sidis::coefficients::_impl` namespace. 

All coefficient functions have the same function signature:
```cpp
vectorization::distribution_product(const sidis::evaluation_parameters &p)
```
`p` contains all values needed for computing the coefficients.

#### `hw` coefficients

The implementation of `hw` coefficients is straightforward – just converting the analytical expressions into numerical code.

#### `gc` coefficients

There are two parts to the `gc` coefficients: generation and implementation in the codebase.

##### Generation

The code generation is done with Mathematica. The zero-mass coefficient functions up to NNLO are given as a [single FORM file](https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2FJHEP10%282025%29016/MediaObjects/13130_2025_27309_MOESM1_ESM.zip). This file is first preprocessed to remove the comments at the start, so that the file begins with the $F_T$ LO coefficient, as in
```
id FTC0q2qM      = (  + Dd([1-x])*Dd([1-z]) * ( 1 ) );
```
This preprocessed file is then imported into Mathematica as a string and processed further by 
- removing all instances of "id"
- renaming `x` to `xi` and `z` to `zeta` (to match the notation)
- replacing occurrences of bracketed expressions like `[1-xi]` with a variable like `onemxi`
- rewriting transcendental functions into their Mathematica formats by e.g. replacing the pattern `Li2(...)` with `PolyLog[2, ...]`

At this point, the string contains a valid Mathematica expression, starting with
```
FTC0q2qM      = (  + Dd[onemxi]*Dd[onemzeta] * ( 1 ) );
```
This expression thus defines the coefficients as variables like `FTC0q2qM`. This expression string can therefore be evaluated with the Mathematica kernel with `ToExpression`, after which e.g. `FTC0q2qM` is a valid Mathematica variable.

Indeed, the coefficients are given as the structure functions $F_T, F_L, F_3$ instead of $F_1, F_2, F_3$. These former are converted to the latter with the simple formulas
```math
\begin{equation}
F_1(F_T)=F_T, \quad F_2(F_T, F_L)=2(F_T+F_L), \quad\text{and}\quad F_3\to 2F_3.
\end{equation}
```
*Side note: there might be an optimization opportunity here, as $F_L$ is typically much simpler in comparison to $F_T$, while the mixed term $F_2$ has high complexity but a lot of shared terms (although this is at least partially alleviated by CSE – see below).*

The coefficients at this point contain delta and plus distributions. For that, a given coefficient is multiplied by `pdf[x/xi]*ff[z/zeta]` and the resulting expression is expanded, yielding e.g.
```mma
Dd[onemxi] Dd[onemzeta] ff[z/zeta] pdf[x/xi]
```
The distributions can effectively be calculated under the integral sign with
```math
\begin{equation}
\delta(1-\xi)f(x/\xi)\to \frac{f(x)}{1-x}
\end{equation}
```
and
```math
\begin{equation}
\left[\frac{\log^n(1-\xi)}{1-\xi}\right]_+f(x/\xi)\to (f(x/\xi)-f(x))\frac{\log^n(1-\xi)}{1-\xi}+f(x)\frac{\log^{n+1}(1-x)}{(1-x)(n+1)}
\end{equation}
```
and similarly for the FF with $\zeta$. The division by $1-x$ is a division by the Jacobian, since e.g.
```math
\begin{equation}
\int_x^1 \mathrm{d}\xi \, \delta(1-\xi)f(x/\xi)=f(x)=\int_x^1 \mathrm{d}\xi \, \frac{f(x)}{1-x}.
\end{equation}
```
By doing this, there's no need to separately keep track of terms inside and outside the integration. 

Note that integrating out the distributions like this is only possible when no terms have any other dependence on the variable of the distribution. For example, applying the simple replacement to a term like
```math
\begin{equation}
\frac{\xi}{(1-\xi)_+}f(x/\xi)
\end{equation}
```
would fail because of the $\xi$ in the numerator. The replacement would simply give
```math
\begin{equation}
\frac{\xi}{(1-\xi)_+}f(x/\xi)\to \xi \frac{f(x/\xi)-f(x)}{1-\xi}+\xi f(x)\frac{\log(1-x)}{1-x}
\end{equation}
```
whereas the correct result is given by
```math
\begin{equation}
\int_x^1 \mathrm{d}\xi \, \frac{\xi}{(1-\xi)_+}f(x/\xi)=\int_x^1 \mathrm{d}\xi \, \frac{\textcolor{red}{\xi} f(x/\xi)-f(x)}{1-\xi}+f(x)\log(1-x)=\int_x^1 \mathrm{d}\xi \left[\frac{\textcolor{red}{\xi} f(x/\xi)-f(x)}{1-\xi}+f(x)\frac{\log(1-x)}{1-x} \right].
\end{equation}
```
It turns out that the coefficient functions are already written in such a way that the above cannot happen, although it is always possible to write them in such a form if needed because of the relation [[PhD](https://jyx.jyu.fi/jyx/Record/jyx_123456789_113028?sid=851737269), §A]
```math
\begin{equation}
[g(x)]_+ h(x)=[g(x)]_+ h(1) +g(x)(h(x)-h(1)).
\end{equation}
```

After the distributions have been dealt with, we are left with an expression for $f_\alpha(\chi/\xi) D_\beta^h(z/\zeta) C_k^{(p),\beta\alpha,s}(\xi, \zeta)$ (see eq. (11)) under the integral sign. The resulting expression can be written in a basis built up of products of the 6 PDF-side terms
- $\frac{f_\alpha(\chi)}{1-\chi}$
- $f_\alpha(\chi/\xi)$
- $(f_\alpha(\chi/\xi)-f_\alpha(\chi))\frac{\log^n(1-\xi)}{1-\xi}+\frac{f_\alpha(\chi)}{1-\chi}\frac{\log^{n+1}(1-\chi)}{n+1}$ for $n=0, 1, 2, 3$

and the 6 analogous FF-side terms. Let's denote these products by $\mathcal{B}_l$ (where e.g. $\mathcal{B}_1=\frac{f_\alpha(\chi)D_\beta^h(z)}{(1-\chi)(1-z)}$), so that
```math
\begin{equation}
f_\alpha(\chi/\xi) D_\beta^h(z/\zeta) C_k^{(p),\beta\alpha,s}(\xi, \zeta)=\sum_l c_l \mathcal{B}_l.
\end{equation}
```
This basis is chosen such that it has the following properties:
- the coefficients $c_l$ depend only on $\xi$ and $\zeta$; the $\chi$- and $z$-dependence, which arises from the distributions, is contained entirely in the basis functions $\mathcal{B}_l$
- the terms arising from the plus distributions stay together (idea being to keep numerical stability, although this is an untested assumption and could be unnecessary)

After the change of basis (`splitByCoefficients`), the coefficients are a list of expressions like
```mma
coefficientExpression[
	(c1 ff[z] pdf[x])/((1 - x) (1 - z)), 
	{coefficient[c1, 1]}, {}
]
```
where 
- the first argument of `coefficientExpression` is the actual sum (eq. (34)),
- the second argument is a list of `coefficient`s that specify the $c_l$ variable name and the associated expression,
- the third argument is a list of region expressions (applicable at NNLO, see below).

The coefficients then go through a series of rewrites and optimizations.

1. Regions are converted into a standard format with `extractRegions`. The NNLO coefficients include $\theta$ functions, which select certain branches in certain regions of the integration space. These $\theta$ functions are collected such that all terms belonging to a given region are together and not scattered throughout the expression. The regions are defined in [[2201.06982](https://inspirehep.net/literature/2012936): eqs. (5.9), (5.12), (5.16)]. For example, the expression 
	```mma
	coefficientExpression[
		(c1 ff[z] pdf[x])/((1 - x) (1 - z)), 
		{coefficient[c1, 1 + 2 xi T[r1]]}, {}
	]
	```
	gets turned into
	```mma
	coefficientExpression[
		(c1 ff[z] pdf[x])/((1 - x) (1 - z)), 
		{coefficient[c1, 1 + Tr1]}, {region[Tr1, 2 xi]}
	]
	```
	by `extractRegions`.

2. Common transcendentals (logs, polylogs, arc tangents, and inverse tangent integrals) are collected by `collapseCommonTranscendentals`. This means that every instance of those functions gets replaced by a variable, which is then stored as a rule. For example, an expression like
	```mma
	Log[x] pdf[x] + Log[x] pdf[x/xi] + PolyLog[2, x] pdf[x/xi]
	```
	gets converted into
	```mma
	logv1 pdf[x] + logv1 pdf[x/xi] + polylogv1 pdf[x/xi]
	```
	with the associated rules
	```mma
	{
		Log[x] -> logv1, 
		PolyLog[2, x] -> polylogv1
	}
	```
	In actuality, `collapseCommonTranscendentals` operates on a list of `coefficientExpressions`. Thus a complete version of the above example is that the expression
	```mma
	{
		coefficientExpression[
			(c1 ff[z] pdf[x])/((1 - x) (1 - z)), 
			{coefficient[c1, Log[x] pdf[x] + PolyLog[2, x] pdf[x/xi] + Log[x] pdf[x/xi]]}, {}
		]
	}
	```
	is turned into
	```mma
	computation[
		{
			coefficientExpression[
				(c1 ff[z] pdf[x])/((1 - x) (1 - z)), 
				{coefficient[c1, logv1 pdf[x] + logv1 pdf[x/xi] + polylogv1 pdf[x/xi]]}, {}
			]
		},
		{Log[x] -> logv1, PolyLog[2, x] -> polylogv1}, {}
	]
	```
	where `computation` now models the entire coefficient evaluation process, where
	1. the first argument is a list of `coefficientExpressions`
	2. the second argument is a list of transcendental replacement rules
	3. the third argument is a list of CSE replacement rules (filled at the next step)

3. Common subexpression elimination (CSE) is performed by `cse`. It goes through all subexpressions and selects subexpressions which
	- appear more than once
	- have a `LeafCount` of at least `minLeafCount` (argument to `cse`; this parameter controls how small the common subexpressions are allowed to be, and thus the amount of them)
	- satisfy the general `acceptance` function (argument to `cse`; used to stop `cse` from picking up e.g. entire `coefficient` expressions with the head included)

4. The `coefficient` and `region` expressions are lifted from individual `coefficientExpressions` to the global `computation` scope, since they are independent of the channel and can therefore be computed once per integration point.

Next, the intermediate representation (IR) is generated. This uses pattern matching to walk through the expression trees and tag them. Examples of these tags are
```mma
intermediateRepresentation[x_?NumberQ] := IRnumber[x]
intermediateRepresentation[x] := IRspecial[x]
intermediateRepresentation[expr_] := applyReplace[IRproduct, intermediateRepresentation, expr] /; Head[expr] == Times
intermediateRepresentation[Rational[a_, b_]] := IRfraction[intermediateRepresentation[a], intermediateRepresentation[b]]
intermediateRepresentation[Power[x_, n_?IntegerQ]] := IRpow[intermediateRepresentation[x], IRinteger[n]]
```
These examples work as follows:
1. Matches numbers and tags them with `IRnumber`. Example: `intermediateRepresentation[2]` $\to$ `IRnumber[2]`
2. Matches the variable `x` (momentum fraction) and tags it with `IRspecial[x]`. Example: `intermediateRepresentation[x]` $\to$ `IRspecial[x]`
3. Matches products, tags them with `IRproduct`, and then recurses over the factors. Example: `intermediateRepresentation[2 x]` $\to$ `IRproduct[IRnumber[2], IRspecial[x]]`
4. Matches rational numbers $a/b$, tags them with `IRfraction`, and recurses over the numerator and denominator. Example: `x/2` $\to$ `IRproduct[IRfraction[IRnumber[1], IRnumber[2]], IRspecial[x]]`
5. Matches powers $x^n$, tags them with `IRpow`, tags the exponent with `IRinteger`, and recurses over the base. Example: `intermediateRepresentation[x^2]` $\to$ `IRpow[IRspecial[x], IRinteger[2]]`

There are many more rules, which tag e.g. the known variables and functions. As a more complicated example,
```mma
intermediateRepresentation[Log[x]^2/(1 - x)]
```
evaluates to
```mma
IRproduct[
	IRpow[
		IRsum[
			IRnumber[1], 
			IRproduct[
				IRnumber[-1], IRspecial[x]
			]
		], 
		IRinteger[-1]
	], 
	IRpow[
		IRlog[IRspecial[x]], IRinteger[2]
	]
]
```

After the IR has been created, the code can be generated. The code generator reads the IR tags and generates code accordingly. The code generation is the reason for making the IR, since the generation becomes more robust when all the necessary information has already been tagged. The generator again uses pattern matching to walk through the IR expression tree. Example rules are:
```mma
codegen[IRvariable[x_]] := SymbolName[x]
codegen[IRnumber[x_]] := ToString[N[x]]
codegen[IRspecial[x]] := "p.pdf_momentum_fraction"
codegen[IRpow[x_, n_]] := templateFunction["math::power", n, x]
```
These work as follows:
1. Matches variables and simply prints the variable name
2. Matches numbers and prints the number as a string
3. Matches the variable `x` specifically and prints `p.pdf_momentum_fraction`, since this is how the code actually stores $x$
4. Matches exponentiation and invokes a general function intended for writing function templates. The code already contains a function template `template <int n> power(auto x)`, where the base (`x`) is taken as a normal argument, but the integer exponent as a template parameter (this is to sidestep performance issues with `std::pow` when the exponent is an integer). The code generator writes out the templated function call as `math::power<n>(x)`.

The variables (of types `v`, `c`, `T`) that are produced throughout the process are used in multiple files handling the caching. First, the coefficient files (one file per coefficient) are generated from `coefficientExpressions` and have the form
```cpp
const double c1t0 = p.c.{pto}.c1t0;
const double c2t0 = p.c.{pto}.c2t0;
// ...
return /* ... */;
```
where the first lines bring the $c_l$ variables to local scope from the per-perturbative-order `pto` cache (`p.c.{pto}`). The $c_l$ variable names are unique across all coefficients through the `t` suffix. Finally, the return statement calculates the sum in eq. (34). As an example, the full `F1C0q2qM.hpp` file reads
```cpp
const double c1t0 = p.c.lo.c1t0;
return (c1t0*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf);
```

Each perturbative order has three common files: `{pto}_cache.hpp`, `{pto}_decl.hpp`, and `{pto}_vars.hpp`. The cache files contain lists of variable assignments for all the common cached variables. The first 5 lines of `nlo_cache.hpp` read
```cpp
logv1 = std::log((1.+(-1.*p.zeta)));
logv2 = std::log(p.zeta);
logv3 = std::log((1.+(-1.*p.xi)));
logv4 = std::log(p.xi);
v0 = (1.+(-2.*p.xi)+(2.*math::power<2>(p.xi)));
// ...
```
The `decl` files contain lists of variable assignments that bring the appropriate variables to scope, e.g.
```cpp
[[maybe_unused]] const double logv1 = p.c.nlo.logv1;
[[maybe_unused]] const double logv2 = p.c.nlo.logv2;
[[maybe_unused]] const double logv3 = p.c.nlo.logv3;
[[maybe_unused]] const double logv4 = p.c.nlo.logv4;
[[maybe_unused]] const double v0 = p.c.nlo.v0;
// ...
```
Finally, the `vars` files contain lists of variable declarations, such as
```cpp
double logv1;
double logv2;
double logv3;
double logv4;
double v0;
// ...
```

##### Implementation

For implementation, the fact that `#include` is a preprocessor directive that literally inserts the contents of the included file at the specified position is (ab)used to wire the incomplete `hpp` files generated by the code generator.

The coefficient files are wired using a lambda trick. As an example, the LO $qq$ coefficient for $F_1$ is written like this:
```cpp
vectorization::distribution_product quark_to_quark(const sidis::evaluation_parameters &p) {
	return [&p]() -> vectorization::distribution_product {
		#include "dimuon/sidis/coefficients/gc/lo_decl.hpp"
		#include "dimuon/sidis/coefficients/gc/F1C0q2qM.hpp"
	}() / p.z / p.pdf_momentum_fraction;
}
```
The function consists of an immediately-invoked lambda expression. The lambda, which just contains `#include` statements, will get filled by the preprocessor with the contents of the `decl` file of the perturbative order (`lo_decl.hpp`, which is empty at LO) and of the `F1C0q2qM.hpp` file (which is shown above). The bound $\zeta\geq z$ isn't checked here, since `evaluate_observable` already enforces it. The lambda is invoked immediately (hence the `()` right after its declaration), which outputs the value of the coefficient function. 

The point of this lambda trick is to keep the code generator as simple as possible, and disconnected from the code itself so that changes to the coefficient evaluation code will hopefully not require any changes to the code generator.

Now, the coefficient files refer to `p.c.{pto}`. Here, `p` is the `sidis::evaluation_parameters` instance that's available in the scope of the lambda body. The `evaluation_parameters` contain an `evaluation_cache_gc` field with the name `c`, which is responsible for handling the caching of variables. 

Inside `evaluation_cache_gc`, there are specific cache variables for each perturbative order, each with its own corresponding cache type `{pto}_cache` within `evaluation_cache_gc`. Each `{pto}_cache` object then contains the necessary variables. The variables are read and updated using the `decl` and `vars` files, e.g.
```cpp
struct nlo_cache {
	void compute(const inputs &p, const evaluation_cache_gc *c) {
		const double LMUF = c->LMUF;
		const double LMUA = c->LMUA;

		#include "dimuon/sidis/coefficients/gc/nlo_cache.hpp"
	}

	#include "dimuon/sidis/coefficients/gc/nlo_vars.hpp"
};
```

### `c++` implementation of math functions

The coefficients (particularly at NNLO) contain some math functions that warrant a mention.

#### Exponentiation

`math::power` implements $x^n$ for $x\in \mathbb{R}$ and $n\in\mathbb{Z}$. The standard function `std::pow` implements the more general case $x^y$ for $x,y\in \mathbb{R}$ and as a consequence can be slower than `math::power`. Whether `std::pow` gets optimized away by the compiler for an integer exponent depends on the compiler, the exponent, and the optimization flags. For these reasons, it's safer to hand-roll the implementation. 

`math::power` is a function template, where the exponentiation is unfolded into repeated multiplication using exponentiation-by-squaring at compile time. That is, the string of multiplications is generated at compile time, and thus the exponent has to be known at compile time. In the case of the coefficient functions, the exponents are indeed known at compile time.

#### Polylogarithms

In general, the polylogarithm $\operatorname{Li}_s$ can be defined as the series
```math
\begin{equation}
    \operatorname{Li}_s(z)\equiv\sum_{k=1}^\infty \frac{z^k}{k^s}
\end{equation}
```
or recursively as an integral with
```math
\begin{equation}
    \operatorname{Li}_1(z)\equiv -\log(1-z), \quad \operatorname{Li}_{s+1}(z)\equiv \int_0^z \mathrm{d} t \, \frac{\operatorname{Li}_s(t)}{t}.
\end{equation}
```
In the coefficients the di- and trilogarithms need to be evaluated on the real axis only. In fact, using a number of polylogarithm identities, the range of any input $x\in\mathbb{R}$ can be reduced to a union of finite intervals. For example, for the dilogarithm ($\operatorname{Li}_2$) it's enough to be able to evaluate it in $[0, 1/2]$. The numerical implementation provided by [[2201.01678]](https://inspirehep.net/literature/2001826), which uses a rational approximation of the form
```math
\begin{equation}
    \hat{\operatorname{Li}}_s(x)=\frac{x\sum_{i=0}^n a_i x^i}{\sum_{j=0}^m b_j x^j}.
\end{equation}
```
The coefficients $\{a_i\}_{i=0}^n$ and $\{b_j\}_{j=0}^m$, which are listed in refs. [[2201.01678]](https://inspirehep.net/literature/2001826), [[2308.11619]](https://inspirehep.net/literature/2690384), can be obtained with the minimax algorithm. 

##### The minimax algorithm

Given any sufficiently well-behaved function $f\colon I\to \mathbb{R}$, the minimax algorithm attempts to find a simpler (in terms of numerical evaluation) function $\hat f\colon I\to \mathbb{R}$ that approximates $f$ to some accuracy $\epsilon$ over the interval $I$. More precisely, such an $\hat f$ has to satisfy $\|f-\hat f\|_\infty <\epsilon$, where $\|\cdot\|_\infty$ is the supremum (or $L^\infty$) norm,
```math
\begin{equation}
    \|f-\hat f\|_\infty \equiv \sup_{x\in I}\left|f(x)-\hat f(x)\right|.
\end{equation}
```

The Stone-Weierstrass theorem says that $\hat f$ can be as simple as a polynomial. However, a rational approximation is sometimes more practical as it can yield a better approximation with the same amount of coefficients as a polynomial approximation. Wolfram Mathematica implements the rational version of this algorithm in the function ``FunctionApproximations`MiniMaxApproximation``. `MiniMaxApproximation` minimizes the relative error (unlike what was written in eq. (38)), which means that $f$ cannot have a zero in the interval $I$. If $f$ has a zero of order $k$ at $a$, one can instead consider approximating the function $x\mapsto f(x)/(x-a)^k$ such that this new function does not have a zero in $I$.

#### Inverse tangent integral

The inverse tangent integral $\operatorname{Ti}_2$ is closely related to the dilogarithm, and can be defined as the integral
```math
\begin{equation}
    \operatorname{Ti}_2(z)\equiv \int_0^z \mathrm{d} t \, \frac{\arctan t}{t}.
\end{equation}
```
Like the polylogarithm, the inverse tangent integral can be generalized with the recursive definition $\operatorname{Ti}_{s+1}(z)\equiv \int_0^z \mathrm{d} t \, \operatorname{Ti}_s(t)/t$, although this is not needed for the present purposes. From eqs. (36) and (39) one can see that
```math
\begin{equation}
    \operatorname{Ti}_2(x)=\operatorname{Im}(\operatorname{Li}_2(ix))
\end{equation}
```
for all $x\in \mathbb{R}$. Again, $\operatorname{Ti}_2$ needs to be evaluated only on the real axis. Using the identities
```math
\begin{align}
    \operatorname{Ti}_2(-x)&=-\operatorname{Ti}_2(x) & \text{for all }x\in \mathbb{R},& \\
    \operatorname{Ti}_2(x)-\operatorname{Ti}_2(1/x)&=\frac{\pi}{2}\log x & \text{for all }x>0,&
\end{align}
```
one can again reduce the range of any argument $x\in \mathbb{R}$ to the unit interval $\left[0, 1\right]$ as follows:
```math
\begin{equation}
    \operatorname{Ti}_2(x)=\begin{cases}
        -\operatorname{Ti}_2(-1/x)-\frac{\pi}{2}\log(-x), &x<-1, \\
        -\operatorname{Ti}_2(-x), &-1\leq x<0, \\
        \operatorname{Ti}_2(x), &0\leq x\leq 1, \\
        \operatorname{Ti}_2(1/x)+\frac{\pi}{2}\log x, &x>1.
    \end{cases}
\end{equation}
```

*Side note: eq. (41) follows from the integral definition in eq. (39) by a change of variables, while eq. (42) can be obtained by using $\operatorname{Ti}_2'(x)=\arctan(x)/x$, the identity $\arctan(x)+\arctan(1/x)=\pi/2$ for $x>0$, and the integral representation of the logarithm $\log x=\int_1^x \mathrm{d} t/t$.*

Let's take the same approach to approximating $\operatorname{Ti}_2$ as [[2201.01678]](https://inspirehep.net/literature/2001826), [[2308.11619]](https://inspirehep.net/literature/2690384) took for the polylogarithm and approximate the inverse tangent integral as the rational function
```math
\begin{equation}
    \hat{\operatorname{Ti}}_2(x)=\frac{x\sum_{i=0}^8 a_i x^i}{\sum_{j=0}^8 b_j x^j}
\end{equation}
```
over the unit interval $\left[0, 1\right]$, where the coefficients $\{a_i\}_{i=0}^8$ and $\{b_j\}_{j=0}^8$ are listed below. The relative error of the approximation is below $1.5\times 10^{-17}$ on the unit interval. This is for exact arithmetic, with non-truncated values for $a_n$ and $b_n$. The accuracy in practice will be lower. Note that achieving a similar accuracy with Chebyshev polynomials requires $\mathcal{O}(20)$ terms.

| $n$	| $a_n$						| $b_n$						|			
| ----- | ------------------------- | ------------------------- |
| 0 	| `1.000000000000000E+00` 	| `1.000000000000000E+00`	|
| 1 	| `1.359606750543622E+00` 	| `1.359606750543613E+00`	|
| 2 	| `2.294603108673590E+00` 	| `2.405714219785837E+00`	|
| 3 	| `1.965431829017859E+00` 	| `2.116499245692448E+00`	|
| 4 	| `1.479377395374053E+00` 	| `1.706678976619235E+00`	|
| 5 	| `7.303876012387028E-01` 	| `9.111698955486797E-01`	|
| 6 	| `2.481653896797473E-01` 	| `3.619761589327320E-01`	|
| 7 	| `4.618705978228984E-02` 	| `9.051409249480067E-02`	|
| 8 	| `6.316042475722170E-04` 	| `9.340084117634812E-03`	|

# Direct mode

In direct mode, a PDF and FFs are supplied and the cross section value (with channel decomposition) is computed normally for that PDF/FF. The PDF/FF values are computed and stored in `sidis::parton_values`. PDFs and FFs are handled by the `lha_interface` class, which uses `LHAPDF` under the hood.

The output from `dimuon` in direct mode has the format
```
===
pdf name, pdf mem, x, y, E, Q2, ren scale, fact scale, frag scale, alpha_s, lo, nlo, nnlo, nlo corr, nnlo corr
===
initial flavor, fragmentation flavor, lo, nlo corr, nnlo corr
...
===
lo sum
nlo sum
nnlo sum
===
```
The first block has one line, containing
- `pdf name`: name of the PDF set
- `pdf mem`: member of the PDF set
- `x`, `y`, `E`, `Q2`: kinematical variables
- `ren scale`: renormalization scale
- `fact scale`: factorization scale
- `frag scale`: fragmentation scale
- `alpha_s`: value of $\alpha_s$ at the renormalization scale
- `lo`, `nlo`, `nnlo`: LO, NLO, NNLO cross section
- `nlo corr`: the NLO correction multiplied by $\alpha_s/2\pi$
- `nnlo corr`: the NNLO correction multiplied by $(\alpha_s/2\pi)^2$

The next block has a line for each physical channel $\alpha\to\beta$, containing
- `initial flavor`: $\alpha$
- `fragmentation flavor`: $\beta$
- `lo`: LO contribution to the cross section
- `nlo corr`, `nnlo corr`: NLO/NNLO contribution to the cross section (i.e. the NLO/NNLO correction *without* the $\alpha_s$ factor)

The final block contains the sum of `lo`, `nlo corr`, `nnlo corr` over the physical channels, and is meant as a sanity check. Since the channel lines don't include the $\alpha_s$ factors, the sums give the `lo`, `nlo corr`, and `nnlo corr` of the first block, after multiplying them by $1$, $\alpha_s/2\pi$, and $(\alpha_s/2\pi)^2$, respectively.

The DIS observables (`dis.differential.charm` and `dis.differential.total`) have the same blocks, without the fragmentation scale and flavor, and only up to NLO:
```
===
pdf name, pdf mem, x, y, E, Q2, ren scale, fact scale, alpha_s, lo, nlo, nlo corr
===
initial flavor, lo, nlo corr
...
===
lo sum
nlo sum
===
```

# Kernel mode

## Interpolation kernels

The PDF and FF can be written as
```math
\begin{align}
f_\alpha(x)&\to \hat f_\alpha(x)= \sum_i \overbrace{f_\alpha(x_i)}^{\equiv f_{\alpha, i}} \phi_i(x) \\
D_\beta^h(z)&\to \hat D_\beta^h(z)= \sum_j \overbrace{D_\beta^h(z_j)}^{\equiv D_{\beta, j}^h} \phi_j(z)
\end{align}
```
where 
- $\{x_i\}_i$ and $\{z_j\}_j$ are grid points contained within $[0, 1]$
	- the grids are fixed
- $\phi$ is an interpolation kernel
	- $\phi_i(x_j)=\delta_{ij}$
	- $\sum_i \phi_i(x)=1$
	- with these properties,
		- $\hat f(x_k)=\sum_i f(x_i)\phi_i(x_k)=\sum_i f(x_i)\delta_{ik}=f(x_k)$
		- $\hat 1(x)=\sum_i\phi_i(x)=1$
- the kernels are fixed
- the interpolation-kernel expansion is approximate for a finite grid
- the grids still depend on parameters like scales, charm mass, and quark-mass thresholds (for NNLO)
- remember again that e.g. the PDF is really the "$xf(x)$" expression, so that $f_{\alpha, i}=x_i f(x_i)$

Suppressing the $\mu\mu$, the per-hadron cross section (for which the tables are produced) then reads
```math
\begin{equation}
\sigma^h=\sum_{\alpha\beta}\sum_{ij} \sum_p \left(\frac{\alpha_s}{2\pi}\right)^p f_{\alpha,i} D_{\beta,j}^h \sigma_{ij,\alpha\beta}^{(p),h}
\end{equation}
```
with
```math
\begin{equation}
\sigma_{ij,\alpha\beta}^{(p),h}=\int \mathrm{d}z \, \mathrm{d}\xi \, \mathrm{d}\zeta \sum_k \sum_{\mathfrak a \mathfrak b} \sum_s \phi_{i,\mathfrak a}^s \phi_{j,\mathfrak b} C_{k, \mathfrak a\mathfrak b}^{(p),\beta\alpha,s}(\xi, \zeta, \chi, z) \omega_k(\mathcal K)B_h(\mathcal K, z)
\end{equation}
```
where $\mathfrak a$ and $\mathfrak b$ denote, as before, at which value the kernels are evaluated. This is now the interpolation table: once the tensor $\sigma_{ij,\alpha\beta}^{(p),h}$ is calculated (which is independent of the PDF and FF), the full cross section becomes a simple tensor contraction.

All the kernel values are computed and contained in `sidis::parton_values`. 

## Choice of kernels and knots

The simplest interpolation is linear interpolation, where the kernel is a triangular hat function. Linear interpolation is sufficient for the task at hand, so that is the kernel used and implemented in `linear_interpolation_basis`. Linear interpolation kernels have the advantage that at any given point, only a maximum of two kernels contribute to that point:
```
																				
                              │                                                   
 \     /\     /\     /\     /\│    /\     /\     /\     /\     /\     /\     /    
  \   /  \   /  \   /  \   /  │   /  \   /  \   /  \   /  \   /  \   /  \   /     
   \ /    \ /    \ /    \ /   │\ /    \ /    \ /    \ /    \ /    \ /    \ /      
    /      \      /      \    │ /      \      /      \      /      \      /       
   / \    / \    / \    / \   │/ \    / \    / \    / \    / \    / \    / \      
  /   \  /   \  /   \  /   \  │   \  /   \  /   \  /   \  /   \  /   \  /   \     
 /     \/     \/     \/     \/│    \/     \/     \/     \/     \/     \/     \    
──────────────────────────────┼──────────────────────────────────────────────────►
                              │                                                   
																				
```
A single kernel is defined by three points: left point $x_l$, middle point $x_m$, and right point $x_r$ with $x_l\leq x_m\leq x_r$. The kernel obtains the value 1 at $x_m$, goes linearly to zero on both sides at $x_l$ and $x_r$, and is zero outside the interval $[x_l, x_r]$. At the edges of the interpolation grid, half-kernels are used where either $x_l=x_m$ or $x_r=x_m$, and there is only one slope.

The interpolation grid (i.e. the knots) is defined by the [Lambert grid](https://eko.readthedocs.io/en/latest/modules/eko/eko.html#eko.interpolation.lambertgrid), which transitions smoothly from logarithmic spacing near 0 to linear spacing near 1. This grid is used for PDFs. For FFs, a modified version where the roles of 0 and 1 are switched is used instead. This is because the behaviour of FFs near the endpoints differs qualitatively from the behaviour of PDFs. 

## Output format

The output from `dimuon` in kernel mode has the format
```
===
x, y, E, Q2, ren scale, fact scale, frag scale
===
pdf kernel idx, ff kernel idx, initial flavor, fragmentation flavor, lo, nlo corr, nnlo corr
...
===
```
The first block has one line, containing
- `x`, `y`, `E`, `Q2`: kinematical variables
- `ren scale`: renormalization scale
- `fact scale`: factorization scale
- `frag scale`: fragmentation scale

The last block has a line for each combination of PDF kernel, FF kernel, and physical channel $\alpha\to\beta$, containing
- `pdf kernel idx`: index of the PDF kernel ($i$)
- `ff kernel idx`: index of the FF kernel ($j$)
- `initial flavor`: $\alpha$
- `fragmentation flavor`: $\beta$
- `lo`: the LO value
- `nlo corr`: the NLO correction (*without* the $\alpha_s$ factor)
- `nnlo corr`: the NNLO correction (*without* the $\alpha_s$ factor)