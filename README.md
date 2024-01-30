# Anisotropic Extensions of Nonlocal Interaction Models

**Research Funding and Support Information:**

-   This research was conducted during the summer of 2021 under the leadership of Dr. Ihsan Topaloglu

-   This research was funded by a VCU UROP Fellowship. Information about this fellowship [can be found
    here](https://academics.provost.vcu.edu/transformative-learning/undergraduate-research/urop-fellowship/)

-   This research was presented at the [VCU Poster Symposium for
    Undergraduate Research and Creativity](https://academics.provost.vcu.edu/transformative-learning/undergraduate-research/poster/) in the spring of 2022. The poster from this presentation is `UROP Poster.pdf`

-   Some of the numerical results for this research were obtained using
    the Teal Cluster from VCU’s High Performance Computing Facility.
    Information about the facility [can be found
    here](https://research.vcu.edu/cores/hprc/facilities/)

-   The code and data in this repository is the work
    of Christopher Flippen unless otherwise stated

**Poster Abstract:**

In this project, we investigated pattern formations which arise from
self-organizing mathematical models significant to physical, biological,
and engineering systems. The mathematical problems examined in our
research have applications in models of materials science and in
collective behavior of multi-agent systems such as superconductor
vortices, robotic swarms, and biological aggregations. These patterns
can be described as minimal energy configurations $N$ interacting
particles. Most nonlocal interaction models in the literature consider
only isotropic interaction energies where the interaction kernel is
radially symmetric. In this project, we first examined pattern
formations induced by several smooth and crystalline anisotropies. The
results of our numerical experiments showed significant quantitative
differences in the properties of the energy-optimal configurations
between smooth and crystalline anisotropies.

**Theory Explanation:**

The Hamiltonian of a particle interaction system with $N$ particles is
given by
$$E(\boldsymbol{x_1}, \ldots, \boldsymbol{x_N}) = \sum\_{\substack{i,j = 1 \\\ i \neq j}}^N K(\boldsymbol{x_i}-\boldsymbol{x_j})$$
where $\boldsymbol{x_1}, \ldots, \boldsymbol{x_N}$ are two-dimensional vectors and $K : \mathbb{R}^2 \to \mathbb{R}$ is the interaction kernel. To determine the ground state of this system, we solve the
following ODE system
$$\frac{d\boldsymbol{x_i}}{dt} = -\frac{1}{N}\sum\_{\substack{i,j = 1 \\\ i \neq j}}^N \nabla K(\boldsymbol{x_i}-\boldsymbol{x_j}).$$
The two types of initial conditions we considered were:

-   particles randomly placed in a rectangle in $\mathbb{R}^2$ $\big([-1,1]\times[-0.5,0.5]\big)$

    -   the files `normalizedGenModel.m` and `parameterizedGenModel.m`
        use the rectangle initial condition

-   particles randomly placed in a ball in $\mathbb{R}^2$ with radius 0.5

    -   the file `parameterizedGenModelBall.m` uses the ball initial
        condition

The types of kernels we considered are as follows.

-   the $\text{“spherical kernel"}$ which uses parameters $p$ and $q$

    -  this kernel is the same as the $c$-norm elliptical kernel where $a = 1$, $b = 1$, and $c = 2$

-   the $\text{“elliptical kernel"}$ which uses parameters $p$, $q$, $a$, $b$

    -   this kernel is the same as the $c$-norm elliptical kernel where $c = 2$

-   the $\text{“elliptical }c\text{-norm kernel"}$ which uses parameters $p$, $q$, $a$, $b$, and $c$

    -   the files `parameterizedSystemGrad.m`, `paramGeneralLcGrad.m`,
        and `paramL1Grad.m` are used to compute the gradient of this
        kernel

-   the $\text{“}\infty\text{-norm kernel"}$ which uses parameters $p$ and $q$

    -   the file `LinfODEsolverSystemGrad.m` is used to compute the
        gradient of this kernel

**Spherical Kernel**

Consider real numbers $p$ and $q$ with $q > p > -2$. The
$\text{“spherical kernel"}$ is defined as
$$K(x,y) = \frac{(x^2 + y^2)^{q/2}}{q} - \frac{(x^2+y^2)^{p/2}}{p}.$$
The gradient of this kernel is

$$
\nabla K(x,y) = \begin{bmatrix} x\left((x^2+y^2)^{q/2-1} - (x^2+y^2)^{p/2-1}\right) \\\ y\left((x^2+y^2)^{q/2-1} - (x^2+y^2)^{p/2-1}\right) \end{bmatrix}.
$$

If $p = 0$ or $q = 0$, we replace 
$$\frac{(x^2+y^2)^{p/2}}{p}\hspace{0.5em}\text{or}\hspace{0.5em}\frac{(x^2+y^2)^{q/2}}{q}\hspace{0.5em}\text{with}\hspace{0.5em} \log\left(\sqrt{x^2+y^2}\right) = \frac{1}{2}\log(x^2+y^2).$$
This gives
$$K\_{p = 0}(x,y) = \frac{(x^2 + y^2)^{q/2}}{q} - \frac{1}{2}\log(x^2+y^2)$$
and
$$K\_{q = 0}(x,y) = \frac{1}{2}\log(x^2+y^2) -  \frac{(x^2+y^2)^{p/2}}{p}.$$
These kernels have gradients

$$
\nabla K\_{p = 0}(x,y) = \begin{bmatrix} x\left(x^2+y^2)^{q/2-1} - \dfrac{1\mathstrut}{x^2+y^2}\right) \\\ y\left(x^2+y^2)^{q/2-1} - \dfrac{1\mathstrut}{x^2+y^2}\right) \end{bmatrix}
$$

and

$$
\nabla K\_{q = 0}(x,y) = \begin{bmatrix} x\left(\dfrac{1\mathstrut}{x^2+y^2} - (x^2+y^2)^{p/2-1}\right) \\\ y\left(\dfrac{1\mathstrut}{x^2+y^2} - (x^2+y^2)^{p/2-1}\right) \end{bmatrix}.
$$

**Elliptical Kernel**

Let $a$ and $b$ be positive real numbers. The $\text{“elliptical kernel"}$ is
defined as
$$K(x,y) = \frac{(a^2x^2+b^2y^2)^{q/2}}{q} - \frac{(a^2x^2+b^2y^2)^{p/2}}{p}.$$
The gradient of this kernel is

$$
\nabla K(x,y) = \begin{bmatrix} a^2x\left((a^2x^2+b^2y^2)^{q/2-1} - (a^2x^2+b^2y^2)^{p/2-1}\right) \\\ b^2y\left((a^2x^2+b^2y^2)^{q/2-1} - (a^2x^2+b^2y^2)^{p/2-1}\right) \end{bmatrix}.
$$

If $p = 0$ or $q = 0$, we replace
$$\frac{(a^2x^2+b^2y^2)^{p/2}}{p}\hspace{0.5em}\text{or}\hspace{0.5em}\frac{(a^2x^2+b^2y^2)^{q/2}}{q}\hspace{0.5em}\text{with}\hspace{0.5em} \log\left(\sqrt{a^2x^2+b^2y^2}\right) = \frac{1}{2}\log(a^2x^2+b^2y^2).$$
This gives
$$K\_{p = 0}(x,y) = \frac{1}{2}\log(a^2x^2+b^2y^2) - \frac{(a^2x^2+b^2y^2)^{p/2}}{p}$$
and
$$K\_{q = 0}(x,y) = \frac{(a^2x^2+b^2y^2)^{q/2}}{q} - \frac{1}{2}\log(a^2x^2+b^2y^2).$$
These kernels have gradients

$$
\nabla K\_{p = 0} = \begin{bmatrix} a^2x\left((a^2x^2+b^2y^2)^{q/2-1} - \dfrac{1\mathstrut}{a^2x^2+b^2y^2}\right) \\\ b^2y\left((a^2x^2+b^2y^2)^{q/2-1} - \dfrac{1\mathstrut}{a^2x^2+b^2y^2}\right) \end{bmatrix}
$$

and

$$
\nabla K\_{q = 0} = \begin{bmatrix} a^2x\left(\dfrac{1\mathstrut}{a^2x^2+b^2y^2} - (a^2x^2+b^2y^2)^{p/2-1}\right) \\\ b^2y\left(\dfrac{1\mathstrut}{a^2x^2+b^2y^2} - (a^2x^2+b^2y^2)^{p/2-1}\right) \end{bmatrix}.
$$

**Elliptical $c$-Norm Kernel**

Let $c$ be a real number with $c \geq 1$, we define the $\text{“elliptical }c\text{-norm kernel"}$ as
$$K(x,y) = \frac{(a^c\|x^c\|+b^c\|y^c\|)^{q/c}}{q} - \frac{(a^c\|x^c\|+b^c\|y^c\|)^{p/c}}{p}.$$
The gradient of this kernel is

$$
\nabla K(x,y) = \begin{bmatrix} \dfrac{a^cx^{2c-1}\mathstrut}{\|x^c\|}\left(\left(a^c\|x^c\|+b^c\|y^c\|\right)^{q/c-1} - \left(a^c\|x^c\| + b^c\|y^c\|\right)^{p/c-1}\right) \\\ \dfrac{b^cy^{2c-1}\mathstrut}{\|y^c\|}\left(\left(a^c\|x^c\|+b^c\|y^c\|\right)^{q/c-1} - \left(a^c\|x^c\| + b^c\|y^c\|\right)^{p/c-1}\right) \end{bmatrix}.
$$

If $p = 0$ or $q = 0$, we replace
$$\frac{(a^c\|x\|^c+b^c\|y\|^c)^{p/c}}{p}\hspace{0.5em}\text{or}\hspace{0.5em}\frac{(a^c\|x\|^c+b^c\|y\|^c)^{q/c}}{q}\hspace{0.5em}\text{with}\hspace{0.5em} \log\left((a^c\|x\|^c+b^c\|y\|^c)^{\frac{1}{c}}\right) = \frac{1}{c}\log(a^c\|x\|^c+b^c\|y\|^c).$$
This gives
$$K\_{q = 0}(x,y) = \frac{1}{c}\log(a^cx^c+b^cy^c) - \frac{(a^c\|x^c\|+b^c\|y^c\|)^{p/c}}{p}$$
and
$$K\_{p = 0}(x,y) = \frac{(a^c\|x^c\|+b^c\|y^c\|)^{q/c}}{q} - \frac{1}{c}\log(a^c\|x\|^c+b^c\|y\|^c).$$
These kernels have gradients

$$
\nabla K\_{q = 0}(x,y) = \begin{bmatrix} \dfrac{a^cx^{2c-1}\mathstrut}{\|x^c\|}\left(\dfrac{1\mathstrut}{a^c\|x^c\|+b^c\|y^c\|} - \left(a^c\|x^c\| + b^c\|y^c\|\right)^{p/c-1}\right) \\\ \dfrac{b^cy^{2c-1}\mathstrut}{\|y^c\|}\left(\dfrac{1\mathstrut}{a^c\|x^c\|+b^c\|y^c\|} - \left(a^c\|x^c\| + b^c\|y^c\|\right)^{p/c-1}\right) \end{bmatrix}
$$

and

$$
\nabla K\_{p = 0}(x,y) =  \begin{bmatrix} \dfrac{a^cx^{2c-1}\mathstrut}{\|x^c\|}\left(\left(a^c\|x^c\|+b^c\|y^c\|\right)^{q/c-1} - \dfrac{1\mathstrut}{a^c\|x^c\|+b^c\|y^c\|}\right) \\\ \dfrac{b^cy^{2c-1}\mathstrut}{\|y^c\|}\left(\left(a^c\|x^c\|+b^c\|y^c\|\right)^{q/c-1} - \dfrac{1\mathstrut}{a^c\|x^c\|+b^c\|y^c\|}\right) \end{bmatrix}.
$$

**$\infty$-Norm Kernel**

The $\text{“}\infty\text{-norm kernel"}$ is
$$K(x,y) = \frac{\max(\|x\|,\|y\|)^q}{q} - \frac{\max(\|x\|,\|y\|)^p}{p}.$$
Note that we do not use the *a* and *b* parameters in this version of
the kernel. Using the fact that
$$\max(x,y) = \frac{1}{2}(x+y) + \frac{1}{2}\|x-y\|,$$
we can write the kernel as
$$K(x,y) = \frac{\displaystyle \left(\frac{1}{2}\bigl(\|x\|+\|y\|\bigr)+\frac{1}{2}\Bigl(\bigl\|\|x\|-\|y\|\bigr\|\Bigr)\right)^q}{\displaystyle q} - \frac{\displaystyle \left(\frac{1}{2}\bigl(\|x\|+\|y\|\bigr)+\frac{1}{2}\Bigl(\bigl\|\|x\|-\|y\|\bigr\|\Bigr)\right)^p}{\displaystyle p}.$$
Using this representation of the kernel, we obtain the gradient

$$
\hspace{-0.34in} \nabla K(x,y) = \begin{bmatrix} \dfrac{x\mathstrut}{2\|x\|}\left(\dfrac{\|x\|-\|y\|\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|}+1\right)\left(\left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{q-1} - \left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{p-1}\right) \\\ \dfrac{y\mathstrut}{2\|y\|}\left(\dfrac{\|x\|-\|y\|\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|}+1\right)\left(\left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{q-1} - \left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{p-1}\right) \end{bmatrix}.
$$

If $p = 0$ or $q = 0$, we replace
$$\frac{\displaystyle \left(\frac{1}{2}\bigl(\|x\|+\|y\|\bigr)+\frac{1}{2}\Bigl(\bigl\|\|x\|-\|y\|\bigr\|\Bigr)\right)^p}{\displaystyle p}\hspace{0.5em}\text{or}\hspace{0.5em}\frac{\displaystyle \left(\frac{1}{2}\bigl(\|x\|+\|y\|\bigr)+\frac{1}{2}\Bigl(\bigl\|\|x\|-\|y\|\bigr\|\Bigr)\right)^q}{\displaystyle q}$$
$$\text{with}\hspace{0.5em}\log\left(\frac{1}{2}(\|x\|+\|y\|)+\frac{1}{2}\left(\bigl\|\|x\|-\|y\|\bigr\|\right)\right).$$
This gives us
$$K\_{p = 0}(x,y) = \frac{\displaystyle \left(\frac{1}{2}\bigl(\|x\|+\|y\|\bigr)+\frac{1}{2}\Bigl(\bigl\|\|x\|-\|y\|\bigr\|\Bigr)\right)^q}{\displaystyle q} - \log\left(\frac{1}{2}(\|x\|+\|y\|)+\frac{1}{2}\left(\bigl\|\|x\|-\|y\|\bigr\|\right)\right)$$
and
$$K\_{q = 0}(x,y) = \log\left(\frac{1}{2}(\|x\|+\|y\|)+\frac{1}{2}\left(\bigl\|\|x\|-\|y\|\bigr\|\right)\right) - \frac{\displaystyle \left(\frac{1}{2}\bigl(\|x\|+\|y\|\bigr)+\frac{1}{2}\Bigl(\bigl\|\|x\|-\|y\|\bigr\|\Bigr)\right)^p}{\displaystyle p}.$$
These kernels have gradients

$$
\nabla K\_{p = 0}(x,y) = \begin{bmatrix} \dfrac{x\mathstrut}{2\|x\|}\left(\dfrac{\|x\|-\|y\|\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|}+1\right)\left(\left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{q-1} - \dfrac{2\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|+\|x\|+\|y\|}\right) \\\ \dfrac{y\mathstrut}{2\|y\|}\left(\dfrac{\|x\|-\|y\|\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|}+1\right)\left(\left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{q-1} - \dfrac{2\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|+\|x\|+\|y\|}\right) \end{bmatrix}
$$

and

$$
\nabla K\_{q = 0}(x,y) = \begin{bmatrix} \dfrac{x\mathstrut}{2\|x\|}\left(\dfrac{\|x\|-\|y\|\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|}+1\right)\left(\dfrac{2\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|+\|x\|+\|y\|} - \left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{p-1}\right) \\\ \dfrac{y\mathstrut}{2\|y\|}\left(\dfrac{\|x\|-\|y\|\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|}+1\right)\left(\dfrac{2\mathstrut}{\bigl\|\|x\|-\|y\|\bigr\|+\|x\|+\|y\|} - \left(\dfrac{\bigl\|\|x\|-\|y\|\bigr\|\mathstrut}{2}+\dfrac{\|x\|+\|y\|\mathstrut}{2}\right)^{p-1}\right) \end{bmatrix}.
$$
