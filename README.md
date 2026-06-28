# BurgersFiniteDifferences
Finite difference solver for the 1D viscous Burgers equation in the Wolfram Language (Mathematica), validated against the exact solution obtained using the Cole-Hopf transformation.

# The finite difference method for the viscous Burgers equation

Explicit finite differences on a periodic, validated against the exact solution (using Cole-Hopf transformation)

Ruben Ranval 

### Introduction

The Burgers equation is a simple canonical one-dimensional model in which a convective nonlinearity competes with viscous diffusion. The reason, I'm using it is that it's the simplest equation that keeps the quadratic nonlinearity of Navier Stokes while remaining exactly solvable (which makes it an ideal benchmark for numerical methods).

In this notebook, I will solve the Burgers equation on a periodic domain with a fully explicit FD scheme, and thn validate the result againt the exact solution obtained through the Cole-Hopf transformation.

### The Mathematical problem

#### Governing equation

We want to solve

![0px5eb32p0vj7](img/0px5eb32p0vj7.png)

with periodic BCs and initial conditions ![1r7lcxwu0cc87](img/1r7lcxwu0cc87.png). The advection term is taken in conservative form, ![14cc91jiguwbl](img/14cc91jiguwbl.png), which is the form used for the discretization below. The constant ![1w0psst4ebg8s](img/1w0psst4ebg8s.png) is the kinematic viscosity.

#### Physical behavior

Two mechanisms compete: the nonlinear term that steepens the profile, and for small \nu tends to form a shock, while diffusion smooths gradients and spreads the solution. Their balance is of course measured by the Reynolds number. With \nu = 0.5 (the default value of this notebook) the run is strongly diffusive: the initinal sine flattens and decays without forming a sharp front (as the results of the simulation will confirm).

### Numerical scheme

#### Spatial discretization

Space is discretized on a uniform periodic grid of N₀ points with spacing ![1jyovcrya1reu](img/1jyovcrya1reu.png). Writing uⱼ for the value at node j, both spatial terms use second - order central differences :

- for the advection term: ![0og0c4iyi2f8y](img/0og0c4iyi2f8y.png)

- for the diffusion term: ![169fnrymdm32k](img/169fnrymdm32k.png)

#### Time integration and stability

Time advances with explicit (forward) Euler, giving the classic foward time central space update: ![1rczsp3amwkss](img/1rczsp3amwkss.png).
FTCS is only conditionally stable. A Von Neumann analysis od the linearized equation gives the following amplification factor: ![1sjkjbcljqavx](img/1sjkjbcljqavx.png), with ![0q7vuzurssr0r](img/0q7vuzurssr0r.png), and ![1e467v96t3mov](img/1e467v96t3mov.png) with ![1ico4erijjl0y](img/1ico4erijjl0y.png) the local wave speed.

Stability (![1jj4rvwsih1d7](img/1jj4rvwsih1d7.png) for every mode ![0djlz41mwuw0e](img/0djlz41mwuw0e.png)) requires both ![108dmqf638bll](img/108dmqf638bll.png) and ![0um0wgju6oed4](img/0um0wgju6oed4.png). The step is therefore chosen as ![125o0h2jo5oz8](img/125o0h2jo5oz8.png) : the first term enforces the diffusion limit, the second the advective CFL limit (both with a safety factor).

### The actual implementation!

#### Setup: parameters and grid

We fix the domain, the physical parameters and the grid . N₀ is entered as a real number (256.0), this forces every downstream quantity to machine-precision arithmetic and prevents an exact symbolic ![0dmxeo4otrn8c](img/0dmxeo4otrn8c.png) from being dragged through the time loop, which would make the run explode in symbolic size.

![089oromsdxcmc](img/089oromsdxcmc.png)

Now the two CFL limits we discussed in the previous section (with their safety factors):

![1fg56an2guyvo](img/1fg56an2guyvo.png)

Now let's generate the timesteps and the initial conditions:

![0fl04prze5pd8](img/0fl04prze5pd8.png)

#### The time-stepping "kernel"

*stepBurgers* advances the whole field by one time step. It is written as a pure function of the state vector so it can be iterated cleanly. *uR* and *uL* represent the right/left neighbors via array rotation.

![13zt7tzp7hrgj](img/13zt7tzp7hrgj.png)

#### Time marching

We'll use *NestList* to build each time step iteratively using *stepBurgers* as the function we're applying, with *initU* as the initial conditions and *instants* the total number of timesteps.

```wl
In[]:= sol = NestList[stepBurgers[#, \[Nu], dt, dx] &, initU, instants];
```

The way we've structured our lists, we have *sol[[1, All]]* representing the solution at all 256 points of the 1D grid at the first instant, and *sol[[instants+1, All]]* representing the solution at all 256 points of the 1D grid at the last instant of the simulation.

```wl
In[]:= sol[[1, All]] == initU
```

```wl
Out[]= True
```

Likewise, *sol[[All, 42]]* represents the position of 42nd point of the grid at all timesteps. 

Now, let's look at the results!

### Results

Snapshots of the numerical solution at t = 0, 0.5, 1 and 2. The initial sine smooths and its amplitude decays monotonically: exactly what is expected in this diffusion dominated regime (here ν=0.5).

```wl
In[]:= GraphicsRow[
    Table[
       ListLinePlot[
     Transpose@{xCoords, sol[[Floor[t/dt] + 1]]}, PlotLabel -> "t = " <> ToString[t], 
     PlotRange -> {-1.2, 1.2}, 
     Frame -> True], 
       {t, {0, 0.5, 1.0, 2.0}}], 
   ImageSize -> 800]
```

![0cvn8akchznn9](img/0cvn8akchznn9.png)

The substitution ![0w8mkx5nrp5ub](img/0w8mkx5nrp5ub.png) turns the nonlinear Burgers equation into the linear heat equation ![0wzcr639vx65o](img/0wzcr639vx65o.png). For ![0dhadalq95xuf](img/0dhadalq95xuf.png) the matching heat-equation data is ![120qbjt8tez7c](img/120qbjt8tez7c.png). Expanding this with the generating function of the modified Bessel functions gives Fourier coefficients ![0cd1nzdwuc72p](img/0cd1nzdwuc72p.png), each mode then decays as ![0i0xn68xl3h9n](img/0i0xn68xl3h9n.png). Reconstructing ![0bhktuqlqgvu1](img/0bhktuqlqgvu1.png) from ![0l5luo7eymijo](img/0l5luo7eymijo.png) yields the reference solution below.

![1huf9ffqyjvhs](img/1huf9ffqyjvhs.png)

Subtracting the exact solution from the numerical one gives the pointwise error. It stays at the order of and below across the domain and shrinks in time as the solution smooths.

![1mmf2d0q3cghb](img/1mmf2d0q3cghb.png)

![1p59o838n0x8y](img/1p59o838n0x8y.png)

The exact solution is an infinite series truncated at nMax modes. The plot (shown at ν = 0.4) compares nMax = 1, 2, 3, 10, 20. A single mode is hopeless at early times, when the profile is steepest and richest in high harmonics, but the series converges fast: a handful of modes already coincide visually, and fewer are needed at later times once diffusion has erased the fine structure.

![05t5907sz0fv8](img/05t5907sz0fv8.png)

![1t08xqfk5htcx](img/1t08xqfk5htcx.png)
