# mcguppy

## Theory

mcguppy is a non-physical Poisson equation solver. It is based on the research by [West (2024)](http://cv.rexwe.st/pdf/srfoe.pdf) and [Sawhney (2020)](http://www.rohansawhney.io/mcgp.pdf).

The *rendering equation* 
$$L(\mathbf{x},\mathbf{y}) = L_e(\mathbf{x},\mathbf{y}) + \int_\mathcal{V} f_r(\mathbf{x},\mathbf{y},\mathbf{z}) G(\mathbf{y},\mathbf{z}) L(\mathbf{y},\mathbf{z}) d\mathbf{z}$$
is a single equation that models global illumination and reproduces photorealistic 3D imagery. 

In 2024, West and Mukherjee proposed a *stylized rendering equation* 
$$L(\mathbf{x},\mathbf{y}) = g_{\theta}\left(L_e(\mathbf{x},\mathbf{y}) + \int_\mathcal{V} f_r(\mathbf{x},\mathbf{y},\mathbf{z}) G(\mathbf{y},\mathbf{z}) L(\mathbf{y},\mathbf{z}) d\mathbf{z}\right)$$
where a *stylization function* $g_\theta$ alters the intensity of light recieved by and reflected from a point on a surface. It generalizes various methods for non-photorealistic rendering into one equation.

Models the potential field of an electric charge, mass density, gravity, fluid pressure and stationary state heat conduction share something in common, which is that they can be modeled by Poisson's equation:
$$u(\mathbf{x}) = \frac{1}{|\partial B(\mathbf{x})|}\int_{\partial B(\mathbf{x})} u(\mathbf{y}) d\mathbf{y} - \int_{B(\mathbf{x})} f(\mathbf{y})G(\mathbf{x}, \mathbf{y}) d\mathbf{y}.$$

The Poisson equation bears some commonalities to the rendering equation in which they are both integral equations which reference themselves in the integrand. Therefore, they can be solved via recursive Monte Carlo methods. 

The purpose of this repository is to inject *modifier functions* into the Poisson equation, similar to how stylization functions alter the rendering equation, to realize "stylized" solutions to the Poisson equation.

$$u(\mathbf{x}) = g_\theta\left(\frac{1}{|\partial B(\mathbf{x})|}\int_{\partial B(\mathbf{x})} u(\mathbf{y}) d\mathbf{y} - \int_{B(\mathbf{x})} f(\mathbf{y})G(\mathbf{x}, \mathbf{y}) d\mathbf{y}\right).$$

In this repository, we have opted for $u$ to be a color field in a 2D domain.

## Results and Discussion
Without any modification to the original Poisson equation solver i.e. $g_\theta(c) = c$, the solutions look nearly harmonic except at solid boundaries. However, the solutions to the modified Poisson equation has some clear discontinuities within the domain. 

The discontinuities arise due to the nearest neighbour query of the Walk on Spheres algorithm, where the hard ridges on the color field define the midpoint of two solid boundaries.

A next step for this project is to replace the nearest neighbour query with $k$-nearest neighbours such that each solid boundary has a weighted contribution to the evaluation of the colour field.

## Implementation guide

### Create solid boundaries and define their boundary conditions, and style.
```
circle₁ = Circle(Point(-2.5,  2.0), 1.2)
bc₁(x) = MIXER_FAKE_NOISE(x, GOLD, 5)
modifier₁ = MODIFIER_identity

circle₂ = Circle(Point(2.5, -1.8), 0.9)
bc₂(x) = MIXER_CHECKER(x, PURPLE, OLIVE, 2)
modifier₂ = MODIFIER_clamp([0.2, 0.2, 0.2], [1.0, 8.0, 1.0]) ∘ MODIFIER_complement

C₁ = boundary₁(circle₁, bc₁, modifier=modifier₁)
C₂ = boundary₂(circle₂, bc₂, modifier=modifier₂)
```

Sample scalar fields and color palettes are described in `BoundaryConditions.jl`. 

Modifiers have to properties of being able to sum and multiply with each other among other operations outlined in `Modifier.jl`.

### Define a scene
Scenes are the union of boundary objects.
```
∂𝕊 = C₁ ∪ C₂
```

### Define the rendering domain
A continuous domain is to be discretised into a lattice
```
Ω = makeDomain(10.0, 10.0) # (x, y) ∈ [-5, 5] × [-5, 5]
resolution = (256 * 2, 256 * 2)
Nx, Ny = resolution
♯Ω = discretize(Ω, Nx, Ny)
```
### Define forcing function
The same color field templates in `BoundaryConditions.jl` can be used to define sources and sinks.
```
A = Point(1.5, 1.0)
B = Point(-1.8, -1.5)
f(x::Point) = 0.5 * MIXER_GAUSSIAN_2D(x, TEAL, A, 1) + 0.5 * MIXER_GAUSSIAN_2D(x, MAGENTA, B, 1.5)
```

### Solve the Poisson problem and the modified Poisson problem
Define the parameters for the solver
```
WoS_depth = 30
num_Samples = 200
ϵ = 0.004
```

Define the color field whose solution solves the modified Poisson equation.
```
u_poisson_modified(p::Point) = solvePoisson(p, ∂𝕊, f, WoS_depth, num_Samples, ϵ)
```

Evaluate the values of the color field at every point in the discretised domain.
```
image_poisson_modified = render(u_poisson_modified, ♯Ω)
```
### Visualize results
Display the colors of the color field in the defined domain.
```
viewImage(image_poisson_modified)
```
