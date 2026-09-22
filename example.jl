"""
Example implementation of the WoS solvers with boundary modifiers.
"""

include("./src/Render.jl")
include("./src/PoissonSolver.jl")

# Set up boundary conditions
circle₁ = Circle(Point(-2.5,  2.0), 1.2)
bc₁(x) = MIXER_FAKE_NOISE(x, GOLD, 5)
modifier₁ = MODIFIER_complement

circle₂ = Circle(Point(2.0, 2.2), 1.0)
bc₂(x) = MIXER_FAKE_NOISE(x, GREEN, 0.2)
modifier₂ = MODIFIER_brightness(1.15) ∘ MODIFIER_channelPermutation([2, 3, 1])

circle₃ = Circle(Point(-2.0, -2.0), 1.4)
bc₃(x) = MIXER_GAUSSIAN_2D(x, RED, Point(0, 0), 2) + 0.5 * MIXER_CONSTANT(x, BLUE)
modifier₃ = MODIFIER_brightness(0.9) ∘ MODIFIER_channelPermutation([1, 3, 2])

circle₄ = Circle(Point(2.5, -1.8), 0.9)
bc₄(x) = MIXER_CHECKER(x, PURPLE, OLIVE, 2)
modifier₄ = MODIFIER_clamp([0.2, 0.2, 0.2], [1.0, 8.0, 1.0]) ∘ MODIFIER_complement

line₁ = line(Point(-3.5, 2.0), Point(3.0, 2.2))
bc_line(x) =MIXER_GAUSSIAN_BLEND(x, ORANGE, LIME, Point(0, 2), 0.3)
modifier_line = MODIFIER_identity

# Make scene
∂𝕊 = boundary(circle₁, bc₁, modifier=MODIFIER_identity) ∪ 
     boundary(circle₂, bc₂, modifier=MODIFIER_identity) ∪ 
     boundary(circle₃, bc₃, modifier=MODIFIER_identity) ∪ 
     boundary(circle₄, bc₄, modifier=MODIFIER_identity) ∪ 
     boundary(line₁, bc_line, modifier=MODIFIER_identity)

∂𝕊_modified = boundary(circle₁, bc₁, modifier=modifier₁) ∪ 
              boundary(circle₂, bc₂, modifier=modifier₂) ∪ 
              boundary(circle₃, bc₃, modifier=modifier₃) ∪ 
              boundary(circle₄, bc₄, modifier=modifier₄) ∪
              boundary(line₁, bc_line, modifier=modifier_line)

# WoS parameters
WoS_depth = 30
num_Samples = 200
ϵ = 0.004

# Rendering domain
Ω = makeDomain(10.0, 10.0)

resolution = (256 * 2, 256 * 2)
Nx, Ny = resolution
♯Ω = discretize(Ω, Nx, Ny)

# Boundary preview
function previewBoundary(∂𝕊::Scene, domain::RenderDomain, resolution::Tuple, ϵ::Float64)
    Nx, Ny = resolution
    ♯domain = discretize(domain, Nx, Ny)
    image = zeros(3, Ny, Nx)
    for i ∈ 1:Nx, j ∈ 1:Ny
        x = ♯domain.grid[i, j]
        candidates = [begin x_nearest = nearestPoint(x, bndry)
                            bndry(x_nearest)
                      end for bndry ∈ ∂𝕊.boundaries if distance(x, bndry) ≤ ϵ]
        if !isempty(candidates)
            image[:, j, i] = candidates[argmax(sum.(candidates))]
        end
    end
    return colorview(RGB, image)
end

previewBoundary(∂𝕊, Ω, resolution, 0.01)

# Localised sources and sinks
A = Point(1.5, 1.0)
B = Point(-1.8, -1.5)
f(x::Point) = 0.5 * MIXER_GAUSSIAN_2D(x, TEAL, A, 1) + 0.5 * MIXER_GAUSSIAN_2D(x, MAGENTA, B, 1.5)

# Source preview
function previewSource(source, ∂𝕊::Scene, domain::RenderDomain, resolution::Tuple, ϵ::Float64)
    f = source
    Nx, Ny = resolution
    ♯domain = discretize(domain, Nx, Ny)
    image = zeros(3, Ny, Nx)
    for i ∈ 1:Nx, j ∈ 1:Ny
        x = ♯domain.grid[i, j]
        candidates = [begin x_nearest = nearestPoint(x, bndry)
                            bndry(x_nearest)
                      end for bndry ∈ ∂𝕊.boundaries if distance(x, bndry) ≤ ϵ]
        if !isempty(candidates)
            image[:, j, i] = candidates[argmax(sum.(candidates))]
        else
            image[:, j, i] = f(x)
        end
    end
    return colorview(RGB, image)
end

previewSource(f, ∂𝕊, Ω, resolution, 0.01)

# Solve the Poisson problem and the modified Poisson problem
u_poisson_identity(p::Point) = solvePoisson(p, ∂𝕊, f, WoS_depth, num_Samples, ϵ)
u_poisson_modified(p::Point) = solvePoisson(p, ∂𝕊_modified, f, WoS_depth, num_Samples, ϵ)

image_poisson_identity = render(u_poisson_identity, ♯Ω)
image_poisson_modified = render(u_poisson_modified, ♯Ω)

# Visualize results
viewImage(image_poisson_identity)
viewImage(image_poisson_modified)