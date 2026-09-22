include("Geometry.jl")
include("Modifier.jl")

# Boundary-condition objects
struct Boundary
    geometry::Geometry
    u::Function # color field
    modifier::Modifier
end

boundary(geometry::Geometry, u::Function; modifier = MODIFIER_identity) = Boundary(geometry, u, modifier)
boundary(geometry::Geometry, color; modifier = MODIFIER_identity) = Boundary(geometry, x -> color, modifier)

(o::Boundary)(x::Point) = o.u(x) # Evaluate the boundary condition.

# Convenience constructors
Dirichlet(geometry::Geometry, color; modifier = MODIFIER_identity) = boundary(geometry, color, modifier)
Dirichlet(geometry::Geometry, u::Function; modifier = MODIFIER_identity) = boundary(geometry, u, modifier)

# Scene
struct Scene
    boundaries::Vector{Boundary}
end

size(𝕊::Scene) = length(𝕊.boundaries)

scene2D() = Scene(Boundary[])
scene(o::Boundary...) = Scene(Boundary[o...])
scene(boundaries::Tuple...) = scene((boundary(geometry, u, modifier) for (geometry, u, modifier) in boundaries)...)

scene(mesh::AbstractVector{<:Geometry}, u::Function; modifier = MODIFIER_identity) = scene((boundary(geometry, u, modifier) for geometry ∈ mesh)...)
scene(mesh::AbstractVector{<:Geometry}, color; modifier = MODIFIER_identity) = scene(mesh, x -> color, modifier)

# Adding boundaries to a scene
Base.:∪(𝕊::Scene, o::Boundary) = Scene([𝕊.boundaries..., o])
Base.:∪(𝕊::Scene, boundaries::Vector{Boundary}) = Scene([𝕊.boundaries... boundaries...])
Base.:∪(𝕊₁::Scene, 𝕊₂::Scene) = Scene([𝕊₁.boundaries..., 𝕊₂.boundaries...])
Base.:∪(o::Boundary, 𝕊::Scene) = 𝕊 ∪ o
Base.:∪(o₁::Boundary, o₂::Boundary) = scene(o₁, o₂)

# Distance to boundaries
distance(x::Point, bndry::Boundary) = distance(x, bndry.geometry)
nearestPoint(x::Point, bndry::Boundary) = nearestPoint(x, bndry.geometry)

# Nearest boundary and point
function nearest(x::Point, 𝕊::Scene)
    isempty(𝕊.boundaries) && throw(ArgumentError("Cannot find nearest boundary in an empty Scene"))
    distances = [distance(x, bndry) for bndry ∈ 𝕊.boundaries]
    bndry_nearest = 𝕊.boundaries[argmin(distances)]
    x_nearest = nearestPoint(x, bndry_nearest)
    return x_nearest, bndry_nearest
end

# Color palette
RED     = [1.0, 0.0, 0.0]
GREEN   = [0.0, 1.0, 0.0]
BLUE    = [0.0, 0.0, 1.0]
CYAN    = [0.0, 1.0, 1.0]
MAGENTA = [1.0, 0.0, 1.0]
YELLOW  = [1.0, 1.0, 0.0]
ORANGE  = [1.0, 0.5, 0.0]
PURPLE  = [0.5, 0.0, 0.5]
LIME    = [0.7, 1.0, 0.0]
TEAL    = [0.0, 0.5, 0.5]
PINK    = [1.0, 0.4, 0.7]
BROWN   = [0.5, 0.25, 0.1]
OLIVE   = [0.4, 0.5, 0.0]
NAVY    = [0.0, 0.0, 0.4]
GOLD    = [1.0, 0.84, 0.0]
SKYBLUE = [0.4, 0.7, 1.0]
VIOLET  = [0.6, 0.3, 1.0]
TURQUOISE = [0.0, 0.9, 0.7]
MAROON  = [0.5, 0.0, 0.0]

# Colour field samples
smoothstep(x) = x*x*(3 - 2x)

# Patterned
MIXER_CONSTANT(p::Point, c) = c
MIXER_STRIPES_x(p::Point, c₁, c₂, λ) = mod(p.y * λ, 1) < 0.5 ? c₁ : c₂
MIXER_STRIPES_Y(p::Point, c₁, c₂, λ) = mod(p.x * λ, 1) < 0.5 ? c₁ : c₂
MIXER_GRID(p::Point, c₁, c₂, λ) = ((mod(p.x*λ, 1) < 0.5) ⊻ (mod(p.y*λ, 1) < 0.5)) ? c₁ : c₂
MIXER_CHECKER(p::Point, c₁, c₂, λ) = (floor(p.x*λ) + floor(p.y*λ)) % 2 == 0 ? c₁ : c₂
MIXER_CELL_NOISE(p::Point, c₁, c₂, λ) = (floor(p.x * λ) + floor(p.y * λ)) % 3 == 0 ? c₁ : c₂
MIXER_FAKE_NOISE(p::Point, c, λ) = c .* (sin(λ * p.x) * sin(λ * p.y))

# Radial color fields
MIXER_GAUSSIAN_2D(p::Point, c, μ::Point, σ) = c .* exp(-norm²(p → μ) / (2σ^2))
MIXER_LORENTZIAN(p::Point, c, μ::Point, γ) = c ./ (1 + norm²(p → μ) / γ^2)
MIXER_MEXICAN_HAT(p::Point, c, μ::Point, σ) = c .* (1 - norm²(p → μ)) .* exp(-norm²(p → μ) / (2σ^2))
MIXER_RADIAL_SINE(p::Point, c, μ::Point, λ) = c .* sin(λ * norm²(p → μ))
function MIXER_GAUSSIAN_RING(p::Point, c, μ::Point, σ, r₀)
    r = norm(p → μ)
    return c .* exp(-(r - r₀)^2 / (2σ^2))
end

# Radial blends
function MIXER_GAUSSIAN_BLEND(p::Point, c₁, c₂, μ::Point, σ)
    t = exp(-norm²(p → μ) / (2σ^2))
    return (1 - t) .* c₁ .+ t .* c₂
end

function MIXER_RADIAL_SMOOTH(p::Point, c₁, c₂, μ::Point, λ)
    r = norm(p → μ)
    t = smoothstep(sin(λ * r))
    return (1 - t) .* c₁ .+ t .* c₂
end

# Directional Gaussians
MIXER_GAUSSIAN_DIRECTIONAL(p::Point, c, μ::Point, σ, v) = c .* exp(-(((p → μ) ⋅ v)^2) / (2σ^2))
MIXER_GAUSSIAN_BAND(p::Point, c, μ::Point, σ, v) = c .* exp(-(((p → μ) ⋅ perp(v))^2) / (2σ^2))

# Gabor
function MIXER_GABOR(p::Point, c, μ::Point, σ, λ, θ)
    d = p → μ
    env = exp(-norm²(d) / (2σ^2))
    proj = d.x * cos(θ) + d.y * sin(θ)
    return c .* env .* cos(λ * proj)
end

function MIXER_GABOR_SINE(p::Point, c, μ::Point, σ, λ, θ)
    d = p → μ
    env = exp(-norm²(d) / (2σ^2))
    proj = d.x * cos(θ) + d.y * sin(θ)
    return (env * sin(λ * proj)) .* c
end

function MIXER_GABOR_SWIRL(p::Point, c, μ::Point, σ, λ)
    d = p → μ
    θ = atan(d.y, d.x)
    env = exp(-norm²(d) / (2σ^2))
    return (env * sin(λ * θ)) .* c
end

# Polar
function MIXER_POLAR_CHECKER(p::Point, c₁, c₂, μ::Point, λ₁, λ₂)
    d = p → μ
    r = norm(d)
    θ = atan(d.y, d.x)
    return ((mod(r * λ₁, 1) < 0.5) ⊻ (mod(θ * λ₂, 1) < 0.5)) ? c₁ : c₂
end

# Miscellaneous
MIXER_INTERFERENCE(p::Point, c₁, c₂, λ₁, λ₂) = c₁ .* (sin(λ₁ * p.x) + sin(λ₂ * p.y)) .+ c₂ .* cos(λ₂ * p.y)

# RGB / channel mixer
MIXER_HUE_SHIFT(p::Point, c, λ) = [sin(λ * p.x) * c[1], sin(λ * p.y) * c[2], sin(λ * (p.x + p.y)) * c[3]]
MIXER_CHANNEL_SWAP(p::Point, c, λ) = [c[2] * sin(λ * p.x), c[3] * sin(λ * p.y), c[1] * sin(λ * (p.x + p.y))]
MIXER_RGB_WARP(p::Point,c, λ₁, λ₂) = [sin(λ₁ * p.x) * c[1], cos(λ₂ * p.y) * c[2], sin(λ₁ * p.y + λ₂ * p.x) * c[3]]

# HASH NOISE
function MIXER_HASH(c::Point, μ::Point, λ)
    d = p → μ
    h = sin(λ * (d.x * 12.9898 + d.y * 78.233)) * 43758.5453
    return c .* (h - floor(h))
end

# THIN-PLATE RBF
function MIXER_RBF_THIN_PLATE(p::Point, c, μ::Point)
    r² = norm²(p → μ)
    return r² == 0 ? zero(c) : c .* (r² * log(sqrt(r²)))
end
