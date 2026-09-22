include("./Vectors.jl")
include("./Geometry.jl")
include("./BoundaryCondition.jl")

# Solver utility functions
mean(x) = sum(x) / length(x)

# Sampling
macro sample(expr)
    """
    Allows one to replace 
        `x = U([0,1])`
            where U([0,1]) = rand()
     with
        `@sample x ~ U([0,1])`
     which is more faithful to notation in probability theory
    """
    @assert expr.head == :call && expr.args[1] == :~
    var  = expr.args[2]
    dist = expr.args[3]
    return :($(esc(var)) = $(esc(dist)))
end

Unif_S²() = unit(vector(randn(), randn())) # Unif(∂B) or Unif(S¹)

function Unif_S²(x::Point, R::Float64)
    dx = R * Unif_S²()
    return x ⊕ dx                # y ~ Unif(B(x, R))
end

function Unif_B²(x::Point, R::Float64)
    ρ = R * sqrt(rand())              # ρ = R√U where U ~ Unif([0,1])
    d = Unif_S²()      # d ~ Unif(S¹)
    return x ⊕ ρ * d                 # y ~ Unif(B(x, R))
end

# Convolution kernels
# Harmonic Green's function of the Laplace operator on a ball or radius R
G(x::Point, y::Point, R::Float64) = 1/(2π) * log(R/distance(x,y))
∇G(x::Point, y::Point, R::Float64) = (x → y)/(2π) * (1/distance(x,y)^2 - 1/R^2)

"""
Poisson Solver





"""

function walkOnSphere(xᵢ::Point, ∂Ω::Scene, f::Function, WoS_depth::Integer, ϵ::Float64)
    WoS_depth ≥ 0 || throw(ArgumentError("WoS_depth must be non-negative"))
    x_nearest, color = nearest(xᵢ, ∂Ω)
    r = distance(xᵢ, x_nearest)
    if WoS_depth == 0 || r ≤ ϵ
        return color(x_nearest)
    else
        @sample xᵢ₊₁ ~ Unif_S²(xᵢ, r) # Unif(∂B²(xᵢ, r))
        boundary_contribution = walkOnSphere(xᵢ₊₁, ∂Ω, f, WoS_depth - 1, ϵ)
        
        @sample y ~ Unif_B²(xᵢ, r)
        μ∂B = 2π * r
        source_contribution = -μ∂B * f(y) * G(xᵢ, y, r)

        WoS_color = boundary_contribution + source_contribution

        g_θ = color.modifier.modificationFunction
        nearestNeighborData = NearestNeighborData(x_nearest, xᵢ, xᵢ₊₁, color.geometry, y)
        
        return g_θ(WoS_color, nearestNeighborData) 
    end
end

function solvePoisson(x::Point, Ω::Scene, f::Function, WoS_depth::Integer, num_Samples::Integer, ϵ::Float64)
    num_Samples > 0 || throw(ArgumentError("num_Samples must be positive"))
    return mean([walkOnSphere(x, Ω, f, WoS_depth, ϵ) for _ in 1:num_Samples])
end
