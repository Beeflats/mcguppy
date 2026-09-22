include("Geometry.jl")

struct Modifier
    modificationFunction
end

struct NearestNeighborData
    x_nearest::Point
    xᵢ::Point
    xᵢ₊₁::Point
    geo::Geometry
    yᵢ::Point
end

get_normal(c::Circle, x::Point, x∂C::Point) = x → x∂C
get_normal(l::LineSegment, x::Point, x∂l::Point) = perp()

# Geometric helper functions
dx(Θ::NearestNeighborData) = Θ.xᵢ → Θ.xᵢ₊₁
ν(Θ::NearestNeighborData)  = Θ.xᵢ → Θ.yᵢ
n̂(Θ::NearestNeighborData) = get_normal(Θ.geo, Θ.xᵢ, Θ.x_nearest)
boundaryDistance(Θ::NearestNeighborData) = distance(Θ.xᵢ, Θ.x_nearest)
stepDistance(Θ::NearestNeighborData) = norm(dx(Θ))
sourceDistance(Θ::NearestNeighborData) = norm(ν(Θ))

# Colour helper function
function getIntensity(color)
    return 0.299 * color[1] + 0.587 * color[2] + 0.114 * color[3]
end

# Modifier algebra
function Base.:+(θ₁::Modifier, θ₂::Modifier)
    g₁ = θ₁.modificationFunction
    g₂ = θ₂.modificationFunction
    function g₁plusg₂(color, data::NearestNeighborData)
        return g₁(color, data) + g₂(color, data)
    end
    θ₁plusθ₂ = Modifier(g₁plusg₂)
    return θ₁plusθ₂
end

function Base.:-(θ₁::Modifier, θ₂::Modifier)
    g₁ = θ₁.modificationFunction
    g₂ = θ₂.modificationFunction
    function g₁minusg₂(color, data::NearestNeighborData)
        return g₁(color, data) - g₂(color, data)
    end
    θ₁minusθ₂ = Modifier(g₁minusg₂)
    return θ₁minusθ₂
end

function Base.:-(θ::Modifier)
    g = θ.modificationFunction
    function g_neg(color, data::NearestNeighborData)
        return -g(color, data)
    end
    negΘ = Modifier(g_neg)
    return negΘ
end

function ∘(θ₁::Modifier, θ₂::Modifier)
    g₁ = θ₁.modificationFunction
    g₂ = θ₂.modificationFunction
    function g₁composeg₂(color, data::NearestNeighborData)
        return g₁(g₂(color, data), data)
    end
    Θ₁composeΘ₂ = Modifier(g₁composeg₂)
    return Θ₁composeΘ₂
end

function Base.:*(θ₁::Modifier, θ₂::Modifier)
    g₁ = θ₁.modificationFunction
    g₂ = θ₂.modificationFunction
    function g₁timesg₂(color, data::NearestNeighborData)
        return g₁(color, data) .* g₂(color, data)
    end
    θ₁timesθ₂ = Modifier(g₁timesg₂)
    return θ₁timesθ₂
end

function Base.:*(c, θ::Modifier)
    g = θ.modificationFunction
    function g_scale(color, data::NearestNeighborData)
        return c .* g(color, data)
    end
    Θ_scale = Modifier(g_scale)
    return Θ_scale
end

# Basic color modifiers
function g_identity(color, data::NearestNeighborData)
    return color
end
MODIFIER_identity = Modifier(g_identity)

function g_complement(color, data::NearestNeighborData)
    return [1, 1, 1] .- color
end
MODIFIER_complement = Modifier(g_complement)

function g_greyscale(color, data::NearestNeighborData)
    intensity = getIntensity(color)
    return intensity .* [1, 1, 1]
end
MODIFIER_greyscale = Modifier(g_greyscale)

function make_g_fill(fillColor)
    function g_fill(color, data::NearestNeighborData)
        return fillColor
    end
    return g_fill
end
MODIFIER_fill(fillColor) = Modifier(make_g_fill(fillColor))

# data-free transformations
function make_g_brightness(k)
    function g_brightness(color, data::NearestNeighborData)
        return k .* color
    end
    return g_brightness
end
MODIFIER_brightness(k) = Modifier(make_g_brightness(k))

function make_g_addColor(offset)
    function g_addColor(color, data::NearestNeighborData)
        return color .+ offset
    end
    return g_addColor
end
MODIFIER_addColor(offset) = Modifier(make_g_addColor(offset))

function make_g_clamp(minColor, maxColor)
    function g_clamp(color, data::NearestNeighborData)
        return clamp.(color, minColor, maxColor)
    end
    return g_clamp
end
MODIFIER_clamp(minColor, maxColor) = Modifier(make_g_clamp(minColor, maxColor))

# Channel modifiers
function g_swapRG(color, data::NearestNeighborData)
    return [color[2], color[1], color[3]]
end
MODIFIER_swapRG = Modifier(g_swapRG)

function g_swapRB(color, data::NearestNeighborData)
    return [color[3], color[2], color[1]]
end
MODIFIER_swapRB = Modifier(g_swapRB)

function g_swapGB(color, data::NearestNeighborData)
    return [color[1], color[3], color[2]]
end
MODIFIER_swapGB = Modifier(g_swapGB)

# Channel extraction
function g_red(color, data::NearestNeighborData)
    return [color[1], 0, 0]
end
MODIFIER_red = Modifier(g_red)

function g_green(color, data::NearestNeighborData)
    return [0, color[2], 0]
end
MODIFIER_green = Modifier(g_green)

function g_blue(color, data::NearestNeighborData)
    return [0, 0, color[3]]
end
MODIFIER_blue = Modifier(g_blue)

function make_g_channelPermutation(order)
    function g_channelPermutation(color, data::NearestNeighborData)
        return color[order]
    end
    return g_channelPermutation
end
MODIFIER_channelPermutation(order) = Modifier(make_g_channelPermutation(order))

# Threshold modifier
function make_g_threshold(threshold, lowColor, highColor)
    function g_threshold(color, data::NearestNeighborData)
        intensity = getIntensity(color)
        if intensity < threshold
            return lowColor
        else
            return highColor
        end
    end
    return g_threshold
end
MODIFIER_threshold(threshold, lowColor, highColor) = Modifier(make_g_threshold(threshold, lowColor, highColor))

function make_g_toonShade(intensityIntervals, toonColors)
    length(toonColors) == length(intensityIntervals) + 1 || throw(ArgumentError("toonColors must contain one more colour than intensityIntervals"))
    function g_toonShade(color, data::NearestNeighborData)
        intensity = getIntensity(color)
        for colorIndex in eachindex(intensityIntervals)
            if intensity < intensityIntervals[colorIndex]
                return toonColors[colorIndex]
            end
        end
        return toonColors[end]
    end
    return g_toonShade
end
MODIFIER_toonShade(intensityIntervals, toonColors) = Modifier(make_g_toonShade(intensityIntervals, toonColors))

# Spatial modifiers
function make_g_boundaryDistance(f::Function)
    function g_boundaryDistance(color, data::NearestNeighborData)
        d = boundaryDistance(data)
        return f(color, d, data)
    end
    return g_boundaryDistance
end
MODIFIER_boundaryDistance(f::Function) = Modifier(make_g_boundaryDistance(f))

function make_g_stepDistance(f::Function)
    function g_stepDistance(color, data::NearestNeighborData)
        d = stepDistance(data)
        return f(color, d, data)
    end
    return g_stepDistance
end
MODIFIER_stepDistance(f::Function) = Modifier(make_g_stepDistance(f))

# Walk-dependent modifiers
function make_g_direction(f::Function)
    function g_direction(color, data::NearestNeighborData)
        direction = dx(data)
        return f(color, direction, data)
    end
    return g_direction
end
MODIFIER_direction(f::Function) = Modifier(make_g_direction(f))

function make_g_directionWeight(v)
    v̂ = normalize(v)
    function g_directionWeight(color, data::NearestNeighborData)
        d̂ = normalize(dx(data))
        weight = dot(d̂, v̂)
        return weight .* color
    end
    return g_directionWeight
end
MODIFIER_directionWeight(v) = Modifier(make_g_directionWeight(v))

function make_g_directionWeightPositive(v)
    v̂ = normalize(v)
    function g_directionWeightPositive(color, data::NearestNeighborData)
        d̂ = normalize(dx(data))
        weight = max(0, dot(d̂, v̂))
        return weight .* color
    end
    return g_directionWeightPositive
end
MODIFIER_directionWeightPositive(v) = Modifier(make_g_directionWeightPositive(v))

function g_normalMap(color, data::NearestNeighborData)
    normal = n̂(data)
    return [normal.x,normal.y, 0]
end
MODIFIER_normalMap = Modifier(g_normalMap)

function g_normalMap01(color, data::NearestNeighborData)
    normal = n̂(data)
    return [(normal.x + 1) / 2, (normal.y + 1) / 2, 0]
end
MODIFIER_normalMap01 = Modifier(g_normalMap01)

function make_g_normalWeight(v)
    v̂ = normalize(v)
    function g_normalWeight(color, data::NearestNeighborData)
        normal = n̂(data)
        weight = dot(normal, v̂)
        return weight .* color
    end
    return g_normalWeight
end
MODIFIER_normalWeight(v) = Modifier(make_g_normalWeight(v))

function make_g_normalWeightPositive(v)
    v̂ = normalize(v)
    function g_normalWeightPositive(color, data::NearestNeighborData)
        normal = n̂(data)
        weight = max(0, dot(normal, v̂))
        return weight .* color
    end
    return g_normalWeightPositive
end
MODIFIER_normalWeightPositive(v) = Modifier(make_g_normalWeightPositive(v))

# Source-sample modifiers
function make_g_sourceDistance(f::Function)
    function g_sourceDistance(color, data::NearestNeighborData)
        d = sourceDistance(data)
        return f(color, d, data)
    end
    return g_sourceDistance
end
MODIFIER_sourceDistance(f::Function) = Modifier(make_g_sourceDistance(f))

function make_g_sourceDirection(f::Function)
    function g_sourceDirection(color, data::NearestNeighborData)
        direction = ν(data)
        return f(color, direction, data)
    end
    return g_sourceDirection
end
MODIFIER_sourceDirection(f::Function) = Modifier(make_g_sourceDirection(f))

function make_g_stripes_x(Δy, stripeColors)
    numColors = length(stripeColors)
    function g_stripes(color, data::NearestNeighborData)
        y = data.xᵢ.y
        colorIndex = mod1(Int(ceil(y / Δy)), numColors)
        return stripeColors[colorIndex]
    end
    return g_stripes
end
MODIFIER_stripes_x(Δy, stripeColors) = Modifier(make_g_stripes_x(Δy, stripeColors))

function make_g_stripes_y(Δx, stripeColors)
    numColors = length(stripeColors)
    function g_stripes(color, data::NearestNeighborData)
        x = data.xᵢ.x
        colorIndex = mod1(Int(ceil(x / Δx)), numColors)
        return stripeColors[colorIndex]
    end
    return g_stripes
end
MODIFIER_stripes_y(Δx, stripeColors) = Modifier(make_g_stripes_y(Δx, stripeColors))

