include("Vectors.jl")

∅ = nothing
∞ = Inf

"""
Geometry
"""
abstract type Geometry end

struct Ray
    #r(t) = O ⊕ t d
    origin::Point
    direction::Vec
end

ray(o::Point, d::Vec) = Ray(o, d)
(r::Ray)(t) = r.origin ⊕ t * r.direction

struct Line <: Geometry
    """
    A line characterized by a point p and direction vector d.
    The line is the set of points
        L(p,d) = { x ∈ ℝ² | x = p ⊕ t d }
    """
    point::Point
    direction::Vec
end

line(p::Point, d::Vec) = Line(p, d)
line(p₁::Point, p₂::Point) = Line(p₁, p₁ → p₂)
lineByNormal(p::Point, n::Vec) = Line(p, perp(n)) # L = { x ∈ ℝ² | (x → p) ⋅ n = 0 }

struct Circle <: Geometry
    """
    A circle is a line (not an area for our purposes) characterized by a center C
    and radius r.
        S¹(C, r) = { x ∈ ℝ² | (C → x) ⋅ (C → x) = r² }
    """
    center::Point
    radius::Float64
end

circle(C::Point, r::Float64) = Circle(C, r)

struct LineSegment <: Geometry
    """
    A line segment is the set of points between p₁ and p₂.
        L[p₁, p₂] = { q ∈ ℝ² | q = p₁ ⊕ λ(p₁→p₂)} where λ ∈ [0,1]
    It can be considered the the convex hull of two points in ℝ²
    """
    p₁::Point
    p₂::Point
end

lineSegment(p₁::Point, p₂::Point) = LineSegment(p₁, p₂)

# Nearest point queries
function nearestPoint(x::Point, S::Circle)
    C = S.center
    R = S.radius
    if C == x                                             # if x is the centre point...
        return C ⊕ (r * unit(Vec(randn(), randn()))) # perturb x randomly and pick the new nearest point 
    end
    return C ⊕ R * (C →ᵘ x)
end

function nearestPoint(x::Point, L::Line)
    p = L.point
    d = L.direction
    return p ⊕ ((p → x) ∥ d)
end

function nearestPoint(x::Point, L::LineSegment)
    P = L.p₁; Q = L.p₂
    if P == Q
        return P
    end
    PQ = P → Q
    Px = P → x
    λ = (Px ⋅ PQ) / norm²(PQ) # scalar multiple coefficient of Px ∥ PQ
    # clamp to segment
    t = clamp(λ, 0, 1)
    return P ⊕ t * PQ
end

# Distance queries
distance(x::Point, g::Geometry) = distance(x, nearestPoint(x, g))