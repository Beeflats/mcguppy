"""
Points
"""

struct Point
    # A Point is an element of Euclidean metric space 𝔼ⁿ
    x::Float64
    y::Float64
end

point(x, y) = Point(x, y)

# Euclidean metric axiom
distance²(P::Point, Q::Point) = (P.x-Q.x)^2 + (P.y-Q.y)^2
distance(P::Point, Q::Point) = √distance²(P, Q)

"""
Vectors
"""
struct Vec
    # A Vector (Vec) is an element of a vector space ℝⁿ
    x::Float64
    y::Float64
end

vector(x, y) = Vec(x, y)

# Vector space axioms
# Additive identity
const ZEROVECTOR = Vec(0, 0)
# Closure, commutativity and associativity of addition
Base.:+(v₁::Vec, v₂::Vec) = Vec(v₁.x + v₂.x, v₁.y + v₂.y)
Base.:-(v₁::Vec, v₂::Vec) = Vec(v₁.x - v₂.x, v₁.y - v₂.y)
# Negation
Base.:-(v::Vec) = Vec(-v.x, -v.y)
# Closure, distributivity and associativity of scalar multiplication
Base.:*(c::Number, v::Vec) = Vec(c * v.x, c * v.y)
Base.:*(v::Vec, c::Number) = Vec(c * v.x, c * v.y)
Base.:/(v::Vec, c::Number) = Vec(v.x / c, v.y / c)

# Displacement axioms:
 # Displacement space is the space (𝔼ⁿ, ℝⁿ, ⊕, →) in which
 # ⊕ [\oplus] and → [\to] are binary operators which satisfying the following
 #   P ⊕ v ∈ 𝔼ⁿ
 #   P → Q ∈ ℝⁿ
 #   P ⊕ (P → Q) = Q
displace(P::Point, v::Vec) = Point(P.x + v.x, P.y + v.y)
⊕(P::Point, v::Vec) = displace(P, v)
# The third axiom can be reformulated as P → (P ⊕ v) = v, therefore
→(P::Point, Q::Point) = Vec(Q.x - P.x, Q.y - P.y)

# Dot product space axioms
⋅(v₁::Vec, v₂::Vec) = v₁.x * v₂.x + v₁.y * v₂.y

# Normed vector space axiom:
#   |P → Q| = d(P, Q)
norm²(v::Vec) = v ⋅ v
length²(v::Vec) = norm²(v) # norm(v) was defined to be equivalent to distance(O, O ⊕ v) where O is an arbitrary point in 𝔼ⁿ 
norm(v::Vec) = sqrt(norm²(v))
unit(v::Vec) = v / norm(v)
normalize(v::Vec) = unit(v)
→ᵘ(P::Point, Q::Point) = unit(P → Q)

# Cross product
×(v₁::Vec, v₂::Vec) = v₁.x * v₂.y - v₁.y * v₂.x

# Orthonormal basis Vectors
const î = Vec(1, 0)
const ĵ = Vec(0, 1)

# Arbitrary perpendicular vector
perp(v::Vec) = Vec(v.y, -v.x)

# Projection
proj(a::Vec, b::Vec) = (a ⋅ b) / length²(b) * b
oproj(a::Vec, b::Vec) = a - proj(a, b)
∥(v₁::Vec, v₂::Vec) = proj(v₁, v₂)
⟂(v₁::Vec, v₂::Vec) = oproj(v₁, v₂)
