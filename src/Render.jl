using Images
using ImageView
using ColorSchemes
include("Vectors.jl")

# Render domains
struct RenderDomain
    p₁::Point               # bottom-left corner
    width::Float64          # horizontal extent
    length::Float64         # vertical extent
    horizontalBasis::Vec   # unit vector along the horizontal direction
    verticalBasis::Vec     # unit vector along the vertical direction
end

struct DiscretizedRenderDomain
    domain::RenderDomain
    Nx::Int
    Ny::Int
    grid::Matrix{Point}     # grid[i,j] is the point at horizontal index i and vertical index j.
end

# Domain construction
function makeDomain(p₁::Point, p₂::Point)
    """
    Construct an axis-aligned rectangular domain from its
    bottom-left and top-right corners.
    """
    width = p₂.x - p₁.x
    length = p₂.y - p₁.y
    width > 0 || throw(ArgumentError("p₂ must lie to the right of p₁"))
    length > 0 || throw(ArgumentError("p₂ must lie above p₁"))
    horizontalBasis = î
    verticalBasis = ĵ
    return RenderDomain(p₁, width, length, horizontalBasis, verticalBasis)
end

function makeDomain(p₁::Point, width::Float64, aspectRatio::Float64)
    """
    Construct a rectangular domain from its bottom-left corner,
    width, and aspect ratio.
    The aspect ratio is width / length.
    """
    width > 0 || throw(ArgumentError("width must be positive"))
    aspectRatio > 0 || throw(ArgumentError("aspectRatio must be positive"))
    length = width / aspectRatio
    return RenderDomain(p₁, Float64(width), Float64(length), î, ĵ)
end

function makeDomain(center::Point, width::Float64, aspectRatio::Float64, orientation::Vec)
    """
    Construct a rectangular domain from its centre, width,
    aspect ratio, and orientation.

    orientation is the horizontal direction of the rectangle.
    """
    width > 0 || throw(ArgumentError("width must be positive"))
    aspectRatio > 0 || throw(ArgumentError("aspectRatio must be positive"))
    norm(orientation) > 0 || throw(ArgumentError("orientation must be non-zero"))
    horizontalBasis = unit(orientation)
    verticalBasis = perp(horizontalBasis)
    length = width / aspectRatio
    # Move from the centre to the bottom-left corner.
    p₁ = center ⊕ (-width / 2) * horizontalBasis ⊕ (-length / 2) * verticalBasis
    return RenderDomain(p₁, Float64(width), Float64(length), horizontalBasis, verticalBasis)
end

function makeDomain(width::Float64, length::Float64)
    """
    Construct an axis-aligned rectangular domain centred at the origin.
    """
    width > 0 || throw(ArgumentError("width must be positive"))
    length > 0 || throw(ArgumentError("length must be positive"))
    return makeDomain(Point(0.0, 0.0), Float64(width), Float64(width / length), î)
end

function makeDomain(p₁::Point, p₂::Point, aspectRatio::Float64)
    """
    Construct a rectangular domain from two opposite corners and
    an aspect ratio.

    p₁ is the bottom-left corner.

    The rectangle is oriented so that the vector p₁ → p₂ is
    the diagonal of the rectangle.
    """
    aspectRatio > 0 || throw(ArgumentError("aspectRatio must be positive"))
    diagonal = p₁ → p₂
    diagonalLength = norm(diagonal)
    diagonalLength > 0 || throw(ArgumentError("p₁ and p₂ must be distinct"))

    # If p₁ and p₂ are intended to be opposite corners, aspect ratio determines the side lengths.
    # For a rectangle with width W and length L:
    #     W / L = aspectRatio
    #     W² + L² = diagonal²

    length = diagonalLength / sqrt(aspectRatio^2 + 1)
    width = aspectRatio * length

    horizontalBasis = unit(diagonal)
    verticalBasis = perp(horizontalBasis)

    # p₂ is the opposite corner, so the diagonal must equal
    #     width * horizontalBasis + length * verticalBasis.

    # Choose the orientation whose positive vertical direction points toward p₂.
    candidate = width * horizontalBasis + length * verticalBasis

    if norm(candidate - diagonal) > 1e-10
        verticalBasis = -verticalBasis
    end

    return RenderDomain(p₁, width, length, horizontalBasis, verticalBasis)
end

# Discretisation
"""
Discretise a rendering domain into an Nx × Ny lattice.
"""
function discretize(domain::RenderDomain, Nx::Integer, Ny::Integer)
    Nx ≥ 2 || throw(ArgumentError("Nx must be at least 2"))
    Ny ≥ 2 || throw(ArgumentError("Ny must be at least 2"))
    Δx = domain.width / (Nx - 1)
    Δy = domain.length / (Ny - 1)
    grid = [domain.p₁ ⊕ (i - 1) * Δx * domain.horizontalBasis ⊕ (j - 1) * Δy * domain.verticalBasis
                for i ∈ 1:Nx, j ∈ 1:Ny]
    return DiscretizedRenderDomain(domain, Int(Nx), Int(Ny), grid)
end

# Rendering
function render(u::Function, ♯domain::DiscretizedRenderDomain)
    """
    Evaluate u at every point in the rendering grid.
    """
    image = zeros(3, ♯domain.Ny, ♯domain.Nx) # The grid is indexed as (x, y), while Julia images are indexed as (row, column) = (y, x). Therefore the image has dimensions (3, Ny, Nx)
    for i ∈ 1:♯domain.Nx, j ∈ 1:♯domain.Ny
        image[:, j, i] = u(♯domain.grid[i, j])
    end
    return image
end

function render(u::Function, domain::RenderDomain, Nx::Integer, Ny::Integer)
    ♯domain = discretize(domain, Nx, Ny)
    return render(u, ♯domain)
end

# Visualisation
function viewImage(rendering)
    colorview(RGB, rendering)
end
