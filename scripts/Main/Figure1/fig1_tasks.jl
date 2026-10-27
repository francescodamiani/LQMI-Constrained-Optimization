# Random square LQMI tasks. Adapted from the submission dimension benchmark.
using LinearAlgebra, Random

function generate_stable_matrix(n; ρ_max=0.9, rng=Random.default_rng())
    A = randn(rng, n, n)
    ρ = maximum(abs, eigvals(A))
    return A * (ρ_max / ρ)
end

function build_random_square_task(n::Int; T::Int=6, ρ_max=0.9,
                                  σ_ξ=0.2, σ_ω=0.2, σ_ρ=0.2, σ_ϵ=0.2, σ_η=0.2,
                                  rng=Random.default_rng())

    dim_x = n
    dim_u = n
    dim_y = n

    A = generate_stable_matrix(n; ρ_max, rng=rng)
    B = generate_stable_matrix(n; ρ_max, rng=rng)
    H = diagm(ones(n))

    D_tmp = abs.(generate_stable_matrix(n; ρ_max, rng=rng))
    D = σ_ρ .* D_tmp
    C = σ_ϵ .* B

    Ω_ω = Diagonal(fill(σ_ω^2, dim_y))
    Ω_ξ = Diagonal(fill(σ_ξ^2, dim_x))
    Ω_η = Diagonal(fill(σ_η^2, dim_x))

    x0 = ones(dim_x)
    z0 = copy(x0)
    x_1_mean = copy(x0)
    z_1_mean = copy(z0)

    Σ_1_x = zeros(dim_x, dim_x)
    Σ_1_z = zeros(dim_x, dim_x)

    Q = zeros(dim_x, dim_x, T)
    R = zeros(dim_u, dim_u, T)
    for t in 1:T
        Q[:, :, t] .= I(dim_x)
        R[:, :, t] .= I(dim_u)
    end
    R[:, :, T] .= zeros(dim_u, dim_u)
    Q[:, :, 1] .= zeros(dim_x, dim_x)

    return (dim_x, dim_u, dim_y, T,
            A, B, H, C, D, Ω_ξ, Ω_ω, Ω_η,
            x0, z0, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q, R)
end

# Helpers for GD baseline (Damiani-style)
function pack_LK(L::Array{Float64,3}, K::Array{Float64,3})
    return vcat(vec(L), vec(K))
end

function unpack_LK(x::Vector{Float64}, dim_u::Int, dim_x::Int, dim_y::Int, T::Int)
    nL = dim_u * dim_x * (T - 1)
    L = reshape(x[1:nL], dim_u, dim_x, T - 1)
    K = reshape(x[nL+1:end], dim_x, dim_y, T - 1)
    return L, K
end


function build_M(A, B, H, L, K, T)
    M = zeros(size(A, 1), size(A, 1), T-1)
    for t in 1:T-1
        M[:, :, t] = A + B * L[:, :, t] - K[:, :, t] * H
    end
    return M
end
