module Scattering_models

using Pkg
Pkg.add(["Revise", "Cubature", "HCubature"])
using Revise
using Cubature: hquadrature
using HCubature: hcubature
# NOTE: Unoptimal default settings:
# using Roots: find_zeros
(@isdefined Root_finding_utils) == false ? include("root_finding_utils.jl") : nothing
using .Root_finding_utils: multiple_root_bisection
using Printf


"""
Ray path equation (17) from Clegg+98 for points source (β_s = 0.0).

# Arguments
- `u::AbstractFloat`: dimensionless observer plane coordinate.
- `α::AbstractFloat`: lens strength.
- `u_::AbstractFloat`: dimensionless lens plane coordinate.
"""
function ray_path(u::T, α::T, u_::T) where {T<:AbstractFloat}
    return u*(1.0 + α*exp(-u^2)) - u_
end


"""
Ray path equation (17) from Clegg+98 for extended source (β_s > 0).

# Arguments
- `u::AbstractFloat`: dimensionless observer plane coordinate.
- `α::AbstractFloat`: lens strength.
- `u_::AbstractFloat`: dimensionless lens plane coordinate.
- `β_s_i::AbstractFloat`: ratio of the angular separation of the point source to angular width of the lens.
"""
function ray_path(u::T, α::T, u_::T, β_s_i::T) where {T<:AbstractFloat}
    return u*(1.0 + α*exp(-u^2)) - u_ - β_s_i 
end


"""
Ray path equation (A11) from Dong+2018.

# Arguments
- `u::AbstractFloat`: dimensionless observer plane coordinate.
- `α::AbstractFloat`: lens strength.
- `u_::AbstractFloat`: dimensionless lens plane coordinate.
- `β_s_i::AbstractFloat`: ratio of the angular separation of the point source to angular width of the lens.
- `ϕ::AbstractFloat`: azimuthal deflection angle.
"""
function ray_path(u::T, α::T, u_::T, β_s_i::T, ϕ::T) where {T<:AbstractFloat}
    γ = sqrt(u_^2 + 2.0*cos(ϕ)*u_*β_s_i + β_s_i^2)
    return u*(1.0 + α*exp(-u^2)) - γ
end


"""
Sum in (20) from Cleagg+98 for point source of unit intensity, i.e. just magnifications.

# Agruments
- `u_::AbstractFloat`: coordinate in the lens plane.
- `α::AbstractFloat`: lens strength.
- `u_min::AbstractFloat`: mininal value of lens plane coordinates (u_) to search for the rays.
- `u_max::AbstractFloat`: maximal value of lens plane coordinates (u_) to search for the rays.
"""
function sum_of_magnifications(u_::T, α::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    I_m = 0.0
    # For given coordinate in the observer plane find all rays in the lens plane
    # u_all = find_zeros(u -> ray_path(u, α, u_), u_min, u_max, no_pts=ceil(Int, (u_max - u_min)/1e-2)) #, xatol=1e-2)
    u_all = multiple_root_bisection(u -> ray_path(u, α, u_), u_min, u_max, 1e-2)
    G_k = 0.0
    G_k_sum = 0.0
    for u_k in u_all 
        # Gain factor for given u_k: (18) in Clegg+98
        G_k = 1.0 /(1.0 + (1.0 - 2u_k^2)*α*exp(-u_k^2))
        G_k_sum += abs(G_k)
    end
	I_m += G_k_sum
    return I_m
end


"""
Expression inside integral in (19) from Clegg+98 for extended source.

# Agruments
- `u_::AbstractFloat`: coordinate in the lens plane.
- `β_s_i::AbstractFloat`: sub-angle of β_s.
- `σ::AbstractFloat`: Width of the source.
- `α::AbstractFloat`: lens strength.
- `u_min::AbstractFloat`: mininal value of lens plane coordinates (u_) to search for the rays.
- `u_max::AbstractFloat`: maximal value of lens plane coordinates (u_) to search for the rays.
"""
function integrand(u_::T, β_s_i::T, α::T, σ::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    I_m = 0.0
    # For given coordinate in the observer plane find all rays in the lens plane
    # u_all = find_zeros(u -> ray_path(u, α, u_, β_s_i), u_min, u_max, no_pts=ceil(Int, (u_max - u_min)/1e-2)) #, xatol=1e-2)
    u_all = multiple_root_bisection(u -> ray_path(u, α, u_, β_s_i), u_min, u_max, 1e-2)
    G_k_sum = 0.0
    # For each ray coming to the specified observer plane position from the
    # sub-angle ``β_s_i`` find the corresponding gain
    for u_k in u_all 
        # Gain factor for given u_k: (18) in Clegg+98
        G_k_sum += abs(1.0 /(1.0 + (1.0 - 2u_k^2)*α*exp(-u_k^2)))
    end
    # Brightness distribution
    B = (1.0/(sqrt(2π)*σ))exp(-β_s_i^2/(2.0*σ^2));
	I_m += B*G_k_sum;
    return I_m
end


"""
Expression inside integral in (A13) from Dong+2018.

# Agruments
- `u_::AbstractFloat`: coordinate in the lens plane.
- `β_s_i::AbstractFloat`: sub-angle of β_s.
- `σ::AbstractFloat`: Width of the source.
- `α::AbstractFloat`: lens strength.
- `u_min::AbstractFloat`: mininal value of lens plane coordinates (u_) to search for the rays.
- `u_max::AbstractFloat`: maximal value of lens plane coordinates (u_) to search for the rays.
"""
function integrand(u_::T, β_s_i::T, ϕ::T, α::T, σ::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    I_m = 0.0
    # For given coordinate in the observer plane find all rays in the lens plane
    # u_all = find_zeros(u -> ray_path(u, α, u_, β_s_i, ϕ), u_min, u_max, no_pts=ceil(Int, (u_max - u_min)/1e-2)) #, xatol=1e-2)
    u_all = multiple_root_bisection(u -> ray_path(u, α, u_, β_s_i, ϕ), u_min, u_max, 1e-2)
    G_k_sum = 0.0
    # For each ray coming to the specified observer plane position from the
    # sub-angle ``β_s_i`` find the corresponding gain
    for u_k in u_all 
        # Radial magnification. (A14)
        G_k_r = abs(1.0 / (1.0 + (1.0 - 2.0 * u_k^2)*α*exp(-u_k^2)))
        # Tangential magnefications. (A15)
        G_k_t = abs(1.0 / (1.0 + α*exp(-u_k^2)))
        # Total magnification
        G_k_sum += G_k_r*G_k_t
    end
    # Brightness distribution
    B = (1.0/(2π*σ^2)) * exp(-β_s_i^2/(2.0*σ^2));
	I_m += B*G_k_sum;
    return I_m
end



"""
(19) from Clegg+98. Integrate over sub-angles of ``β_s``.
"""
function integrate_numerically_clegg(α::T, σ::T, β_s_min::T, β_s_max::T, u_::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    # Cubature
    (val, err) = hquadrature(β_s_i -> integrand(u_, β_s_i, α, σ, u_min, u_max),
                             β_s_min, β_s_max; reltol=1e-7, abstol=1e-03,
                             maxevals=0)
    return val, err
end


"""
(A13) from Dong+2018. Integrate over sub-angles of ``β_s`` and sub-anlges of ``ϕ``.
"""
function integrate_numerically_dong(α::T, σ::T, β_s_min::T, β_s_max::T, u_::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    # HCubature
    (val, err) = hcubature(x -> integrand(u_, x[1], x[2], α, σ, u_min, u_max),
                           [β_s_min, 0.0], [β_s_max, 2*π]; rtol=1e-4,
                           atol=1e-2, maxevals=0, initdiv=4)
    return val, err
end


"""
Scaterring of the point source from Clegg+1998.

# Arguments
- `u__array::Vector{AbstractFloat}` - sequence of the observer plane coordinates where we
the scaterred intensity is calculated.
- `α::AbstractFloat` - strength of the lens.
"""
function scattering_point_clegg1998(u__array::Vector{T}, α::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    I_result = Vector{T}(undef, size(u__array, 1))
    Threads.@threads for i in 1:size(u__array, 1)
    # for i in 1:size(u__array, 1)
        I_result[i] = sum_of_magnifications(u__array[i], α, u_min, u_max)
    end
    return I_result
end


"""
Scaterring of the extended source from Clegg+1998.

# Arguments
- `u__array::Vector{AbstractFloat}` - sequence of the observer plane coordinates where we
the scaterred intensity is calculated.
- `α::AbstractFloat` - strength of the lens.
- `β_s::AbstractFloat` - angular size of the source in units of angular size of the
lens.
"""
function scattering_clegg1998(u__array::Vector{T}, β_s::T, α::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    # Angular size (std) of the Gaussian source
    σ = β_s/2.355
    I_result = Vector{T}(undef, size(u__array, 1))
    I_err = Vector{T}(undef, size(u__array, 1))
    Threads.@threads for i in 1:size(u__array, 1)
    # for i in 1:size(u__array, 1)
        I_result[i], I_err[i] = integrate_numerically_clegg(α, σ, -4σ, 4σ, u__array[i], u_min, u_max) 
    end
    return I_result, I_err
end


"""
Scaterring of the extended source on the axisymmetric lens from Dong+2018.

# Arguments
- ``u__array::Vector{AbstractFloat}`` - sequence of the observer plane coordinates where we
the scaterred intensity is calculated.
- ``α::AbstractFloat`` - strength of the lens.
- ``β_s::AbstractFloat`` - angular size of the source in units of angular size of the
the lens.
- ``b::AbstractFloat`` - Impact parameter (perpendicular distance from the symmetry
axis, in units of lens size.
"""
function scattering_dong2018(u__array::Vector{T}, β_s::T, α::T, b::T, u_min::T, u_max::T) where {T<:AbstractFloat}
    # Angular size (std) of the Gaussian source
    σ = β_s/2.355
    I_result = Vector{T}(undef, size(u__array, 1))
    I_err = Vector{T}(undef, size(u__array, 1))
    # for i in 1:size(u__array, 1)
    Threads.@threads for i in 1:size(u__array, 1)
        # (A13) from Dong+2018. Note that integration borders for ``phi`` is ``[0.0, 2.0*π]``.
        # Note that for negative ``u_`` we pass non-negative to the function
        # (using hypot). It is OK, because the picture is symmetric.
        I_result[i], I_err[i] = integrate_numerically_dong(α, σ, 0.0, 4.0*σ, hypot(u__array[i], b), u_min, u_max) 
    end
    return I_result, I_err
end

end
