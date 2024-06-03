module Root_finding_utils


# Home made function for searching the multiple roots
function rootsearch(f::Function, a::T, b::T, dx::T) where(T<:AbstractFloat)
	x1 = a
	f1 = f(a)
	x2 = a + dx
	f2 = f(x2)
	
	while f1*f2 > 0
		# No such interval!
		if x1 >= b
			return (a, b)
        end
		x1 = x2
		f1 = f2
		x2 = x1 + dx
		f2 = f(x2)
    end
	return (x1, x2)
end

# https://mmas.github.io/bisection-method-julia
function bisection(f::Function, a::AbstractFloat, b::AbstractFloat;
                   tol::AbstractFloat=1e-5, maxiter::Integer=100)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
        i += 1
        i != maxiter || error("Max iteration exceeded")
        c = (a+b)/2
        fc = f(c)
        if fc == 0
            break
        elseif fa*fc > 0
            a = c  # Root is in the right half of [a,b].
            fa = fc
        else
            b = c  # Root is in the left half of [a,b].
        end
    end
    return c
end


function multiple_root_bisection(f::Function, a::T, b::T, eps::T) where(T<:AbstractFloat)
    roots = Vector{T}()
	
	while true
		# Find interval ``[x1, x2]`` of width ``eps``, where ``f`` changes sign
		x1x2 = rootsearch(f, a, b, eps)
        if x1x2[1] !== a || x1x2[2] !== b
			# Increase the low interval end
            a = x1x2[2]
			# Find root of ``f`` inside the found interval
            root = try bisection(f, x1x2[1], x1x2[2])
            catch y
                # If failed to find - just pass
                if isa(y, ErrorException)
                    println("Failed to find root using bisect!")
                    continue
                end
            end
			append!(roots, root)
		else
			break
        end
    end
	
	return roots
end

end
