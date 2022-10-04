using LinearAlgebra

function main()

    # Set PDE boundaries
    L::Float64 = 2π
    Tf::Float64 = 5

    # PDE parameters
    nx::Int64 = 5
    nt::Int64 = 5000
    Δx::Float64 = L/nx
    Δt::Float64 = Tf/nt
    v::Float64 = 0.01

    x = [Δx*i for i = 1:nx]

    u = [sin(x[i])*exp(-(x[i]-π)^2) for i = 1:nx]

    T = upwind(u,nx,Δx,Δt,v)

    # explicit(u,nx,nt,Δx,Δt,v)

    return x,u,T

end

function explicit(u::Array{Float64}, nx::Int64, nt::Int64, Δx::Float64,Δt::Float64, v::Float64)

    # Main loop
    for i=1:nt
        T = upwind(u,nx,Δx,Δt,v)
        u += Δt.*(T*u)
    end

    return 0

end

function upwind(u::Array{Float64}, nx::Int64, Δx::Float64, Δt::Float64, v::Float64)

    # 1st Order Finite Difference
    # Main diagonal of the operator matrix
    # d1 = [if u[i] > 0 1 else (-1) end for i = 1:nx]
    # l1 = [if u[i] > 0 (-1) else 0 end for i = 2:nx]
    # u1 = [if u[i] > 0 0 else 1 end for i = 1:nx-1]


    # Pure Advection Diffusion
    d1 = [1. for i = 1:nx]
    l1 = [-1. for i = 1:nx-1]
    u1 = [0.0 for i = 1:nx-1]

    # Periodic Boundary Conditions
    T1 = [0 for i = 1:nx, j = 1:nx] + Tridiagonal(l1,d1,u1)
    T1[1,nx] = u[1] > 0 ? (-1) : 0
    T1[nx,1] = u[nx] > 0 ? 0 : 1

    # 2nd Order Finite Difference
    d2 = [-2 for i = 1:nx]
    ul2 = [1 for i = 1:nx-1]
    T2 = [0 for i = 1:nx, j = 1:nx] + Tridiagonal(ul2,d2,ul2)
    T2[1,nx] = 1
    T2[nx,1] = 1

    T = (T2.*(v/Δx^2)-T1./Δx)

    return T1

end

main()

