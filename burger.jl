function main()
    L::Float64 = 2π
    T::Float64 = 5

    nx::Int16 = 100
    nsteps::Int16 = 5000

    dx::Float64 = L/nx
    dt::Float64 = T/nsteps

    u = Array{Float64}(undef, nx)
    x = [dx*i for i=1:nx]
    init(u,nx,dx)
    😿(u,nx,nsteps,dx,dt)

    return x,u

end

function init(u::Array{Float64},nx::Int16,dx::Float64)
    for i = 1:nx
        x = dx*i
        u[i] = sin(x)*exp(-(x-π)^2)
    end

    return 0
end

function 😿(u::Array{Float64},nx::Int16,nt::Int16,dx::Float64,dt::Float64)
    v::Float64 = 0.01
    f::Float64 = 0
    coeff::Int8 = 0

    for i = 1:nt
        u[1] = u[nx-1]
        u[2] = u[nx]
        for j = 2:nx-1
            coeff = Int(u[j]>0)
            f = v*(dx^2)*(u[j+1]-2*u[j]+u[j-1])-u[j]/dx*(u[j-coeff+1]-u[j-coeff])
            u[j] += dt*f
        end
    end
    return 0
end

main()

