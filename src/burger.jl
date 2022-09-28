using LinearAlgebra

function main()
    L::Float64 = 2π
    ⏳::Float64 = 5

    nx::Int16 = 100
    nwg::Int16 = nx + 2 # 2 ghost points
    nsteps::Int16 = 5000

    Δx::Float64 = L/nx
    Δ⏳::Float64 = ⏳/nsteps

    u = Array{Float64}(undef, nwg)
    x = [Δx*i for i=1:nwg]
    init(u,nwg,Δx)
    expl(u,nx,nsteps,Δx,Δ⏳)
    #implicit(u,nx,nsteps,Δx,Δ⏳)

    return x,u

end

function init(u::Array{Float64},nx::Int16,Δx::Float64)
    for i = 1:nx
        x = Δx*i
        u[i] = sin(x)*exp(-(x-π)^2)
    end

    return 0
end

function expl(u::Array{Float64},nx::Int16,nt::Int16,Δx::Float64,Δ⏳::Float64)
    v::Float64 = 0.01
    f::Float64 = 0
    coeff::Int8 = 0

    for i = 1:nt
        u[1] = u[nx-1]
        u[2] = u[nx]
        u[nx+1] = u[3]
        u[nx+2] = u[4]
        for j = 2:nx+1
            coeff = Int(u[j]>0)
            f = v*(Δx^2)*(u[j+1]-2*u[j]+u[j-1])-u[j]/Δx*(u[j-coeff+1]-u[j-coeff])
            u[j] += Δ⏳*f
        end
    end
    return 0
end

function implicit(u::Array{Float64},nx::Int16,nt::Int16,Δx::Float64,Δ⏳::Float64)
    # general flow: make first guess using forward euler, then do Newton iterations

    uk = Array{Float64}(undef, nx+2)
    a = Array{Float64}(undef, nx+2)


    for t = 1:nt
        uk = copy(u)
        uk[1] = u[nx-1]
    end


    return 0
end

main()

