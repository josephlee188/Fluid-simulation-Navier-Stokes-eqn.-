mutable struct Fluid
    size::Int
    dt::Float64
    diff::Float64
    visc::Float64

    src_v::Tuple{Float64,Float64}

    # s::Array
    # density::Array
    # Vx::Array
    # Vy::Array
    # Vx0::Array
    # Vy0::Array
    # aux::Array

    s::CuArray
    density::CuArray
    Vx::CuArray
    Vy::CuArray
    Vx0::CuArray
    Vy0::CuArray
    aux::CuArray

  
end

function Fluid(dt, diffusion, viscosity)
    # GPU arrays 
    # s =         CUDA.fill(0.0f0, (N,N))
    # density =   CUDA.fill(0.0f0, (N,N))
    # Vx =        CUDA.fill(0.0f0, (N,N))
    # Vy =        CUDA.fill(0.0f0, (N,N))
    # Vx0 =       CUDA.fill(0.0f0, (N,N))
    # Vy0 =       CUDA.fill(0.0f0, (N,N))
    # aux =       CUDA.fill(0.0f0, (N-2,N-2,12))
    # aux[:,:,11] .= CuArray([i for i in 2:N-1, j in 2:N-1])
    # aux[:,:,12] .= CuArray([j for i in 2:N-1, j in 2:N-1])

    src_v = (0.0, 80.0)

    s =         zeros(N,N)
    density =   zeros(N,N)
    Vx =        zeros(N,N)
    Vy =        zeros(N,N)
    Vx0 =       zeros(N,N)
    Vy0 =       zeros(N,N)
    aux =       zeros(N-2,N-2,12)   # Container for preallocation
    aux[:,:,11] .= [i for i in 2:N-1, j in 2:N-1]
    aux[:,:,12] .= [j for i in 2:N-1, j in 2:N-1]

    fluid = Fluid(N, dt, diffusion, viscosity,
    src_v, s, density, Vx, Vy, Vx0, Vy0, aux)
    return fluid
end

function simul_step!(fl::Fluid)
    visc =  fl.visc
    dt =    fl.dt
    diff =  fl.diff

    diffuse!(1, fl.Vx0, fl.Vx, visc, dt)
    diffuse!(2, fl.Vy0, fl.Vy, visc, dt)
    
    project!(fl.Vx0, fl.Vy0, fl.Vx, fl.Vy)
    
    advect!(1, fl.Vx, fl.Vx0, fl.Vx0, fl.Vy0, dt, fl.aux)
    advect!(2, fl.Vy, fl.Vy0, fl.Vx0, fl.Vy0, dt, fl.aux)

    project!(fl.Vx, fl.Vy, fl.Vx0, fl.Vy0)

    diffuse!(0, fl.s, fl.density, diff, dt)
    advect!(0, fl.density, fl.s, fl.Vx, fl.Vy, dt, fl.aux)

end

function add_density!(fluid::Fluid, x, y, amount)
    fluid.density[[CartesianIndex(x,y),]] .+= cu([amount,])
end

# External force
function add_velocity!(fluid::Fluid, x, y, amountx, amounty)
    fluid.Vx[[CartesianIndex(x,y),]] .+= cu([amountx,])
    fluid.Vy[[CartesianIndex(x,y),]] .+= cu([amounty,])
end

function add_density_circle!(fluid::Fluid, r, amount)
    for θ in LinRange(0, 2π*(1-1/n_θ), n_θ)
        x = Int(round(r * cos(θ) + N/2))
        y = Int(round(r * sin(θ) + N/2))
        add_density!(fluid, x, y, amount)
    end
end

function add_velocity_circle!(fluid::Fluid, r, amount)
    for θ in LinRange(0, 2π*(1-1/n_θ), n_θ)
        x = Int(round(r * cos(θ) + N/2))
        y = Int(round(r * sin(θ) + N/2))
        vx = -amount * cos(θ) + η * randn()
        vy = -amount * sin(θ) + η * randn()
        add_velocity!(fluid, x, y, vx, vy)
    end
end

get_rand_range(a,b) = (2*rand() - 1)*(a-b) + b  

normalize(a,b) = (a, b)./√(a^2+b^2)

get_angle(a,b) = begin
    ua, ub = normalize(a,b)
    ub >= 0 ? acos(ua) : 2π - acos(ua)
end

function simul_step2!(fl)
    # 'simul_step!' with sources added, and update the Node for 
    # visualization.

    #####################  Sources  #######################
    # for i = -2:2
    #     for j = -2:2
    #         add_density!(fl, cx+i,cy+j, get_rand_range(200,300))
    #     end
    # end

    # vx, vy = fl.src_v

    # add_velocity!(fl, cx, cy, vx, vy)

    # # Update velocity source direction
    # θ = get_angle(vx, vy)

    # θ += randn() * η

    # fl.src_v = V.*(cos(θ), sin(θ))
    
    add_density_circle!(fl, r, amount_d)
    add_velocity_circle!(fl, r, amount_v)

    ########################################################

    simul_step!(fl)

    # anode[] = Array(fluid.density)
end

function simul_live!(scene, dnode, pnode, fl, model, 
    nframe, framerate)
    for t = 1:nframe
        simul_step2!(fl)

        scene[2].colorrange = color_range

        sleep(1/framerate)
    end
end


