
using Agents

mutable struct Bead <: AbstractAgent
    id::Int
    pos::Tuple{Float64,Float64}
    vel::Tuple{Float64,Float64}
    mass::Float64
end

function initialize(;
    N_beads=100,
    extend=(10,10),
    dt=0.1, 
    v0 = 5.0,
    friction = 1.0)

    fluid = Fluid(0.2, 1e-7, 2e-5)

    space2d = ContinuousSpace(2; periodic=false, extend=extend)
    model = ABM(Bead, space2d, 
    properties = 
    Dict(
        :dt => dt, 
        :fluid_obj => fluid,
        :friction => friction)
    )

    for ind in 1:N_beads
        pos = Tuple(rand(2)) .* extend[1]
        vel = Tuple(randn(2)) .* v0
        add_agent!(pos, model, vel, 3.0)
    end

    index!(model)
    return model
end

# agent_step!(agent, model) = move_agent!(agent, model, model.dt)
function agent_step!(agent, model)
    Vx = Array(model.fluid_obj.Vx)
    Vy = Array(model.fluid_obj.Vy)
    fric = model.friction
    posx, posy = Int.(round.(agent.pos))
    if posx >= 1 && posx <= N && posy >= 1 && posy <= N
        vx = agent.vel[1] + Vx[posx, posy] * model.friction
        vy = agent.vel[2] + Vy[posx, posy] * model.friction
    else
        vx = agent.vel[1]
        vy = agent.vel[2]
    end
    agent.vel = (vx, vy)
    move_agent!(agent, model, model.dt)
end

function model_step!(model)
    for (a1, a2) in interacting_pairs(model, 20.0, :nearest)
        elastic_collision!(a1, a2, :mass)
    end
    simul_step2!(model.fluid_obj)
end

get_position_points(model) = 
    [Point2f0(ag.pos[1], ag.pos[2]) for (_, ag) in model.agents]

function bead_model_step!(pnode, dnode, model)

    Agents.step!(model, agent_step!, model_step!, 1)

    pnode[] = get_position_points(model)

    dnode[] = Array(model.fluid_obj.density)
end

function bead_model_live!(scene, pnode, dnode, model, nframe, 
    framerate)
    for _ in 1:nframe
        bead_model_step!(pnode, dnode, model)

        scene[2].colorrange = color_range

        sleep(1/framerate)
    end
end




