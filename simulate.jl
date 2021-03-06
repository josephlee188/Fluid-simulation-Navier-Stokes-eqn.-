"""
Simulation of Navier Stokes equation (2D) with CUDA

This is written solely as practice for Julia language.

The algorithms presented here are almost direct copy of the following sources;

https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
https://www.youtube.com/watch?v=alhpH6ECFvQ&t (YouTube Channel The Coding Train)
"Real-Time Fluid Dynamics for Games" by Joe Stam

Advices and comments are always welcome :)
"""

using CUDA
using GLMakie, AbstractPlotting
include("Fluid.jl")
include("Fluid_Motion.jl")
include("bead_agents.jl")

CUDA.allowscalar(false)
N = 500
iter = 16
SCALE = 4
t = 0

bnd_ind1 = [CartesianIndex(1,1),
            CartesianIndex(1,N),
            CartesianIndex(N,1),
            CartesianIndex(N,N)]
bnd_ind2 = [CartesianIndex(2,1),
            CartesianIndex(2,N),
            CartesianIndex(N-1,1),
            CartesianIndex(N-1,N)]
bnd_ind3 = [CartesianIndex(1,2),
            CartesianIndex(1,N-1),
            CartesianIndex(N,2),
            CartesianIndex(N,N-1)]

fluid = 0.0; # Free the previous gpu array
fluid = Fluid(0.2, 1e-7, 2e-5)  # dt, 

model = initialize(;extend=(N, N), fluid_obj=fluid)

points = get_position_points(model)

scene = Scene(resolution=(1000,1000))

fluid_d = Array(fluid.density)

d_Node = Node(reshape(fluid_d, (N,N)))

p_Node = Node(points)

heatmap!(scene, d_Node)

scatter!(scene, p_Node)

scatter_obj = scene[end]
scatter_obj.attributes.markersize[] = 5
scatter_obj.attributes.strokewidth[] = 3

scene[Axis].attributes.names.textsize[] = (2,2)
scene[Axis].attributes.ticks.textsize[] = (2,2)

scene

## Live
# cx = Int(N/2)   # Source position x
# cy = Int(N/2)   # Source position y
# η = 0.1        # Rate of change of source angle
# V = 20.0

amount_d = 30.0
amount_v = 1.0
r = 100.0
n_θ = 400
η = 1.0


color_range = (0.0, 30.0)
nframe = 100

simul_live!(scene, d_Node, p_Node, fluid, model,
            nframe, 60)

## Save video
n_frames = 500
framerate = 60
t_iterator = 1:n_frames
cx = Int(N/2)
cy = Int(N/2)
ax = 0.0
ay = 100.0
vrate = 100.0

GLMakie.record(scene, "fluid_simulation.gif", t_iterator; framerate=framerate) do t 

    simul_step2!(d_Node, fluid)

    scene[end].colorrange = (0.0, 15.0)
    
end

