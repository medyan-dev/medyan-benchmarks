import MEDYANSimRunner
using StaticArrays
using LinearAlgebra
using Random
using Setfield
using OrderedCollections: OrderedDict
using SmallZarrGroups
include("myofilament-model.jl")

const Δt = 0.001 #s

const height = 7500.0 # nm
const radius = 1000.0 # nm
const chem_voxel_size = 500.0 # nm
const grid = CubicGrid((cld(2*radius,chem_voxel_size),cld(2*radius,chem_voxel_size),cld(height,chem_voxel_size)),chem_voxel_size)

const chemboundingcapsule = MEDYAN.boundary_capsule(;
    axis=SA[0.0,0.0,height-2radius],
    radius,
)
# offset chemboundary by 10.0 nm
const mechboundoffset = 10.0
const mechboundingcapsule = MEDYAN.boundary_capsule(;
    axis=SA[0.0,0.0,height-2radius],
    radius=radius - mechboundoffset,
    stiffness=100.0,
)

const myo_params = MyosinParameters(;
    number_of_monomers_per_side=8,
    load_force=0.0,
    on_rate=0.5,
    off_rate=0.35,
    step_distance=6.0,
    off_bell_distance=20.0,
    step_bell_distance=-10.0,
    weak_off_rate=100.0,
    weak_on_off_ratio=1/10,
    myo_motor_mech_params=(;
        k= 0.04,# pN/nm
        maxdist= 30.0, # nm
    ),
)
(;c, used_ATP) = make_context(myo_params;
    grid,
    context_kwargs = (;
        g_tol=5.0,
        max_cylinder_force=6000.0,
        nthreads=4,
        cylinder_skin_radius=7.0
    ),
    treadmilling_rate_factor=1.0,
)

group = SmallZarrGroups.load_zip(joinpath(@__DIR__,"traj", MEDYANSimRunner.step_path(10000)))["snap"]

function bench(c; seed=1234)
    Random.seed!(seed)
    MEDYAN.load_snapshot!(c, group["medyan"]::ZGroup)
    for step in 1:5
        @info "$step mechanics"
        MEDYAN.minimize_energy!(c; brownian_motion_time=Δt)
        @info "$step chemistry"
        MEDYAN.run_chemistry!(c, Δt)
    end
    c
end

# @time bench(c)
