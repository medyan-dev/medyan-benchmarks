using ColorSchemes
using JLD2
using LinearAlgebra
using MEDYAN
using Plots
using StaticArrays
using Statistics
using StatsBase

using Random
Random.seed!(2)


function is_new_medyan()
    hasfield(MEDYAN.MembraneMechParams, :edgelength_mem3dg_min)
end


function runwith(;
        speciespotentialenergy=MEDYAN.default_membranespeciespotentialenergy,
        init_diffusion=0,
        init_diffusion_actin = 30000,
        init_diffusion_linker = 1000,
        time::Real = 40,
        report_frame_interval::Integer = 100,
        minimize_time_interval::Real = 0.01,
        compartment_size=500,
        grid_size=(3,3,3),
        membranemechparams = [MEDYAN.MembraneMechParams()],
        membrane_species_params = SA[MEDYAN.MembraneSpeciesParams()],
        meshsize=60,
        vis=nothing,
    )
    @info "Start running"

    statedef = MEDYAN.StateDef(
        diffusingspeciesnames = [:a, :diffusing_actin, :ld],
        membranediffusingspeciesnames = [:ma],
        filamentnames = [
            (:actin, [:plusend, :minusend, :middle, :bound, :restrained]),
        ],
        link_2mon_names = [:crosslinker, :bead_fixer],
    )
    s = MEDYAN.SysDef(statedef)

    #----------------------------------
    # Some constants.
    #----------------------------------
    r_on = 120.0
    r_off = 0.1
    monomerspacing = 2.7
    numpercylinder = 40
    diffusion_coeffs = [2e6, 20e6, 2e6]
    num_fil = 10
    fil_bot_z = 50
    fil_init_center_xy = SA[grid_size[1] * compartment_size / 2, grid_size[2] * compartment_size / 2]
    fil_init_spacing = 30

    filamentmechparams = [MEDYAN.FilamentMechParams(;
        radius= 3.0,
        spacing= monomerspacing,
        klength= 40*100.0,
        kangle= 40*672.0,
        numpercylinder,
    )]

    grid = MEDYAN.CubicGrid(grid_size, compartment_size)
    mesh_adapt_params = MEDYAN.MeshAdaptParams(
        max_size = meshsize,
        min_size = meshsize / 9,
    )


    # All reactions.
    MEDYAN.add_membranediffusion_asbulkreaction!(s)
    MEDYAN.add_membranesitereaction!(;
        s,
        name_newmembranesite = :ms_desorption,
        membranediffusingreactants = [:ma],
        membranediffusingproducts = Symbol[],
        reactionexpr_extra = "--> diffusing.a",
        rate = r_off,
        canchangerate_bypotentialenergy = true,
        invvolumepower = 0
    )
    MEDYAN.add_membranesitereaction!(;
        s,
        name_newmembranesite = :ms_adsorption,
        membranediffusingreactants = Symbol[],
        membranediffusingproducts = [:ma],
        reactionexpr_extra = "diffusing.a -->",
        rate = r_on,
        canchangerate_bypotentialenergy = false,
        invvolumepower = 1
    )

    add_link_2mon!(s,
        :crosslinker,
        Link2MonState((;), (L0=NaN,)),
        MEDYAN.DistanceRestraintMechParams(k=8.0),
    )
    add_link_2mon!(s,
        :bead_fixer,
        Link2MonState((;), (mr0=SA[NaN,NaN,NaN], mv̂0=SA[NaN,NaN,NaN],)),
        MEDYAN.RestraintMechParams(kr=200.0, kv̂=5000.0),
    )
    #minus end polymerization
    addfilamentend_reaction!(s,
        :actin, #filament type name
        :pm, #new site name
        true, #is minus end
        [:minusend]=>[:minusend,:middle], #change in monomer states
        monomerspacing, #free space needed for reaction (nm)
        "diffusing.diffusing_actin -->", #reaction expression
        0.0173 * 500^3, #rate of reaction ((nm^3)^(invvolumepower)/s)
        1, #inverse volume power
    )
    #plus end polymerization
    addfilamentend_reaction!(s, :actin, :pp, false,
        [:plusend]=>[:middle,:plusend], monomerspacing,
        "diffusing.diffusing_actin -->", 0.154 * 500^3, 1,
    )
    #minus end depolymerization
    addfilamentend_reaction!(s, :actin, :dpm, true,
        [:minusend,:middle]=>[:minusend], 0.0,
        "--> diffusing.diffusing_actin", 0.8, 0,
    )
    #plus end depolymerization
    addfilamentend_reaction!(s, :actin, :dpp, false,
        [:middle,:plusend]=>[:plusend], 0.0,
        "--> diffusing.diffusing_actin", 1.4, 0,
    )
    #crosslinker binding site
    let
        site = MEDYAN.LinkableSiteMinAngleRange(
            s.filament.actin,
            s.filament.actin,
            10,
            10,
            s.state.actin.middle,
            s.state.actin.middle,
            30.0,
            40.0,
            cos(5*π/180)
        )
        addlinkablesite!(s,:crosslinkbinding,site)
        sitecallback = MEDYAN.SimpleCrosslinkBindCallback(
            s.linkablesite.crosslinkbinding.id,
            s.link_2mon.crosslinker,
            s.state.actin.bound,
            [s.diffusing.ld=>-1],
        )
        addreactioncallback!(s,
            "linkablesite.crosslinkbinding + diffusing.ld",
            0.01*500^3/2,
            1,
            sitecallback,
        )
    end
	#crosslinker unbinding
    let
    	site = MEDYAN.Link2MonSiteSlipBond(f0 = inv(0.24*MEDYAN.default_β) , k0 = 0.3)
	    addunbindinglink_2mon_site!(s, 
            :crosslinker, :unbinding, site,
            :actin, :middle, :actin, :middle,
            "--> diffusing.ld", 1.0, 0, 
        )
    end


    #----------------------------------
    # Create context.
    #----------------------------------
    c = MEDYAN.Context(s,grid;
        diffusion_coeffs,
        g_tol = 0.1,
        shake_before_minimization = false,
        membrane_species_params,
        membranemechparams,
        filamentmechparams,
        func_membranespeciespotentialenergy = speciespotentialenergy,
        sharedtypedconfigs = MEDYAN.SharedTypedConfigs(
            mesh_boundary_pinning_mode = Val(:border2),
            bending_mode = Val(:bashkirov),
        ),
    )

    #----------------------------------
    # Initializing.
    #----------------------------------
    # Initialize membrane.
    m1 = let
        newmembrane!(c; type=1, meshinit=MEDYAN.MeshInitPlane(
            boxorigin = SA[0,0,0],
            boxwidths = SVector(grid_size) * compartment_size,
            normal = SA[0,0,1],
            normal_multiplier_from_origin = 300,
            samples = (grid_size[1] * compartment_size ÷ meshsize, grid_size[2] * compartment_size ÷ meshsize, 2),
        ))
    end

    # Adapt membranes.
    MEDYAN.adapt_membranes!(c; params = mesh_adapt_params)
    @info "Initial mech..."
    MEDYAN.minimize_energy!(c)
    MEDYAN.adapt_membranes!(c; params = mesh_adapt_params)
    compute_all_membrane_geometry!_system(c; include_ff=true)
    MEDYAN.compartmentalize!(c)
    if !isnothing(vis)
        MEDYAN.drawcontext!(vis, c, s)
    end


    # Initialize filaments.
    function add_one_fil(xy)
        local eqlen = numpercylinder * monomerspacing
        local numcyl = 2
        local nummon = numpercylinder * numcyl
        local monomerstates = fill(s.state.actin.middle, nummon)
        monomerstates[begin] = s.state.actin.restrained
        monomerstates[end] = s.state.actin.plusend

        local nodepositions = [SA[xy[1], xy[2], fil_bot_z + k * eqlen] for k ∈ 0:numcyl]
        local node_mids = [numpercylinder * k for k ∈ 0:numcyl]
        local fid = MEDYAN.chem_newfilament!(c;
            ftid = s.filament.actin,
            monomerstates,
            nodepositions,
            node_mids,
        )
        local mon1 = MonomerName(s.filament.actin, fid, 0)
        MEDYAN.chem_newlink_2mon!(
            c,
            s.link_2mon.bead_fixer,
            mon1 => mon1;
            changedmechstate = (; mr0 = MEDYAN.mon_position(c, mon1), mv̂0 = SA[0,0,1], ),
        )
    end
    let
        num_fil_rem = num_fil
        if num_fil_rem > 0
            add_one_fil(fil_init_center_xy)
            num_fil_rem -= 1
        end
        dir_in_fan(fan) = ((1 + fan) % 6) * (π / 3)
        starting_ang_in_fan(fan) = ((fan + 6 - 1) % 6) * (π / 3)

        ring = 1
        fan = 1
        fan_serial = 1
        while num_fil_rem > 0
            st = starting_ang_in_fan(fan)
            dir = dir_in_fan(fan)
            new_xy = fil_init_center_xy + fil_init_spacing * (ring * SA[cos(st), sin(st)] + fan_serial * SA[cos(dir), sin(dir)])
            add_one_fil(new_xy)
            num_fil_rem -= 1
            # Next
            if fan_serial == ring
                if fan == 6
                    ring += 1
                    fan = 1
                else
                    fan += 1
                end
                fan_serial = 1
            else
                fan_serial += 1
            end
        end
    end


    # Find left part of the membrane and send some copies to it.
    @info "Adding for $init_diffusion diffusing molecules."
    adddiffusingcount_rand!(c, s.diffusing.a, init_diffusion)
    adddiffusingcount_rand!(c, s.diffusing.diffusing_actin, init_diffusion_actin)
    adddiffusingcount_rand!(c, s.diffusing.ld, init_diffusion_linker)

    function getconc(m)
        3 .* m.vertices.attr.copynumbers.ma ./ m.vertices.attr.astar
    end

    coords = Vector{Vector{SVector{3,Float64}}}() # Coordinate history, indexed by [frame index, vertex index, dimension]
    curvs = Vector{Vector{Float64}}() # Curvatures, indexed by [frame index, vertex index]
    concs = Vector{Vector{Float64}}() # Indexed by [frame index, vertex index]
    counts = Vector{Vector{Int}}() # Indexed by [frame index, vertex index]

    #----------------------------------
    # Main loop.
    #----------------------------------
    for f ∈ 1:round(Int, time / minimize_time_interval)
        if f % 1 == 0
            @info "f=$f t=$(f * minimize_time_interval)"
            flush(stderr)
        end
        begin
            # Chemistry.
            compute_all_membrane_geometry!_system(c; include_ff=true)
            MEDYAN.run_chemistry!(c, minimize_time_interval)

            # Resolve membrane.
            compute_all_membrane_geometry!_system(c)
            let tree = MEDYAN.build_aabbtree_frommeshes(c.membranes)
                MEDYAN.resolve_all_filament_mesh_crossing!(c, tree)
            end
    
            # Mechanics.
            @info "Doing mech"
            @time MEDYAN.minimize_energy!(c)
            # Check for edge lengths.
            mindist2::Float64 = Inf
            for eindex in eachindex(m1.edges)
                e = MEDYAN.IE(eindex)
                if !MEDYAN.onborder(m1, e)
                    h = MEDYAN.halfedge(m1, e)
                    ho = MEDYAN.oppo(m1, h)
                    v1 = MEDYAN.target(m1, h)
                    v2 = MEDYAN.target(m1, ho)
                    c1 = m1.vertices.attr.coord[v1.value]
                    c2 = m1.vertices.attr.coord[v2.value]

                    u = c2 - c1
                    dist2 = dot(u, u)
                    mindist2 = min(mindist2, dist2)
                end
            end
            @info "Min edge length is $(sqrt(mindist2))"
            compute_all_membrane_geometry!_system(c; include_ff=true)
            MEDYAN.adapt_membranes!(c; params = mesh_adapt_params)

            # Misc.
            yield()
        end

        if f % report_frame_interval == 0
            push!(counts, copy(m1.vertices.attr.copynumbers.ma))
            push!(concs, m1.vertices.attr.copynumbers.ma ./ m1.vertices.attr.astar .* 3)
            push!(curvs, copy(m1.vertices.attr.curv))
            push!(coords, copy(m1.vertices.attr.coord))
        end
        if !isnothing(vis)
            MEDYAN.drawcontext!(vis, c, s)
        end
    end

    # Return values.
    (; coords, curvs, concs, counts, )
end



function filopodia_config1()
    kbend_lipid = 100
    if is_new_medyan()
        @info "MEDYAN: new"
        membranemechparams = [
            MEDYAN.MembraneMechParams(;
                kbend=kbend_lipid, tension=0.02, kpinning=500,
                edgelength_mem3dg_min = 6.0,
                edgelength_mem3dg_k = 80.0,
            ),
        ]
    else
        @info "MEDYAN: old"
        membranemechparams = [
            MEDYAN.MembraneMechParams(; kbend=kbend_lipid, tension=0.02, kpinning=500, ),
        ]
    end
    kbends_prot = SA[100]
    eqcurvs = SA[1/1]
    diffusion_init = SA[0]

    (; kbend_lipid, membranemechparams, kbends_prot, eqcurvs, diffusion_init)
end
function filopodia_eachrun1()
    (; membranemechparams, kbends_prot, eqcurvs, diffusion_init,) = filopodia_config1()
    all_indices = collect(Iterators.product(eachindex(kbends_prot), eachindex(eqcurvs), eachindex(diffusion_init)))
    index = all_indices[1]
    (ikbend_prot, ieqcurv, idiffusion_init) = index

    membrane_species_params = SA[
        MEDYAN.MembraneSpeciesParams(;
            diffusion_coeff = 2e4,
            area = 40,
            kbend = kbends_prot[ikbend_prot],
            eqcurv = eqcurvs[ieqcurv],
        ),
    ]

    runwith(;
        time = 0.4,
        init_diffusion = diffusion_init[idiffusion_init],
        membranemechparams, membrane_species_params,
        speciespotentialenergy = MEDYAN.MembraneSpeciesPotentialEnergy_BendingBashkirov(;
            params_mech = membranemechparams,
            params_species = membrane_species_params,
        ),
    )
end


res = filopodia_eachrun1();
