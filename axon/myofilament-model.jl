# using Myosin model from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768639/

using MEDYAN
using StaticArrays
using LinearAlgebra
using Setfield
using Random
using OrderedCollections: OrderedDict
using MEDYAN:
    add_link_type!,
    BondConfig,
    bond_energy_force,
    Link


# First create a agent names with actin and myosin filaments
# Also include a myosin end link2mon and a myosin end to actin link
const agent_names = MEDYAN.AgentNames(;
    diffusingspeciesnames=[
        :A, # diffusing actin
        :B, # diffusing active brancher
        :L, # diffusing linker
    ],
    filamentnames=[
        (:actin, [
            :middle,
            :weak_bound,
            :strong_bound,
            :plusend,
            :minusend,
            :boundARP23, #arp23 on minus end of attached daughter filament
            :freeARP23, #arp23 on minus end of attached daughter filament
            :bound, # general bound actin
        ]),
        (:myo, [:myo,]),
    ],
)


const monomerspacing = 2.7 # nm

const myo_spacing = 20.0

"""
Parameters for a model of non muscle myosin 2.
"""
@kwdef struct MyosinParameters
    number_of_monomers_per_side::Int
    load_force::Float64 # pN
    on_rate::Float64 # 1/s
    off_rate::Float64 # 1/s
    weak_off_rate::Float64 # 1/s
    weak_on_off_ratio::Float64 # unitless
    step_distance::Float64 # nm
    off_bell_distance::Float64 # nm
    step_bell_distance::Float64 # nm
    myo_motor_mech_params=(;
        k= 0.05,# pN/nm
        maxdist= 50.0, # nm
    )
end

function make_context(p::MyosinParameters; grid = CubicGrid((12,2,2),500.0), context_kwargs, treadmilling_rate_factor= 1.0)

    # Create a SysDef
    s = MEDYAN.SysDef(agent_names)

    "net number of ATP used by myosin, just for tracking purposes"
    used_ATP::Base.RefValue{Int} = Ref(0)

    add_diffusion_coeff!(s,
        :A,
        20E6, # nm^2/s
    )

    add_diffusion_coeff!(s,
        :B,
        0.2E6, # nm^2/s
    )

    add_diffusion_coeff!(s,
        :L,
        2.0E6, # nm^2/s
    )

    # #define restraints
    # add_link_type!(s;
    #     name=:restraint,
    #     description="harmonic position restraint",
    #     places=[FilaMonoIdx(),],
    #     bonds=[BondConfig(;
    #         bond=MEDYAN.PositionRestraint(),
    #         input=(1,),
    #         param=(;k=0.2),
    #         state=(;r0=SA[NaN,NaN,NaN]),
    #     )]
    # )
    # add_link_type!(s;
    #     name=:constforce,
    #     description="constant force",
    #     places=[FilaMonoIdx(),],
    #     bonds=[
    #         BondConfig(;
    #             bond=MEDYAN.ConstantForce(),
    #             input=(1,),
    #             param=(;),
    #             state=(;f=SA[NaN,NaN,NaN]),
    #         ),
    #     ]
    # )

    # Add Actin filament parameters
    add_filament_params!(s, :actin, MEDYAN.ACTIN_FIL_PARAMS)

    actin_middle_state::UInt8 = s.state.actin.middle
    actin_weak_bound_state::UInt8 = s.state.actin.weak_bound
    actin_strong_bound_state::UInt8 = s.state.actin.strong_bound

    # Add Myosin filament parameters
    add_filament_params!(s, :myo, MEDYAN.FilamentMechParams(
        radius= 15.0,
        spacing= myo_spacing,
        klength= 10*100.0,
        kangle= 1,
        numpercylinder= 11,
        max_num_unmin_end= 1,
    ))

    # Define a callback for the myo motor binding reaction
    function bind_motor(c::MEDYAN.Context; link, place, kwargs...)
        if get_chem_state(c, place).monomer_state != actin_middle_state
            return 0 # the monomer is already bound
        end
        local link_state = get_state(c, link)
        local numUnbound::Int = link_state.numUnbound
        local numBound::Int = link_state.numBound

        local new_translation = -p.step_distance/2
        local new_unbound_state = (numBound=numBound+1, numUnbound=numUnbound-1)

        local myo_end_tag = only(link2tags(c,link))
        local m_pos = get_position(c, myo_end_tag)
        local a_pos = get_position(c, place)
        local a_plusvec = get_directions(c, place)[1]
        # get the forces and energies if the motor were to bind.
        E, junk = bond_energy_force(
            MEDYAN.MaxDistanceRestraint(),
            ((m_pos,), (a_pos, a_plusvec),),
            merge(p.myo_motor_mech_params, (translated=new_translation,)),
        )
        # TODO add catch-slip dynamics here. Or not.
        boltzmann_factor = exp(-c.β*E)
        if rand() < boltzmann_factor
            # update the myo_fil_end link number of bound and unbound motors
            update_link!(c, link; state = new_unbound_state)
            # update the state of the bound monomer, setting it to bound to the myo_fil_end link
            update_fila_mono_state!(c, place, actin_weak_bound_state)
            # add a myo_motor link between the monomer and the myo_fil_end link myo filament monomer.
            # This will be used to apply a mechanical force to the myo filament monomer
            make_link!(c;
                type = :myo_motor,
                places = (myo_end_tag, place),
                state = (;parent_tag=link,),
                bond_states = ((;translated=new_translation,),),
            )
            return 1
        else
            return 0
        end
    end

    weak_on_rate::Float64 = p.weak_off_rate*p.weak_on_off_ratio
    add_link_type!(s;
        name=:myo_fil_end,
        description="This has no mechanical force field but keeps track of the number of bound and unbound motors at the end of a myosin filament",
        places=[FilaMonoIdx(),],
        state= (;
            numBound = 0,
            numUnbound = 0,
        ),
        reactions=[[
            (; # binding to an actin filament
                affect! = bind_motor,
                rate = ((c::MEDYAN.Context; link_state, kwargs...) -> link_state.numUnbound*weak_on_rate),
                fila_cutoff = (:actin, 150.0), # this should be significantly larger than the myo_fil_end maxdist
            )
        ]]
    )

    # Add myo_motor parameters
    # This has a mechanical force field that is used to apply a force to the myo filament monomer
    add_link_type!(s;
        name = :myo_motor,
        description = "myosin motor",
        places = [FilaMonoIdx(), FilaMonoIdx()],
        state = (;
            parent_tag = Link(),
        ),
        bonds = [
            BondConfig(;
                bond = MEDYAN.MaxDistanceRestraint(),
                input = (1,2,),
                param = p.myo_motor_mech_params,
                state = (translated=NaN,),
            ),
        ],
        reactions = [[
            (; # ADP unbinds, ATP binds, actin unbinds, and ATP hydrolyses
                # A⋅M⋅D
                # Assume Pi rebinding is impossible
                # Assume after ADP unbinds, ATP binding, actin unbinding, and ATP hydrolysis are instant and irreversible.
                # Force dependent unbinding rate.
                affect! = (c::MEDYAN.Context; link, kwargs...) -> let
                    used_ATP[] += 1
                    local link_state = get_state(c, link)
                    local parent_link = link_state.parent_tag
                    local parent_link_state = get_state(c, parent_link)
                    @assert parent_link_state.numBound > 0
                    update_link!(c, parent_link; state = (numBound=parent_link_state.numBound-1, numUnbound=parent_link_state.numUnbound+1,))
                    local myo_end_tag, actin_tag = link2tags(c, link)
                    update_fila_mono_state!(c, actin_tag, actin_middle_state)
                    remove_link!(c, link)
                    1
                end,
                rate = (c::MEDYAN.Context; link_state, link, link_data, kwargs...) -> let
                    # get force on the actin filament and the myosin projected on the positive actin filament direction.
                    local out = get_link_mechanics(c, link, link_data)
                    local F = out.inputs[2][2] ⋅ out.forces[2]
                    p.off_rate*exp(c.β*F*p.off_bell_distance)
                end,
                enabled = false,
            ),
            (; # weak unbinding
                # A-M⋅D⋅Pi
                affect! = (c::MEDYAN.Context; link, kwargs...) -> let
                    local link_state = get_state(c, link)
                    local parent_link = link_state.parent_tag
                    local parent_link_state = get_state(c, parent_link)
                    @assert parent_link_state.numBound > 0
                    update_link!(c, parent_link; state = (numBound=parent_link_state.numBound-1, numUnbound=parent_link_state.numUnbound+1,))
                    local myo_end_tag, actin_tag = link2tags(c, link)
                    update_fila_mono_state!(c, actin_tag, actin_middle_state)
                    remove_link!(c, link)
                    1
                end,
                rate = Returns(p.weak_off_rate),
                enabled = true,
            ),
            (; # power stroke
                # A-M⋅D⋅Pi
                affect! = (c::MEDYAN.Context; link, kwargs...) -> let
                    update_link!(c, link;
                        reaction_enabled = ((true, false, false,),),
                        bond_states = ((translated=p.step_distance/2,),),
                    )
                    update_fila_mono_state!(c, link2tags(c, link)[2], actin_strong_bound_state)
                    1
                end,
                rate = (c::MEDYAN.Context; link_state, link, link_data, kwargs...) -> let
                    local out = get_link_mechanics(c, link, link_data)
                    local F = out.inputs[2][2] ⋅ out.forces[2]
                    p.on_rate*exp(c.β*F*p.step_bell_distance)
                end,
                enabled = true,
            ),
        ]]
    )

    MEDYAN.add_link_type!(s;
        name=:crosslinker,
        description="actin crosslinker",
        places=[MEDYAN.FilaMonoIdx(), MEDYAN.FilaMonoIdx()],
        bonds=[
            (;
                bond=MEDYAN.DistanceRestraint(),
                input=(1,2),
                state=(;L0=NaN,),
                param=(;k=8.0,),
            ),
        ],
        reactions=[
            [
                (; # unbinding
                    affect! = (c::MEDYAN.Context; link, chem_voxel, kwargs...) -> let
                        local mt, pt = link2tags(c, link)
                        remove_link!(c, link)
                        #set the new monomer states
                        update_fila_mono_state!(c, mt, :middle)
                        update_fila_mono_state!(c, pt, :middle)
                        add_diffusing_count!(c; species=:L, chem_voxel, amount=+1)
                        1
                    end,
                    rate= MEDYAN.LinkRateSlipBond(f0 = inv(0.24*MEDYAN.default_β), k0 = 0.3*treadmilling_rate_factor, kmax=1E4),
                ),
            ],
        ],
    )

    MEDYAN.add_link_type!(s;
        name=:brancher,
        description="actin brancher",
        places=[MEDYAN.FilaMonoIdx(), MEDYAN.FilaMonoIdx()],
        bonds=[
            (;
                bond=MEDYAN.BranchBendingCosine(),
                input=(1,2),
                param=(;kr=100.0, kbend=10.0, cos0=cos(1.22), sin0=sin(1.22)),
                no_collide=true,
            ),
        ],
        reactions=[
            [
                (; # unbinding
                    affect! = (c::MEDYAN.Context; link, chem_voxel, kwargs...) -> let
                        local mt, pt = link2tags(c, link)
                        #set the new monomer states
                        update_fila_mono_state!(c, mt, :middle)
                        update_fila_mono_state!(c, pt, :freeARP23)
                        remove_link!(c, link)
                        1
                    end,
                    rate= MEDYAN.LinkRateSlipBond(f0 = 6.0, k0 = 0.02*treadmilling_rate_factor, kmax=1E4),
                ),
            ],
        ],
    )

    #plus end polymerization
    addfilamentend_reaction!(s, :actin, :pp, false,
        [:plusend]=>[:middle,:plusend], monomerspacing,
        "diffusing.A -->", 0.154*500^3*treadmilling_rate_factor, 1,
    )
    #plus end depolymerization
    addfilamentend_reaction!(s, :actin, :dpp, false,
        [:middle,:plusend]=>[:plusend], 0.0,
        "--> diffusing.A", 1.4*treadmilling_rate_factor, 0,
    )

    #minus end polymerization
    addfilamentend_reaction!(s, :actin, :mp, true,
        [:minusend]=>[:minusend,:middle], monomerspacing,
        "diffusing.A -->", 0.0173*500^3*treadmilling_rate_factor, 1,
    )
    #minus end depolymerization
    addfilamentend_reaction!(s, :actin, :dmp, true,
        [:minusend,:middle]=>[:minusend], 0.0,
        "--> diffusing.A", 0.8*treadmilling_rate_factor, 0,
    )

    #Destruction
    addfilamentendsite!(s,:actin,:destroy_actin_fil,
        MEDYAN.FilamentEndSiteGeneral(false,[s.state.actin.minusend,s.state.actin.plusend],0.0)
    )
    MEDYAN.addreactioncallback!(s, "filamentendsite.actin.destroy_actin_fil", 1.0*treadmilling_rate_factor, 0,
        MEDYAN.FilamentDestructionCallback(s.filament.actin, s.filamentendsite.actin.destroy_actin_fil.id, [s.diffusing.A=>2])
    )

    addfilamentendsite!(s,:actin,:destroy_arp23_fil,
        MEDYAN.FilamentEndSiteGeneral(false,[s.state.actin.freeARP23,s.state.actin.plusend],0.0)
    )
    #note ARP23 takes up two monomer spaces, so no actin is created.
    #Assume very fast.
    MEDYAN.addreactioncallback!(s, "filamentendsite.actin.destroy_arp23_fil", 100.0, 0,
        MEDYAN.FilamentDestructionCallback(s.filament.actin, s.filamentendsite.actin.destroy_arp23_fil.id, [s.diffusing.B=>1])
    )


    #branching site
    site = MEDYAN.FilamentSiteGeneral(2,fill(s.state.actin.middle,3))
    addfilamentsite!(s, :actin, :branch, site)
    branchcallback = MEDYAN.FilamentSiteBranchingCallback(
        MEDYAN.GeneralFilamentSiteCallback(
            s.filament.actin,
            s.filamentsite.actin.branch.id,
            1,
            [s.state.actin.bound],
            [],
        ),
        s.link.brancher.id,
        s.filament.actin,
        true,
        true, 
        [s.state.actin.boundARP23,s.state.actin.plusend], 
        # note, this plusend isn't actually an actin, its just part of the ARP23 that acts like an actin.
        # That is why no diffusing actin is used up.
        [s.diffusing.B=>-1],
    )
    addreactioncallback!(s,
        "filamentsite.actin.branch + diffusing.B",
        0.154*500^3,
        1,
        branchcallback,
    )

    #minus end free ARP23 depolymerization, note ARP23 takes up two monomer spaces.
    addfilamentend_reaction!(s, :actin, :dpolymbr, true,
        [:freeARP23,:middle,:middle]=>[:minusend], 0.0,
        "--> diffusing.B", 0.8*treadmilling_rate_factor, 0,
    )

    #crosslinker binding
    site = MEDYAN.Decimated2MonSiteMinAngleRange(
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
    add_decimated_2mon_site!(s,:crosslinkbinding,site)
    sitecallback = MEDYAN.SimpleCrosslinkBindCallback(
        s.decimated_2mon_site.crosslinkbinding.id,
        s.link.crosslinker.id,
        s.state.actin.bound,
        [s.diffusing.L=>-1],
    )
    addreactioncallback!(s,
        "decimated_2mon_site.crosslinkbinding + diffusing.L",
        0.01*500^3/2*treadmilling_rate_factor,
        1,
        sitecallback,
    )

    used_ATP[] = 0
    c = MEDYAN.Context(s, grid;
        cylinder_skin_radius=7.0,
        context_kwargs...,
    )
    (;c, used_ATP)
end
