module CanopyFluxesMultilayerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates canopy fluxes
  !
  ! !USES:
  use shr_kind_mod               , only : r8 => shr_kind_r8
  use abortutils                 , only : endrun
  use decompMod                  , only : bounds_type
  use pftconMod                  , only : pftcon
  use PatchType                  , only : patch
  use atm2lndType                , only : atm2lnd_type
  use TemperatureType            , only : temperature_type
  use WaterStateType             , only : waterstate_type
  use WaterFluxType              , only : waterflux_type
  use EnergyFluxType             , only : energyflux_type
  use FrictionVelocityMod        , only : frictionvel_type
  use SoilStateType              , only : soilstate_type
  use SurfaceAlbedoType          , only : surfalb_type
  use CanopyFluxesMultilayerType , only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxesMultilayer  ! Compute canopy fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CanopyFluxesSum        ! Sum leaf and soil fluxes
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine CanopyFluxesMultilayer (bounds, num_exposedvegp, filter_exposedvegp, &
  num_nolakec, filter_nolakec,declin, atm2lnd_inst, temperature_inst, waterstate_inst, &
  waterflux_inst, energyflux_inst, frictionvel_inst, soilstate_inst, surfalb_inst, &
  mlcanopy_inst,urbanparams_inst,soil_water_retention_curve,canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Compute fluxes for sunlit and shaded leaves at each level
    ! and for soil surface.
    !
    ! !USES:
    use clm_varpar, only : ivis, inir, isun, isha, nlevsno, nlevgrnd, nlevcanml, nleaf
    use clm_varcon, only : mmh2o, mmdry, cpd, cpw, rgasc, denice, denh2o
    use clm_varctl, only : gstyp, turb_type, dtime_sub,iulog
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use SolarRadiationMod, only: SolarRadiation
    use SoilTemperatureMod, only : SoilThermProp
    use CanopyWaterMod, only : CanopyInterception, CanopyEvaporation
    use PlantHydraulicsMod, only: SoilResistance, PlantResistance
    use LeafPhotosynthesisMod, only : PhotosynthesisParam
    use CanopyNitrogenProfileMod, only : CanopyNitrogenProfile
    use LongwaveRadiationMod, only : LongwaveRadiation
    use LeafTemperatureMod, only : LeafHeatCapacity
    use LeafBoundaryLayerMod, only : LeafBoundaryLayer
    use LeafFluxesMod, only : LeafFluxes
    use SoilFluxesMultilayerMod, only : SoilFluxesMultilayer
    use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq ! jymody
    use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt ! jymody
    !use clm_time_manager,  only : get_curr_calday ! jymody
    use CanopyTurbulenceMod, only : CanopyTurbulence, CanopyTurbulenceDummy
    use UrbanParamsType , only : urbanparams_type
	use MathToolsMod    , only : beta_function
	use GridcellType        , only : grc
	use shr_orb_mod   , only : shr_orb_decl, shr_orb_cosz
	use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
	use CanopyStateType     , only : canopystate_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_exposedvegp           ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:)     ! CLM patch filter for non-snow-covered vegetation
    integer, intent(in) :: num_nolakec               ! Number of non-lake points in CLM column filter
    integer, intent(in) :: filter_nolakec(:)         ! CLM column filter for non-lake points
    real(r8), intent(in) :: declin

    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    type(waterflux_type)   , intent(inout) :: waterflux_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(frictionvel_type) , intent(inout) :: frictionvel_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    type(surfalb_type)     , intent(inout) :: surfalb_inst
    type(mlcanopy_type)    , intent(inout) :: mlcanopy_inst
    type(urbanparams_type)    , intent(inout) :: urbanparams_inst
    class(soil_water_retention_curve_type) , intent(in)            :: soil_water_retention_curve
    type(canopystate_type) , intent(inout)    :: canopystate_inst

    !
    ! !LOCAL VARIABLES:
    integer  :: i
    integer  :: f                                    ! Filter index
    integer  :: p                                    ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                    ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                    ! Gridcell index for CLM g/l/c/p hierarchy
    integer  :: ic                                   ! Aboveground layer index
    integer  :: il                                   ! Sunlit (1) or shaded (2) leaf index
    integer  :: nstep                                ! Current model time step number
    real(r8) :: dtime                                ! Model time step (s)
    integer  :: num_sub_steps                        ! Number of sub-time steps
    integer  :: niter                                ! Current sub-time step
    integer  :: fn                                   ! number of values in vegetated patch filter
    logical  :: perchroot     = .false.
    integer  :: filterc_tmp(1:num_exposedvegp)

    ! true  => btran is based on active layer (defined over two years);
    ! false => btran is based on currently unfrozen levels
    logical  :: perchroot_alt = .false.
    integer  :: jtop(bounds%begc:bounds%endc)            ! lbning

    real(r8) :: cv (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! CLM: soil heat capacity (J/m2/K)
    real(r8) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! CLM: soil thermal conductivity at layer interface (W/m/K)
    real(r8) :: tk_h2osfc(bounds%begc:bounds%endc)                 ! CLM: thermal conductivity of h2osfc (W/m/K)

    real(r8) :: qflx_prec_intr_sum(bounds%begp:bounds%endp)
    real(r8) :: ircan_sum(bounds%begp:bounds%endp)
    real(r8) :: ustar_sum(bounds%begp:bounds%endp)
    real(r8) :: irsoi_sum(bounds%begp:bounds%endp)
    real(r8) :: rnsoi_sum(bounds%begp:bounds%endp)
    real(r8) :: shsoi_sum(bounds%begp:bounds%endp)
    real(r8) :: lhsoi_sum(bounds%begp:bounds%endp)
    real(r8) :: etsoi_sum(bounds%begp:bounds%endp)
    real(r8) :: gsoi_sum(bounds%begp:bounds%endp)
    real(r8) :: ga_sum(bounds%begp:bounds%endp,0:nlevcanml)
    real(r8) :: shair_sum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: etair_sum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: cfair_sum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: stair_sum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: irleaf_sum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: rnleaf_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: stleaf_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: shleaf_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: lhleaf_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: trleaf_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: evleaf_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: an_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: ag_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: gs_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: gbv_sum(bounds%begp:bounds%endp,1:nlevcanml,1:nleaf)
    real(r8) :: gs_veg(bounds%begp)
    real(r8) :: gbv_veg(bounds%begp)
    real(r8) :: ga_veg(bounds%begp)
    real(r8) :: minlwp
    real(r8) :: zwinput(1:nlevcanml)
    real(r8) :: zsinput(1:nlevcanml)
    real(r8) :: betap
    real(r8) :: betaq
    real(r8) :: betav
    real(r8) :: betaf
    real(r8) :: dlaiinput(1:nlevcanml),dlaiinput1(1:nlevcanml),dlaiinput2(1:nlevcanml)
    real(r8) :: dpaiinput(1:nlevcanml),dpaiinput1(1:nlevcanml),dpaiinput2(1:nlevcanml)
    real(r8) :: fracsunsum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: fracshasum(bounds%begp:bounds%endp,1:nlevcanml)
    real(r8) :: coszen_patch(1:num_exposedvegp)
    !real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8) :: dayl_factor(bounds%begp:bounds%endp) ! scalar (0-1) for daylength effect on Vcmax
    integer  :: bbf
    integer  :: LADtype = .true.  ! true = double canopy & false = single canopy
    real(r8) :: totalsum,totalsums, betax, betax1, betax2
    real(r8) :: zuref = 33._r8
    real(r8) :: r = 0.65_r8
    real(r8) :: H = 13.5_r8
    !---------------------------------------------------------------------

    associate ( &
         dayl                 => grc%dayl                                  , & ! Input:  [real(r8) (:)   ]  daylength (s)
         max_dayl             => grc%max_dayl                              , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)
                                                    ! *** Miscellaneous input ***
    begc        => bounds%begc                 , &  ! Column index
    endc        => bounds%endc                 , &  ! Column index
    !minlwp      => pftcon%minlwp               , &  ! Minimum leaf water potential (MPa)
    lai         => mlcanopy_inst%lai           , &  ! Leaf area index of canopy (m2/m2)
    sai         => mlcanopy_inst%sai           , &  ! Stem area index of canopy (m2/m2)
    ncan        => mlcanopy_inst%ncan          , &  ! Number of aboveground layers
    nbot        => mlcanopy_inst%nbot          , &  ! Index for bottom leaf layer
    ntop        => mlcanopy_inst%ntop          , &  ! Index for top leaf layer
    dpai        => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    tacclim     => mlcanopy_inst%tacclim       , &  ! Average air temperature for acclimation (K)
    btran       => mlcanopy_inst%btran         , &  ! Ball-Berry soil wetness factor (-)
    swleaf      => mlcanopy_inst%swleaf        , &  ! Leaf absorbed solar radiation (W/m2 leaf)
                                                    ! *** Atmospheric forcing data ***
    zref        => mlcanopy_inst%zref          , &  ! Reference height (m)
    zref_old    => mlcanopy_inst%zref_old      , &  ! Reference height for previous timestep (m)
    tref        => mlcanopy_inst%tref          , &  ! Air temperature at reference height (K)
    qref        => mlcanopy_inst%qref          , &  ! Specific humidity at reference height (kg/kg)
    rhref       => mlcanopy_inst%rhref         , &  ! Relative humidity at reference height (%) 
    uref        => mlcanopy_inst%uref          , &  ! Wind speed at reference height (m/s)
    uref33        => mlcanopy_inst%uref33          , &  ! Wind speed at measurement height (m/s)
    pref        => mlcanopy_inst%pref          , &  ! Air pressure at reference height (Pa)
    swskyb      => mlcanopy_inst%swskyb        , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd        , &  ! Atmospheric diffuse solar radiation (W/m2)
    irsky       => mlcanopy_inst%irsky         , &  ! Atmospheric longwave radiation (W/m2)
    qflx_rain   => mlcanopy_inst%qflx_rain     , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow   => mlcanopy_inst%qflx_snow     , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
    co2ref      => mlcanopy_inst%co2ref        , &  ! Atmospheric CO2 at reference height (umol/mol)
    o2ref       => mlcanopy_inst%o2ref         , &  ! Atmospheric O2 at reference height (mmol/mol)
                                                    ! *** Derived atmospheric data ***
    thref       => mlcanopy_inst%thref         , &  ! Atmospheric potential temperature (K)
    thvref      => mlcanopy_inst%thvref        , &  ! Atmospheric virtual potential temperature (K)
    eref        => mlcanopy_inst%eref          , &  ! Vapor pressure at reference height (Pa)
    rhoair      => mlcanopy_inst%rhoair        , &  ! Air density at reference height (kg/m3)
    rhomol      => mlcanopy_inst%rhomol        , &  ! Molar density at reference height (mol/m3)
    mmair       => mlcanopy_inst%mmair         , &  ! Molecular mass of air at reference height (kg/mol)
    cpair       => mlcanopy_inst%cpair         , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
                                                    ! *** Calculated ***
    sumpai      => mlcanopy_inst%sumpai        , &  ! Cumulative plant area index (m2/m2)
    tair        => mlcanopy_inst%tair          , &  ! Air temperature profile (K)
    eair        => mlcanopy_inst%eair          , &  ! Vapor pressure profile (Pa)
    cair        => mlcanopy_inst%cair          , &  ! Atmospheric CO2 profile (umol/mol)
    tveg        => mlcanopy_inst%tveg          , &  ! Vegetation temperature profile (K)
    tvegsun     => mlcanopy_inst%tvegsun         , &  ! Sensible heat flux, vegetation (W/m2)
    tvegsha     => mlcanopy_inst%tvegsha         , &  ! Sensible heat flux, vegetation (W/m2)
    tair_old    => mlcanopy_inst%tair_old      , &  ! Air temperature profile for previous timestep (K)
    eair_old    => mlcanopy_inst%eair_old      , &  ! Vapor pressure profile for previous timestep (Pa)
    cair_old    => mlcanopy_inst%cair_old      , &  ! Atmospheric CO2 profile for previous timestep (umol/mol)
    tveg_old    => mlcanopy_inst%tveg_old      , &  ! Vegetation temperature profile for previous timestep (K)
    tleaf       => mlcanopy_inst%tleaf         , &  ! Leaf temperature (K)
    tleaf_old   => mlcanopy_inst%tleaf_old     , &  ! Leaf temperature for previous timestep (K)
    shair       => mlcanopy_inst%shair         , &  ! Canopy air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair         , &  ! Canopy air water vapor flux (mol H2O/m2/s)
    cfair       => mlcanopy_inst%cfair         , &  ! Canopy air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair         , &  ! Canopy air storage heat flux (W/m2)
    irleaf      => mlcanopy_inst%irleaf        , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf        , &  ! Leaf net radiation (W/m2 leaf)
    stleaf      => mlcanopy_inst%stleaf        , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf        , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf        , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf        , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf        , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    wind_1st    => mlcanopy_inst%wind_1st      , &
    wind_1stold => mlcanopy_inst%wind_1stold   , &
    ag          => mlcanopy_inst%ag            , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an          => mlcanopy_inst%an            , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    gs          => mlcanopy_inst%gs            , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    gbv         => mlcanopy_inst%gbv           , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    fracsun     => mlcanopy_inst%fracsun       , &  ! Sunlit fraction of canopy layer
    fracsha     => mlcanopy_inst%fracsha       , &  ! Shaded fraction of canopy layer
    psil        => mlcanopy_inst%psil          , &  ! Leaf water potential (MPa)
    lwp         => mlcanopy_inst%lwp           , &  ! Leaf water potential of canopy layer (MPa)
    fracminlwp  => mlcanopy_inst%fracminlwp    , &  ! Fraction of canopy with lwp < minlwp
    irsoi       => mlcanopy_inst%irsoi         , &  ! Absorbed longwave radiation, ground (W/m2)
    rnsoi       => mlcanopy_inst%rnsoi         , &  ! Net radiation, ground (W/m2)
    shsoi       => mlcanopy_inst%shsoi         , &  ! Sensible heat flux, ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi         , &  ! Latent heat flux, ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi          , &  ! Soil heat flux (W/m2)
    etsoi       => mlcanopy_inst%etsoi         , &  ! Water vapor flux, ground (mol H2O/m2/s)
    ircan       => mlcanopy_inst%ircan         , &  ! Upward longwave radiation above canopy (W/m2)
    ustar       => mlcanopy_inst%ustar         , &  ! Friction velocity (m/s)
    ga_prof     => mlcanopy_inst%ga_prof       , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    qflx_prec_intr => mlcanopy_inst%qflx_prec_intr , &  ! Intercepted precipitation (kg H2O/m2/s)

    ztop39           =>  mlcanopy_inst%ztop39         , &
	ztop           =>  mlcanopy_inst%ztop         , &  ! Canopy height (m)
	nuforce           =>  mlcanopy_inst%nuforce         , &
	zuforce           =>  mlcanopy_inst%zuforce         , &
	root_biomass   =>  mlcanopy_inst%root_biomass , & 
	dlai           =>  mlcanopy_inst%dlai         , & 
	dsai           =>  mlcanopy_inst%dsai         , & 
	zs             =>  mlcanopy_inst%zs           , & 
	zw             =>  mlcanopy_inst%zw           , & 

	solar_zen      =>  mlcanopy_inst%solar_zen    , & 
	uforc      	  =>  mlcanopy_inst%uforc      	  , & 
                                                               ! *** CLM variables ***
    forc_hgt    => atm2lnd_inst%forc_hgt_grc              , &  ! CLM: Atmospheric reference height (m)
    forc_u      => atm2lnd_inst%forc_u_grc                , &  ! CLM: Atmospheric wind speed in east direction (m/s)
    forc_v      => atm2lnd_inst%forc_v_grc                , &  ! CLM: Atmospheric wind speed in north direction (m/s)
    forc_rh     => atm2lnd_inst%forc_rh_grc               , &  ! CLM: Atmospheric relative humidity (%)
    forc_pco2   => atm2lnd_inst%forc_pco2_grc             , &  ! CLM: Atmospheric CO2 partial pressure (Pa)
    forc_po2    => atm2lnd_inst%forc_po2_grc              , &  ! CLM: Atmospheric O2 partial pressure (Pa)
    forc_solad  => atm2lnd_inst%forc_solad_grc            , &  ! CLM: Atmospheric direct beam radiation (W/m2)
    forc_solai  => atm2lnd_inst%forc_solai_grc            , &  ! CLM: Atmospheric diffuse radiation (W/m2)
    forc_t      => atm2lnd_inst%forc_t_downscaled_col     , &  ! CLM: Atmospheric temperature (K)
    forc_q      => atm2lnd_inst%forc_q_downscaled_col     , &  ! CLM: Atmospheric specific humidity (kg/kg)
    forc_pbot   => atm2lnd_inst%forc_pbot_downscaled_col  , &  ! CLM: Atmospheric pressure (Pa)
    forc_lwrad  => atm2lnd_inst%forc_lwrad_downscaled_col , &  ! CLM: Atmospheric longwave radiation (W/m2)
    forc_rain   => atm2lnd_inst%forc_rain_downscaled_col  , &  ! CLM: Rainfall rate (mm/s)
    forc_snow   => atm2lnd_inst%forc_snow_downscaled_col  , &  ! CLM: Snowfall rate (mm/s)
    btran_patch => energyflux_inst%btran_patch            , &  ! CLM: Transpiration wetness factor (0 to 1)
    t_a10_patch => temperature_inst%t_a10_patch           , &  ! CLM: 10-day running mean of the 2-m temperature (K)
    htop_patch  => canopystate_inst%htop_patch            , &
    tlai_patch  => canopystate_inst%tlai_patch            , &
    tsai_patch  => canopystate_inst%tsai_patch            , &

         watsat                 => soilstate_inst%watsat_col                    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)
         watdry                 => soilstate_inst%watdry_col                    , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=0                      (constant)
         watopt                 => soilstate_inst%watopt_col                    , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=1                      (constant)
         eff_porosity           => soilstate_inst%eff_porosity_col              , & ! Output: [real(r8) (:,:) ]  effective soil porosity
         soilbeta               => soilstate_inst%soilbeta_col                  , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity
         rootr                  => soilstate_inst%rootr_patch                   , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         h2osoi_ice             => waterstate_inst%h2osoi_ice_col               , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_vol             => waterstate_inst%h2osoi_vol_col               , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
         h2osoi_liq             => waterstate_inst%h2osoi_liq_col               , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_liqvol          => waterstate_inst%h2osoi_liqvol_col            , & ! Output: [real(r8) (:,:) ]  volumetric liquid water (v/v)


         smpso         => pftcon%smpso                      , & ! Input:  soil water potential at full stomatal opening (mm)
         smpsc         => pftcon%smpsc                      , & ! Input:  soil water potential at full stomatal closure (mm)

         t_soisno      => temperature_inst%t_soisno_col     , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         !watsat        => soilstate_inst%watsat_col         , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)
         sucsat        => soilstate_inst%sucsat_col         , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                        (constant)
         bsw           => soilstate_inst%bsw_col            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                         (constant)
         !eff_porosity  => soilstate_inst%eff_porosity_col   , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice
         rootfr        => soilstate_inst%rootfr_patch       , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer
         !rootr         => soilstate_inst%rootr_patch        , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer

         !btrane         => energyflux_inst%btran_patch       , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1) (integrated soil water stress)
         btran2        => energyflux_inst%btran2_patch      , & ! Output: [real(r8) (:)   ]  integrated soil water stress square

         !cpleaf    => mlcanopy_inst%cpleaf  , &

         rresis        => energyflux_inst%rresis_patch      & ! Output: [real(r8) (:,:) ]  root soil water stress (resistance) by layer (0-1)  (nlevgrnd)

         !h2osoi_vol    => waterstate_inst%h2osoi_vol_col    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         !h2osoi_liqvol => waterstate_inst%h2osoi_liqvol_col   & ! Output: [real(r8) (:,:) ]  liquid volumetric moisture, will be used for BeTR+
         !thk          =>    soilstate_inst%thk_col
    )


!write (iulog,*) 'num_exposedvegp=',num_exposedvegp, bounds%endp, bounds%endc, bounds%endg

    do f = 1, num_exposedvegp!p = bounds%begp,bounds%endp
       p = filter_exposedvegp(f)
       c = patch%column(p)
       g = patch%gridcell(p)

       !write (iulog,*) 'num_exposedvegp=',num_exposedvegp, 'p=',p,'c=',c,'g=',g

 !       btran(p)  = 0.0_r8
 !       dayl_factor(p)=min(1._r8,max(0.01_r8,(dayl(g)*dayl(g))/(max_dayl(g)*max_dayl(g))))
	    zref(p) = forc_hgt(g)!-forc_hgt(g)/float(nlevcanml)
		!ztop		(p) =htop_patch(p)
		lai			(p) =tlai_patch(p)
		sai			(p) =tsai_patch(p)
		root_biomass(p) = 500._r8 ! bonan 2014        
		ncan        (p) = nlevcanml
		nbot        (p) = 1

        htop_patch(p) = 43.725_r8 ! htop+H*(1-r) = 39+13.5*(1-0.65)

		!bbf =      int(ceiling(nlevcanml-ztop(p)/(zref(p)/nlevcanml))) ! nint
		!ntop        (p) = nlevcanml-bbf
		if (zref(p)+0.001_r8 < htop_patch(p)) then
		    write(iulog,*) zref(p) , htop_patch(p)
		    call endrun (msg=' zref(p) < htop_patch(p)')
		end if
		if (abs(zref(p)-htop_patch(p))<nlevcanml/zref(p)) then
		    ntop(p) = nlevcanml - 1
		else
		    ntop(p) = nint(htop_patch(p)*nlevcanml/zref(p)+0.4999_r8)
		end if
        nuforce(p) = nint(zuref*nlevcanml/zref(p)+0.4999_r8)
        write(iulog,*) 'zref(p)=',zref(p) , 'htop_patch(p)=',htop_patch(p)
		
		zwinput(1:nlevcanml) = (/ ((float(i)/float(nlevcanml)*zref(p)), i = 1, nlevcanml) /) ! interface
		zsinput(1:nlevcanml) = (/ (((float(i)-0.5_r8)/float(nlevcanml)*zref(p)), i = 1, nlevcanml) /) ! mid

        ztop(p) = zwinput(ntop(p))
        zuforce(p) = zwinput(nuforce(p))

        ztop39(p) = 39._r8

        ! shape parameters

		betap = 7.5_r8 !11.5_r8
		betaq = 2.75_r8 !3.5_r8
		betav = beta_function (betap, betaq)

        dlaiinput(:) = 0.0_r8 ! mid
        dpaiinput(:) = 0.0_r8 ! mid
        dlaiinput1(:) = 0.0_r8 ! mid
        dpaiinput1(:) = 0.0_r8 ! mid
        dlaiinput2(:) = 0.0_r8 ! mid
        dpaiinput2(:) = 0.0_r8 ! mid

        if (sai(p)/zwinput(ntop(p)) < 0.001_r8) then
             call endrun (msg=' ERROR: slai each layer: less than 0.001')
        end if

        totalsum = 0._r8
        totalsums = 0._r8

		do ic = nbot(p), ntop(p)

            if (LADtype) then
                betax1 = (zsinput(ic)+(r-1._r8)*H)/zwinput(ntop(p))
                betax =  min(max(betax1,0._r8),1._r8);

                dlaiinput1(ic) =  lai(p)/zwinput(ntop(p))*betax**(betap-1._r8) * &
                                (1._r8-betax)**(betaq-1._r8)/betav

                betax2 = (zsinput(ic)+r*H)/zwinput(ntop(p))
                betax =  min(max(betax2,0._r8),1._r8);
                dlaiinput2(ic) =  lai(p)/zwinput(ntop(p))*betax**(betap-1._r8) * &
                                (1._r8-betax)**(betaq-1._r8)/betav

                dlaiinput(ic) = dlaiinput1(ic)*r+dlaiinput2(ic)*(1._r8-r)

                if ((0._r8<betax1) .and. (betax1<1._r8)) then
                    dpaiinput1(ic) = sai(p)/zwinput(ntop(p))
                end if
                if ((0._r8<betax2) .and. (betax2<1._r8)) then
                    dpaiinput2(ic) = sai(p)/zwinput(ntop(p))
                end if
                dpaiinput(ic) = dpaiinput1(ic)*r+dpaiinput2(ic)*(1._r8-r)
            else
                dlaiinput(ic) =  lai(p)/zwinput(ntop(p))*(zsinput(ic)/zwinput(ntop(p)))**(betap-1._r8) * &
                                (1._r8-zsinput(ic)/zwinput(ntop(p)))**(betaq-1._r8)/betav
                dpaiinput(ic) = sai(p)/zwinput(ntop(p))
            end if



		    if (dlaiinput(ic)  < 0.0001_r8) then
		        dlaiinput(ic) = 0._r8
		    end if
		    totalsum = totalsum + dlaiinput(ic)
		    totalsums = totalsums + dpaiinput(ic)
		end do
		


		zs  (p,0) = (0._r8-0.5_r8)/float(nlevcanml)*zref(p)
		zw  (p,0) = (0._r8/float(nlevcanml)*zref(p))
		do ic = 1, nlevcanml
			dlai(p,ic) = dlaiinput(ic)*lai(p)/totalsum
			if (totalsums > 0._r8) then
			    dsai(p,ic) = dpaiinput(ic)*sai(p)/totalsums
			else
			    dsai(p,ic) = 0._r8
			end if
			dpai(p,ic) = dlai(p,ic)+dsai(p,ic)
			zs  (p,ic) = zsinput(ic)
			zw  (p,ic) = zwinput(ic)
		end do
		
       coszen_patch(p) = shr_orb_cosz (get_curr_calday(), grc%lat(g), grc%lon(g), declin)
       !write(iulog,*) 'coszen_patch(p)=' ,coszen_patch(p),acos(coszen_patch(p))
		solar_zen(p)                   = acos(min(max(coszen_patch(p), 0._r8) ,1._r8))

		!tair_old (p,0:nlevcanml)
		!eair_old (p,0:nlevcanml)
		!cair_old (p,0:nlevcanml)
		!tveg_old (p,0:nlevcanml,1:nleaf)
		!tleaf_old(p,1:nlevcanml,1:nleaf)
		uforc    (p)                   = forc_u(g)
	end do


    minlwp = -2._r8
    ! Get current step (counter) and step size (seconds)

    nstep = get_nstep()
    dtime = get_step_size()
    if (nstep < 2) then
        write(iulog,*) 'dpai =', dpai
    end if
    ! Set time step for sub-stepping of flux calculations and number of sub-steps

    !dtime_sub = 5._r8 * 60._r8
    num_sub_steps = int(dtime / dtime_sub)

    !---------------------------------------------------------------------
    ! Copy CLM variables to multilayer canopy variables
    !---------------------------------------------------------------------

     write(iulog,*) 'MLstart -> ---------------------------------------  '
    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       c = patch%column(p)
       g = patch%gridcell(p)
		
		!tacclim  (p)                   = btran_patch(p) !jymody
       ! Atmospheric forcing: CLM grid cell (g) variables -> patch (p) variables

       !zref(p) = forc_hgt(g)
       uref33(p) = max (0.001_r8, sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       rhref(p) = forc_rh(g)
       swskyb(p,ivis) = forc_solad(g,ivis)
       swskyb(p,inir) = forc_solad(g,inir)
       swskyd(p,ivis) = forc_solai(g,ivis)
       swskyd(p,inir) = forc_solai(g,inir)

       ! Atmospheric forcing: CLM column (c) variables -> patch (p) variables

       tref(p) = forc_t(c)
       qref(p) = forc_q(c)
       pref(p) = forc_pbot(c)
       irsky(p) = forc_lwrad(c)
       qflx_rain(p) = forc_rain(c)
       qflx_snow(p) = forc_snow(c)

       ! CO2 and O2: note unit conversion

       co2ref(p) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8  ! Pa -> umol/mol
       o2ref(p)  = forc_po2(g) / forc_pbot(c) * 1.e03_r8   ! Pa -> mmol/mol

       ! Miscellaneous

       btran(p) = btran_patch(p)
       tacclim(p) = t_a10_patch(p)

       ! Check to see if forcing height has changed

       if (zref(p) /= zref_old(p)) then
             write(iulog,*) zref(p), zref_old(p)
          call endrun (msg=' ERROR: CanopyFluxesMultilayer: forcing height is not constant')
       end if
       zref_old(p) = zref(p)
       !write(iulog,*) 'Input(0-1)  swskyb(p,ivis)', swskyb(p,ivis), 'swskyb(p,inir)=',swskyb(p,inir), 'swskyd(p,ivis)=',swskyd(p,ivis), 'swskyd(p,inir)=', swskyd(p,inir)
       !write(iulog,*) 'Input(0-2)  tref(p)', tref(p), 'rhref(p)=',rhref(p), 'uref(p)=',uref(p), 'qref(p)=', qref(p),'pref(p)=', pref(p),'irsky(p)=',irsky(p)

    end do

    !---------------------------------------------------------------------
    ! Derived atmospheric input
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       eref(p) = qref(p) * pref(p) / (mmh2o / mmdry + (1._r8 - mmh2o / mmdry) * qref(p))  ! Vapor pressure at reference height (Pa)
       rhomol(p) = pref(p) / (rgasc * tref(p))                                            !Molar density at reference height (mol/m3)
       rhoair(p) = rhomol(p) * mmdry * (1._r8 - (1._r8 - mmh2o/mmdry) * eref(p) / pref(p))! Air density at reference height (kg/m3)
       mmair(p) = rhoair(p) / rhomol(p)  !Molecular mass of air at reference height (kg/mol)
       cpair(p) = cpd * (1._r8 + (cpw/cpd - 1._r8) * qref(p)) * mmair(p)
       thref(p) = tref(p) + 0.0098_r8 * zref(p)
       thvref(p) = thref(p) * (1._r8 + 0.61_r8 * qref(p))

       !write(iulog,*) 'Input(1)  eref(p)=', eref(p), 'rhomol(p)=',rhomol(p),'rhoair(p)=',rhoair(p), 'mmair(p)=',mmair(p), 'cpair(p)=', cpair(p),'thref(p)=', thref(p),'thvref(p)=',thvref(p)
    end do

    !---------------------------------------------------------------------
    ! Cumulative plant area index (lai+sai)
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Layers above the canopy have no vegetation

       do ic = ntop(p)+1, ncan(p)
          sumpai(p,ic) = 0._r8
       end do

       ! Fill in canopy layers (at the midpoint), starting from the top

       do ic = ntop(p), 1, -1
          if (ic == ntop(p)) then
             sumpai(p,ic) = 0.5_r8 * dpai(p,ic)
          else
             sumpai(p,ic) = sumpai(p,ic+1) + 0.5_r8 * (dpai(p,ic+1) + dpai(p,ic))
          end if
       end do
      ! write(iulog,*) 'sumpai=', sumpai(p,:)
    end do



    !---------------------------------------------------------------------
    ! Soil thermal conductivity and heat capacity - Do not need the returned
    ! values. Only need thk(c,snl(c)+1), which is the thermal conductivity of
    ! the first snow/soil layer.
    !---------------------------------------------------------------------

    !call SoilThermProp (bounds, num_nolakec, filter_nolakec, tk(begc:endc,:), &
    !cv(begc:endc,:), tk_h2osfc(begc:endc), waterstate_inst, temperature_inst, &
    !soilstate_inst)

    call SoilThermProp (bounds, num_nolakec, filter_nolakec, tk(begc:endc,:), &
    cv(begc:endc,:), tk_h2osfc(begc:endc), urbanparams_inst, temperature_inst, waterstate_inst, soilstate_inst)

    !write(iulog,*) 'tk=',tk, 'cv=',cv,'tk_h2osfc=',tk_h2osfc
    !---------------------------------------------------------------------
    ! Solar radiation transfer through the canopy
    !---------------------------------------------------------------------

    call SolarRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
    surfalb_inst, mlcanopy_inst)

    !---------------------------------------------------------------------
    ! Plant hydraulics
    !---------------------------------------------------------------------

    call SoilResistance (num_exposedvegp, filter_exposedvegp, &
    soilstate_inst, waterstate_inst, mlcanopy_inst)

    call PlantResistance (num_exposedvegp, filter_exposedvegp, mlcanopy_inst)

    !---------------------------------------------------------------------
    ! Canopy profiles of photosynthetic capacity
    !---------------------------------------------------------------------

    call CanopyNitrogenProfile (num_exposedvegp, filter_exposedvegp, mlcanopy_inst)

    !---------------------------------------------------------------------
    ! Use sub-stepping to calculate fluxes over the full time step
    !---------------------------------------------------------------------

    ! Initialize fluxes that are summed over sub-time steps

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       qflx_prec_intr_sum(p) = 0._r8
       ircan_sum(p) = 0._r8
       ustar_sum(p) = 0._r8
       irsoi_sum(p) = 0._r8
       rnsoi_sum(p) = 0._r8
       shsoi_sum(p) = 0._r8
       lhsoi_sum(p) = 0._r8
       etsoi_sum(p) = 0._r8
       gsoi_sum(p) = 0._r8
       ga_sum(p,0) = 0._r8
       do ic = 1, ncan(p)
          irleaf_sum(p,ic) = 0._r8
          shair_sum(p,ic) = 0._r8
          etair_sum(p,ic) = 0._r8
          cfair_sum(p,ic) = 0._r8
          stair_sum(p,ic) = 0._r8
          ga_sum(p,ic) = 0._r8
             fracsunsum(p,ic) = 0._r8
             fracshasum(p,ic) =  0._r8
          do il = 1, nleaf
             rnleaf_sum(p,ic,il) = 0._r8
             stleaf_sum(p,ic,il) = 0._r8
             shleaf_sum(p,ic,il) = 0._r8
             lhleaf_sum(p,ic,il) = 0._r8
             trleaf_sum(p,ic,il) = 0._r8
             evleaf_sum(p,ic,il) = 0._r8
             an_sum(p,ic,il) = 0._r8
             ag_sum(p,ic,il) = 0._r8
             gs_sum(p,ic,il) = 0._r8
             gbv_sum(p,ic,il) = 0._r8
          end do
       end do
    end do

    ! Sub-stepping loop

    do niter = 1, num_sub_steps

       ! Save values for previous timestep

       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          do ic = 1, ncan(p)
             tleaf_old(p,ic,isun) = tleaf(p,ic,isun)
             tleaf_old(p,ic,isha) = tleaf(p,ic,isha)
          end do

          do ic = 0, ncan(p)
             tair_old(p,ic) = tair(p,ic)
             !wind_1stold(p,ic) = wind_1st(p,ic)
             eair_old(p,ic) = eair(p,ic)
             cair_old(p,ic) = cair(p,ic)
             tveg_old(p,ic,isun) = tveg(p,ic,isun)
             tveg_old(p,ic,isha) = tveg(p,ic,isha)
          end do
       end do

       ! Canopy interception

       call CanopyInterception (num_exposedvegp, filter_exposedvegp, mlcanopy_inst)

       ! Longwave radiation transfer through canopy

       call LongwaveRadiation (bounds, num_exposedvegp, filter_exposedvegp, mlcanopy_inst)

       ! Loop through all patches to calculate fluxes

       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)

          ! Photosynthesis parameters
 
          call PhotosynthesisParam (p, mlcanopy_inst)

          ! Leaf fluxes for each canopy layer

          do ic = 1, ncan(p)

             ! Leaf heat capacity

             call LeafHeatCapacity (p, ic, mlcanopy_inst)

             ! Sunlit leaves: net radiation, boundary layer conductances and fluxes

             rnleaf(p,ic,isun) = swleaf(p,ic,isun,ivis) + swleaf(p,ic,isun,inir) + irleaf(p,ic)
             call LeafBoundaryLayer (p, ic, isun, mlcanopy_inst)
             call LeafFluxes (p, ic, isun, mlcanopy_inst)

             ! Shaded leaves: net radiation, boundary layer conductances and fluxes

             rnleaf(p,ic,isha) = swleaf(p,ic,isha,ivis) + swleaf(p,ic,isha,inir) + irleaf(p,ic)
             call LeafBoundaryLayer (p, ic, isha, mlcanopy_inst)
             call LeafFluxes (p, ic, isha, mlcanopy_inst)

          end do

         ! write(iulog,*) 'cpleaf',cpleaf(p,:)  1861.3333333333330

          ! Soil fluxes

          call SoilFluxesMultilayer (p, soilstate_inst, temperature_inst, &
          energyflux_inst, waterflux_inst, mlcanopy_inst)

          ! Canopy turbulence, aeorodynamic conductances, and scalar profiles using above-
          ! and within-canopy coupling with a roughness sublayer (RSL) parameterization

          select case (turb_type)
             case (0)
                call CanopyTurbulenceDummy (p, mlcanopy_inst)
             case (1:4)
                call CanopyTurbulence (p, niter, soilstate_inst, temperature_inst, &
                energyflux_inst, waterflux_inst, mlcanopy_inst)
             case default
                call endrun (msg=' ERROR: CanopyFluxesMultilayer: turbulence type not valid')
          end select

          ! Update canopy intercepted water for evaporation and dew

          call CanopyEvaporation (p, mlcanopy_inst)

          ! Canopy water and leaf energy fluxes need to be accumulated over all
          ! sub-time steps because canopy water is updated every time step.
          ! All other variables are instantaneous for the final time step.

          qflx_prec_intr_sum(p) = qflx_prec_intr_sum(p) + qflx_prec_intr(p)
          ircan_sum(p) = ircan_sum(p) + ircan(p)
          ustar_sum(p) = ustar_sum(p) + ustar(p)
          irsoi_sum(p) = irsoi_sum(p) + irsoi(p)
          rnsoi_sum(p) = rnsoi_sum(p) + rnsoi(p)
          shsoi_sum(p) = shsoi_sum(p) + shsoi(p)
          lhsoi_sum(p) = lhsoi_sum(p) + lhsoi(p)
          etsoi_sum(p) = etsoi_sum(p) + etsoi(p)
          gsoi_sum(p) = gsoi_sum(p) + gsoi(p)
          ga_sum(p,0) = ga_sum(p,0) + ga_prof(p,0)

          do ic = 1, ncan(p)
             irleaf_sum(p,ic) = irleaf_sum(p,ic) + irleaf(p,ic)
             shair_sum(p,ic) = shair_sum(p,ic) + shair(p,ic)
             etair_sum(p,ic) = etair_sum(p,ic) + etair(p,ic)
             cfair_sum(p,ic) = cfair_sum(p,ic) + cfair(p,ic)
             stair_sum(p,ic) = stair_sum(p,ic) + stair(p,ic)
             ga_sum(p,ic) = ga_sum(p,ic) + ga_prof(p,ic)
             fracsunsum(p,ic) = fracsunsum(p,ic) + fracsun(p,ic)
             fracshasum(p,ic) = fracshasum(p,ic) + fracsha(p,ic)
             do il = 1, nleaf
                rnleaf_sum(p,ic,il) = rnleaf_sum(p,ic,il) + rnleaf(p,ic,il)
                stleaf_sum(p,ic,il) = stleaf_sum(p,ic,il) + stleaf(p,ic,il)
                shleaf_sum(p,ic,il) = shleaf_sum(p,ic,il) + shleaf(p,ic,il)
                lhleaf_sum(p,ic,il) = lhleaf_sum(p,ic,il) + lhleaf(p,ic,il)
                trleaf_sum(p,ic,il) = trleaf_sum(p,ic,il) + trleaf(p,ic,il)
                evleaf_sum(p,ic,il) = evleaf_sum(p,ic,il) + evleaf(p,ic,il)
                an_sum(p,ic,il) = an_sum(p,ic,il) + an(p,ic,il)
                ag_sum(p,ic,il) = ag_sum(p,ic,il) + ag(p,ic,il)
                gs_sum(p,ic,il) = gs_sum(p,ic,il) + gs(p,ic,il)
                gbv_sum(p,ic,il) = gbv_sum(p,ic,il) + gbv(p,ic,il)
             end do
          end do



       end do ! End patch loop
    end do    ! End sub-stepping loop

    !---------------------------------------------------------------------
    ! Sum leaf and soil fluxes
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       qflx_prec_intr(p) = qflx_prec_intr_sum(p) / float(num_sub_steps)
       ircan(p) = ircan_sum(p) / float(num_sub_steps)
       ustar(p) = ustar_sum(p) / float(num_sub_steps)
       irsoi(p) = irsoi_sum(p) / float(num_sub_steps)
       rnsoi(p) = rnsoi_sum(p) / float(num_sub_steps)
       shsoi(p) = shsoi_sum(p) / float(num_sub_steps)
       lhsoi(p) = lhsoi_sum(p) / float(num_sub_steps)
       etsoi(p) = etsoi_sum(p) / float(num_sub_steps)
       gsoi(p) = gsoi_sum(p) / float(num_sub_steps)
       ga_prof(p,0) = ga_sum(p,0) / float(num_sub_steps)


       ga_veg(p) = ga_prof(p,0)
       gs_veg(p) = 0._r8
       gbv_veg(p) = 0._r8

       do ic = 1, ncan(p)
          irleaf(p,ic) = irleaf_sum(p,ic) / float(num_sub_steps)
          shair(p,ic) = shair_sum(p,ic) / float(num_sub_steps)
          etair(p,ic) = etair_sum(p,ic) / float(num_sub_steps)
          cfair(p,ic) = cfair_sum(p,ic) / float(num_sub_steps)
          stair(p,ic) = stair_sum(p,ic) / float(num_sub_steps)
          ga_prof(p,ic) = ga_sum(p,ic) / float(num_sub_steps)
          fracsun(p,ic) = fracsunsum(p,ic) / float(num_sub_steps)
          fracsha(p,ic) = fracshasum(p,ic) / float(num_sub_steps)
          ga_veg(p) = ga_veg(p) + ga_prof(p,ic)
          do il = 1, nleaf
             rnleaf(p,ic,il) = rnleaf_sum(p,ic,il) / float(num_sub_steps)
             stleaf(p,ic,il) = stleaf_sum(p,ic,il) / float(num_sub_steps)
             shleaf(p,ic,il) = shleaf_sum(p,ic,il) / float(num_sub_steps)
             lhleaf(p,ic,il) = lhleaf_sum(p,ic,il) / float(num_sub_steps)
             trleaf(p,ic,il) = trleaf_sum(p,ic,il) / float(num_sub_steps)
             evleaf(p,ic,il) = evleaf_sum(p,ic,il) / float(num_sub_steps)
             an(p,ic,il) = an_sum(p,ic,il) / float(num_sub_steps)
             ag(p,ic,il) = ag_sum(p,ic,il) / float(num_sub_steps)
             gs(p,ic,il) = gs_sum(p,ic,il) / float(num_sub_steps)
             gbv(p,ic,il) = gbv_sum(p,ic,il) / float(num_sub_steps)
          end do
          gs_veg(p) = gs_veg(p) +   ( gs(p,ic,isun)*fracsun(p,ic)+gs(p,ic,isha)*fracsha(p,ic))*dpai(p,ic)
          gbv_veg(p) = gbv_veg(p) + (gbv(p,ic,isun)*fracsun(p,ic)+gbv(p,ic,isha)*fracsha(p,ic))*dpai(p,ic)
       end do
       ga_veg(p) = ga_veg(p)/float(ncan(p)+1)
       gs_veg(p) = gs_veg(p)/float(ncan(p))
       gbv_veg(p) = gbv_veg(p)/float(ncan(p))
       write(iulog,*) 'ga_veg(p)=',ga_veg(p),'gs_veg(p)=',gs_veg(p),'gbv_veg(p)=',gbv_veg(p)
       write(iulog,*) 'ga_prof(p,0)=',ga_prof(p,0),'rhomol(p)=',rhomol(p)



       call CanopyFluxesSum (p, mlcanopy_inst)

       !------------------------------------------------------------------
       ! Need to merge temperature and leaf water potential for sunlit and
       ! shaded leaves because sun/shade fractions change over time
       !------------------------------------------------------------------
       do ic  = 0, ncan(p)
          tvegsun(p,ic)=tveg(p,ic,isun)
          tvegsha(p,ic)=tveg(p,ic,isha)
       end do
       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             tveg(p,ic,isun) = tveg(p,ic,isun) * fracsun(p,ic) + tveg(p,ic,isha) * fracsha(p,ic)
             tveg(p,ic,isha) = tveg(p,ic,isun)
          end if
       end do
       do ic = nbot(p), ntop(p)
          tleaf(p,ic,isun) = tleaf(p,ic,isun) * fracsun(p,ic) + tleaf(p,ic,isha) * fracsha(p,ic)
          tleaf(p,ic,isha) = tleaf(p,ic,isun)
       end do

       ! Leaf water potential for each canopy layer and fraction of the canopy
       ! that is water stressed

       select case (gstyp)

       case (0, 1)

       fracminlwp(p) = 0._r8
       do ic = 1, ncan(p)
          lwp(p,ic) = 0.
       end do

       case (2)

       fracminlwp(p) = 0._r8
       do ic = nbot(p), ntop(p)
          lwp(p,ic) = psil(p,ic,isun) * fracsun(p,ic) + psil(p,ic,isha) * fracsha(p,ic)
          if (lwp(p,ic) <= minlwp) then
             fracminlwp(p) = fracminlwp(p) + dpai(p,ic)
          end if
       end do

       if (lai(p) > 0._r8) then
          fracminlwp(p) = fracminlwp(p) / (lai(p) + sai(p))
       end if

       end select

    end do     ! End patch loop

    !---------------------------------------------------------------------
    ! Copy multilayer canopy variables to CLM variables
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

!      temperature_inst%t_ref2m_patch (p)
!      waterstate_inst%q_ref2m_patch (p)
!      energyflux_inst%taux_patch (p)
!      energyflux_inst%tauy_patch (p)
!      energyflux_inst%eflx_sh_tot_patch (p)
!      energyflux_inst%eflx_lh_tot_patch (p)
!      energyflux_inst%eflx_lwrad_out_patch (p)
!      waterflux_inst%qflx_evap_tot_patch (p)
!      frictionvel_inst%u10_clm_patch (p)
!      frictionvel_inst%fv_patch (p)
!      frictionvel_inst%ram1_patch (p)

    end do

    end associate
  end subroutine CanopyFluxesMultilayer

  !-----------------------------------------------------------------------
  subroutine CanopyFluxesSum (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Sum leaf and soil fluxes
    !
    ! !USES:
	use clm_time_manager, only : get_nstep
     use clm_varctl, only : iulog
    use clm_varctl, only : turb_type
    use clm_varpar, only : ivis, inir, isun, isha
    use WaterVaporMod, only : LatVap
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p        ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                  ! Aboveground layer index
	integer  :: nstep
    real(r8) :: err                 ! Energy imbalance
    real(r8) :: radin               ! Incoming radiation
    real(r8) :: radout              ! Outgoing radiation
    real(r8) :: avail               ! Available energy
    real(r8) :: flux                ! Turbulent fluxes + storage
    real(r8) :: fracgreen           ! Green fraction of plant area index: lai/(lai+sai)
    real(r8) :: stleafcan
    real(r8) :: staircan
    !---------------------------------------------------------------------

    associate ( &
    ncan        => mlcanopy_inst%ncan          , &  ! Number of aboveground layers
    nbot        => mlcanopy_inst%nbot          , &  ! Index for bottom leaf layer
    ntop        => mlcanopy_inst%ntop          , &  ! Index for top leaf layer
    dpai        => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    tref        => mlcanopy_inst%tref          , &  ! Air temperature at reference height (K)
    fwet        => mlcanopy_inst%fwet          , &  ! Fraction of plant area index that is wet
    fdry        => mlcanopy_inst%fdry          , &  ! Fraction of plant area index that is green and dry
    shair       => mlcanopy_inst%shair         , &  ! Canopy air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair         , &  ! Canopy air water vapor flux (mol H2O/m2/s)
    cfair       => mlcanopy_inst%cfair         , &  ! Canopy air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair         , &  ! Canopy air storage heat flux (W/m2)
    irleaf      => mlcanopy_inst%irleaf        , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf        , &  ! Leaf net radiation (W/m2 leaf)
    stleaf      => mlcanopy_inst%stleaf        , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf        , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf        , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf        , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf        , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)

    trleafsun   => mlcanopy_inst%trleafsun     , &  ! Leaf ...
    evleafsun   => mlcanopy_inst%evleafsun     , &  ! Leaf ...
    trleafsha   => mlcanopy_inst%trleafsha     , &  ! Leaf ...
    evleafsha   => mlcanopy_inst%evleafsha     , &  ! Leaf ...
    swleaf      => mlcanopy_inst%swleaf        , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    an          => mlcanopy_inst%an            , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    ag          => mlcanopy_inst%ag            , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    fracsun     => mlcanopy_inst%fracsun       , &  ! Sunlit fraction of canopy layer
    fracsha     => mlcanopy_inst%fracsha       , &  ! Shaded fraction of canopy layer
    swveg       => mlcanopy_inst%swveg         , &  ! Absorbed solar radiation, vegetation (W/m2)
    rnsoi       => mlcanopy_inst%rnsoi         , &  ! Net radiation, ground (W/m2)
    shsoi       => mlcanopy_inst%shsoi         , &  ! Sensible heat flux, ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi         , &  ! Latent heat flux, ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi          , &  ! Soil heat flux (W/m2)
    swsoi       => mlcanopy_inst%swsoi         , &  ! Absorbed solar radiation, ground (W/m2)
    irsoi       => mlcanopy_inst%irsoi         , &  ! Absorbed longwave radiation, ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi         , &  ! Water vapor flux, ground (mol H2O/m2/s)
    swskyb      => mlcanopy_inst%swskyb        , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd        , &  ! Atmospheric diffuse solar radiation (W/m2)
    irsky       => mlcanopy_inst%irsky         , &  ! Atmospheric longwave radiation (W/m2)
    albcan      => mlcanopy_inst%albcan        , &  ! Albedo above canopy
    ircan       => mlcanopy_inst%ircan         , &  ! Upward longwave radiation above canopy (W/m2)
                                                    ! *** Output ***
    sw_prof     => mlcanopy_inst%sw_prof       , &  ! Canopy layer absorbed solar radiation (W/m2)
    ir_prof     => mlcanopy_inst%ir_prof       , &  ! Canopy layer absorbed longwave radiation (W/m2)
    rn_prof     => mlcanopy_inst%rn_prof       , &  ! Canopy layer net radiation (W/m2)
    st_prof     => mlcanopy_inst%st_prof       , &  ! Canopy layer storage heat flux (W/m2)
    sh_prof     => mlcanopy_inst%sh_prof       , &  ! Canopy layer sensible heat flux (W/m2)
    lh_prof     => mlcanopy_inst%lh_prof       , &  ! Canopy layer latent heat flux (W/m2)
    et_prof     => mlcanopy_inst%et_prof       , &  ! Canopy layer water vapor flux (mol H2O/m2/s)
    tr_prof     => mlcanopy_inst%tr_prof       , &  ! Canopy layer water vapor flux (mol H2O/m2/s)
    ev_prof     => mlcanopy_inst%ev_prof       , &  ! Canopy layer water vapor flux (mol H2O/m2/s)
    fc_prof     => mlcanopy_inst%fc_prof       , &  ! Canopy layer CO2 flux (umol CO2/m2/s)
    rnet        => mlcanopy_inst%rnet          , &  ! Net radiation (W/m2)
    stflx       => mlcanopy_inst%stflx         , &  ! Canopy storage heat flux (W/m2)
    shflx       => mlcanopy_inst%shflx         , &  ! Sensible heat flux (W/m2)
    lhflx       => mlcanopy_inst%lhflx         , &  ! Latent heat flux (W/m2)
    etflx       => mlcanopy_inst%etflx         , &  ! Water vapor flux (mol H2O/m2/s)
    irveg       => mlcanopy_inst%irveg         , &  ! Absorbed longwave radiation, vegetation (W/m2)
    irvegsun    => mlcanopy_inst%irvegsun      , &  ! Absorbed longwave radiation, sunlit canopy (W/m2)
    irvegsha    => mlcanopy_inst%irvegsha      , &  ! Absorbed longwave radiation, shaded canopy (W/m2)
    !tveg        => mlcanopy_inst%tveg         , &  ! Sensible heat flux, vegetation (W/m2)
    !tvegsun     => mlcanopy_inst%tvegsun         , &  ! Sensible heat flux, vegetation (W/m2)
    !tvegsha     => mlcanopy_inst%tvegsha         , &  ! Sensible heat flux, vegetation (W/m2)
    shveg       => mlcanopy_inst%shveg         , &  ! Sensible heat flux, vegetation (W/m2)
    shvegsun    => mlcanopy_inst%shvegsun      , &  ! Sensible heat flux, sunlit canopy (W/m2)
    shvegsha    => mlcanopy_inst%shvegsha      , &  ! Sensible heat flux, shaded canopy (W/m2)
    lhveg       => mlcanopy_inst%lhveg         , &  ! Latent heat flux, vegetation (W/m2)
    lhvegsun    => mlcanopy_inst%lhvegsun      , &  ! Latent heat flux, sunlit canopy (W/m2)
    lhvegsha    => mlcanopy_inst%lhvegsha      , &  ! Latent heat flux, shaded canopy (W/m2)
    etveg       => mlcanopy_inst%etveg         , &  ! Water vapor flux, vegetation (mol H2O/m2/s)
    evveg       => mlcanopy_inst%evveg         , &  ! Water vapor flux, vegetation (mol H2O/m2/s)
    trveg       => mlcanopy_inst%trveg         , &  ! Water vapor flux, vegetation (mol H2O/m2/s)
    etvegsun    => mlcanopy_inst%etvegsun      , &  ! Water vapor flux, sunlit canopy (mol H2O/m2/s)
    etvegsha    => mlcanopy_inst%etvegsha      , &  ! Water vapor flux, shaded canopy (mol H2O/m2/s)
    gppveg      => mlcanopy_inst%gppveg        , &  ! Gross primary production (umol CO2/m2/s)
    gppvegsun   => mlcanopy_inst%gppvegsun     , &  ! Gross primary production, sunlit canopy (umol CO2/m2/s)
    gppvegsha   => mlcanopy_inst%gppvegsha       &  ! Gross primary production, shaded canopy (umol CO2/m2/s)
    )

    ! Leaf flux profiles
	nstep = get_nstep()
    do ic = 1, ncan(p)
       if (dpai(p,ic) > 0._r8) then
          ir_prof(p,ic) = irleaf(p,ic) * dpai(p,ic)
          sw_prof(p,ic,ivis) = (swleaf(p,ic,isun,ivis)*fracsun(p,ic) + swleaf(p,ic,isha,ivis)*fracsha(p,ic)) * dpai(p,ic)
          sw_prof(p,ic,inir) = (swleaf(p,ic,isun,inir)*fracsun(p,ic) + swleaf(p,ic,isha,inir)*fracsha(p,ic)) * dpai(p,ic)
          rn_prof(p,ic) = (rnleaf(p,ic,isun)*fracsun(p,ic) + rnleaf(p,ic,isha)*fracsha(p,ic)) * dpai(p,ic)
          st_prof(p,ic) = (stleaf(p,ic,isun)*fracsun(p,ic) + stleaf(p,ic,isha)*fracsha(p,ic)) * dpai(p,ic)
          sh_prof(p,ic) = (shleaf(p,ic,isun)*fracsun(p,ic) + shleaf(p,ic,isha)*fracsha(p,ic)) * dpai(p,ic)
          lh_prof(p,ic) = (lhleaf(p,ic,isun)*fracsun(p,ic) + lhleaf(p,ic,isha)*fracsha(p,ic)) * dpai(p,ic)
          et_prof(p,ic) = (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) &  ! Canopy layer water vapor flux (mol H2O/m2/s)
                        + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * fracsha(p,ic) * dpai(p,ic)
          fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
          fc_prof(p,ic) = (an(p,ic,isun)*fracsun(p,ic) + an(p,ic,isha)*fracsha(p,ic)) * dpai(p,ic) * fracgreen

          trleafsun(p,ic) = trleaf(p,ic,isun) ! Leaf transpiration flux (mol H2O/m2 leaf/s)
          evleafsun(p,ic) = evleaf(p,ic,isun) ! Leaf evaporation flux (mol H2O/m2 leaf/s)
          trleafsha(p,ic) = trleaf(p,ic,isha)
          evleafsha(p,ic) = evleaf(p,ic,isha)

          tr_prof(p,ic) = trleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) &  ! Canopy layer water vapor flux (mol H2O/m2/s)
                        + trleaf(p,ic,isha) * fracsha(p,ic) * dpai(p,ic)
          ev_prof(p,ic) = evleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) &  ! Canopy layer water vapor flux (mol H2O/m2/s)
                        + evleaf(p,ic,isha) * fracsha(p,ic) * dpai(p,ic)
       else
          ir_prof(p,ic) = 0._r8
          sw_prof(p,ic,ivis) = 0._r8
          sw_prof(p,ic,inir) = 0._r8
          rn_prof(p,ic) = 0._r8
          st_prof(p,ic) = 0._r8
          sh_prof(p,ic) = 0._r8
          lh_prof(p,ic) = 0._r8
          et_prof(p,ic) = 0._r8
          fc_prof(p,ic) = 0._r8

          trleafsun(p,ic) = 0._r8
          evleafsun(p,ic) = 0._r8
          trleafsha(p,ic) = 0._r8
          evleafsha(p,ic) = 0._r8

          tr_prof(p,ic) = 0._r8
          ev_prof(p,ic) = 0._r8
       end if
    end do

    ! Sum leaf fluxes

    irveg(p) = 0._r8
    stflx(p) = 0._r8
    shveg(p) = 0._r8
    lhveg(p) = 0._r8
    etveg(p) = 0._r8
    gppveg(p) = 0._r8
    trveg(p) = 0._r8
    evveg(p) = 0._r8

    do ic = 1, ncan(p)
       irveg(p) = irveg(p) + ir_prof(p,ic)
       stflx(p) = stflx(p) + st_prof(p,ic)
       shveg(p) = shveg(p) + sh_prof(p,ic)
       lhveg(p) = lhveg(p) + lh_prof(p,ic)
       etveg(p) = etveg(p) + et_prof(p,ic)
       trveg(p) = trveg(p) + tr_prof(p,ic)
       evveg(p) = evveg(p) + ev_prof(p,ic)
       if (dpai(p,ic) > 0._r8) then
          fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
          gppveg(p) = gppveg(p) + (ag(p,ic,isun)*fracsun(p,ic) + ag(p,ic,isha)*fracsha(p,ic)) * dpai(p,ic) * fracgreen
       end if
    end do
    stleafcan = stflx(p)
    ! Check energy balance for conservation

    err = swveg(p,ivis) + swveg(p,inir) + irveg(p) - shveg(p) - lhveg(p) - stflx(p)
    if ((abs(err) > 0.03_r8) .and. (nstep>1))  then
       call endrun (msg=' ERROR: CanopyFluxesSum: energy conservation error (1)')
    end if

    ! Soil fluxes

    ic = 0
    sw_prof(p,ic,ivis) = swsoi(p,ivis)
    sw_prof(p,ic,inir) = swsoi(p,inir)
    rn_prof(p,ic) = rnsoi(p)
    st_prof(p,ic) = 0._r8
    sh_prof(p,ic) = shsoi(p)
    lh_prof(p,ic) = lhsoi(p)
    et_prof(p,ic) = etsoi(p)
    fc_prof(p,ic) = 0._r8
          tr_prof(p,ic) = 0._r8
          ev_prof(p,ic) = 0._r8
    ! Turbulent fluxes

    select case (turb_type)
       case (0)
          ! Sum of layer fluxes
          shflx(p) = shveg(p) + shsoi(p)
          lhflx(p) = lhveg(p) + lhsoi(p)
          etflx(p) = etveg(p) + etsoi(p)
       case (1:4)
          ! Turbulent fluxes are at the top of the canopy
          shflx(p) = shair(p,ntop(p))
          etflx(p) = etair(p,ntop(p))
          lhflx(p) = etair(p,ntop(p)) * LatVap(tref(p))

          write(iulog,*) 'Sum of layer fluxes => :: shflx =', shveg(p) + shsoi(p),'lhflx =',lhveg(p) + lhsoi(p),'etflx=',etveg(p) + etsoi(p)
       case default
          call endrun (msg=' ERROR: CanopyFluxesMultilayer: turbulence type not valid')
    end select


    ! Add canopy air heat storage to storage flux
    !stleafcan = stflx(p)
    staircan = 0._r8
    do ic = 1, ntop(p)
       staircan = staircan + stair(p,ic)
       !stflx(p) = stflx(p) + stair(p,ic)
    end do
    stflx(p) = staircan + stflx(p)
    ! Overall energy balance check:
    ! radiation in - radiation out - soil heat = available energy = turbulent flux + canopy storage flux

    rnet(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir) + irveg(p) + irsoi(p)
    radin = swskyb(p,ivis) + swskyd(p,ivis) + swskyb(p,inir) + swskyd(p,inir) + irsky(p)
    radout = albcan(p,ivis)*(swskyb(p,ivis)+swskyd(p,ivis)) + albcan(p,inir)*(swskyb(p,inir)+swskyd(p,inir)) + ircan(p)

    err = rnet(p) - (radin - radout)
    if (abs(err) > 0.01_r8) then
       call endrun (msg=' ERROR: CanopyFluxesSum: energy conservation error (2)')
    end if

    avail = radin - radout - gsoi(p)
    flux = shflx(p) + lhflx(p) + stflx(p)
    err = avail - flux
    if ((abs(err) > 0.1_r8) .and. (nstep>2))  then
        write(iulog,*) '**-----------------------energy conservation error (3) err=', err
        write(iulog,*) 'rnet(p)=',rnet(p),'staircan=',staircan,'stleafcan=',stleafcan
        write(iulog,*) 'radin =', radin,'radout =', radout,'gsoi(p) =', gsoi(p)
        write(iulog,*) 'shflx =', shflx(p),'lhflx =', lhflx(p),'stflx(p) =', stflx(p)
        write(iulog,*) 'radin - radout - gsoi(p) = ', avail
        write(iulog,*) 'shflx(p) + lhflx(p) + stflx(p) = ', flux
       call endrun (msg=' ERROR: CanopyFluxesSum: energy conservation error (3)')
    end if



    ! Sunlit and shaded canopy fluxes

    irvegsun(p) = 0._r8; irvegsha(p) = 0._r8
    shvegsun(p) = 0._r8; shvegsha(p) = 0._r8
    lhvegsun(p) = 0._r8; lhvegsha(p) = 0._r8
    etvegsun(p) = 0._r8; etvegsha(p) = 0._r8
    gppvegsun(p) = 0._r8; gppvegsha(p) = 0._r8

    do ic = 1, ncan(p)
       if (dpai(p,ic) > 0._r8) then
          irvegsun(p) = irvegsun(p) + irleaf(p,ic) * fracsun(p,ic) * dpai(p,ic)
          irvegsha(p) = irvegsha(p) + irleaf(p,ic) * fracsha(p,ic) * dpai(p,ic)
          shvegsun(p) = shvegsun(p) + shleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
          shvegsha(p) = shvegsha(p) + shleaf(p,ic,isha) * fracsha(p,ic) * dpai(p,ic)
          lhvegsun(p) = lhvegsun(p) + lhleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
          lhvegsha(p) = lhvegsha(p) + lhleaf(p,ic,isha) * fracsha(p,ic) * dpai(p,ic)
          etvegsun(p) = etvegsun(p) + (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic)
          etvegsha(p) = etvegsha(p) + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * fracsha(p,ic) * dpai(p,ic)
          fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
          gppvegsun(p) = gppvegsun(p) + ag(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) * fracgreen
          gppvegsha(p) = gppvegsha(p) + ag(p,ic,isha) * fracsha(p,ic) * dpai(p,ic) * fracgreen
       end if
    end do
write(iulog,*) 'rnet    =',   rnet(p),'stflx    =',   stflx(p),'shflx    =',   shflx(p)
write(iulog,*) 'lhflx=',   lhflx(p),'etveg=', etveg(p),'trveg=',trveg(p),'evveg', evveg(p)
!write(iulog,*) 'stflx    =',   stflx     (p)    ! Canopy storage heat flux (W/m2)
!write(iulog,*) 'shflx    =',   shflx     (p)    ! Sensible heat flux (W/m2)
!write(iulog,*) 'lhflx    =',   lhflx     (p)    ! Latent heat flux (W/m2)
!write(iulog,*) 'etflx    =',   etflx     (p)    ! Water vapor flux (mol H2O/m2/s)
!write(iulog,*) 'irveg    =',   irveg     (p)    ! Absorbed longwave radiation, vegetation (W/m2)
!!write(iulog,*) 'irvegsun =',   irvegsun  (p)    ! Absorbed longwave radiation, sunlit canopy (W/m2)
!!write(iulog,*) 'irvegsha =',   irvegsha  (p)    ! Absorbed longwave radiation, shaded canopy (W/m2)
!write(iulog,*) 'shveg    =',   shveg     (p)    ! Sensible heat flux, vegetation (W/m2)
!!write(iulog,*) 'shvegsun =',   shvegsun  (p)    ! Sensible heat flux, sunlit canopy (W/m2)
!!write(iulog,*) 'shvegsha =',   shvegsha  (p)    ! Sensible heat flux, shaded canopy (W/m2)
!write(iulog,*) 'lhveg    =',   lhveg     (p)    ! Latent heat flux, vegetation (W/m2)
!!write(iulog,*) 'lhvegsun =',   lhvegsun  (p)    ! Latent heat flux, sunlit canopy (W/m2)
!!write(iulog,*) 'lhvegsha =',   lhvegsha  (p)    ! Latent heat flux, shaded canopy (W/m2)
!write(iulog,*) 'etveg    =',   etveg     (p)    ! Water vapor flux, vegetation (mol H2O/m2/s)
!!write(iulog,*) 'etvegsun =',   etvegsun  (p)    ! Water vapor flux, sunlit canopy (mol H2O/m2/s)
!!write(iulog,*) 'etvegsha =',   etvegsha  (p)    ! Water vapor flux, shaded canopy (mol H2O/m2/s)
!write(iulog,*) 'gppveg   =',   gppveg    (p)    ! Gross primary production (umol CO2/m2/s)
!!write(iulog,*) 'gppvegsun=',   gppvegsun (p)    ! Gross primary production, sunlit canopy (umol CO2/m2/s)
!!write(iulog,*) 'gppvegsha=',   gppvegsha (p)    ! Gross primary production, shaded canopy (umol CO2/m2/s)

    end associate
  end subroutine CanopyFluxesSum

end module CanopyFluxesMultilayerMod
