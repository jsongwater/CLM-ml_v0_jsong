module CanopyFluxesMultilayerType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Multilayer canopy module data structure
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : nlevgrnd, numrad, nlevcanml, nleaf
  use clm_varcon, only : ispval, spval, nan => spval
  use abortutils, only : endrun
  use decompMod , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA TYPES:

  type, public :: mlcanopy_type

    ! Vegetation input variables

    real(r8), pointer :: ztop(:)            ! Canopy height (m)
    real(r8), pointer :: ztop39(:)            ! Canopy height (m)
    real(r8), pointer :: lai(:)             ! Leaf area index of canopy (m2/m2)
    real(r8), pointer :: sai(:)             ! Stem area index of canopy (m2/m2)
    real(r8), pointer :: root_biomass(:)    ! Fine root biomass (g biomass / m2)
    integer , pointer :: ncan(:)            ! Number of aboveground layers
    integer , pointer :: nbot(:)            ! Index for bottom leaf layer
    integer , pointer :: ntop(:)            ! Index for top leaf layer
    integer , pointer :: nuforce(:)
    real(r8), pointer :: zuforce(:)
    real(r8), pointer :: dlai(:,:)          ! Layer leaf area index (m2/m2)
    real(r8), pointer :: dsai(:,:)          ! Layer stem area index (m2/m2)
    real(r8), pointer :: dpai(:,:)          ! Layer plant area index (m2/m2)
    real(r8), pointer :: sumpai(:,:)        ! Cumulative plant area index (m2/m2) [for nlevcanml layers]
    real(r8), pointer :: zs(:,:)            ! Canopy height for scalar concentration and source (m)
    real(r8), pointer :: zw(:,:)            ! Canopy heights at layer interfaces (m)

    ! Atmospheric input variables

    real(r8), pointer :: zref(:)            ! Reference height (m)
    real(r8), pointer :: zref_old(:)        ! Reference height for previous timestep (m)
    real(r8), pointer :: tref(:)            ! Air temperature at reference height (K)
    real(r8), pointer :: uref(:)            ! Wind speed at reference height (m/s)
    real(r8), pointer :: uref33(:)            ! Wind speed at measurement height (m/s)
    real(r8), pointer :: rhref(:)           ! Relative humidity at reference height (%)
    real(r8), pointer :: pref(:)            ! Air pressure at reference height (Pa)
    real(r8), pointer :: co2ref(:)          ! Atmospheric CO2 at reference height (umol/mol)
    real(r8), pointer :: o2ref(:)           ! Atmospheric O2 at reference height (mmol/mol)
    real(r8), pointer :: solar_zen(:)       ! Solar zenith angle (radians)
    real(r8), pointer :: swskyb(:,:)        ! Atmospheric direct beam solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: swskyd(:,:)        ! Atmospheric diffuse solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: irsky(:)           ! Atmospheric longwave radiation (W/m2)
    real(r8), pointer :: qflx_rain(:)       ! Rainfall (mm H2O/s = kg H2O/m2/s)
    real(r8), pointer :: qflx_snow(:)       ! Snowfall (mm H2O/s = kg H2O/m2/s)
    real(r8), pointer :: tacclim(:)         ! Average air temperature for acclimation (K)

    ! Additional derived input variables

    real(r8), pointer :: eref(:)            ! Vapor pressure at reference height (Pa)
    real(r8), pointer :: qref(:)            ! Specific humidity at reference height (kg/kg)
    real(r8), pointer :: rhoair(:)          ! Air density at reference height (kg/m3)
    real(r8), pointer :: rhomol(:)          ! Molar density at reference height (mol/m3)
    real(r8), pointer :: mmair(:)           ! Molecular mass of air at reference height (kg/mol)
    real(r8), pointer :: cpair(:)           ! Specific heat of air at constant pressure, at reference height (J/mol/K)

    ! Canopy layer variables [for nlevcanml layers]

    real(r8), pointer :: vcmax25(:,:)       ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    real(r8), pointer :: jmax25(:,:)        ! C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    real(r8), pointer :: kp25(:,:)          ! C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)
    real(r8), pointer :: rd25(:,:)          ! Leaf respiration rate at 25C for canopy layer (umol/m2/s)
    real(r8), pointer :: wind(:,:)          ! Wind speed profile (m/s)
    real(r8), pointer :: wind_1stoldH(:,:)
    real(r8), pointer :: wind_1st(:,:)          ! Wind speed profile (m/s)
    real(r8), pointer :: wind2LB(:,:)          ! Wind speed profile (m/s)
    real(r8), pointer :: wind_1stold(:,:)          ! Wind speed profile (m/s)
    real(r8), pointer :: wind_most(:,:)     ! Wind speed profile from MOST (m/s)
    real(r8), pointer :: tair(:,:)          ! Air temperature profile (K)
    real(r8), pointer :: tair_most(:,:)     ! Air temperature profile from MOST (K)
    real(r8), pointer :: eair(:,:)          ! Vapor pressure profile (Pa)
    real(r8), pointer :: cair(:,:)          ! Atmospheric CO2 profile (umol/mol)
    real(r8), pointer :: tveg(:,:,:)        ! Vegetation temperature profile (K)
    real(r8), pointer :: tvegsun(:,:)        ! Vegetation temperature profile (K)
    real(r8), pointer :: tvegsha(:,:)        ! Vegetation temperature profile (K)
    real(r8), pointer :: tair_old(:,:)      ! Air temperature profile for previous timestep (K)
    real(r8), pointer :: eair_old(:,:)      ! Vapor pressure profile for previous timestep (Pa)
    real(r8), pointer :: cair_old(:,:)      ! Atmospheric CO2 profile for previous timestep (umol/mol)
    real(r8), pointer :: tveg_old(:,:,:)    ! Vegetation temperature profile for previous timestep (K)
    real(r8), pointer :: fracsun(:,:)       ! Sunlit fraction of canopy layer
    real(r8), pointer :: fracsha(:,:)       ! Shaded fraction of canopy layer
    real(r8), pointer :: irleaf(:,:)        ! Leaf absorbed longwave radiation for canopy layer(W/m2 leaf)
    real(r8), pointer :: lwp(:,:)           ! Leaf water potential of canopy layer (MPa)
    real(r8), pointer :: lsc(:,:)           ! Leaf-specific conductance of canopy layer (mmol H2O/m2 leaf/s/MPa)
    real(r8), pointer :: h2ocan(:,:)        ! Canopy layer intercepted water (kg H2O/m2)
    real(r8), pointer :: fwet(:,:)          ! Fraction of plant area index that is wet
    real(r8), pointer :: fdry(:,:)          ! Fraction of plant area index that is green and dry
    real(r8), pointer :: shair(:,:)         ! Canopy air sensible heat flux (W/m2)
    real(r8), pointer :: etair(:,:)         ! Canopy air water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: cfair(:,:)
    real(r8), pointer :: stair(:,:)         ! Canopy air storage heat flux (W/m2)
    real(r8), pointer :: sw_prof(:,:,:)     ! Canopy layer absorbed solar radiation (W/m2)
    real(r8), pointer :: ir_prof(:,:)       ! Canopy layer absorbed longwave radiation (W/m2)
    real(r8), pointer :: rn_prof(:,:)       ! Canopy layer net radiation (W/m2)
    real(r8), pointer :: st_prof(:,:)       ! Canopy layer storage heat flux (W/m2)
    real(r8), pointer :: sh_prof(:,:)       ! Canopy layer sensible heat flux (W/m2)
    real(r8), pointer :: lh_prof(:,:)       ! Canopy layer latent heat flux (W/m2)
    real(r8), pointer :: et_prof(:,:)       ! Canopy layer water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: tr_prof(:,:)       ! Canopy layer water vapor flux (mol H2O/m2/s
    real(r8), pointer :: ev_prof(:,:)       ! Canopy layer water vapor flux (mol H2O/m2/s
    real(r8), pointer :: fc_prof(:,:)       ! Canopy layer CO2 flux (umol CO2/m2/s)
    real(r8), pointer :: ga_prof(:,:)       ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)

    ! Leaf variables [for nlevcanml layers and nleaf leaves (sunlit or shaded)]

    real(r8), pointer :: tleaf(:,:,:)       ! Leaf temperature (K)
    real(r8), pointer :: tleaf_old(:,:,:)   ! Leaf temperature for previous timestep (K)
    real(r8), pointer :: rnleaf(:,:,:)      ! Leaf net radiation (W/m2 leaf)
    real(r8), pointer :: stleaf(:,:,:)      ! Leaf storage heat flux (W/m2 leaf)
    real(r8), pointer :: shleaf(:,:,:)      ! Leaf sensible heat flux (W/m2 leaf)
    real(r8), pointer :: lhleaf(:,:,:)      ! Leaf latent heat flux (W/m2 leaf)
    real(r8), pointer :: swleaf(:,:,:,:)    ! Leaf absorbed solar radiation (W/m2 leaf) [for numrad wavebands]
    real(r8), pointer :: trleaf(:,:,:)      ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: trleafsun(:,:)      ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: trleafsha(:,:)      ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: evleaf(:,:,:)      ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: evleafsun(:,:)      ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: evleafsha(:,:)      ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: psil(:,:,:)        ! Leaf water potential (MPa)

    real(r8), pointer :: gbh(:,:,:)         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    real(r8), pointer :: gbv(:,:,:)         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    real(r8), pointer :: gbc(:,:,:)         ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)

    real(r8), pointer :: apar(:,:,:)        ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    real(r8), pointer :: aparsun(:,:)        ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    real(r8), pointer :: aparsha(:,:)        ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    real(r8), pointer :: ac(:,:,:)          ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: aj(:,:,:)          ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: ap(:,:,:)          ! Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: ag(:,:,:)          ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: an(:,:,:)          ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: rd(:,:,:)          ! Leaf respiration rate (umol CO2/m2 leaf/s)
    real(r8), pointer :: ci(:,:,:)          ! Leaf intercellular CO2 (umol/mol)
    real(r8), pointer :: cs(:,:,:)          ! Leaf surface CO2 (umol/mol)
    real(r8), pointer :: hs(:,:,:)          ! Leaf fractional humidity at leaf surface (-)
    real(r8), pointer :: vpd(:,:,:)         ! Leaf vapor pressure deficit (Pa)
    real(r8), pointer :: gs(:,:,:)          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    real(r8), pointer :: alphapsn(:,:,:)    ! Leaf 13C fractionation factor for photosynthesis (-)

    real(r8), pointer :: cpleaf(:,:)        ! Leaf heat capacity (J/m2 leaf/K)

    ! Canopy variables (fluxes are per m2 ground area)

    real(r8), pointer :: swveg(:,:)         ! Absorbed solar radiation, vegetation (W/m2) [for numrad wavebands]
    real(r8), pointer :: swvegsun(:,:)      ! Absorbed solar radiation, sunlit canopy (W/m2) [for numrad wavebands]
    real(r8), pointer :: swvegsha(:,:)      ! Absorbed solar radiation, shaded canopy (W/m2) [for numrad wavebands]

    real(r8), pointer :: irveg(:)           ! Absorbed longwave radiation, vegetation (W/m2)
    real(r8), pointer :: irvegsun(:)        ! Absorbed longwave radiation, sunlit canopy (W/m2)
    real(r8), pointer :: irvegsha(:)        ! Absorbed longwave radiation, shaded canopy (W/m2)

    real(r8), pointer :: shveg(:)           ! Sensible heat flux, vegetation (W/m2)
    real(r8), pointer :: shvegsun(:)        ! Sensible heat flux, sunlit canopy (W/m2)
    real(r8), pointer :: shvegsha(:)        ! Sensible heat flux, shaded canopy (W/m2)

    real(r8), pointer :: lhveg(:)           ! Latent heat flux, vegetation (W/m2)
    real(r8), pointer :: lhvegsun(:)        ! Latent heat flux, sunlit canopy (W/m2)
    real(r8), pointer :: lhvegsha(:)        ! Latent heat flux, shaded canopy (W/m2)

    real(r8), pointer :: etveg(:)           ! Water vapor flux, vegetation (mol H2O/m2/s)
    real(r8), pointer :: trveg(:)           ! Water vapor flux, vegetation (mol H2O/m2/s)
    real(r8), pointer :: evveg(:)           ! Water vapor flux, vegetation (mol H2O/m2/s)
    real(r8), pointer :: etvegsun(:)        ! Water vapor flux, sunlit canopy (mol H2O/m2/s)
    real(r8), pointer :: etvegsha(:)        ! Water vapor flux, shaded canopy (mol H2O/m2/s)

    real(r8), pointer :: gppveg(:)          ! Gross primary production (umol CO2/m2/s)
    real(r8), pointer :: gppvegsun(:)       ! Gross primary production, sunlit canopy (umol CO2/m2/s)
    real(r8), pointer :: gppvegsha(:)       ! Gross primary production, shaded canopy (umol CO2/m2/s)

    real(r8), pointer :: albcan(:,:)        ! Albedo above canopy [for numrad wavebands]
    real(r8), pointer :: ircan(:)           ! Upward longwave radiation above canopy (W/m2)
    real(r8), pointer :: rnet(:)            ! Net radiation (W/m2)
    real(r8), pointer :: stflx(:)           ! Canopy storage heat flux (W/m2)
    real(r8), pointer :: shflx(:)           ! Sensible heat flux (W/m2)
    real(r8), pointer :: lhflx(:)           ! Latent heat flux (W/m2)
    real(r8), pointer :: etflx(:)           ! Water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: fracminlwp(:)      ! Fraction of canopy with lwp < minlwp

    real(r8), pointer :: ustar(:)           ! Friction velocity (m/s)
    real(r8), pointer :: uforc(:)           ! Wind speed at reference height including stability effect (m/s)
    real(r8), pointer :: uaf(:)             ! Wind speed at canopy top (m/s)
    real(r8), pointer :: taf(:)             ! Air temperature at canopy top (K)
    real(r8), pointer :: qaf(:)             ! Specific humidity at canopy top (kg/kg)
    real(r8), pointer :: eaf(:)             ! Vapor pressure at canopy top (Pa)
    real(r8), pointer :: obu(:)             ! Obukhov length (m)
    real(r8), pointer :: obu_gah(:)         ! Obukhov length used for gah (m)
    real(r8), pointer :: obuold(:)          ! Obukhov length from previous iteration
    integer,  pointer :: nmozsgn(:)         ! Number of times stability changes sign during iteration

    real(r8), pointer :: z0mg(:)            ! Roughness length of ground (m)
    real(r8), pointer :: thref(:)           ! Atmospheric potential temperature (K)
    real(r8), pointer :: thvref(:)          ! Atmospheric virtual potential temperature (K)
    real(r8), pointer :: gah(:)             ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    real(r8), pointer :: PrSc(:)            ! Prandtl (Schmidt) number at canopy top
    real(r8), pointer :: Lc(:)              ! Canopy density length scale (m)
    real(r8), pointer :: zdisp(:)           ! Displacement height (m)
    real(r8), pointer :: tstar(:)           ! Temperature scale (K)
    real(r8), pointer :: qstar(:)           ! Water vapor scale (kg/kg)

    real(r8), pointer :: td(:,:)            ! Exponential transmittance of diffuse radiation through a single leaf layer

    ! Soil energy balance

    real(r8), pointer :: rnsoi(:)           ! Net radiation, ground (W/m2)
    real(r8), pointer :: shsoi(:)           ! Sensible heat flux, ground (W/m2)
    real(r8), pointer :: lhsoi(:)           ! Latent heat flux, ground (W/m2)
    real(r8), pointer :: gsoi(:)            ! Soil heat flux (W/m2)
    real(r8), pointer :: swsoi(:,:)         ! Absorbed solar radiation, ground (W/m2) [for numrad wavebands]
    real(r8), pointer :: irsoi(:)           ! Absorbed longwave radiation, ground (W/m2)
    real(r8), pointer :: etsoi(:)           ! Water vapor flux, ground (mol H2O/m2/s)
    real(r8), pointer :: tg(:)              ! Soil surface temperature (K)

    ! Soil moisture variables

    real(r8), pointer :: btran(:)           ! Ball-Berry soil wetness factor (-)
    real(r8), pointer :: psis(:)            ! Weighted soil water potential (MPa)
    real(r8), pointer :: rsoil(:)           ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    real(r8), pointer :: soil_et_loss(:,:)  ! Fraction of total transpiration from each soil layer (-)
    real(r8), pointer :: eg(:)              ! Soil surface vapor pressure (Pa)
    real(r8), pointer :: rhg(:)             ! Relative humidity of airspace at soil surface (fraction)

    ! Water flux variables

    real(r8), pointer :: qflx_prec_intr(:)  ! Intercepted precipitation (kg H2O/m2/s)

  contains

    procedure, public  :: Init              ! CLM initialization of data type
    procedure, private :: InitAllocate      ! CLM initialization: allocate module data structure
    procedure, private :: InitHistory       ! CLM initialization: setup history file variables
    procedure, private :: InitCold          ! CLM initialization: cold-start initialization
    procedure, public  :: Restart           ! CLM restart file

  end type mlcanopy_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)
    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start
    !
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate module data structure
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp      ! Beginning patch index for CLM g/l/c/p hierarchy
    integer :: endp      ! Ending patch index for CLM g/l/c/p hierarchy
    !---------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate (this%ztop           (begp:endp))                            ; this%ztop           (:)       = nan
    allocate (this%ztop39           (begp:endp))                            ; this%ztop39           (:)       = nan
    allocate (this%lai            (begp:endp))                            ; this%lai            (:)       = nan
    allocate (this%sai            (begp:endp))                            ; this%sai            (:)       = nan
    allocate (this%root_biomass   (begp:endp))                            ; this%root_biomass   (:)       = nan
    allocate (this%ncan           (begp:endp))                            ; this%ncan           (:)       = ispval
    allocate (this%nbot           (begp:endp))                            ; this%nbot           (:)       = ispval
    allocate (this%ntop           (begp:endp))                            ; this%ntop           (:)       = ispval
    allocate (this%nuforce           (begp:endp))                            ; this%nuforce           (:)       = ispval
    allocate (this%zuforce           (begp:endp))                            ; this%zuforce           (:)       = ispval
    allocate (this%dlai           (begp:endp,1:nlevcanml))                  ; this%dlai           (:,:)     = nan
    allocate (this%dsai           (begp:endp,1:nlevcanml))                  ; this%dsai           (:,:)     = nan
    allocate (this%dpai           (begp:endp,1:nlevcanml))                  ; this%dpai           (:,:)     = nan
    allocate (this%sumpai         (begp:endp,1:nlevcanml))                  ; this%sumpai         (:,:)     = nan
    allocate (this%zs             (begp:endp,0:nlevcanml))                  ; this%zs             (:,:)     = nan
    allocate (this%zw             (begp:endp,0:nlevcanml))                  ; this%zw             (:,:)     = nan
    allocate (this%zref           (begp:endp))                            ; this%zref           (:)       = nan
    allocate (this%zref_old       (begp:endp))                            ; this%zref_old       (:)       = nan
    allocate (this%tref           (begp:endp))                            ; this%tref           (:)       = nan
    allocate (this%uref           (begp:endp))                            ; this%uref           (:)       = nan
    allocate (this%uref33           (begp:endp))                            ; this%uref33           (:)       = nan
    allocate (this%rhref          (begp:endp))                            ; this%rhref          (:)       = nan
    allocate (this%pref           (begp:endp))                            ; this%pref           (:)       = nan
    allocate (this%co2ref         (begp:endp))                            ; this%co2ref         (:)       = nan
    allocate (this%o2ref          (begp:endp))                            ; this%o2ref          (:)       = nan
    allocate (this%solar_zen      (begp:endp))                            ; this%solar_zen      (:)       = nan
    allocate (this%swskyb         (begp:endp,1:numrad))                   ; this%swskyb         (:,:)     = nan
    allocate (this%swskyd         (begp:endp,1:numrad))                   ; this%swskyd         (:,:)     = nan
    allocate (this%irsky          (begp:endp))                            ; this%irsky          (:)       = nan
    allocate (this%qflx_rain      (begp:endp))                            ; this%qflx_rain      (:)       = nan
    allocate (this%qflx_snow      (begp:endp))                            ; this%qflx_snow      (:)       = nan
    allocate (this%tacclim        (begp:endp))                            ; this%tacclim        (:)       = nan
    allocate (this%eref           (begp:endp))                            ; this%eref           (:)       = nan
    allocate (this%qref           (begp:endp))                            ; this%qref           (:)       = nan
    allocate (this%rhoair         (begp:endp))                            ; this%rhoair         (:)       = nan
    allocate (this%rhomol         (begp:endp))                            ; this%rhomol         (:)       = nan
    allocate (this%mmair          (begp:endp))                            ; this%mmair          (:)       = nan
    allocate (this%cpair          (begp:endp))                            ; this%cpair          (:)       = nan
    allocate (this%vcmax25        (begp:endp,1:nlevcanml))                  ; this%vcmax25        (:,:)     = nan
    allocate (this%jmax25         (begp:endp,1:nlevcanml))                  ; this%jmax25         (:,:)     = nan
    allocate (this%kp25           (begp:endp,1:nlevcanml))                  ; this%kp25           (:,:)     = nan
    allocate (this%rd25           (begp:endp,1:nlevcanml))                  ; this%rd25           (:,:)     = nan
    allocate (this%wind           (begp:endp,0:nlevcanml))                  ; this%wind           (:,:)     = nan
    allocate (this%wind_1stoldH           (begp:endp,0:nlevcanml))                  ; this%wind_1stoldH           (:,:)     = nan
    allocate (this%wind_1st           (begp:endp,0:nlevcanml))              ; this%wind_1st           (:,:)     = nan
    allocate (this%wind2LB           (begp:endp,0:nlevcanml))              ; this%wind2LB           (:,:)     = nan
    allocate (this%wind_1stold           (begp:endp,0:nlevcanml))              ; this%wind_1stold           (:,:)     = nan
    allocate (this%wind_most      (begp:endp,0:nlevcanml))                  ; this%wind_most      (:,:)     = nan
    allocate (this%tair           (begp:endp,0:nlevcanml))                  ; this%tair           (:,:)     = nan
    allocate (this%tair_most      (begp:endp,0:nlevcanml))                  ; this%tair_most      (:,:)     = nan
    allocate (this%eair           (begp:endp,0:nlevcanml))                  ; this%eair           (:,:)     = nan
    allocate (this%cair           (begp:endp,0:nlevcanml))                  ; this%cair           (:,:)     = nan
    allocate (this%tveg           (begp:endp,0:nlevcanml,1:nleaf))          ; this%tveg           (:,:,:)   = nan
    allocate (this%tvegsun           (begp:endp,0:nlevcanml))          ; this%tvegsun           (:,:)   = nan
    allocate (this%tvegsha           (begp:endp,0:nlevcanml))          ; this%tvegsha           (:,:)   = nan
    allocate (this%tair_old       (begp:endp,0:nlevcanml))                  ; this%tair_old       (:,:)     = nan
    allocate (this%eair_old       (begp:endp,0:nlevcanml))                  ; this%eair_old       (:,:)     = nan
    allocate (this%cair_old       (begp:endp,0:nlevcanml))                  ; this%cair_old       (:,:)     = nan
    allocate (this%tveg_old       (begp:endp,0:nlevcanml,1:nleaf))          ; this%tveg_old       (:,:,:)   = nan
    allocate (this%fracsun        (begp:endp,1:nlevcanml))                  ; this%fracsun        (:,:)     = nan
    allocate (this%fracsha        (begp:endp,1:nlevcanml))                  ; this%fracsha        (:,:)     = nan
    allocate (this%irleaf         (begp:endp,1:nlevcanml))                  ; this%irleaf         (:,:)     = nan
    allocate (this%lwp            (begp:endp,1:nlevcanml))                  ; this%lwp            (:,:)     = nan
    allocate (this%lsc            (begp:endp,1:nlevcanml))                  ; this%lsc            (:,:)     = nan
    allocate (this%h2ocan         (begp:endp,1:nlevcanml))                  ; this%h2ocan         (:,:)     = nan
    allocate (this%fwet           (begp:endp,1:nlevcanml))                  ; this%fwet           (:,:)     = nan
    allocate (this%fdry           (begp:endp,1:nlevcanml))                  ; this%fdry           (:,:)     = nan
    allocate (this%shair          (begp:endp,1:nlevcanml))                  ; this%shair          (:,:)     = nan
    allocate (this%etair          (begp:endp,1:nlevcanml))                  ; this%etair          (:,:)     = nan
    allocate (this%cfair          (begp:endp,1:nlevcanml))                  ; this%cfair          (:,:)     = nan
    allocate (this%stair          (begp:endp,1:nlevcanml))                  ; this%stair          (:,:)     = nan
    allocate (this%sw_prof        (begp:endp,0:nlevcanml,1:numrad))         ; this%sw_prof        (:,:,:)   = nan
    allocate (this%ir_prof        (begp:endp,0:nlevcanml))                  ; this%ir_prof        (:,:)     = nan
    allocate (this%rn_prof        (begp:endp,0:nlevcanml))                  ; this%rn_prof        (:,:)     = nan
    allocate (this%st_prof        (begp:endp,0:nlevcanml))                  ; this%st_prof        (:,:)     = nan
    allocate (this%sh_prof        (begp:endp,0:nlevcanml))                  ; this%sh_prof        (:,:)     = nan
    allocate (this%lh_prof        (begp:endp,0:nlevcanml))                  ; this%lh_prof        (:,:)     = nan
    allocate (this%et_prof        (begp:endp,0:nlevcanml))                  ; this%et_prof        (:,:)     = nan
    allocate (this%tr_prof        (begp:endp,0:nlevcanml))                  ; this%tr_prof        (:,:)     = nan
    allocate (this%ev_prof        (begp:endp,0:nlevcanml))                  ; this%ev_prof        (:,:)     = nan
    allocate (this%fc_prof        (begp:endp,0:nlevcanml))                  ; this%fc_prof        (:,:)     = nan
    allocate (this%ga_prof        (begp:endp,0:nlevcanml))                  ; this%ga_prof        (:,:)     = nan
    allocate (this%tleaf          (begp:endp,1:nlevcanml,1:nleaf))          ; this%tleaf          (:,:,:)   = nan
    allocate (this%tleaf_old      (begp:endp,1:nlevcanml,1:nleaf))          ; this%tleaf_old      (:,:,:)   = nan
    allocate (this%rnleaf         (begp:endp,1:nlevcanml,1:nleaf))          ; this%rnleaf         (:,:,:)   = nan
    allocate (this%stleaf         (begp:endp,1:nlevcanml,1:nleaf))          ; this%stleaf         (:,:,:)   = nan
    allocate (this%shleaf         (begp:endp,1:nlevcanml,1:nleaf))          ; this%shleaf         (:,:,:)   = nan
    allocate (this%lhleaf         (begp:endp,1:nlevcanml,1:nleaf))          ; this%lhleaf         (:,:,:)   = nan
    allocate (this%swleaf         (begp:endp,1:nlevcanml,1:nleaf,1:numrad)) ; this%swleaf         (:,:,:,:) = nan
    allocate (this%trleaf         (begp:endp,1:nlevcanml,1:nleaf))          ; this%trleaf         (:,:,:)   = nan
    allocate (this%evleaf         (begp:endp,1:nlevcanml,1:nleaf))          ; this%evleaf         (:,:,:)   = nan
    allocate (this%trleafsun      (begp:endp,1:nlevcanml))                  ; this%trleafsun      (:,:)   = nan
    allocate (this%evleafsun      (begp:endp,1:nlevcanml))                  ; this%evleafsun      (:,:)   = nan
    allocate (this%trleafsha      (begp:endp,1:nlevcanml))                  ; this%trleafsha      (:,:)   = nan
    allocate (this%evleafsha      (begp:endp,1:nlevcanml))                  ; this%evleafsha      (:,:)   = nan
    allocate (this%psil           (begp:endp,1:nlevcanml,1:nleaf))          ; this%psil           (:,:,:)   = nan
    allocate (this%gbh            (begp:endp,1:nlevcanml,1:nleaf))          ; this%gbh            (:,:,:)   = nan
    allocate (this%gbv            (begp:endp,1:nlevcanml,1:nleaf))          ; this%gbv            (:,:,:)   = nan
    allocate (this%gbc            (begp:endp,1:nlevcanml,1:nleaf))          ; this%gbc            (:,:,:)   = nan
    allocate (this%apar           (begp:endp,1:nlevcanml,1:nleaf))          ; this%apar           (:,:,:)   = nan
    allocate (this%aparsun           (begp:endp,1:nlevcanml))          ; this%aparsun           (:,:)   = nan
    allocate (this%aparsha           (begp:endp,1:nlevcanml))          ; this%aparsha           (:,:)   = nan
    allocate (this%ac             (begp:endp,1:nlevcanml,1:nleaf))          ; this%ac             (:,:,:)   = nan
    allocate (this%aj             (begp:endp,1:nlevcanml,1:nleaf))          ; this%aj             (:,:,:)   = nan
    allocate (this%ap             (begp:endp,1:nlevcanml,1:nleaf))          ; this%ap             (:,:,:)   = nan
    allocate (this%ag             (begp:endp,1:nlevcanml,1:nleaf))          ; this%ag             (:,:,:)   = nan
    allocate (this%an             (begp:endp,1:nlevcanml,1:nleaf))          ; this%an             (:,:,:)   = nan
    allocate (this%rd             (begp:endp,1:nlevcanml,1:nleaf))          ; this%rd             (:,:,:)   = nan
    allocate (this%ci             (begp:endp,1:nlevcanml,1:nleaf))          ; this%ci             (:,:,:)   = nan
    allocate (this%cs             (begp:endp,1:nlevcanml,1:nleaf))          ; this%cs             (:,:,:)   = nan
    allocate (this%hs             (begp:endp,1:nlevcanml,1:nleaf))          ; this%hs             (:,:,:)   = nan
    allocate (this%vpd            (begp:endp,1:nlevcanml,1:nleaf))          ; this%vpd            (:,:,:)   = nan
    allocate (this%gs             (begp:endp,1:nlevcanml,1:nleaf))          ; this%gs             (:,:,:)   = nan
    allocate (this%alphapsn       (begp:endp,1:nlevcanml,1:nleaf))          ; this%alphapsn       (:,:,:)   = nan
    allocate (this%cpleaf         (begp:endp,1:nlevcanml))                  ; this%cpleaf         (:,:)     = nan
    allocate (this%swveg          (begp:endp,1:numrad))                   ; this%swveg          (:,:)     = nan
    allocate (this%swvegsun       (begp:endp,1:numrad))                   ; this%swvegsun       (:,:)     = nan
    allocate (this%swvegsha       (begp:endp,1:numrad))                   ; this%swvegsha       (:,:)     = nan
    allocate (this%irveg          (begp:endp))                            ; this%irveg          (:)       = nan
    allocate (this%irvegsun       (begp:endp))                            ; this%irvegsun       (:)       = nan
    allocate (this%irvegsha       (begp:endp))                            ; this%irvegsha       (:)       = nan
    allocate (this%shveg          (begp:endp))                            ; this%shveg          (:)       = nan
    allocate (this%shvegsun       (begp:endp))                            ; this%shvegsun       (:)       = nan
    allocate (this%shvegsha       (begp:endp))                            ; this%shvegsha       (:)       = nan
    allocate (this%lhveg          (begp:endp))                            ; this%lhveg          (:)       = nan
    allocate (this%lhvegsun       (begp:endp))                            ; this%lhvegsun       (:)       = nan
    allocate (this%lhvegsha       (begp:endp))                            ; this%lhvegsha       (:)       = nan
    allocate (this%etveg          (begp:endp))                            ; this%etveg          (:)       = nan
    allocate (this%trveg          (begp:endp))                            ; this%trveg          (:)       = nan
    allocate (this%evveg          (begp:endp))                            ; this%evveg          (:)       = nan
    allocate (this%etvegsun       (begp:endp))                            ; this%etvegsun       (:)       = nan
    allocate (this%etvegsha       (begp:endp))                            ; this%etvegsha       (:)       = nan
    allocate (this%gppveg         (begp:endp))                            ; this%gppveg         (:)       = nan
    allocate (this%gppvegsun      (begp:endp))                            ; this%gppvegsun      (:)       = nan
    allocate (this%gppvegsha      (begp:endp))                            ; this%gppvegsha      (:)       = nan
    allocate (this%albcan         (begp:endp,1:numrad))                   ; this%albcan         (:,:)     = nan
    allocate (this%ircan          (begp:endp))                            ; this%ircan          (:)       = nan
    allocate (this%rnet           (begp:endp))                            ; this%rnet           (:)       = nan
    allocate (this%stflx          (begp:endp))                            ; this%stflx          (:)       = nan
    allocate (this%shflx          (begp:endp))                            ; this%shflx          (:)       = nan
    allocate (this%lhflx          (begp:endp))                            ; this%lhflx          (:)       = nan
    allocate (this%etflx          (begp:endp))                            ; this%etflx          (:)       = nan
    allocate (this%fracminlwp     (begp:endp))                            ; this%fracminlwp     (:)       = nan
    allocate (this%ustar          (begp:endp))                            ; this%ustar          (:)       = nan
    allocate (this%uforc          (begp:endp))                            ; this%uforc          (:)       = nan
    allocate (this%uaf            (begp:endp))                            ; this%uaf            (:)       = nan
    allocate (this%taf            (begp:endp))                            ; this%taf            (:)       = nan
    allocate (this%qaf            (begp:endp))                            ; this%qaf            (:)       = nan
    allocate (this%eaf            (begp:endp))                            ; this%eaf            (:)       = nan
    allocate (this%obu            (begp:endp))                            ; this%obu            (:)       = nan
    allocate (this%obu_gah        (begp:endp))                            ; this%obu_gah        (:)       = nan
    allocate (this%obuold         (begp:endp))                            ; this%obuold         (:)       = nan
    allocate (this%nmozsgn        (begp:endp))                            ; this%nmozsgn        (:)       = ispval
    allocate (this%z0mg           (begp:endp))                            ; this%z0mg           (:)       = nan
    allocate (this%thref          (begp:endp))                            ; this%thref          (:)       = nan
    allocate (this%thvref         (begp:endp))                            ; this%thvref         (:)       = nan
    allocate (this%gah            (begp:endp))                            ; this%gah            (:)       = nan
    allocate (this%PrSc           (begp:endp))                            ; this%PrSc           (:)       = nan
    allocate (this%Lc             (begp:endp))                            ; this%Lc             (:)       = nan
    allocate (this%zdisp          (begp:endp))                            ; this%zdisp          (:)       = nan
    allocate (this%tstar          (begp:endp))                            ; this%tstar          (:)       = nan
    allocate (this%qstar          (begp:endp))                            ; this%qstar          (:)       = nan
    allocate (this%td             (begp:endp,1:nlevcanml))                  ; this%td             (:,:)     = nan
    allocate (this%rnsoi          (begp:endp))                            ; this%rnsoi          (:)       = nan
    allocate (this%shsoi          (begp:endp))                            ; this%shsoi          (:)       = nan
    allocate (this%lhsoi          (begp:endp))                            ; this%lhsoi          (:)       = nan
    allocate (this%gsoi           (begp:endp))                            ; this%gsoi           (:)       = nan
    allocate (this%swsoi          (begp:endp,1:numrad))                   ; this%swsoi          (:,:)     = nan
    allocate (this%irsoi          (begp:endp))                            ; this%irsoi          (:)       = nan
    allocate (this%etsoi          (begp:endp))                            ; this%etsoi          (:)       = nan
    allocate (this%tg             (begp:endp))                            ; this%tg             (:)       = nan
    allocate (this%btran          (begp:endp))                            ; this%btran          (:)       = nan
    allocate (this%psis           (begp:endp))                            ; this%psis           (:)       = nan
    allocate (this%rsoil          (begp:endp))                            ; this%rsoil          (:)       = nan
    allocate (this%soil_et_loss   (begp:endp,1:nlevgrnd))                 ; this%soil_et_loss   (:,:)     = nan
    allocate (this%eg             (begp:endp))                            ; this%eg             (:)       = nan
    allocate (this%rhg            (begp:endp))                            ; this%rhg            (:)       = nan
    allocate (this%qflx_prec_intr (begp:endp))                            ; this%qflx_prec_intr (:)       = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory (this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup the fields that can be output on history files
    !
    ! !USES:
    use histFileMod, only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%dlai(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_dlai', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml dlai', &
         ptr_patch=this%dlai, default='inactive')

    this%dsai(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_dsai', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml dsai', &
         ptr_patch=this%dsai, default='inactive')

    this%sumpai(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_sumpai', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml sumpai', &
         ptr_patch=this%sumpai, default='inactive')

    this%zs(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_zs', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml zs', &
         ptr_patch=this%zs, default='inactive')
!!!! ---'ml_aparsun','ml_aparsha'
    this%aparsun(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_aparsun', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml aparsun', &
         ptr_patch=this%aparsun, default='inactive')
    this%aparsha(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_aparsha', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml aparsha', &
         ptr_patch=this%aparsha, default='inactive')

!!!!!---'ml_tvegsun','ml_tvegsha'


    this%tvegsun(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_tvegsun', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml tvegsun', &
         ptr_patch=this%tvegsun, default='inactive')

    this%tvegsha(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_tvegsha', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml tvegsha', &
         ptr_patch=this%tvegsha, default='inactive')

!!!!!!!---'ml_tr_prof','ml_ev_prof','ml_trveg','ml_evveg'

    this%tr_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_tr_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml tr_prof', &
         ptr_patch=this%tr_prof, default='inactive')
    this%ev_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_ev_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml ev_prof', &
         ptr_patch=this%ev_prof, default='inactive')

    this%trveg(begp:endp) = nan
    call hist_addfld1d (fname='ml_trveg', units='', &
         avgflag='A', long_name='ml trveg', &
         ptr_patch=this%trveg, set_lake=0._r8, set_urb=0._r8, default='inactive')
    this%evveg(begp:endp) = nan
    call hist_addfld1d (fname='ml_evveg', units='', &
         avgflag='A', long_name='ml evveg', &
         ptr_patch=this%evveg, set_lake=0._r8, set_urb=0._r8, default='inactive')
!!!!!!!!----          'ml_trleafsun','ml_evleafsun','ml_trleafsha','ml_evleafsha'


    this%trleafsun(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_trleafsun', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml trleafsun', &
         ptr_patch=this%trleafsun, default='inactive')

    this%evleafsun(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_evleafsun', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml evleafsun', &
         ptr_patch=this%evleafsun, default='inactive')

     this%trleafsha(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_trleafsha', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml trleafsha', &
         ptr_patch=this%trleafsha, default='inactive')

    this%evleafsha(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_evleafsha', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml evleafsha', &
         ptr_patch=this%evleafsha, default='inactive')


!----------------------------------------------------
!'solar_zen','ml_swskyb','ml_swskyd','ml_irsky','ml_qflx_rain','ml_qflx_snow','ml_vcmax25','ml_jmax25','ml_kp25','ml_rd25','ml_wind','ml_wind_most','ml_tair','ml_tair_most','ml_eair','ml_cair','ml_fracsun','ml_irleaf','ml_lwp','ml_lsc','ml_h2ocan','ml_fwet','ml_shair','ml_etair','ml_stair','ml_ir_prof','ml_rn_prof','ml_st_prof','ml_lh_prof','ml_et_prof','ml_fc_prof','ml_ga_prof','ml_td', 'ml_cpleaf','ml_swsoi','ml_swveg','ml_swvegsun','ml_swvegsha','ml_albcan','ml_soil_et_loss'
    this%solar_zen(begp:endp) = nan
    call hist_addfld1d (fname='solar_zen', units='', &
         avgflag='A', long_name='solar_zen', &
         ptr_patch=this%solar_zen, set_lake=0._r8, set_urb=0._r8, default='inactive')

    this%swskyb(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_swskyb', units='', type2d='numrad', &
         avgflag='A', long_name='ml swskyb', &
         ptr_patch=this%swskyb, default='inactive')

    this%swskyd(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_swskyd', units='', type2d='numrad', &
         avgflag='A', long_name='ml swskyd', &
         ptr_patch=this%swskyd, default='inactive')

    this%irsky(begp:endp) = nan
    call hist_addfld1d (fname='ml_irsky', units='', &
         avgflag='A', long_name='ml irsky', &
         ptr_patch=this%irsky, set_lake=0._r8, set_urb=0._r8, default='inactive')

    this%qflx_rain(begp:endp) = nan
    call hist_addfld1d (fname='ml_qflx_rain', units='', &
         avgflag='A', long_name='ml qflx_rain', &
         ptr_patch=this%qflx_rain, set_lake=0._r8, set_urb=0._r8, default='inactive')

    this%qflx_snow(begp:endp) = nan
    call hist_addfld1d (fname='ml_qflx_snow', units='', &
         avgflag='A', long_name='ml qflx_snow', &
         ptr_patch=this%qflx_snow, set_lake=0._r8, set_urb=0._r8, default='inactive')



    this%vcmax25(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_vcmax25', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml vcmax25', &
         ptr_patch=this%vcmax25, default='inactive')
    this%jmax25(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_jmax25', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml jmax25', &
         ptr_patch=this%jmax25, default='inactive')
    this%kp25(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_kp25', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml kp25', &
         ptr_patch=this%kp25, default='inactive')
    this%rd25(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_rd25', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml rd25', &
         ptr_patch=this%rd25, default='inactive')


    this%wind(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_wind', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml wind', &
         ptr_patch=this%wind, default='inactive')

    this%wind_1st(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_wind_1st', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml wind_1st', &
         ptr_patch=this%wind_1st, default='inactive')

    this%wind_most(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_wind_most', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml wind_most', &
         ptr_patch=this%wind_most, default='inactive')

    this%tair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_tair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml tair', &
         ptr_patch=this%tair, default='inactive')
    this%tair_most(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_tair_most', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml tair_most', &
         ptr_patch=this%tair_most, default='inactive')

    this%eair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_eair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml eair', &
         ptr_patch=this%eair, default='inactive')

    this%cair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_cair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml cair', &
         ptr_patch=this%cair, default='inactive')
    this%fracsun(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_fracsun', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml fracsun', &
         ptr_patch=this%fracsun, default='inactive')

    this%irleaf(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_irleaf', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml irleaf', &
         ptr_patch=this%irleaf, default='inactive')
    this%lwp(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_lwp', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml lwp', &
         ptr_patch=this%lwp, default='inactive')

    this%lsc(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_lsc', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml lsc', &
         ptr_patch=this%lsc, default='inactive')


    this%h2ocan (begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_h2ocan', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml h2ocan ', &
         ptr_patch=this%h2ocan , default='inactive')

    this%fwet(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_fwet', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml fwet', &
         ptr_patch=this%fwet, default='inactive')


    this%shair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_shair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml shair', &
         ptr_patch=this%shair, default='inactive')

    this%etair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_etair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml etair', &
         ptr_patch=this%etair, default='inactive')
    this%cfair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_cfair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml cfair', &
         ptr_patch=this%cfair, default='inactive')
    this%stair(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_stair', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml stair', &
         ptr_patch=this%stair, default='inactive')

    this%ir_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_ir_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml ir_prof', &
         ptr_patch=this%ir_prof, default='inactive')
    this%rn_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_rn_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml rn_prof', &
         ptr_patch=this%rn_prof, default='inactive')

    this%st_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_st_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml st_prof', &
         ptr_patch=this%st_prof, default='inactive')
    this%lh_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_lh_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml lh_prof', &
         ptr_patch=this%lh_prof, default='inactive')

    this%et_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_et_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml et_prof', &
         ptr_patch=this%et_prof, default='inactive')

    this%fc_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_fc_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml fc_prof', &
         ptr_patch=this%fc_prof, default='inactive')
    this%ga_prof(begp:endp,0:nlevcanml) = nan
    call hist_addfld2d (fname='ml_ga_prof', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml ga_prof', &
         ptr_patch=this%ga_prof, default='inactive')

    this%td(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_td', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml td', &
         ptr_patch=this%td, default='inactive')

    this%cpleaf(begp:endp,1:nlevcanml) = nan
    call hist_addfld2d (fname='ml_cpleaf', units='', type2d='nlevcanml', &
         avgflag='A', long_name='ml cpleaf', &
         ptr_patch=this%cpleaf, default='inactive')


    this%swsoi(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_swsoi', units='', type2d='numrad', &
         avgflag='A', long_name='ml swsoi', &
         ptr_patch=this%swsoi, default='inactive')

    this%swveg(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_swveg', units='', type2d='numrad', &
         avgflag='A', long_name='ml swveg', &
         ptr_patch=this%swveg, default='inactive')

    this%swvegsun(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_swvegsun', units='', type2d='numrad', &
         avgflag='A', long_name='ml swvegsun', &
         ptr_patch=this%swvegsun, default='inactive')

    this%swvegsha(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_swvegsha', units='', type2d='numrad', &
         avgflag='A', long_name='ml swvegsha', &
         ptr_patch=this%swvegsha, default='inactive')

    this%albcan(begp:endp,1:numrad) = nan
    call hist_addfld2d (fname='ml_albcan', units='', type2d='numrad', &
         avgflag='A', long_name='ml albcan', &
         ptr_patch=this%albcan, default='inactive')

    this%soil_et_loss(begp:endp,1:nlevgrnd) = nan
    call hist_addfld2d (fname='ml_soil_et_loss', units='', type2d='nlevgrnd', &
         avgflag='A', long_name='ml soil_et_loss', &
         ptr_patch=this%soil_et_loss, default='inactive')





  !  this%dvdv(begp:endp,1:nlevcanml) = nan
  !  call hist_addfld2d (fname='ml_dvdv', units='', type2d='nlevcanml', &
  !       avgflag='A', long_name='ml dvdv', &
  !       ptr_patch=this%dvdv, default='inactive')
  !  this%dvdv(begp:endp,1:nlevcanml) = nan
  !  call hist_addfld2d (fname='ml_dvdv', units='', type2d='nlevcanml', &
  !       avgflag='A', long_name='ml dvdv', &
  !       ptr_patch=this%dvdv, default='inactive')
  !
  !  this%dvdv(begp:endp,1:nlevcanml) = nan
  !  call hist_addfld2d (fname='ml_dvdv', units='', type2d='nlevcanml', &
  !       avgflag='A', long_name='ml dvdv', &
  !       ptr_patch=this%dvdv, default='inactive')
  !  this%dvdv(begp:endp,1:nlevcanml) = nan
  !  call hist_addfld2d (fname='ml_dvdv', units='', type2d='nlevcanml', &
  !       avgflag='A', long_name='ml dvdv', &
  !       ptr_patch=this%dvdv, default='inactive')
  !
  !  this%dvdv(begp:endp,1:nlevcanml) = nan
  !  call hist_addfld2d (fname='ml_dvdv', units='', type2d='nlevcanml', &
  !       avgflag='A', long_name='ml dvdv', &
  !       ptr_patch=this%dvdv, default='inactive')
  !

    this%irveg          (begp:endp)       = nan
    this%irvegsun       (begp:endp)       = nan
    this%irvegsha       (begp:endp)       = nan
    this%shveg          (begp:endp)       = nan
    this%shvegsun       (begp:endp)       = nan
    this%shvegsha       (begp:endp)       = nan
    this%lhveg          (begp:endp)       = nan
    this%lhvegsun       (begp:endp)       = nan
    this%lhvegsha       (begp:endp)       = nan
    this%etveg          (begp:endp)       = nan
    this%etvegsun       (begp:endp)       = nan
    this%etvegsha       (begp:endp)       = nan
    this%gppveg         (begp:endp)       = nan
    this%gppvegsun      (begp:endp)       = nan
    this%gppvegsha      (begp:endp)       = nan
    this%ircan          (begp:endp)       = nan
    this%rnet           (begp:endp)       = nan
    this%stflx          (begp:endp)       = nan
    this%shflx          (begp:endp)       = nan
    this%lhflx          (begp:endp)       = nan
    this%etflx          (begp:endp)       = nan
    this%fracminlwp     (begp:endp)       = nan
    this%ustar          (begp:endp)       = nan
    this%uforc          (begp:endp)       = nan
    this%uaf            (begp:endp)       = nan
    this%taf            (begp:endp)       = nan
    this%qaf            (begp:endp)       = nan
    this%eaf            (begp:endp)       = nan
    this%obu            (begp:endp)       = nan
    this%obu_gah        (begp:endp)       = nan
    this%obuold         (begp:endp)       = nan
    this%nmozsgn        (begp:endp)       = ispval
    this%z0mg           (begp:endp)       = nan
    this%thref          (begp:endp)       = nan
    this%thvref         (begp:endp)       = nan
    this%gah            (begp:endp)       = nan
    this%PrSc           (begp:endp)       = nan
    this%Lc             (begp:endp)       = nan
    this%zdisp          (begp:endp)       = nan
    this%tstar          (begp:endp)       = nan
    this%qstar          (begp:endp)       = nan
    this%rnsoi          (begp:endp)       = nan
    this%shsoi          (begp:endp)       = nan
    this%lhsoi          (begp:endp)       = nan
    this%gsoi           (begp:endp)       = nan
    this%irsoi          (begp:endp)       = nan
    this%etsoi          (begp:endp)       = nan
    this%tg             (begp:endp)       = nan
    this%btran          (begp:endp)       = nan
    this%psis           (begp:endp)       = nan
    this%rsoil          (begp:endp)       = nan
    this%eg             (begp:endp)       = nan
    this%rhg            (begp:endp)       = nan
    this%qflx_prec_intr (begp:endp)       = nan
! ----------------------------
!'ml_irveg','ml_irvegsun','ml_irvegsha','ml_shveg','ml_shvegsun','ml_shvegsha','ml_lhveg','ml_lhvegsun','ml_lhvegsha','ml_etveg','ml_etvegsun','ml_etvegsha','ml_gppveg','ml_gppvegsun','ml_gppvegsha','ml_ircan','ml_rnet','ml_stflx','ml_shflx','ml_lhflx','ml_etflx','ml_fracminlwp','ml_ustar','ml_uforc','ml_uaf','ml_taf','ml_qaf','ml_eaf','ml_obu','ml_obu_gah','ml_obuold','ml_nmozsgn','ml_z0mg','ml_thref','ml_thvref','ml_gah','ml_PrSc','ml_Lc','ml_zdisp','ml_tstar','ml_qstar','ml_rnsoi','ml_shsoi','ml_lhsoi','ml_gsoi','ml_irsoi','ml_etsoi','ml_tg','ml_btran','ml_psis','ml_rsoil','ml_eg','ml_rhg','ml_qflx_prec_intr'

    call hist_addfld1d (fname= 'ml_irveg'              , units='', avgflag='A', long_name=  'ml irveg         ', ptr_patch=this%irveg         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_irvegsun'           , units='', avgflag='A', long_name=  'ml irvegsun      ', ptr_patch=this%irvegsun      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_irvegsha'           , units='', avgflag='A', long_name=  'ml irvegsha      ', ptr_patch=this%irvegsha      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_shveg'              , units='', avgflag='A', long_name=  'ml shveg         ', ptr_patch=this%shveg         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_shvegsun'           , units='', avgflag='A', long_name=  'ml shvegsun      ', ptr_patch=this%shvegsun      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_shvegsha'           , units='', avgflag='A', long_name=  'ml shvegsha      ', ptr_patch=this%shvegsha      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_lhveg'              , units='', avgflag='A', long_name=  'ml lhveg         ', ptr_patch=this%lhveg         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_lhvegsun'           , units='', avgflag='A', long_name=  'ml lhvegsun      ', ptr_patch=this%lhvegsun      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_lhvegsha'           , units='', avgflag='A', long_name=  'ml lhvegsha      ', ptr_patch=this%lhvegsha      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_etveg'              , units='', avgflag='A', long_name=  'ml etveg         ', ptr_patch=this%etveg         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_etvegsun'           , units='', avgflag='A', long_name=  'ml etvegsun      ', ptr_patch=this%etvegsun      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_etvegsha'           , units='', avgflag='A', long_name=  'ml etvegsha      ', ptr_patch=this%etvegsha      , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_gppveg'             , units='', avgflag='A', long_name=  'ml gppveg        ', ptr_patch=this%gppveg        , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_gppvegsun'          , units='', avgflag='A', long_name=  'ml gppvegsun     ', ptr_patch=this%gppvegsun     , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_gppvegsha'          , units='', avgflag='A', long_name=  'ml gppvegsha     ', ptr_patch=this%gppvegsha     , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_ircan'              , units='', avgflag='A', long_name=  'ml ircan         ', ptr_patch=this%ircan         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_rnet'               , units='', avgflag='A', long_name=  'ml rnet          ', ptr_patch=this%rnet          , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_stflx'              , units='', avgflag='A', long_name=  'ml stflx         ', ptr_patch=this%stflx         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_shflx'              , units='', avgflag='A', long_name=  'ml shflx         ', ptr_patch=this%shflx         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_lhflx'              , units='', avgflag='A', long_name=  'ml lhflx         ', ptr_patch=this%lhflx         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_etflx'              , units='', avgflag='A', long_name=  'ml etflx         ', ptr_patch=this%etflx         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_fracminlwp'         , units='', avgflag='A', long_name=  'ml fracminlwp    ', ptr_patch=this%fracminlwp    , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_ustar'              , units='', avgflag='A', long_name=  'ml ustar         ', ptr_patch=this%ustar         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_uforc'              , units='', avgflag='A', long_name=  'ml uforc         ', ptr_patch=this%uforc         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_uaf'                , units='', avgflag='A', long_name=  'ml uaf           ', ptr_patch=this%uaf           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_taf'                , units='', avgflag='A', long_name=  'ml taf           ', ptr_patch=this%taf           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_qaf'                , units='', avgflag='A', long_name=  'ml qaf           ', ptr_patch=this%qaf           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_eaf'                , units='', avgflag='A', long_name=  'ml eaf           ', ptr_patch=this%eaf           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_obu'                , units='', avgflag='A', long_name=  'ml obu           ', ptr_patch=this%obu           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_obu_gah'            , units='', avgflag='A', long_name=  'ml obu_gah       ', ptr_patch=this%obu_gah       , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_obuold'             , units='', avgflag='A', long_name=  'ml obuold        ', ptr_patch=this%obuold        , set_lake=0._r8, set_urb=0._r8, default='inactive')
    !call hist_addfld1d (fname= 'ml_nmozsgn'            , units='', avgflag='A', long_name=  'ml nmozsgn       ', ptr_patch=this%nmozsgn       , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_z0mg'               , units='', avgflag='A', long_name=  'ml z0mg          ', ptr_patch=this%z0mg          , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_thref'              , units='', avgflag='A', long_name=  'ml thref         ', ptr_patch=this%thref         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_thvref'             , units='', avgflag='A', long_name=  'ml thvref        ', ptr_patch=this%thvref        , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_gah'                , units='', avgflag='A', long_name=  'ml gah           ', ptr_patch=this%gah           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_PrSc'               , units='', avgflag='A', long_name=  'ml PrSc          ', ptr_patch=this%PrSc          , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_Lc'                 , units='', avgflag='A', long_name=  'ml Lc            ', ptr_patch=this%Lc            , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_zdisp'              , units='', avgflag='A', long_name=  'ml zdisp         ', ptr_patch=this%zdisp         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_tstar'              , units='', avgflag='A', long_name=  'ml tstar         ', ptr_patch=this%tstar         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_qstar'              , units='', avgflag='A', long_name=  'ml qstar         ', ptr_patch=this%qstar         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_rnsoi'              , units='', avgflag='A', long_name=  'ml rnsoi         ', ptr_patch=this%rnsoi         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_shsoi'              , units='', avgflag='A', long_name=  'ml shsoi         ', ptr_patch=this%shsoi         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_lhsoi'              , units='', avgflag='A', long_name=  'ml lhsoi         ', ptr_patch=this%lhsoi         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_gsoi'               , units='', avgflag='A', long_name=  'ml gsoi          ', ptr_patch=this%gsoi          , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_irsoi'              , units='', avgflag='A', long_name=  'ml irsoi         ', ptr_patch=this%irsoi         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_etsoi'              , units='', avgflag='A', long_name=  'ml etsoi         ', ptr_patch=this%etsoi         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_tg'                 , units='', avgflag='A', long_name=  'ml tg            ', ptr_patch=this%tg            , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_btran'              , units='', avgflag='A', long_name=  'ml btran         ', ptr_patch=this%btran         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_psis'               , units='', avgflag='A', long_name=  'ml psis          ', ptr_patch=this%psis          , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_rsoil'              , units='', avgflag='A', long_name=  'ml rsoil         ', ptr_patch=this%rsoil         , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_eg'                 , units='', avgflag='A', long_name=  'ml eg            ', ptr_patch=this%eg            , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_rhg'                , units='', avgflag='A', long_name=  'ml rhg           ', ptr_patch=this%rhg           , set_lake=0._r8, set_urb=0._r8, default='inactive')
    call hist_addfld1d (fname= 'ml_qflx_prec_intr'     , units='', avgflag='A', long_name=  'ml qflx_prec_intr', ptr_patch=this%qflx_prec_intr, set_lake=0._r8, set_urb=0._r8, default='inactive')



!    this%gppveg(begp:endp) = spval
!    call hist_addfld1d (fname='GPPVEG', units='umol/m2s', &
!         avgflag='A', long_name='Gross primary production', &
!         ptr_patch=this%gppveg, set_lake=0._r8, set_urb=0._r8, default='inactive')
!
!    this%lwp(begp:endp,1:nlevcanml) = spval
!    call hist_addfld2d (fname='LWP', units='MPa', type2d='nlevcanml', &
!         avgflag='A', long_name='Leaf water potential of canopy layer', &
!         ptr_patch=this%lwp, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold (this, bounds)
    !
    ! !DESCRIPTION:
    ! Cold-start initialization for multilayer canopy
    !
    ! !USES:
    use clm_varpar, only : nlevcanml,nleaf
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: iv               ! Canopy leaf layer index
    !---------------------------------------------------------------------

    ! Initialize leaf water potential and intercepted water

    do p = bounds%begp, bounds%endp
          this%wind       (p,0)    = 0._r8
          this%wind_1stoldH       (p,0)    = 1._r8/float(nlevcanml)
          this%wind_1st (p,0)    = 0._r8
          this%wind_1stold (p,0)    = 0._r8
          this%wind2LB     (p,0)    = 0._r8
          this%tair       (p,0)    = 25._r8+273.15_r8
          this%eair       (p,0)    = 3395._r8
          this%cair       (p,0)    = 380._r8
          this%ga_prof       (p,0)    = 1._r8

          this%zref_old       (p) = 44._r8
          this%tair_old       (p,0)    = 25._r8+273.15_r8
          this%eair_old       (p,0)    = 3395._r8
          this%cair_old       (p,0)    = 380._r8
          this%tveg_old       (p,0,1)  = 25._r8+273.15_r8
          this%tveg_old       (p,0,2)  = 25._r8+273.15_r8

          this%tveg       (p,0,1)  = 25._r8+273.15_r8
          this%tleaf      (p,0,1)  = 25._r8+273.15_r8
          this%tveg       (p,0,2)  = 25._r8+273.15_r8
          this%tleaf      (p,0,2)  = 25._r8+273.15_r8
          this%tg         (p)      = 25._r8+273.15_r8


       do iv = 1, nlevcanml
          this%lwp(p,iv) = -0.1_r8
          this%h2ocan(p,iv) = 0._r8

          this%wind       (p,iv)    = 1._r8*float(iv)/float(nlevcanml)
          this%wind_1stoldH       (p,iv)    = 1._r8/float(nlevcanml)
          this%wind_1st       (p,iv)    = 1._r8*float(iv)/float(nlevcanml)
          this%wind_1stold       (p,iv)    = 1._r8*float(iv)/float(nlevcanml)
          this%wind2LB       (p,iv)    = 1._r8*float(iv)/float(nlevcanml)
          this%tair       (p,iv)    = 25._r8+273.15_r8
          this%eair       (p,iv)    = 3395._r8
          this%cair       (p,iv)    = 380._r8
          this%ga_prof       (p,iv)    = 1._r8

          this%tair_old       (p,iv)    = 25._r8+273.15_r8
          this%eair_old       (p,iv)    = 3395._r8
          this%cair_old       (p,iv)    = 380._r8
          this%tveg_old       (p,iv,1)  = 25._r8+273.15_r8
          this%tleaf_old      (p,iv,1)  = 25._r8+273.15_r8
          this%tveg_old       (p,iv,2)  = 25._r8+273.15_r8
          this%tleaf_old      (p,iv,2)  = 25._r8+273.15_r8

          this%tveg       (p,iv,1)  = 25._r8+273.15_r8
          this%tleaf      (p,iv,1)  = 25._r8+273.15_r8
          this%tveg       (p,iv,2)  = 25._r8+273.15_r8
          this%tleaf      (p,iv,2)  = 25._r8+273.15_r8
       end do
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod, only : restartvar
    use clm_varpar, only : nlevcanml,nleaf
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !---------------------------------------------------------------------
    call InitCold (this, bounds)
    ! Example for 2-d patch variable

    !call restartvar(ncid=ncid, flag=flag, varname='lwp', xtype=ncd_double,  &
    !   dim1name='pft', dim2name='levcan', switchdim=.true., &
    !   long_name='leaf water potential of canopy layer', units='MPa', &
    !   interpinic_flag='interp', readvar=readvar, data=this%lwp)

  end subroutine Restart

end module CanopyFluxesMultilayerType
