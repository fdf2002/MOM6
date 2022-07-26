!> Implemented non-local mixing due to Ocean Thermal Energy Conversion
module MOM_otec

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : register_static_field, time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_density_derivs
use MOM_EOS,           only : EOS_type

implicit none ; private

#include <MOM_memory.h>

public otec_step, otec_init

!> Control structure for OTEC
type, public :: otec_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: apply_otec !< If true, OTEC will be applied.

  !! OTEC input variables
  real    :: w_cw !< Cold-water pumping rate [Z T-1 ~> m s-1]
  real    :: w_ww !< Warm-water pumping rate [Z T-1 ~> m s-1]

  real    :: depth_cold, depth_warm, depth_out !< Pipe depths [m]
  

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing

end type otec_CS

!> Control structure that stores some variables from thermo_var_ptrs,
!! but only for a single grid cell, so as to make column calculations more efficient.
type, private :: thermo_var_1d
  ! If allocated, the following variables have nz layers.
  real, pointer :: T(:) => NULL() !< Potential temperature [degC].
  real, pointer :: S(:) => NULL() !< Salinity [PSU] or [gSalt/kg], generically [ppt].
  real          :: p_surf = 0.0 !< Ocean surface pressure used in equation of state
                                !! calculations [R L2 T-2 ~> Pa]

  type(EOS_type), pointer :: eqn_of_state => NULL() !< Type that indicates the
                                                    !! equation of state to use.
  real :: C_p            !<   The heat capacity of seawater [Q degC-1 ~> J degC-1 kg-1].
                         !! When conservative temperature is used, this is
                         !! constant and exactly 3991.86795711963 J degC-1 kg-1.

end type thermo_var_1d

contains

!> Given a depth, finds the appropriate layer that contains that depth.
subroutine find_layer(h1d, GV, target_depth, k, layer_depth)
  type(verticalGrid_type),   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)), intent(in)  :: h1d !< Layer thicknesses at the grid cell [H ~> m or kg m-2].
  real,                      intent(in)  :: target_depth !< The depth we are searching for [H ~> m or kg m-2].

  integer, intent(out) :: k !< The layer corresponding to a depth of target_depth.
                            !! If k > GV%ke, then the depth does not exist at this location.
  real,    intent(out) :: layer_depth !< The depth of layer k [H ~> m or kg m-2].

  k = 0
  layer_depth = 0.0

  do while (layer_depth <= target_depth)
    k = k + 1
    if (k > GV%ke) return ! ocean is not deep enough
    layer_depth = layer_depth + h1d(k)
  enddo

end subroutine find_layer

!> Drains the layer at a cell, starting at a minimum depth of z.
!! If this depletes the layer fully, then uses layer k+1 to finish draining.
!! Upon depleting the deepest layer, stops with a warning in stdout.
subroutine mass_sink(h1d, tv1d, GV, sink_depth, dThickness, &
                                    netMassOut, netSaltOut, netHeatOut)
  type(verticalGrid_type),   intent(in)    :: GV !< The ocean's vertical grid structure.

  ! 1-dimensional copies of arrays from a thermo_vars_ptr object.
  real, dimension(SZK_(GV)), intent(inout) :: h1d !< Layer thicknesses at the grid cell [H ~> m or kg m-2]
  type(thermo_var_1d),       intent(in)    :: tv1d !< 1-dimensional copies of S and T

  ! Mass Sink Parameters
  real,                      intent(in)    :: sink_depth !< The depth of this mass sink [H ~> m or kg m-2].
  real,                      intent(inout) :: dThickness !< Amount to change layer thickness [H ~> m or kg m-2]
                                                         !! Must be negative at the start of the subroutine.
                                                         !! Upon returning, dThickness is the actual change in thickness.
  real,                      intent(inout) :: netMassOut !< The total mass being extracted per unit area [kg m-2].
  real,                      intent(inout) :: netSaltOut !< The total amount of salt being extracted
                                                         !! [ppt H ~> ppt m or ppt kg m-2].
  real,                      intent(inout) :: netHeatOut !< The total heat being extracted [degC H ~> degC m or degC kg m-2].

  integer :: k
  real    :: layer_depth, dh, m, maximum_drainage, pressure, rho, dh_total

  call find_layer(h1d, GV, sink_depth, k, layer_depth)
  if (k > GV%ke) then
    call MOM_error(WARNING, "MOM_otec: Ocean floor reached before intake.")
    return
  endif
  pressure = tv1d%p_surf + GV%H_to_RZ * GV%g_Earth * sink_depth

  dh_total = dThickness
  do while (dh_total < 0) ! as long as there is still more to take out

    ! The maximum drainage from this layer is everything below sink_depth.
    ! Ensure the layer thickness is always at least Angstrom.
    maximum_drainage = min(h1d(k) - GV%Angstrom_H, layer_depth - sink_depth)
    ! Drain as much as specified by input, or until the layer vanishes.
    dh = max(dh_total, -maximum_drainage)
    h1d(k) = max(GV%Angstrom_H, h1d(k) + dh)
    dh_total = dh_total - dh
    layer_depth = layer_depth - dh ! Layer bottom has moved up

    call calculate_density(tv1d%T(k), tv1d%S(k), pressure, rho, tv1d%eqn_of_state)
    ! Net out tracers will be per area, i.e. [kg m-2].
    m = -rho*dh

    netMassOut = netMassOut + m
    netSaltOut = netSaltOut + m*tv1d%S(k)
    netHeatOut = netHeatOut + m*tv1d%T(k)

    ! Increment for the next iteration
    k = k + 1

    ! If ocean bottom is reached
    if (k > GV%ke) then
      call MOM_error(WARNING, "MOM_otec: Ocean floor reached during intake.")
      return
    endif

    layer_depth = layer_depth + h1d(k)

  enddo

  ! Update dThickness to reflect actual amount removed
  dThickness = min(dThickness-dh_total, -GV%H_subroundoff)

end subroutine mass_sink



subroutine mass_source(h1d, tv1d, GV, src_depth, dThickness, netMassIn, netSaltIn, netHeatIn)
  type(verticalGrid_type),   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)), intent(inout) :: h1d !< Layer thicknesses at the grid cell [H ~> m or kg m-2]
  type(thermo_var_1d),       intent(in)    :: tv1d !< A structure containing pointers
                                                 !! to any available thermodynamic fields.

  real,                      intent(in)    :: src_depth !< The depth of this mass sink [H ~> m or kg m-2].
  real,                      intent(in)    :: dThickness !< The change in layer thickness [H ~> m or kg m-2].

  real, optional,            intent(in)    :: netMassIn !< The total MASS being added per unit area [kg m-2].
  real, optional,            intent(in)    :: netSaltIn !< The total amount of salt being added with the water
                                                        !! [ppt H ~> ppt m or ppt kg m-2].
  real, optional,            intent(in)    :: netHeatIn !< The total heat content of the water being added
                                                        !! [degC H ~> degC m or degC kg m-2].

  ! Local variables
  integer :: k
  real :: layer_depth, &
          iMassIn, & ! Inverse of the mass being added
          T_add, S_add, p_add, rho_add, iRho_add, & ! Properties of the water being added.
          p_k, rho_k, iRho_k, oldMass, & ! Properties of the layer before water is added.
          iNewMass ! Inverse of total mass after injection [m2 kg-1].

  k = 0; layer_depth = 0.0
  call find_layer(h1d, GV, src_depth, k, layer_depth)
  if (k > GV%ke) then
    call MOM_error(WARNING, "MOM_otec: Ocean floor reached before returning water.")
    return
  endif
  
  iMassIn = 1./netMassIn
  ! Find pressure and density of the water being added.
  T_add = netHeatIn * iMassIn
  S_add = netSaltIn * iMassIn
  p_add = tv1d%p_surf + GV%H_to_RZ * GV%g_Earth * src_depth
  call calculate_density(T_add, S_add, p_add, rho_add, tv1d%eqn_of_state)
  iRho_add = 1./rho_add ! The inverse of density [m3 kg-1]

  ! Now find average pressure & density for the layer as a whole.
  p_k = tv1d%p_surf + GV%H_to_RZ*GV%g_Earth * (layer_depth - 0.5*h1d(k))
  call calculate_density(tv1d%T(k), tv1d%S(k), p_k, rho_k, tv1d%eqn_of_state)
  iRho_k = 1./rho_k
  oldMass = rho_k*h1d(k)
  iNewMass = 1./(oldMass + netMassIn)

  ! Update mass and tracers of the layer.
  h1d(k) = h1d(k) + dThickness
  tv1d%S(k) = (oldMass*tv1d%S(k) + netSaltIn) * iNewMass
  tv1d%T(k) = (oldMass*tv1d%T(k) + netHeatIn) * iNewMass

end subroutine mass_source

!> Applies two mass sinks in every lateral grid cell.
!! Add description here.
subroutine otec_step(h, tv, dt, G, GV, US, CS, halo)
  type(ocean_grid_type),                     intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv !< A structure containing pointers
                                                                 !! to any available thermodynamic fields.
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(otec_CS),                             intent(in)    :: CS !< The control structure returned by
                                                           !! a previous call to
                                                           !! otec_init.
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work
  ! Local variables

  integer :: i, j, k, is, ie, js, je, nz, k2
  real :: dh_cold, dh_warm, dh_mixed
  real :: massTransport, saltTransport, heatTransport
  real, dimension(SZK_(GV)) :: h1d
  real, dimension(SZK_(GV)), target :: T1d, S1d
  type(thermo_var_1d) :: tv1d !< 1-dimensional copy of thermodynamic fields

  tv1d%C_p = tv%C_p
  tv1d%eqn_of_state => tv%eqn_of_state
  tv1d%T => T1d
  tv1d%S => S1d

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_otec: "//&
         "Module must be initialized before it is used.")

  do k=1,GV%ke
    !print *, "thickness of layer",k, "=",h(2,2,k)
  enddo
  
  if (.not.CS%apply_otec) return

  !print *, "OTEC is initialized"

  do j=js,je
    do i=is,ie

      tv1d%p_surf = tv%p_surf(i,j)
      ! Copy this column into a 1D array (for runtime efficiency)
      do k=1,GV%ke
        h1d(k) = h(i,j,k)
        T1d(k) = tv%T(i,j,k)
        S1d(k) = tv%S(i,j,k)
        !print *, "thickness of layer", k, "=", h1d(k)
      enddo

      ! Prepare tracers to be moved between layers
      massTransport = 0.0
      saltTransport = 0.0
      heatTransport = 0.0

      ! Remove warm water from the surface layer
      dh_cold = -CS%w_cw * dt
      dh_warm = -CS%w_ww * dt

      call mass_sink(h1d, tv1d, GV, CS%depth_cold, dh_cold, &
                      massTransport, saltTransport, heatTransport)
      !print *, "m=",massTransport, " salt=",saltTransport, " heat=",heatTransport
      call mass_sink(h1d, tv1d, GV, CS%depth_warm, dh_warm, &
                      massTransport, saltTransport, heatTransport)
      ! We conserve thickness, not mass.
      dh_mixed = -(dh_cold + dh_warm)
      call mass_source(h1d, tv1d, GV, CS%depth_out, dh_mixed, &
                        massTransport, saltTransport, heatTransport)

      ! Copy the 1D working arrays back into the originals
      do k=1,GV%ke
        h(i,j,k) = h1d(k)
        tv%T(i,j,k) = T1d(k)
        tv%S(i,j,k) = S1d(k)
      enddo
    enddo ! i-loop
  enddo ! j-loop

end subroutine otec_step

!> Initialize parameters and allocate memory associated with the OTEC module.
subroutine otec_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure used to regulate diagnostic output.
  type(otec_CS),     intent(inout)       :: CS   !< OTEC heating control struct

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_otec"  ! module name
  character(len=48)  :: thickness_units
  ! Local variables
  character(len=200) :: inputdir, otec_file, filename, otec_var
  real :: w_cw  ! A uniform pumping rate [Z T-1 ~> m s-1]
  real :: gamma ! Ratio of warm- to cold-water pumping rate
  real :: depth_cold, depth_warm, depth_out ! Depths of the pipes [m]
  integer :: i, j, isd, ied, jsd, jed, id
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time

  ! write parameters to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "OTEC_W_CW", w_cw, &
                 "The constant OTEC cold-water pumping rate or 0 to "//&
                 "disable OTEC.", &
                 units="m s-1", default=0.0)
  call get_param(param_file, mdl, "OTEC_GAMMA", gamma, &
                 "The ratio of warm- to cold-water OTEC pumping rates.",&
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "OTEC_COLD_DEPTH", depth_cold, &
                  "The depth of the cold-water intake for OTEC.", &
                  units="m", default=1000.0)
  call get_param(param_file, mdl, "OTEC_WARM_DEPTH", depth_warm, &
                  "The depth of the warm-water intake for OTEC.", &
                  units="m", default=20.0)
  call get_param(param_file, mdl, "OTEC_OUT_DEPTH", depth_out, &
                  "The depth of the mixed-water output for OTEC.", &
                  units="m", default=500.0)

  CS%w_cw = w_cw
  CS%w_ww = gamma*w_cw

  CS%depth_cold = depth_cold
  CS%depth_warm = depth_warm
  CS%depth_out  = depth_out

  CS%apply_otec = .not.(w_cw == 0.0)

end subroutine otec_init

end module MOM_otec
