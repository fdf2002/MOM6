!> Implemented non-local mixing due to Ocean Thermal Energy Conversion
module MOM_otec

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : register_static_field, time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public interior_mass_sink, otec_init

!> Control structure for OTEC
type, public :: otec_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: w_otec !< Pumping rate [Z T-1 ~> m s-1]
  logical :: apply_otec !< If true, OTEC will be applied.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing

end type otec_CS

contains

!> Applies OTEC...
!! Add description here.
subroutine interior_mass_sink(h, tv, dt, G, GV, US, CS, halo)
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
  integer :: isj, iej, num_left, nkmb, k_tgt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_otec: "//&
         "Module must be initialized before it is used.")
  
  if (.not.CS%apply_otec) return

  ! Fixed value of k (for now) !
  ! THIS WILL NOT STAY         !
  k = 2

  do j=js,je
    do i=is,ie
      ! Remove mass from layer k at a rate set by w_otec
      h(i, j, k) = h(i, j, k) - CS%w_otec * dt
    enddo ! i-loop
  enddo ! j-loop

end subroutine interior_mass_sink

!> Initialize parameters and allocate memory associated with the OTEC module.
subroutine otec_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure used to regulate diagnostic output.
  type(otec_CS),     intent(inout) :: CS   !< OTEC heating control struct

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_otec"  ! module name
  character(len=48)  :: thickness_units
  ! Local variables
  character(len=200) :: inputdir, otec_file, filename, otec_var
  real :: w_otec  ! A uniform pumping rate [Z T-1 ~> m s-1]
  integer :: i, j, isd, ied, jsd, jed, id
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time

  ! write parameters to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "W_OTEC", w_otec, &
                 "The constant OTEC pumping rate or a rescaling or "//&
                 "0 to disable the OTEC pumping.", &
                 units="m s-1", default=0.0)
  CS%w_otec = w_otec
  CS%apply_otec = .not.(w_otec == 0.0)
  if (.not.CS%apply_otec) return

end subroutine otec_init

end module MOM_otec
