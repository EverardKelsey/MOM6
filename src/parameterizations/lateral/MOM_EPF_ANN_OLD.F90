!> Calculates horizontal viscosity and viscous stresses
module MOM_EPF_ANN

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_grid,                  only : ocean_grid_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_diag_mediator,         only : diag_ctrl, time_type, post_data, register_diag_field
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_domains,               only : create_group_pass, do_group_pass, group_pass_type, &
                                      start_group_pass, complete_group_pass, &
                                      To_North, To_East, pass_var, CORNER, pass_vector
use MOM_coms,                  only : reproducing_sum
use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                                      CLOCK_MODULE, CLOCK_ROUTINE
use MOM_ANN,                   only : ANN_init, ANN_apply, ANN_end, ANN_CS
use MOM_error_handler,         only : MOM_mesg 
use MOM_spatial_means,         only : global_volume_mean
use MOM_interface_heights,     only : find_eta
use MOM_lateral_mixing_coeffs, only : VarMix_CS, calc_slope_functions
use MOM_variables,             only : thermo_var_ptrs, cont_diag_ptrs
use MOM_open_boundary,         only : ocean_OBC_type
use MOM_isopycnal_slopes,      only : calc_isoneutral_slopes

implicit none ; private

#include <MOM_memory.h>

public EPF_lateral_stress, EPF_init, EPF_end, EPF_copy_gradient_and_thickness

!> Control structure for EPF parameterization.
type, public :: EPF_CS ; private
  ! Parameters
  real      :: amplitude        !< The nondimensional scaling factor in EPF model (tuning parameter),
                                !! typically 0.1 - 10 [nondim].
  real      :: DT               !< Baroclinic timestep
  real      :: H_cutoff         !< cutoff depth from which isopycnal slopes are calculated, should be min ocean depth as default
  real      :: GFS              !< Reduced gravity at free surface
  real      :: GINT             !< reduced gravity at layer interface
  logical   :: fz_negative   !< If true, the vertical convergence of dual form stresses subtracted from diffu and diffv
  logical   :: add_EPF_fxy       !< update the diffusivities with horizontal divergence of EPF? 
  ! logical   :: add_EPF_fy       !< update the diffusivities with horizontal divergence of EPF? 
  ! logical   :: add_EPF_fyx       !< update the diffusivities with horizontal divergence of EPF? 
  ! logical   :: add_EPF_fyy       !< update the diffusivities with horizontal divergence of EPF? 
  logical   :: add_EPF_fz       !< update the diffusivities with vertical divergence of EPF?
  logical   :: add_EPF_subPE    !< update the diffusivities with horizontal divergence of subfilter PE?
  
  ! allocating memory in the heap. Should allocate memory here for things needed in subsequent timesteps
  real, dimension(:,:,:), allocatable :: &
          f_u, &       !< the zonal diffusivity that is updated by EPF divergence
          f_v, &       !< the meridional diffusivity that is updated by EPF divergence
          f_x, &        !< the zonal acceleration due to horiztonal divergence of EPF 
          f_y, &        !< the meridional acceleration due to horiztonal divergence of EPF 
          fx_PE, &        !< the zonal acceleration due to horiztonal divergence of EPF 
          fy_PE, &        !< the meridional acceleration due to horizontal divergence of EPF
          fxz, &      !< the zonal acceleration due to vertical "divergence" of EPF
          fyz         !< the meridional acceleration due to vertical "divergence" of EPF

  real, dimension(:,:,:), allocatable :: &
          sh_xx,   &   !< Horizontal tension (du/dx - dv/dy) in h (CENTER)
                       !! points including metric terms [T-1 ~> s-1]
          sh_xy,   &   !< Horizontal shearing strain (du/dy + dv/dx) in q (CORNER)
                       !! points including metric terms [T-1 ~> s-1]
          vort_xy, &   !< Vertical vorticity (dv/dx - du/dy) in q (CORNER)
                       !! points including metric terms [T-1 ~> s-1]
          hq           !< Thickness in CORNER points [H ~> m or kg m-2]


  real, dimension(:,:,:), allocatable :: &
          Txy_h, &
          Txx,     & !< Subgrid stress xx component in h [L2 T-2 ~> m2 s-2]
          Tyy,     & !< Subgrid stress yy component in h [L2 T-2 ~> m2 s-2]
          Txx_PE,     & !< Subgrid stress xx component in h [L2 T-2 ~> m2 s-2]
          Tyy_PE,     & !< Subgrid stress yy component in h [L2 T-2 ~> m2 s-2]
          Txy,     & !< Subgrid stress xy component in q [L2 T-2 ~> m2 s-2]
          Txz,     & !< Subgrid x component of form stress in h [L2 T-2 ~> m2 s-2]
          Tyz      !& !< Subgrid y component of form stress in h [L2 T-2 ~> m2 s-2]

  real, dimension(:,:,:), allocatable :: &
          mom_norm_h, &
          buoy_norm_h, &
          depth_mask, &      !< the depth mask used when applying ANN at h points
          slope_x, &         !< zonal gradient of interface at u points
          slope_y, &         !< zmeridional gradient of interface at v points
          slope_x_top,     & !< Zonal gradient of top interface at h points
          slope_x_bot,     & !< zonal gradient of bottom interface at h points
          slope_y_top,     & !< meridional gradient of top interface at h points
          slope_y_bot,     & !< meridional gradient of bottom interface at h points
          sh_xy_h,         & !< horizontal shearing strain in h
          vort_xy_h,       & !< vertical vorticity in h
          eta,             & !< interface location as calculated using find_eta
          h_eta              !< the thickness used as input to find_eta function
          
  real, dimension(:,:), allocatable :: &
          Delsq_h, &
          Coriolis, &
          kappa_h, & !< Scaling coefficient in h points [L2 ~> m2]
          kappa_q    !< Scaling coefficient in q points [L2 ~> m2]

  real, allocatable ::    &
        Coriolis_h(:,:)     !< Coriolis parameter at h points [T ~> s]

  !integer :: use_ann  !< 0: ANN is turned off, 1: default ANN for EPF
  integer :: ann_type  !< 0: ANN is layer agnostic, 1: ANN is a column-type model, 2: two layer agnostic ANNs
  integer :: n_inputs !< Number of inputs to the ANN, default is 7
  integer :: n_outputs !< Number of outputs from the ANN, default is 5
  integer :: n_inputs_intfc !< Number of inputs to the interface ANN, default is 5
  integer :: n_outputs_intfc !< Number of outputs from the interface ANN, default is 2
  integer :: n_inputs_centre !< Number of inputs to the centre fluxANN, default is 3
  integer :: n_outputs_centre !< Number of outputs from the centre ANN, default is3
  
  type(ANN_CS) :: ann_instance !< ANN instance
  type(ANN_CS) :: ann_instance_intfc !< ANN instance
  type(ANN_CS) :: ann_instance_centre !< ANN instance
  character(len=200)   :: ann_file
  character(len=200)   :: intfc_ann_file
  character(len=200)   :: centre_ann_file
  !character(len=200) :: ann_file = "/scratch/kae10022/Training_ANNs/trained_ANNs/NET_64_64.nc" !< Default test ANN with EPF model for column model

  real :: subroundoff_shear

  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_Txx = -1
  integer :: id_Tyy = -1
  integer :: id_Txx_PE = -1
  integer :: id_Tyy_PE = -1
  integer :: id_Txy = -1
  integer :: id_Txy_h = -1
  integer :: id_Txz = -1
  integer :: id_Tyz = -1
  integer :: id_slope_x_top = -1, id_slope_y_top = -1
  integer :: id_slope_x_bot = -1, id_slope_y_bot = -1
  integer :: id_slope_x = -1, id_slope_y = -1
  integer :: id_sh_xx = -1, id_sh_xy = -1, id_vort_xy = -1
  integer :: id_sh_xy_h = -1, id_vort_xy_h = -1
  integer :: id_eta = -1, id_h_eta = -1
  integer :: id_depth_mask = -1
  integer :: id_Coriolis_h = -1
  integer :: id_Delsq_h = -1
  integer :: id_mom_norm_h = -1, id_buoy_norm_h = -1
  integer :: id_f_x = -1, id_f_y = -1
  integer :: id_fxz = -1, id_fyz = -1
  integer :: id_f_u = -1, id_f_v = -1
  integer :: id_fx_PE = -1, id_fy_PE = -1
  !>@}

  !>@{ CPU time clock IDs
  integer :: id_clock_module
  integer :: id_clock_copy
  integer :: id_clock_cdiss
  integer :: id_clock_stress
  integer :: id_clock_stress_ANN
  integer :: id_clock_divergence
  integer :: id_clock_mpi
  integer :: id_clock_filter
  integer :: id_clock_post
  integer :: id_clock_source
  !>@}

  !>@{ MPI group passes
  type(group_pass_type) :: &
      pass_Tq, pass_Th, &        !< handles for halo passes of Txy and Txx, Tyy
      pass_xx, pass_xy           !< handles for halo passes of sh_xx and sh_xy, vort_xy
  !integer :: Stress_halo = -1!, & !< The halo size in filter of the stress tensor
            ! HPF_halo = -1       !< The halo size in filter of the velocity gradient
  !>@}

end type EPF_CS

contains
!> Read parameters, allocate and precompute arrays,
!! register diagnosicts used in Zanna_Bolton_2020().
subroutine EPF_init(Time, G, GV, US, param_file, diag, CS, use_EPF_ANN)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(EPF_CS),            intent(inout) :: CS         !< EPF control structure.
  logical,                 intent(out)   :: use_EPF_ANN !< If true, turns on EPF parameterisation.
  

  real :: subroundoff_Cor     ! A negligible parameter which avoids division by zero
                              ! but small compared to Coriolis parameter [T-1 ~> s-1]

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: i, j

  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_EPF_ANN" ! This module's name.

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "USE_EPF_ANN", use_EPF_ANN, &
                 "If true, turns on EPF " //&
                 "subgrid momentum parameterization of mesoscale eddies.", default=.false.)
  if (.not. use_EPF_ANN) return

  call get_param(param_file, mdl, "EPF_ANN_TYPE", CS%ann_type, &
                 "If 0, use single layer agnostic ANN. If 1, use column-model type ANN. If 2, " //&
                 "use two layer agnostic ANNs (one for interface fluxes, one for central).", fail_if_missing=.true.)
  if (CS%ann_type==2) then
    call get_param(param_file, mdl, "INTFC_ANN_FILE", CS%intfc_ann_file, &
                    "pathway to the saved ann to use for parameterisation of interface fluxes", &
                    fail_if_missing=.true.)
    call get_param(param_file, mdl, "N_INPUTS_INTFC", CS%n_inputs_intfc, &
                    "Number of input features to interface flux ANN", default=5)
    call get_param(param_file, mdl, "N_OUTPUTS_INTFC", CS%n_outputs_intfc, &
                    "Number of output features to interface flux ANN", default=2)
    call get_param(param_file, mdl, "CENTRE_ANN_FILE", CS%centre_ann_file, &
                    "pathway to the saved ann to use for parameterisation of layer centre fluxes", &
                    fail_if_missing=.true.)   
    call get_param(param_file, mdl, "N_INPUTS_CENTRE", CS%n_inputs_centre, &
                    "Number of input features to centre flux ANN", default=3)
    call get_param(param_file, mdl, "N_OUTPUTS_CENTRE", CS%n_outputs_centre, &
                    "Number of output features to centre flux ANN", default=3)
    call ANN_init(CS%ann_instance_intfc, CS%intfc_ann_file)
    call ANN_init(CS%ann_instance_centre, CS%centre_ann_file)
  else
    call get_param(param_file, mdl, "ANN_FILE", CS%ann_file, &
                  "Pathway to the saved ann to use in parameterisation", &
                  fail_if_missing=.false.)
    call get_param(param_file, mdl, "N_INPUTS", CS%n_inputs, &
                  "Number of input features to ANN", default=7)
    call get_param(param_file, mdl, "N_OUTPUTS", CS%n_outputs, &
                  "Number of out features from ANN", default=5)
    call ANN_init(CS%ann_instance, CS%ann_file)
  end if
  call get_param(param_file, mdl, "ADD_FXY", CS%add_EPF_fxy, &
                 "If true, updates diffusivities with horizontal divergence of EPF.", &
                  default=.false.)

  call get_param(param_file, mdl, "ADD_FZ", CS%add_EPF_fz, &
                  "If true, updates diffusivities with vertical divergence of EPF.", &
                   default=.false.)  
  
  call get_param(param_file, mdl, "ADD_subPE", CS%add_EPF_subPE, &
                  "If true, add the subfilter PE term to the Txx and Tyy terms", &
                  default = .false.)

  call get_param(param_file, mdl, "EPF_SCALING", CS%amplitude, &
                 "The nondimensional tuning parameter for EPF ANN, " //&
                 "typically 0.5-2.5", units="nondim", default=0.5)
             
  call get_param(param_file, mdl, "DT", CS%dt, &
                 "The (baroclinic) dynamics time step.", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)

  call get_param(param_file, mdl, "GFS", CS%GFS, &
                 "Reduced gravity at free surface", units="m s-2", scale=US%L_T2_to_m_s2, &
                 fail_if_missing=.true.)
                 
  call get_param(param_file, mdl, "GINT", CS%GINT, &
                 "Reduced gravity at interface between layers", units="m s-2", scale=US%L_T2_to_m_s2, &
                 fail_if_missing=.true.)

  call get_param(param_file, mdl, "FZ_NEGATIVE", CS%fz_negative, &
                 "Subtract form stress convergence?", &
                 default=.true.)

  call get_param(param_file, mdl, "H_CUTOFF", CS%H_cutoff, &
                 "The minimum ocean depth for use in calculation of interface slopes", units="m", &
                 scale = US%Z_to_L, default = 1.0)

  CS%diag => diag
  !!type(axes_grp)  :: axesBL, axesTL, axesCuL, axesCvL
  ! L indicates that it is a layer variable
  ! i indicates that it is an interface variable
  ! B indicates a corner point
  ! T indicates centre point
  ! Cu indicates centered at u points
  ! Cv indicates centered at v points
  !!type(axes_grp)  :: axesBi, axesTi, axesCui, axesCvi
  CS%id_Txx = register_diag_field('ocean_model', 'Txx', diag%axesTL, Time, &
      'Diagonal term (Txx) in the EPF stress tensor (Reynolds stress component)', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Tyy = register_diag_field('ocean_model', 'Tyy', diag%axesTL, Time, &
      'Diagonal term (Tyy) in the EPF stress tensor (Reynolds stress component)', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Txx_PE = register_diag_field('ocean_model', 'Txx_PE', diag%axesTL, Time, &
      'Diagonal term (Txx_PE) in the EPF stress tensor corresponding to subfilter potential energy', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Tyy_PE = register_diag_field('ocean_model', 'Tyy_PE', diag%axesTL, Time, &
      'Diagonal term (Tyy_PE) in the EPF stress tensor corresponding to subfilter potential energy', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Txy = register_diag_field('ocean_model', 'Txy', diag%axesBL, Time, &
      'Off-diagonal term (Txy) in the EPF stress tensor', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Txy_h = register_diag_field('ocean_model', 'Txy_h', diag%axesTL, Time, &
      'Off-diagonal term (Txy_h) in the EPF stress tensor at the grid centre', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Txz = register_diag_field('ocean_model', 'Txz', diag%axesTL, Time, &
      'Zonal form stress difference in the EPF stress tensor', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  CS%id_Tyz = register_diag_field('ocean_model', 'Tyz', diag%axesTL, Time, &
      'Meridional form stress difference in the EPF stress tensor', 'm2 s-2', conversion=US%L_T_to_m_s**2)

  !CS%id_fx = register_diag_field('ocean_model', 'fx', diag%axes)

  ! Registering fields for debugging purposes
  CS%id_slope_x = register_diag_field('ocean_model', 'slope_x', diag%axesCui, Time, &
      'zonal slope of the interfaces at u points', 'nondim')

  CS%id_slope_y = register_diag_field('ocean_model', 'slope_y', diag%axesCvi, Time, &
      'zonal slope of the interfaces at v points', 'nondim')

  CS%id_slope_x_top = register_diag_field('ocean_model', 'slope_x_top', diag%axesTL, Time, &
      'zonal slope of the top interfaces at h points', 'nondim')

  CS%id_slope_y_top = register_diag_field('ocean_model', 'slope_y_top', diag%axesTL, Time, &
      'meridional slope of the top interfaces at h points', 'nondim')

  CS%id_slope_x_bot = register_diag_field('ocean_model', 'slope_x_bot', diag%axesTL, Time, &
      'zonal slope of the bottom interfaces at h points', 'nondim')

  CS%id_slope_y_bot = register_diag_field('ocean_model', 'slope_y_bot', diag%axesTL, Time, &
      'meridional slope of the bottom interfaces at h points', 'nondim')

  CS%id_sh_xx = register_diag_field('ocean_model','sh_xx', diag%axesTL, Time, &
      'Horizontal tension (du/dx - dv/dy) in h', 'm s-1',conversion=US%L_T_to_m_s)

  CS%id_sh_xy = register_diag_field('ocean_model','sh_xy', diag%axesBL, Time, &
      'Horizontal shearing strain (du/dy + dv/dx) in q', 'm s-1',conversion=US%L_T_to_m_s)

  CS%id_vort_xy = register_diag_field('ocean_model','vort_xy', diag%axesBL, Time, &
      'Vertical vorticity (dv/dx - du/dy) in q', 'm s-1',conversion=US%L_T_to_m_s)

  CS%id_sh_xy_h = register_diag_field('ocean_model','sh_xy_h', diag%axesTL, Time, &
      'Horizontal shearing strain (du/dy + dv/dx) in q', 'm s-1',conversion=US%L_T_to_m_s)

  CS%id_vort_xy_h = register_diag_field('ocean_model','vort_xy_h', diag%axesTL, Time, &
      'Vertical vorticity (dv/dx - du/dy) in q', 'm s-1',conversion=US%L_T_to_m_s)

  CS%id_eta = register_diag_field('ocean_model','eta',diag%axesTi, Time, &
      'layer interfaces as computed from the find_eta submodule', 'm', conversion=US%Z_to_L)

  CS%id_h_eta = register_diag_field('ocean_model','h_eta',diag%axesTL, Time, &
      'layer thickness when input to the find_eta function', 'm', conversion=US%Z_to_L)

  CS%id_depth_mask = register_diag_field('ocean_model','depth_mask',diag%axesTL, Time, &
      'Mask based off of a cutoff depth for application to ANN input', 'nondim')

  CS%id_Coriolis_h = register_diag_field('ocean_model', 'Coriolis_h', diag%axesTL, Time, &
       'Coriolis parameter', 's-1', conversion=US%T_to_s)

  CS%id_Delsq_h = register_diag_field('ocean_model', 'Delsq_h', diag%axesTL, Time, &
       'Squared Harmonic mean of grid spacing', 'm2', conversion=US%L_to_m**2)

  CS%id_mom_norm_h = register_diag_field('ocean_model', 'mom_norm_h', diag%axesTL, Time, &
       'Norm of velocity gradients (from input features)', 's-1', conversion=US%T_to_s)

  CS%id_buoy_norm_h = register_diag_field('ocean_model', 'buoy_norm_h', diag%axesTL, Time, &
       'Norm of interface slopes (from input features)', 'nondim')
  
  CS%id_f_u = register_diag_field('ocean_model','f_u', diag%axesCuL, Time, &
       'saved diffu variable that is updated by divergence of EPF flux', 'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_f_v = register_diag_field('ocean_model','f_v', diag%axesCvL, Time, &
       'saved diffv variable that is updated by divergence of EPF flux', 'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_f_x = register_diag_field('ocean_model','f_x', diag%axesCuL, Time, &
       'horizontal divergence of zonal Reynolds stress', &
       'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_f_y = register_diag_field('ocean_model','f_y', diag%axesCvL, Time, &
       'horizontal divergence of meridional Reynolds stress', &
       'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_fx_PE = register_diag_field('ocean_model','fx_PE', diag%axesCuL, Time, &
       'horizontal divergence of zonal EPF due to subfilter potential energy', 'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_fxz = register_diag_field('ocean_model','fxz', diag%axesCuL, Time, &
       'vertical divergence of zonal EPF', 'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  ! CS%id_fyx = register_diag_field('ocean_model','fyx', diag%axesCvL, Time, &
  !      'horizontal divergence of meridional irrotational component of Reynolds stress', &
  !      'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  ! CS%id_fyy = register_diag_field('ocean_model','fyy', diag%axesCvL, Time, &
  !      'horizontal divergence of meridional rotational component of Reynolds stress', &
  !      'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_fy_PE = register_diag_field('ocean_model','fy_PE', diag%axesCvL, Time, &
       'horizontal divergence of meridional EPF due to subfilter potential energy', 'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  CS%id_fyz = register_diag_field('ocean_model','fyz', diag%axesCvL, Time, &
       'vertical divergence of meridional EPF', 'm s-2', conversion=(US%L_T_to_m_s**2)/US%L_to_m)

  ! Clock IDs
  ! Only module is measured with syncronization. While smaller
  ! parts are measured without - because these are nested clocks.
  CS%id_clock_module = cpu_clock_id('(Ocean EPF param)', grain=CLOCK_MODULE)
  CS%id_clock_copy = cpu_clock_id('(EPF copy fields)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_stress = cpu_clock_id('(EPF compute stress)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_stress_ANN = cpu_clock_id('(EPF compute stress ANN)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_divergence = cpu_clock_id('(EPF compute divergence)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_mpi = cpu_clock_id('(EPF filter MPI exchanges)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_post = cpu_clock_id('(EPF post data)', grain=CLOCK_ROUTINE, sync=.false.)
  CS%id_clock_source = cpu_clock_id('(EPF compute energy source)', grain=CLOCK_ROUTINE, sync=.false.)

  CS%subroundoff_shear = 1e-30 * US%T_to_s
  
  

  ! Allocate memory
  ! We set the stress tensor and velocity gradient tensor to zero
  ! with full halo because they potentially may be filtered
  ! with marching halo algorithm
  allocate(CS%sh_xx(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%sh_xy(SZIB_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%vort_xy(SZIB_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%hq(SZIB_(G),SZJB_(G),SZK_(GV)))
  allocate(CS%sh_xy_h(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%vort_xy_h(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)

  allocate(CS%eta(SZI_(G), SZJ_(G), SZK_(GV)+1))
  allocate(CS%h_eta(SZI_(G), SZJ_(G), SZK_(GV)))
  allocate(CS%depth_mask(SZI_(G), SZJ_(G), SZK_(GV)), source=1.)
  
  allocate(CS%slope_x(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke+1), source=0.)
  allocate(CS%slope_y(G%isd:G%ied,G%JsdB:G%JedB,GV%ke+1), source=0.)

  allocate(CS%slope_x_top(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%slope_x_bot(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%slope_y_top(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%slope_y_bot(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)

  allocate(CS%Txx(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Tyy(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Txx_PE(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Tyy_PE(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Txy(SZIB_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%Txy_h(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Txz(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%Tyz(SZI_(G),SZJ_(G),SZK_(GV)), source=0.)
  
  allocate(CS%kappa_h(SZI_(G),SZJ_(G)))
  allocate(CS%kappa_q(SZIB_(G),SZJB_(G)))
  
  allocate(CS%Delsq_h(SZI_(G),SZJ_(G)))
  allocate(CS%Coriolis_h(SZI_(G),SZJ_(G)))

  allocate(CS%mom_norm_h(SZI_(G),SZJ_(G),SZK_(GV)))
  allocate(CS%buoy_norm_h(SZI_(G),SZJ_(G),SZK_(GV)))

  allocate(CS%f_u(SZIB_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%f_x(SZIB_(G),SZJ_(G),SZK_(GV)), source=0.)
  ! allocate(CS%fxy(SZIB_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%fx_PE(SZIB_(G),SZJ_(G),SZK_(GV)), source=0.)
  allocate(CS%fxz(SZIB_(G),SZJ_(G),SZK_(GV)), source=0.)

  allocate(CS%f_v(SZI_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%f_y(SZI_(G),SZJB_(G),SZK_(GV)), source=0.)
  ! allocate(CS%fyy(SZI_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%fy_PE(SZI_(G),SZJB_(G),SZK_(GV)), source=0.)
  allocate(CS%fyz(SZI_(G),SZJB_(G),SZK_(GV)), source=0.)

  subroundoff_Cor = 1e-30 * US%T_to_s
  ! Precomputing f
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%Coriolis_h(i,j) = (abs(0.25 * ((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) &
                          + (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1)))) + subroundoff_Cor) 
  enddo; enddo


  do j=js-2, je+2 ; do i=is-2, ie+2
    CS%Delsq_h(i,j) = 2 * (((G%dxT(i,j)**2) * (G%dyT(i,j)**2)) / ((G%dxT(i,j)**2) + (G%dyT(i,j)**2)))
  enddo ;enddo

  ! Precomputing the scaling coefficient
  ! Mask is included to automatically satisfy B.C.
  do j=js-2,je+2 ; do i=is-2,ie+2
    CS%kappa_h(i,j) = -CS%amplitude * G%areaT(i,j) * G%mask2dT(i,j)
  enddo; enddo

  do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
    CS%kappa_q(I,J) = -CS%amplitude * G%areaBu(I,J) * G%mask2dBu(I,J)
  enddo; enddo


end subroutine EPF_init

!> Save precomputed velocity gradients and thickness
!! from the horizontal eddy viscosity module
!! We save as much halo for velocity gradients as possible
!! In symmetric (preferable) memory model: halo 2 for sh_xx
!! and halo 1 for sh_xy and vort_xy
!! We apply zero boundary conditions to velocity gradients
!! which is required for filtering operations
subroutine EPF_copy_gradient_and_thickness(                        &
                                      sh_xx, sh_xy, vort_xy,  &
                                      hq,                             &    
                                      G, GV, CS, k)
  type(ocean_grid_type),         intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(EPF_CS),                  intent(inout) :: CS     !< EPF control structure.

  real, dimension(SZIB_(G),SZJB_(G)), &
    intent(in) :: sh_xy       !< horizontal shearing strain (du/dy + dv/dx)
                              !! including metric terms [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
    intent(in) :: vort_xy     !< Vertical vorticity (dv/dx - du/dy)
                              !! including metric terms [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
    intent(in) :: hq          !< harmonic mean of the harmonic means
                              !! of the u- & v point thicknesses [H ~> m or kg m-2]

  real, dimension(SZI_(G),SZJ_(G)), &
    intent(in) :: sh_xx       !< horizontal tension (du/dx - dv/dy)
                              !! including metric terms [T-1 ~> s-1]

  integer, intent(in) :: k    !< The vertical index of the layer to be passed.

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: i, j

  call cpu_clock_begin(CS%id_clock_copy)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do J=js-1,Jeq ; do I=is-1,Ieq
    CS%hq(I,J,k) = hq(I,J)
  enddo; enddo

  ! No physical B.C. is required for
  ! sh_xx in ZB2020. However, filtering
  ! may require BC
  do j=Jsq-1,je+2 ; do i=Isq-1,ie+2
    CS%sh_xx(i,j,k) = sh_xx(i,j) * G%mask2dT(i,j)
  enddo ; enddo

  ! We multiply by mask to remove
  ! implicit dependence on CS%no_slip
  ! flag in hor_visc module
  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    CS%sh_xy(I,J,k) = sh_xy(I,J) * G%mask2dBu(I,J)
  enddo; enddo

  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    CS%vort_xy(I,J,k) = vort_xy(I,J) * G%mask2dBu(I,J)
  enddo; enddo

  call cpu_clock_end(CS%id_clock_copy)

end subroutine EPF_copy_gradient_and_thickness

!subroutine slope_calc(e, G, GV, slope_x, slope_y)
subroutine slope_calc(G, GV, CS)
  type(ocean_grid_type),                         intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                       intent(in)     :: GV  !< The ocean's vertical grid structure.
  type(EPF_CS),                                  intent(inout)  :: CS   !< EPS control structure.
  
  ! local variables 
  integer i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! loop over u points 
  do I = Isq-1, Ieq+1
    do j = js-2, je+2 
      do k = 1, nz+1
      ! Note the negative sign in front of difference! This is because of sign error (see MOM_lateral_mixing_coeffs.F90 lines 895-904
        CS%slope_x(I,j,k) = (-(CS%eta(i,j,k)-CS%eta(i+1,j,k))*G%IdxCu(I,j)) * G%OBCmaskCu(I,j) 
      enddo

    enddo
  enddo

  ! loop over v points 
  do i = is-2, ie+2
    do J = Jsq-1, Jeq+1 
      do k = 1, nz+1
        CS%slope_y(i,J,k)= (-(CS%eta(i,j,k)-CS%eta(i,j+1,k))*G%IdyCv(i,J)) * G%OBCmaskCv(i,J)
      enddo

    enddo
  enddo
  
end subroutine slope_calc

!> Compute stress tensor T =
!! (Txx, Txy;
!!  Txy, Tyy;
!!  Txz, Tyz)
!!  with ANN 

subroutine make_mask(h, G, GV, CS)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(EPF_CS),              intent(inout) :: CS   !< EPS control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2]

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  ! real, dimension(SZI_(G), SZJ_(G), SZK_(GV)), intent(inout) :: depth_mask

  do k = 1, nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      if (h(i,j,k) < CS%H_cutoff) then
        CS%depth_mask(i,j,k) = 0.0
      else
        CS%depth_mask(i,j,k) = 1.0
      endif 
    enddo ; enddo
  enddo
end subroutine make_mask

subroutine compute_stress_ANN_collocated(h, tv, G, GV, CS, US, VarMix, OBC)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(EPF_CS),              intent(inout) :: CS   !< EPS control structure.
  type(unit_scale_type),     intent(in)    :: US    !< A dimensional unit scaling type
  type(VarMix_CS), target,   intent(inout)    :: VarMix !< Variable mixing coefficients
  type(thermo_var_ptrs),     intent(in)    :: tv  !<thermodynamics structure
  type(ocean_OBC_type),      pointer       :: OBC !< Open boundaries control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2]
  !real,                                      intent(in)    :: dt !< Time increment [T ~> s]
   
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  real :: x(CS%n_inputs), y(CS%n_outputs)
  real :: input_norm
  real :: input_norm_mom
  real :: input_norm_buoy

  logical :: &
        use_VarMix, &
        use_stored_slopes
  
  ! real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
  !       mom_norm_h, &    ! Norm in h points for momentum inputs [T-1 ~ s-1]
  !       buoy_norm_h      ! Norm in h points for interface inputs [T-1 ~ s-1]
        
  real, dimension(SZIB_(G), SZJ_(G), SZK_(GV)+1) :: &
        slope_x, &
        slope_y
        
  real, dimension(SZI_(G),SZJ_(G)) :: &
        sqr_h, & ! Sum of squares in h points
        sqr_eta_h, &
        Txy      ! Predicted Txy in center points to be interpolated to corners

  call cpu_clock_begin(CS%id_clock_stress_ANN)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! mom_norm_h = 0.
  ! buoy_norm_h = 0.

  call pass_var(CS%sh_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%sh_xx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%vort_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  
  call find_eta(h, tv, G, GV, US, CS%eta, halo_size=1) ! calculating the interface height
  call pass_var(CS%eta, G%Domain, clock=CS%id_clock_mpi)
  CS%h_eta = h
  call slope_calc(G, GV, CS) ! calculating the interface slopes instead of relying on VarMix
  call pass_vector(CS%slope_x, CS%slope_y, G%Domain, clock=CS%id_clock_mpi)
  
  call make_mask(h, G, GV, CS) !creating the mask based off of depth
  ! Interpolate input features
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      ! It is assumed that B.C. is applied to sh_xy and vort_xy
      CS%sh_xy_h(i,j,k) = 0.25 * ( (CS%sh_xy(I-1,J-1,k) + CS%sh_xy(I,J,k)) &
                       + (CS%sh_xy(I-1,J,k) + CS%sh_xy(I,J-1,k)) ) * G%mask2dT(i,j)

      CS%vort_xy_h(i,j,k) = 0.25 * ( (CS%vort_xy(I-1,J-1,k) + CS%vort_xy(I,J,k)) &
                         + (CS%vort_xy(I-1,J,k) + CS%vort_xy(I,J-1,k)) ) * G%mask2dT(i,j)
      
      CS%slope_x_top(i,j,k) = 0.5 * (( (CS%slope_x(I-1,j,k) + CS%slope_x(I,j,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j) 
      CS%slope_y_top(i,j,k) = 0.5 * (( (CS%slope_y(i,J-1,k) + CS%slope_y(i,J,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j)
      CS%slope_x_bot(i,j,k) = 0.5 * ( (CS%slope_x(I-1,j,k+1) + CS%slope_x(I,j,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)
      CS%slope_y_bot(i,j,k) = 0.5 * ( (CS%slope_y(i,J-1,k+1) + CS%slope_y(i,J,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)


      sqr_eta_h(i,j) = (CS%slope_x_top(i,j,k)**2) + (CS%slope_y_top(i,j,k)**2) + (CS%slope_x_bot(i,j,k)**2) + (CS%slope_y_bot(i,j,k)**2)
      sqr_h(i,j) = (CS%sh_xx(i,j,k)**2) + (CS%sh_xy_h(i,j,k)**2) + (CS%vort_xy_h(i,j,k)**2)
      
      CS%mom_norm_h(i,j,k) = sqrt(sqr_h(i,j))
      CS%buoy_norm_h(i,j,k) = sqrt(sqr_eta_h(i,j))
      
    enddo; enddo
  enddo

  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      x(1) = CS%vort_xy_h(i,j,k)
      x(2) = CS%sh_xx(i,j,k)
      x(3) = CS%sh_xy_h(i,j,k)

      x(4) = CS%slope_x_top(i,j,k)
      x(5) = CS%slope_y_top(i,j,k)
      x(6) = CS%slope_x_bot(i,j,k)
      x(7) = CS%slope_y_bot(i,j,k)

      input_norm_mom = CS%mom_norm_h(i,j,k)
      input_norm_buoy = CS%buoy_norm_h(i,j,k)

      x(1:3) = x(1:3) / (input_norm_mom + CS%subroundoff_shear) !normalising momentum inputs
      x(4:7) = x(4:7) / (input_norm_buoy + CS%subroundoff_shear) !normalising buoyancy inputs

      call ANN_apply(x, y, CS%ann_instance)

      y(1:3) = y(1:3) * CS%Delsq_h(i,j) * input_norm_mom * input_norm_mom
      y(4:5) = y(4:5) * CS%Delsq_h(i,j) * CS%Coriolis_h(i,j) * input_norm_mom * input_norm_buoy

      CS%Txx(i,j,k) = y(1)
      Txy(i,j)      = y(2)
      CS%Tyy(i,j,k) = y(3)
      CS%Txz(i,j,k) = y(4)
      CS%Tyz(i,j,k) = y(5)
      CS%Txy_h(i,j,k) = y(2)
    enddo ; enddo

    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      CS%Txy(I,J,k) = 0.25 * ( (Txy(i+1,j+1) + Txy(i,j)) &
                             + (Txy(i+1,j)   + Txy(i,j+1))) * G%mask2dBu(I,J)
    enddo; enddo

  enddo ! end of k loop

  call pass_var(CS%Txy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%Txx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Tyy, G%Domain, clock=CS%id_clock_mpi)

  call cpu_clock_end(CS%id_clock_stress_ANN)

end subroutine compute_stress_ANN_collocated

subroutine compute_stress_separate_ANN(h, tv, G, GV, CS, US, VarMix, OBC)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(EPF_CS),              intent(inout) :: CS   !< EPS control structure.
  type(unit_scale_type),     intent(in)    :: US    !< A dimensional unit scaling type
  type(VarMix_CS), target,   intent(inout)    :: VarMix !< Variable mixing coefficients
  type(thermo_var_ptrs),     intent(in)    :: tv  !<thermodynamics structure
  type(ocean_OBC_type),      pointer       :: OBC !< Open boundaries control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2]
  !real,                                      intent(in)    :: dt !< Time increment [T ~> s]
   
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  real :: x_intfc(CS%n_inputs_intfc), y_intfc(CS%n_outputs_intfc)
  real :: x_centre(CS%n_inputs_centre), y_centre(CS%n_outputs_centre)
  
  real :: input_norm_mom
  real :: input_norm_buoy

  logical :: &
        use_VarMix, &
        use_stored_slopes
  
  ! real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
  !       mom_norm_h, &    ! Norm in h points for momentum inputs [T-1 ~ s-1]
  !       buoy_norm_h      ! Norm in h points for interface inputs [T-1 ~ s-1]
        
  real, dimension(SZIB_(G), SZJ_(G), SZK_(GV)+1) :: &
        slope_x, &
        slope_y
  
  real, dimension(SZIB_(G), SZJ_(G), SZK_(GV)-1) :: & ! to keep consistent with Txz and Tyz being the differences between top and bottom stress, we have to intorduce this intermediate variable
        intfc_stress_x, &
        intfc_stress_y

  real, dimension(SZI_(G),SZJ_(G)) :: &
        sqr_h, & ! Sum of squares in h points
        sqr_eta_h, &
        Txy      ! Predicted Txy in center points to be interpolated to corners

  call cpu_clock_begin(CS%id_clock_stress_ANN)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! mom_norm_h = 0.
  ! buoy_norm_h = 0.

  call pass_var(CS%sh_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%sh_xx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%vort_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  
  call find_eta(h, tv, G, GV, US, CS%eta, halo_size=1) ! calculating the interface height
  call pass_var(CS%eta, G%Domain, clock=CS%id_clock_mpi)
  CS%h_eta = h
  call slope_calc(G, GV, CS) ! calculating the interface slopes instead of relying on VarMix
  call pass_vector(CS%slope_x, CS%slope_y, G%Domain, clock=CS%id_clock_mpi)
  
  call make_mask(h, G, GV, CS) !creating the mask based off of depth
  ! Interpolate input features
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      ! It is assumed that B.C. is applied to sh_xy and vort_xy
      CS%sh_xy_h(i,j,k) = 0.25 * ( (CS%sh_xy(I-1,J-1,k) + CS%sh_xy(I,J,k)) &
                       + (CS%sh_xy(I-1,J,k) + CS%sh_xy(I,J-1,k)) ) * G%mask2dT(i,j)

      CS%vort_xy_h(i,j,k) = 0.25 * ( (CS%vort_xy(I-1,J-1,k) + CS%vort_xy(I,J,k)) &
                         + (CS%vort_xy(I-1,J,k) + CS%vort_xy(I,J-1,k)) ) * G%mask2dT(i,j)
      
      CS%slope_x_top(i,j,k) = 0.5 * (( (CS%slope_x(I-1,j,k) + CS%slope_x(I,j,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j) 
      CS%slope_y_top(i,j,k) = 0.5 * (( (CS%slope_y(i,J-1,k) + CS%slope_y(i,J,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j)
      CS%slope_x_bot(i,j,k) = 0.5 * ( (CS%slope_x(I-1,j,k+1) + CS%slope_x(I,j,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)
      CS%slope_y_bot(i,j,k) = 0.5 * ( (CS%slope_y(i,J-1,k+1) + CS%slope_y(i,J,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)


      sqr_eta_h(i,j) = (CS%slope_x_top(i,j,k)**2) + (CS%slope_y_top(i,j,k)**2)
      sqr_h(i,j) = (CS%sh_xx(i,j,k)**2) + (CS%sh_xy_h(i,j,k)**2) + (CS%vort_xy_h(i,j,k)**2)
      
      CS%mom_norm_h(i,j,k) = sqrt(sqr_h(i,j))
      CS%buoy_norm_h(i,j,k) = sqrt(sqr_eta_h(i,j))
      
    enddo; enddo
  enddo


  ! a loop which calculates the interface flux(es)
  do k=1, nz-1 ! for two layers, we will only calculate one interface flux... the lower interface
    do j=js-2,je+2 ; do i=is-2,ie+2
      x_intfc(1) = CS%vort_xy_h(i,j,k)
      x_intfc(2) = CS%sh_xx(i,j,k)
      x_intfc(3) = CS%sh_xy_h(i,j,k)
      x_intfc(4) = CS%slope_x_bot(i,j,k)
      x_intfc(5) = CS%slope_y_bot(i,j,k)

      input_norm_mom = CS%mom_norm_h(i,j,k)
      input_norm_buoy = CS%buoy_norm_h(i,j,k+1) !the interface norm is based off of top interface slopes, we we want norm from lower layer

      x_intfc(1:3) = x_intfc(1:3) / (input_norm_mom + CS%subroundoff_shear) !normalising momentum inputs
      x_intfc(4:5) = x_intfc(4:5) / (input_norm_buoy + CS%subroundoff_shear) !normalising momentum inputs

      call ANN_apply(x_intfc, y_intfc, CS%ann_instance_intfc)
      y_intfc = y_intfc * CS%Delsq_h(i,j) * CS%Coriolis_h(i,j) * input_norm_mom * input_norm_buoy
      intfc_stress_x(i,j,k) = y_intfc(1)
      intfc_stress_y(i,j,k) = y_intfc(2)
    enddo; enddo
  enddo

  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      x_centre(1) = CS%vort_xy_h(i,j,k)
      x_centre(2) = CS%sh_xx(i,j,k)
      x_centre(3) = CS%sh_xy_h(i,j,k)

      input_norm_mom = CS%mom_norm_h(i,j,k)

      x_centre(1:3) = x_centre(1:3) / (input_norm_mom + CS%subroundoff_shear) !normalising momentum inputs
      
      call ANN_apply(x_centre, y_centre, CS%ann_instance_centre)

      y_centre(1:3) = y_centre(1:3) * CS%Delsq_h(i,j) * input_norm_mom * input_norm_mom
      
      CS%Txx(i,j,k) = y_centre(1)
      Txy(i,j)      = y_centre(2)
      CS%Tyy(i,j,k) = y_centre(3)
      CS%Txy_h(i,j,k) = y_centre(2)

      !there are only nz-1 interface stresses due to zero boundary conditions
      if (k==1) then !zero boundary condition minus the interface stress
        CS%Txz(i,j,k) = -intfc_stress_x(i,j,k) 
        CS%Tyz(i,j,k) = -intfc_stress_y(i,j,k)
      elseif (k==nz) then !interface stress minus the zero boundary condition
        CS%Txz(i,j,k) = intfc_stress_x(i,j,k-1) 
        CS%Tyz(i,j,k) = intfc_stress_y(i,j,k-1)
      else !if there are more than 2 layers, then there is more than 1 calculated interface flux, and we need to take difference between these
        CS%Txz(i,j,k) = intfc_stress_x(i,j,k-1) - intfc_stress_x(i,j,k)
        CS%Tyz(i,j,k) = intfc_stress_y(i,j,k-1) - intfc_stress_y(i,j,k)
      end if
    enddo ; enddo
    
    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      CS%Txy(I,J,k) = 0.25 * ( (Txy(i+1,j+1) + Txy(i,j)) &
                             + (Txy(i+1,j)   + Txy(i,j+1))) * G%mask2dBu(I,J)
    enddo; enddo

  enddo ! end of k loop

  call pass_var(CS%Txy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%Txx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Tyy, G%Domain, clock=CS%id_clock_mpi)

  call cpu_clock_end(CS%id_clock_stress_ANN)

end subroutine compute_stress_separate_ANN

subroutine compute_stress_ANN_collocated_column(h, tv, G, GV, CS, US, VarMix, OBC)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(EPF_CS),              intent(inout) :: CS   !< EPS control structure.
  type(unit_scale_type),     intent(in)    :: US    !< A dimensional unit scaling type
  type(VarMix_CS), target,   intent(inout)    :: VarMix !< Variable mixing coefficients
  type(thermo_var_ptrs),     intent(in)    :: tv  !<thermodynamics structure
  type(ocean_OBC_type),      pointer       :: OBC !< Open boundaries control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2]
  !real,                                      intent(in)    :: dt !< Time increment [T ~> s]
   
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  real :: x(CS%n_inputs), y(CS%n_outputs)
  real :: input_norm
  real :: input_norm_mom
  real :: input_norm_buoy

  logical :: &
        use_VarMix, &
        use_stored_slopes
  
  ! real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
  !       mom_norm_h, &    ! Norm in h points for momentum inputs [T-1 ~ s-1]
  !       buoy_norm_h      ! Norm in h points for interface inputs [T-1 ~ s-1]
        
  real, dimension(SZIB_(G), SZJ_(G), SZK_(GV)+1) :: &
        slope_x, &
        slope_y
        
  real, dimension(SZI_(G),SZJ_(G)) :: &
        sqr_h, & ! Sum of squares in h points
        sqr_eta_h
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
        Txy      ! Predicted Txy in center points to be interpolated to corners

  call cpu_clock_begin(CS%id_clock_stress_ANN)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! mom_norm_h = 0.
  ! buoy_norm_h = 0.

  call pass_var(CS%sh_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%sh_xx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%vort_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  
  call find_eta(h, tv, G, GV, US, CS%eta, halo_size=1) ! calculating the interface height
  call pass_var(CS%eta, G%Domain, clock=CS%id_clock_mpi)
  CS%h_eta = h
  call slope_calc(G, GV, CS) ! calculating the interface slopes instead of relying on VarMix
  call pass_vector(CS%slope_x, CS%slope_y, G%Domain, clock=CS%id_clock_mpi)
  
  call make_mask(h, G, GV, CS) !creating the mask based off of depth
  ! Interpolate input features
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      ! It is assumed that B.C. is applied to sh_xy and vort_xy
      CS%sh_xy_h(i,j,k) = 0.25 * ( (CS%sh_xy(I-1,J-1,k) + CS%sh_xy(I,J,k)) &
                       + (CS%sh_xy(I-1,J,k) + CS%sh_xy(I,J-1,k)) ) * G%mask2dT(i,j)

      CS%vort_xy_h(i,j,k) = 0.25 * ( (CS%vort_xy(I-1,J-1,k) + CS%vort_xy(I,J,k)) &
                         + (CS%vort_xy(I-1,J,k) + CS%vort_xy(I,J-1,k)) ) * G%mask2dT(i,j)
      
      CS%slope_x_top(i,j,k) = 0.5 * (( (CS%slope_x(I-1,j,k) + CS%slope_x(I,j,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j) 
      CS%slope_y_top(i,j,k) = 0.5 * (( (CS%slope_y(i,J-1,k) + CS%slope_y(i,J,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j)
      CS%slope_x_bot(i,j,k) = 0.5 * ( (CS%slope_x(I-1,j,k+1) + CS%slope_x(I,j,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)
      CS%slope_y_bot(i,j,k) = 0.5 * ( (CS%slope_y(i,J-1,k+1) + CS%slope_y(i,J,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)


      sqr_eta_h(i,j) = (CS%slope_x_top(i,j,k)**2) + (CS%slope_y_top(i,j,k)**2) 
      sqr_h(i,j) = (CS%sh_xx(i,j,k)**2) + (CS%sh_xy_h(i,j,k)**2) + (CS%vort_xy_h(i,j,k)**2)
      
      CS%mom_norm_h(i,j,k) = sqrt(sqr_h(i,j))
      CS%buoy_norm_h(i,j,k) = sqrt(sqr_eta_h(i,j))
      
    enddo; enddo
  enddo


  do j=js-2,je+2 ; do i=is-2,ie+2
    x(1) = CS%vort_xy_h(i,j,1) / (CS%mom_norm_h(i,j,1) + CS%subroundoff_shear)
    x(2) = CS%sh_xx(i,j,1) / (CS%mom_norm_h(i,j,1) + CS%subroundoff_shear)
    x(3) = CS%sh_xy_h(i,j,1) / (CS%mom_norm_h(i,j,1) + CS%subroundoff_shear)
    x(4) = CS%vort_xy_h(i,j,2) / (CS%mom_norm_h(i,j,2) + CS%subroundoff_shear)
    x(5) = CS%sh_xx(i,j,2) / (CS%mom_norm_h(i,j,2) + CS%subroundoff_shear)
    x(6) = CS%sh_xy_h(i,j,2) / (CS%mom_norm_h(i,j,2) + CS%subroundoff_shear)

    x(7) = CS%slope_x_top(i,j,1) / (CS%buoy_norm_h(i,j,1) + CS%subroundoff_shear)
    x(8) = CS%slope_y_top(i,j,1) / (CS%buoy_norm_h(i,j,1) + CS%subroundoff_shear)
    x(9) = CS%slope_x_top(i,j,2) / (CS%buoy_norm_h(i,j,2) + CS%subroundoff_shear)
    x(10) = CS%slope_y_top(i,j,2) / (CS%buoy_norm_h(i,j,2) + CS%subroundoff_shear)

    call ANN_apply(x, y, CS%ann_instance)

    y(1:3) = y(1:3) * CS%Delsq_h(i,j) * CS%mom_norm_h(i,j,1) * CS%mom_norm_h(i,j,1) !upper layer Reynolds stresses, uu, uv, vv
    y(4:6) = y(4:6) * CS%Delsq_h(i,j) * CS%mom_norm_h(i,j,2) * CS%mom_norm_h(i,j,2) !lower layer Reynolds stresses, uu, uv, vv
    y(7:8) = y(7:8) * CS%Delsq_h(i,j) * CS%Coriolis_h(i,j) * CS%mom_norm_h(i,j,1) * CS%buoy_norm_h(i,j,2) !dual form stresses at interface
    y(9) = y(9) * CS%Delsq_h(i,j) * CS%GFS * CS%buoy_norm_h(i,j,1) * CS%buoy_norm_h(i,j,1) !subfilter PE in upper layer
    y(10) = y(10) * CS%Delsq_h(i,j) * CS%GINT * CS%buoy_norm_h(i,j,2) * CS%buoy_norm_h(i,j,2) !subfilter PE in lower layer


    CS%Txx(i,j,1) = y(1)
    CS%Txx_PE(i,j,1) = (y(9) / h(i,j,1))
    !CS%Txx(i,j,1) = y(1) + (y(9) / h(i,j,1))
    Txy(i,j,1) = y(2) 
    !CS%Tyy(i,j,1) = y(3) + (y(9) / h(i,j,1))
    CS%Tyy(i,j,1) = y(3)
    CS%Tyy_PE(i,j,1) = (y(9) / h(i,j,1))
    !CS%Txx(i,j,2) = y(4) + (y(10) / h(i,j,2))
    CS%Txx(i,j,2) = y(4)
    CS%Txx_PE(i,j,2) = (y(10) / h(i,j,2))
    Txy(i,j,2) = y(5)
    !CS%Tyy(i,j,2) = y(6) + (y(10) / h(i,j,2))
    CS%Tyy(i,j,2) = y(6)
    CS%Tyy_PE(i,j,2) = (y(10) / h(i,j,2))
    CS%Txz(i,j,1) = -y(7) !zero top BC - dual form stress
    CS%Txz(i,j,2) = y(7) !dual form stress - zero bottom BC
    CS%Tyz(i,j,1) = -y(8) !zero top BC - dual form stress
    CS%Tyz(i,j,2) = y(8) !dual form stress - zero bottom BC
  enddo; enddo

  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    CS%Txy(I,J,1) = 0.25 * ( (Txy(i+1,j+1,1) + Txy(i,j,1)) &
                           + (Txy(i+1,j,1)   + Txy(i,j+1,1))) * G%mask2dBu(I,J)
    CS%Txy(I,J,2) = 0.25 * ( (Txy(i+1,j+1,2) + Txy(i,j,2)) &
                           + (Txy(i+1,j,2)   + Txy(i,j+1,2))) * G%mask2dBu(I,J)
  enddo; enddo
  
  call pass_var(CS%Txy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%Txx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Tyy, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Txx_PE, G%Domain, clock = CS%id_clock_mpi)
  call pass_var(CS%Tyy_PE, G%Domain, clock = CS%id_clock_mpi)

  call cpu_clock_end(CS%id_clock_stress_ANN)

end subroutine compute_stress_ANN_collocated_column

!> Compute the divergence of subgrid stress
!! weighted with thickness, i.e.
!! (fx,fy) = 1/h Div(h * [Txx, Txy; Txy, Tyy; Txz, Tyz])
!! and update the acceleration due to eddy viscosity as
!! diffu = diffu + dx; diffv = diffv + dy
!! Optionally, before computing the divergence, we attenuate the stress
!! according to the Klower formula.


subroutine compute_stress_ANN_collocated_column(h, tv, G, GV, CS, US, VarMix, OBC)
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure
  type(EPF_CS),              intent(inout) :: CS   !< EPS control structure.
  type(unit_scale_type),     intent(in)    :: US    !< A dimensional unit scaling type
  type(VarMix_CS), target,   intent(inout)    :: VarMix !< Variable mixing coefficients
  type(thermo_var_ptrs),     intent(in)    :: tv  !<thermodynamics structure
  type(ocean_OBC_type),      pointer       :: OBC !< Open boundaries control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thickness [H ~> m or kg m-2]
  !real,                                      intent(in)    :: dt !< Time increment [T ~> s]
   
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  real :: x(CS%n_inputs), y(CS%n_outputs)
  real :: input_norm
  real :: input_norm_mom
  real :: input_norm_buoy

  logical :: &
        use_VarMix, &
        use_stored_slopes
  
  ! real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
  !       mom_norm_h, &    ! Norm in h points for momentum inputs [T-1 ~ s-1]
  !       buoy_norm_h      ! Norm in h points for interface inputs [T-1 ~ s-1]
        
  real, dimension(SZIB_(G), SZJ_(G), SZK_(GV)+1) :: &
        slope_x, &
        slope_y
        
  real, dimension(SZI_(G),SZJ_(G)) :: &
        sqr_h, & ! Sum of squares in h points
        sqr_eta_h
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
        Txy      ! Predicted Txy in center points to be interpolated to corners

  call cpu_clock_begin(CS%id_clock_stress_ANN)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! mom_norm_h = 0.
  ! buoy_norm_h = 0.

  call pass_var(CS%sh_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%sh_xx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%vort_xy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  
  call find_eta(h, tv, G, GV, US, CS%eta, halo_size=1) ! calculating the interface height
  call pass_var(CS%eta, G%Domain, clock=CS%id_clock_mpi)
  CS%h_eta = h
  call slope_calc(G, GV, CS) ! calculating the interface slopes instead of relying on VarMix
  call pass_vector(CS%slope_x, CS%slope_y, G%Domain, clock=CS%id_clock_mpi)
  
  call make_mask(h, G, GV, CS) !creating the mask based off of depth
  ! Interpolate input features
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      ! It is assumed that B.C. is applied to sh_xy and vort_xy
      CS%sh_xy_h(i,j,k) = 0.25 * ( (CS%sh_xy(I-1,J-1,k) + CS%sh_xy(I,J,k)) &
                       + (CS%sh_xy(I-1,J,k) + CS%sh_xy(I,J-1,k)) ) * G%mask2dT(i,j)

      CS%vort_xy_h(i,j,k) = 0.25 * ( (CS%vort_xy(I-1,J-1,k) + CS%vort_xy(I,J,k)) &
                         + (CS%vort_xy(I-1,J,k) + CS%vort_xy(I,J-1,k)) ) * G%mask2dT(i,j)
      
      CS%slope_x_top(i,j,k) = 0.5 * (( (CS%slope_x(I-1,j,k) + CS%slope_x(I,j,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j) 
      CS%slope_y_top(i,j,k) = 0.5 * (( (CS%slope_y(i,J-1,k) + CS%slope_y(i,J,k)))* CS%depth_mask(i,j,k)) * G%mask2dT(i,j)
      CS%slope_x_bot(i,j,k) = 0.5 * ( (CS%slope_x(I-1,j,k+1) + CS%slope_x(I,j,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)
      CS%slope_y_bot(i,j,k) = 0.5 * ( (CS%slope_y(i,J-1,k+1) + CS%slope_y(i,J,k+1)) ) * G%mask2dT(i,j) * CS%depth_mask(i,j,k)


      sqr_eta_h(i,j) = (CS%slope_x_top(i,j,k)**2) + (CS%slope_y_top(i,j,k)**2) 
      sqr_h(i,j) = (CS%sh_xx(i,j,k)**2) + (CS%sh_xy_h(i,j,k)**2) + (CS%vort_xy_h(i,j,k)**2)
      
      CS%mom_norm_h(i,j,k) = sqrt(sqr_h(i,j))
      CS%buoy_norm_h(i,j,k) = sqrt(sqr_eta_h(i,j))
      
    enddo; enddo
  enddo


  do j=js-2,je+2 ; do i=is-2,ie+2
    x(1) = CS%vort_xy_h(i,j,1) / (CS%mom_norm_h(i,j,1) + CS%subroundoff_shear)
    x(2) = CS%sh_xx(i,j,1) / (CS%mom_norm_h(i,j,1) + CS%subroundoff_shear)
    x(3) = CS%sh_xy_h(i,j,1) / (CS%mom_norm_h(i,j,1) + CS%subroundoff_shear)
    x(4) = CS%vort_xy_h(i,j,2) / (CS%mom_norm_h(i,j,2) + CS%subroundoff_shear)
    x(5) = CS%sh_xx(i,j,2) / (CS%mom_norm_h(i,j,2) + CS%subroundoff_shear)
    x(6) = CS%sh_xy_h(i,j,2) / (CS%mom_norm_h(i,j,2) + CS%subroundoff_shear)

    x(7) = CS%slope_x_top(i,j,1) / (CS%buoy_norm_h(i,j,1) + CS%subroundoff_shear)
    x(8) = CS%slope_y_top(i,j,1) / (CS%buoy_norm_h(i,j,1) + CS%subroundoff_shear)
    x(9) = CS%slope_x_top(i,j,2) / (CS%buoy_norm_h(i,j,2) + CS%subroundoff_shear)
    x(10) = CS%slope_y_top(i,j,2) / (CS%buoy_norm_h(i,j,2) + CS%subroundoff_shear)

    call ANN_apply(x, y, CS%ann_instance)

    y(1:3) = y(1:3) * CS%Delsq_h(i,j) * CS%mom_norm_h(i,j,1) * CS%mom_norm_h(i,j,1) !upper layer Reynolds stresses, uu, uv, vv
    y(4:6) = y(4:6) * CS%Delsq_h(i,j) * CS%mom_norm_h(i,j,2) * CS%mom_norm_h(i,j,2) !lower layer Reynolds stresses, uu, uv, vv
    y(7:8) = y(7:8) * CS%Delsq_h(i,j) * CS%Coriolis_h(i,j) * CS%mom_norm_h(i,j,1) * CS%buoy_norm_h(i,j,2) !dual form stresses at interface
    y(9) = y(9) * CS%Delsq_h(i,j) * CS%GFS * CS%buoy_norm_h(i,j,1) * CS%buoy_norm_h(i,j,1) !subfilter PE in upper layer
    y(10) = y(10) * CS%Delsq_h(i,j) * CS%GINT * CS%buoy_norm_h(i,j,2) * CS%buoy_norm_h(i,j,2) !subfilter PE in lower layer


    CS%Txx(i,j,1) = y(1)
    CS%Txx_PE(i,j,1) = (y(9) / h(i,j,1))
    !CS%Txx(i,j,1) = y(1) + (y(9) / h(i,j,1))
    Txy(i,j,1) = y(2) 
    !CS%Tyy(i,j,1) = y(3) + (y(9) / h(i,j,1))
    CS%Tyy(i,j,1) = y(3)
    CS%Tyy_PE(i,j,1) = (y(9) / h(i,j,1))
    !CS%Txx(i,j,2) = y(4) + (y(10) / h(i,j,2))
    CS%Txx(i,j,2) = y(4)
    CS%Txx_PE(i,j,2) = (y(10) / h(i,j,2))
    Txy(i,j,2) = y(5)
    !CS%Tyy(i,j,2) = y(6) + (y(10) / h(i,j,2))
    CS%Tyy(i,j,2) = y(6)
    CS%Tyy_PE(i,j,2) = (y(10) / h(i,j,2))
    CS%Txz(i,j,1) = -y(7) !zero top BC - dual form stress
    CS%Txz(i,j,2) = y(7) !dual form stress - zero bottom BC
    CS%Tyz(i,j,1) = -y(8) !zero top BC - dual form stress
    CS%Tyz(i,j,2) = y(8) !dual form stress - zero bottom BC
  enddo; enddo

  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    CS%Txy(I,J,1) = 0.25 * ( (Txy(i+1,j+1,1) + Txy(i,j,1)) &
                           + (Txy(i+1,j,1)   + Txy(i,j+1,1))) * G%mask2dBu(I,J)
    CS%Txy(I,J,2) = 0.25 * ( (Txy(i+1,j+1,2) + Txy(i,j,2)) &
                           + (Txy(i+1,j,2)   + Txy(i,j+1,2))) * G%mask2dBu(I,J)
  enddo; enddo
  
  call pass_var(CS%Txy, G%Domain, clock=CS%id_clock_mpi, position=CORNER)
  call pass_var(CS%Txx, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Tyy, G%Domain, clock=CS%id_clock_mpi)
  call pass_var(CS%Txx_PE, G%Domain, clock = CS%id_clock_mpi)
  call pass_var(CS%Tyy_PE, G%Domain, clock = CS%id_clock_mpi)

  call cpu_clock_end(CS%id_clock_stress_ANN)

end subroutine compute_stress_ANN_collocated_column



subroutine compute_stress_divergence(u, v, h, diffu, diffv, dx2h, dy2h, dx2q, dy2q, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(EPF_CS),         intent(inout) :: CS   !< EPF param control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
        intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
        intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
        intent(in) :: h             !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
        intent(out) :: diffu           !< Zonal acceleration due to convergence of
                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
        intent(out) :: diffv           !< Meridional acceleration due to convergence
                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJ_(G)),           &
        intent(in) :: dx2h          !< dx^2 at h points [L2 ~> m2]
  real, dimension(SZI_(G),SZJ_(G)),           &
        intent(in) :: dy2h          !< dy^2 at h points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)),         &
        intent(in) :: dx2q          !< dx^2 at q points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)),         &
        intent(in) :: dy2q          !< dy^2 at q points [L2 ~> m2]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
        Mxx, & ! Subgrid stress Txx multiplied by thickness and dy^2 [H L4 T-2 ~> m5 s-2]
        Myy, & ! Subgrid stress Tyy multiplied by thickness and dx^2 [H L4 T-2 ~> m5 s-2]
        Mxx_PE, & ! Subgrid stress Txx_PE multiplied by thickness and dy^2 [H L4 T-2 ~> m5 s-2]
        Myy_PE ! Subgrid stress Tyy_PE multiplied by thickness and dx^2 [H L4 T-2 ~> m5 s-2]

  real, dimension(SZIB_(G),SZJB_(G)) :: &
        Mxy    ! Subgrid stress Txy multiplied by thickness [H L2 T-2 ~> m3 s-2]

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
        EPFu              !< Zonal acceleration due to convergence of
                          !! along-coordinate stress tensor for EPF model
                          !! [L T-2 ~> m s-2]

  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
        EPFv              !< Meridional acceleration due to convergence
                          !! of along-coordinate stress tensor for EPF model
                          !! [L T-2 ~> m s-2]
        

  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) !< A term in the kinetic energy budget
                                            ! [H L2 T-3 ~> m3 s-3 or W m-2]

  real :: h_u ! Thickness interpolated to u points [H ~> m or kg m-2].
  real :: h_v ! Thickness interpolated to v points [H ~> m or kg m-2].
  real :: f_x  ! Zonal acceleration due to Reynolds stresses [L T-2 ~> m s-2]
  real :: f_y  ! Meridional acceleration due to Reynolds stresses [L T-2 ~> m s-2]
  ! real :: fyx  ! Meridional acceleration due to Reynolds stresses [L T-2 ~> m s-2]
  ! real :: fyy  ! Meridional acceleration due to Reynolds stresses [L T-2 ~> m s-2]
  real :: fx_PE  ! Zonal acceleration due to subfilter PE [L T-2 ~> m s-2]
  real :: fy_PE  ! Meridional acceleration due to subfilter PE [L T-2 ~> m s-2]
  real :: fxz ! Zonal acceleration due to form stress divergence
  real :: fyz ! Meridional acceleration due to form stress divergence

  real :: h_neglect    ! Thickness so small it can be lost in
                       ! roundoff and so neglected [H ~> m or kg m-2]

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  logical :: save_EPFu, save_EPFv ! Save the acceleration due to EPF model
  
  character(len=100) :: message

  call cpu_clock_begin(CS%id_clock_divergence)


  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff

  EPFu = 0.
  EPFv = 0.

  ! multiplying the stresses by thickness in preparation for taking the divergence
  do k=1,nz
    do J=js-1,Jeq ; do I=is-1,Ieq
      Mxy(I,J) = CS%Txy(I,J,k) * CS%hq(I,J,k)
    enddo ; enddo

    do j=js-1,je+1 ; do i=is-1,ie+1
      Mxx(i,j) = ((CS%Txx(i,j,k)) * h(i,j,k)) * dy2h(i,j)
      Myy(i,j) = ((CS%Tyy(i,j,k)) * h(i,j,k)) * dx2h(i,j) 
      ! Mxx_PE(i,j) = ((CS%Txx_PE(i,j,k)) * h(i,j,k)) * dy2h(i,j)
      ! Myy_PE(i,j) = ((CS%Tyy_PE(i,j,k)) * h(i,j,k)) * dx2h(i,j)    
    enddo ; enddo


    ! Evaluate 1/h x.Div(h S) (Line 1495 of MOM_hor_visc.F90)
    ! Minus occurs because in original file (du/dt) = - div(S),
    ! but here is the discretization of div(S)
    do j=js,je ; do I=Isq,Ieq
      h_u = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k)) + h_neglect 
      ! fxx = ((G%IdyCu(I,j)*(Mxx(i,j) - Mxx(i+1,j))) * G%IareaCu(I,j)) / h_u
      ! fxy = ((G%IdxCu(I,j)*(dx2q(I,J-1)*Mxy(I,J-1) - dx2q(I,J)*Mxy(I,J))) * G%IareaCu(I,j)) / h_u
      f_x = ((G%IdyCu(I,j)*(Mxx(i,j)                - &
                            Mxx(i+1,j))             + &
              G%IdxCu(I,j)*(dx2q(I,J-1)*Mxy(I,J-1)  - &
                            dx2q(I,J)  *Mxy(I,J)))  * &
              G%IareaCu(I,j)) / h_u
      ! fx_PE = ((G%IdyCu(I,j)*(Mxx_PE(i,j)                - &
      !                       Mxx_PE(i+1,j)))  * &
      !         G%IareaCu(I,j)) / h_u
      fxz = - 0.5 * (CS%Txz(i,j,k) + CS%Txz(i+1,j,k)) / h_u
      !diffu(I,j,k) = diffu(I,j,k) !+ fx + fx_z
      CS%f_u(I,j,k)  = diffu(I,j,k)
      CS%f_x(I,j,k) = f_x
      ! CS%fxy(I,j,k) = fxy
      ! CS%fx_PE(I,j,k) = fx_PE
      CS%fxz(I,j,k) = fxz
      ! if (CS%add_EPF_fxx) then
      !   diffu(I,j,k) = diffu(I,j,k) + (CS%amplitude * fxx)
      ! endif
      if (CS%add_EPF_fxy) then
        diffu(I,j,k) = diffu(I,j,k) + (CS%amplitude * f_x)
      endif

      if (CS%add_EPF_subPE) then
        diffu(I,j,k) = diffu(I,j,k) + (CS%amplitude * fx_PE)
      endif

      if (CS%add_EPF_fz) then
        if (CS%fz_negative) then
          diffu(I,j,k) = diffu(I,j,k) - (CS%amplitude * fxz)
        else
          diffu(I,j,k) = diffu(I,j,k) + (CS%amplitude * fxz)
        endif
      endif
    enddo ; enddo

    ! Evaluate 1/h y.Div(h S) (Line 1517 of MOM_hor_visc.F90)
    do J=Jsq,Jeq ; do i=is,ie
      h_v = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k)) + h_neglect
      ! fyx = ((G%IdyCv(i,J)*(dy2q(I-1,J)*Mxy(I-1,J) - dy2q(I,J)*Mxy(I,J)))*G%IareaCv(i,J)) / h_v
      ! fyy = ((G%IdxCv(i,J)*(Myy(i,j) - Myy(i,j+1))) * G%IareaCv(i,J)) / h_v
      f_y = ((G%IdyCv(i,J)*(dy2q(I-1,J)*Mxy(I-1,J)  - &
                            dy2q(I,J)  *Mxy(I,J))   + & ! NOTE this plus
              G%IdxCv(i,J)*(Myy(i,j)                - &
                            Myy(i,j+1)))            * &
              G%IareaCv(i,J)) / h_v
      ! fy_PE = ((G%IdxCv(i,J)*(Myy_PE(i,j)                - &
      !                       Myy_PE(i,j+1)))            * &
      !         G%IareaCv(i,J)) / h_v
      fyz = - 0.5 * (CS%Tyz(i,j,k) + CS%Tyz(i,j+1,k)) / h_v
      !diffv(i,J,k) = diffv(i,J,k) !+ fy + fy_z
      CS%f_v(i,J,k) = diffv(I,j,k)
      CS%f_y(i,J,k) = f_y
      ! CS%fyy(i,J,k) = fyy
      ! CS%fy_PE(i,J,k) = fy_PE
      CS%fyz(i,J,k) = fyz

      if (CS%add_EPF_fxy) then
        diffv(i,J,k) = diffv(i,J,k) + (CS%amplitude * f_y)
      endif

      ! if (CS%add_EPF_fyy) then
      !   diffv(i,J,k) = diffv(i,J,k) + (CS%amplitude * fyy)
      ! endif

      if (CS%add_EPF_subPE) then
        diffv(i,J,k) = diffv(i,J,k) + (CS%amplitude * fy_PE)
      endif

      if (CS%add_EPF_fz) then
        if (CS%fz_negative) then
          diffv(i,J,k) = diffv(i,J,k) - (CS%amplitude * fyz)
        else
          diffv(i,J,k) = diffv(i,J,k) + (CS%amplitude * fyz)
        endif
      endif
      
    enddo ; enddo

  enddo ! end of k loop  
  call cpu_clock_end(CS%id_clock_divergence)

end subroutine compute_stress_divergence


!> Baroclinic Zanna-Bolton-2020 parameterization, see
!! eq. 6 in https://laurezanna.github.io/files/Zanna-Bolton-2020.pdf
!! We compute the lateral stress tensor according to ZB2020 model
!! and update the acceleration due to eddy viscosity (diffu, diffv)
!! as follows:
!! diffu = diffu + ZB2020u
!! diffv = diffv + ZB2020v
subroutine EPF_lateral_stress(u, v, h, tv, diffu, diffv, G, GV, CS, &
                             dx2h, dy2h, dx2q, dy2q, VarMix, US, OBC)
  type(ocean_grid_type),         intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV !< The ocean's vertical grid structure.
  type(EPF_CS),                  intent(inout) :: CS !< EPF control structure.
  type(VarMix_CS),               intent(inout)    :: VarMix !< Variable mixing coefficients
  type(thermo_var_ptrs),         intent(in)    :: tv  !<thermodynamic structure
  type(unit_scale_type),         intent(in)    :: US    !< A dimensional unit scaling type
  type(ocean_OBC_type),      pointer           :: OBC !< Open boundaries control structure

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                 intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2].

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                        intent(inout) :: diffu   !< Zonal acceleration due to eddy viscosity.
                                                 !! It is updated with ZB closure [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                        intent(inout) :: diffv   !< Meridional acceleration due to eddy viscosity.
                                                 !! It is updated with ZB closure [L T-2 ~> m s-2]

  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: dx2h    !< dx^2 at h points [L2 ~> m2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: dy2h    !< dy^2 at h points [L2 ~> m2]

  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: dx2q    !< dx^2 at q points [L2 ~> m2]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: dy2q    !< dy^2 at q points [L2 ~> m2]

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  call cpu_clock_begin(CS%id_clock_module)

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Compute the stress tensor given the
  ! (optionally sharpened) velocity gradients
  if (CS%ann_type == 0) then
    call compute_stress_ANN_collocated(h, tv, G, GV, CS, US, VarMix, OBC)
  else if (CS%ann_type == 1) then
    call compute_stress_ANN_collocated_column(h, tv, G, GV, CS, US, VarMix, OBC)
  else if (CS%ann_type == 2) then
    call compute_stress_separate_ANN(h, tv, G, GV, CS, US, VarMix, OBC)
  endif
  !call compute_stress_ANN_collocated(G, GV, CS)

  ! Update the acceleration due to eddy viscosity (diffu, diffv)
  ! with the ZB2020 lateral parameterization
  call compute_stress_divergence(u, v, h, diffu, diffv,    &
                                 dx2h, dy2h, dx2q, dy2q, &
                                 G, GV, CS)

  call cpu_clock_begin(CS%id_clock_post)
  if (CS%id_Txx > 0) call post_data(CS%id_Txx, CS%Txx, CS%diag)
  if (CS%id_Tyy > 0) call post_data(CS%id_Tyy, CS%Tyy, CS%diag)
  if (CS%id_Txx_PE > 0) call post_data(CS%id_Txx_PE, CS%Txx_PE, CS%diag)
  if (CS%id_Tyy_PE > 0) call post_data(CS%id_Tyy_PE, CS%Tyy_PE, CS%diag)
  if (CS%id_Txy > 0) call post_data(CS%id_Txy, CS%Txy, CS%diag)
  if (CS%id_Txz > 0) call post_data(CS%id_Txz, CS%Txz, CS%diag)
  if (CS%id_Tyz > 0) call post_data(CS%id_Tyz, CS%Tyz, CS%diag)
  if (CS%id_Txy_h > 0) call post_data(CS%id_Txy_h, CS%Txy_h, CS%diag)

  if (CS%id_sh_xx > 0) call post_data(CS%id_sh_xx, CS%sh_xx, CS%diag)
  if (CS%id_sh_xy > 0) call post_data(CS%id_sh_xy, CS%sh_xy, CS%diag)
  if (CS%id_vort_xy > 0) call post_data(CS%id_vort_xy, CS%vort_xy, CS%diag)
  if (CS%id_sh_xy_h > 0) call post_data(CS%id_sh_xy_h, CS%sh_xy_h, CS%diag)
  if (CS%id_vort_xy_h > 0) call post_data(CS%id_vort_xy_h, CS%vort_xy_h, CS%diag)
  
  if (CS%id_eta > 0) call post_data(CS%id_eta, CS%eta, CS%diag)
  if (CS%id_h_eta > 0) call post_data(CS%id_h_eta, CS%h_eta, CS%diag)
  if (CS%id_slope_x > 0) call post_data(CS%id_slope_x, CS%slope_x, CS%diag)
  if (CS%id_slope_y > 0) call post_data(CS%id_slope_y, CS%slope_y, CS%diag)
  if (CS%id_depth_mask > 0) call post_data(CS%id_depth_mask, CS%depth_mask, CS%diag)

  if (CS%id_slope_x_top > 0) call post_data(CS%id_slope_x_top, CS%slope_x_top, CS%diag)
  if (CS%id_slope_y_top > 0) call post_data(CS%id_slope_y_top, CS%slope_y_top, CS%diag)
  if (CS%id_slope_x_bot > 0) call post_data(CS%id_slope_x_bot, CS%slope_x_bot, CS%diag)
  if (CS%id_slope_y_bot > 0) call post_data(CS%id_slope_y_bot, CS%slope_y_bot, CS%diag)
  if (CS%id_Coriolis_h > 0) call post_data(CS%id_Coriolis_h, CS%Coriolis_h, CS%diag)
  if (CS%id_Delsq_h > 0) call post_data(CS%id_Delsq_h, CS%Delsq_h, CS%diag)
  if (CS%id_mom_norm_h > 0) call post_data(CS%id_mom_norm_h, CS%mom_norm_h, CS%diag)
  if (CS%id_buoy_norm_h > 0) call post_data(CS%id_buoy_norm_h, CS%buoy_norm_h, CS%diag)

  if (CS%id_f_u > 0) call post_data(CS%id_f_u, CS%f_u, CS%diag)
  if (CS%id_f_x > 0) call post_data(CS%id_f_x, CS%f_x, CS%diag)
  if (CS%id_fx_PE > 0) call post_data(CS%id_fx_PE, CS%fx_PE, CS%diag)
  if (CS%id_fxz > 0) call post_data(CS%id_fxz, CS%fxz, CS%diag)

  if (CS%id_f_v > 0) call post_data(CS%id_f_v, CS%f_v, CS%diag)
  if (CS%id_f_y > 0) call post_data(CS%id_f_y, CS%f_y, CS%diag)
  if (CS%id_fy_PE > 0) call post_data(CS%id_fy_PE, CS%fy_PE, CS%diag)
  if (CS%id_fyz > 0) call post_data(CS%id_fyz, CS%fyz, CS%diag)

  call cpu_clock_end(CS%id_clock_post)

  call cpu_clock_end(CS%id_clock_module)

end subroutine EPF_lateral_stress

!> Deallocate any variables allocated in EPF_init
subroutine EPF_end(CS)
  type(EPF_CS), intent(inout) :: CS  !< EPF control structure.
  deallocate(CS%sh_xx)
  deallocate(CS%sh_xy)
  deallocate(CS%vort_xy)
  deallocate(CS%sh_xy_h)
  deallocate(CS%vort_xy_h)
  deallocate(CS%hq)

  deallocate(CS%Txx)
  deallocate(CS%Tyy)
  deallocate(CS%Txx_PE)
  deallocate(CS%Tyy_PE)
  deallocate(CS%Txy)
  deallocate(CS%Txy_h)
  deallocate(CS%Txz)
  deallocate(CS%Tyz)
  deallocate(CS%eta)
  deallocate(CS%h_eta)
  deallocate(CS%depth_mask)

  deallocate(CS%kappa_h)
  deallocate(CS%kappa_q)

  deallocate(CS%Coriolis_h)

  deallocate(CS%slope_x)
  deallocate(CS%slope_y)
  deallocate(CS%slope_x_top)
  deallocate(CS%slope_y_top)
  deallocate(CS%slope_x_bot)
  deallocate(CS%slope_y_bot)

  deallocate(CS%Delsq_h)
  deallocate(CS%mom_norm_h)
  deallocate(CS%buoy_norm_h)

  deallocate(CS%f_u)
  deallocate(CS%f_x)
  deallocate(CS%fx_PE)
  deallocate(CS%fxz)

  deallocate(CS%f_v)
  deallocate(CS%f_y)
  deallocate(CS%fy_PE)
  deallocate(CS%fyz)
  
end subroutine EPF_end


end module MOM_EPF_ANN