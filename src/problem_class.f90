!define problem

module problem_class

  use fault_stress, only : kernel_type
  use mesh, only : mesh_type, spectral_mesh_type

  implicit none

  public

  ! Placeholder for output quantities
  type optr
    ! Vector quantities
    double precision, pointer :: v(:) => null()
    ! Scalar quantities (double)
    double precision, pointer :: s => null()
    ! Scalar quantities (integer)
    integer, pointer :: i => null()
  end type optr

 ! timeseries outputs: at every time step, but only macroscopic quantities
  type ot_type
    double precision :: v_th
    double precision, dimension(:), allocatable :: xsta, ysta, zsta, v_pre, v_pre2
    integer :: ic=-1, ntout=0, not=-1, not_vmax=-1
    integer, dimension(:), allocatable :: iasp, iot
    character(len=16), dimension(:), allocatable :: fmt, fmt_vmax
    type(optr), dimension(:), allocatable :: objects_ot, objects_vmax
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    integer :: ntout=0, nrup=-1
    integer :: count, dyn_count, dyn_stat, i_ox_dyn, i_ox_seq, seq_count
    integer :: nwout, nwout_dyn, nxout, nxout_dyn
    double precision :: pot_pre
    character(len=256) :: header
    character(len=16), dimension(:), allocatable :: fmt
    type(optr), dimension(:), allocatable :: objects_rup
  end type ox_type

  ! SEISMIC: structure that holds the CNS model parameters
  ! See input.f90 for a description of the parameters
  type cns_type
    double precision, dimension(:), allocatable :: &
      a_tilde, mu_tilde_star, H, L, & ! Granular flow and microstructural params
      y_gr_star, phi_c, phi0, lambda, &
      A, n, m, A_bulk, n_bulk, m_bulk ! Creep mechanism parameters
    integer :: N_creep=-1 ! Number of active creep mechanisms
  end type cns_type
  ! End of the CNS model structure

  ! SEISMIC: structure that holds the thermal pressurisation (TP) model parameters
  ! See input.f90 for a description of the parameters
  ! Spectral mesh parameters (Dlogl, lw_max, Nl) are hard-coded in mesh.f90
  type tp_type
    type (spectral_mesh_type) :: mesh
    double precision, dimension(:), allocatable :: &
      rhoc, inv_rhoc, beta, eta, k_t, k_p, l, w, inv_w, P_a, T_a, &
      Pi, Theta, PiTheta, Omega, dP_dt, Theta_prev, PiTheta_prev, &
      alpha_th, alpha_hy, Lam, Lam_prime, Lam_T, phi_b, dilat_factor
    double precision :: t_prev=0d0
  end type tp_type
  ! End of the TP model structure

  ! Structure that holds the fluid injection model parameters
  ! See input.f90 for a description of the parameters
  type injection_type
    !double precision, dimension(:), allocatable :: k_min, k_max !p_a, DP_DT, dk_dt, FMCk, qBmkh   
    double precision :: p_a = 0d0, dP_dt = 0d0, q = -0.002d0, B0 = 1d0, fi = 0.30d0, &
                        muv = 0.09d-3, ct = 1.45d-10, k = 1d-15, k_min = 1d-18, k_max = 1d-14, &  ! k_min = 1d-18, k_max = 1d-14
                        rw = 0.10d0, re = 500, h = 500, nroots = 10, t0 = 0d0, tf = 0d0, Y0 = 0d0, X0 = 0d0, theta = 30d0  !q = -0.002d0, B0 = 1d0, muv = 0.09d-3, k = 1d-12, ct = 1.45d-10, fi = 0.30d0 &        
    real, dimension(10) :: an = (/0.000766341288526, 0.001403117645524,0.002034694279645,0.002664739504977, &
                                  0.003294127716977, 0.003923174121984, 0.004552020133089, 0.005180738635555, &
                                  0.005809371010047 ,  0.006437942494247/)                    
  end type injection_type
  ! End of the fluid injection structure  
  ! End of the fluid injection structure

  ! SEISMIC: requested features structure
  ! stress_coupling: normal stress variations at subduction interfaces
  ! tp: thermal pressurisation
  ! cohesion: time-dependent cohesion (CNS model)
  ! localisation: localisation of deformation (CNS model)
  type features_type
    integer :: stress_coupling, tp, localisation, injection
  end type features_type
  ! End of features structure

  ! SEISMIC: Runge-Kutta-Fehlberg solver parameters
  type rk45_type
    integer :: iflag
    double precision, dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork
  end type rk45_type
  ! End of features structure

  ! SEISMIC: Runge-Kutta-Fehlberg solver parameters
  type rk45_2_type
    double precision, dimension(5) :: t_coeffs
    double precision, dimension(5) :: error_coeffs
    double precision, dimension(19) :: coeffs
    double precision, dimension(:), allocatable :: k1, k2, k3, k4, k5, k6
    double precision, dimension(:), allocatable :: dy1, dy2, dy3, dy4, dy5, dy6
    double precision, dimension(:), allocatable :: y_new, h, e
  end type rk45_2_type
  ! End of features structure

  ! Structure for (unit)testing
  type test_type
    logical :: test_mode = .false.
    logical :: test_passed = .true.
  end type test_type
  ! End of testing structure

  type problem_type
    type (mesh_type) :: mesh
    type (kernel_type) :: kernel

    ! Containers for local and global quantities
    type(optr), dimension(:), allocatable :: objects_glob, objects_loc
    ! Number of objects in the containers
    integer :: nobj=9

    ! Basic variables
    ! NOTE: if theta2 needs to be in output sequence, add here
    double precision, pointer ::  tau(:) => null(), dtau_dt(:) => null(), &
                                  sigma(:) => null(), slip(:) => null(), &
                                  v(:) => null(), theta(:) => null(), &
                                  P(:) => null(), T(:) => null()
    double precision, dimension(:), allocatable :: dtheta_dt, dtheta2_dt
    double precision, dimension(:), allocatable :: theta2, alpha
    double precision, pointer ::  tau_max(:) => null(), t_rup(:) => null(), &
                                  v_max(:) => null(), t_vmax(:) => null()
    double precision, pointer :: pot => null(), pot_rate => null()
   ! Boundary conditions
    integer :: i_sigma_cpl=0, finite=0
   ! Friction properties
    double precision, dimension(:), allocatable ::  a, b, dc, v1, v2, &
                                                    mu_star, v_star, &
                                                    theta_star, coh, v_pl
    integer :: itheta_law=1, i_rns_law=1, neqs=3
   ! Elastic properties
    double precision :: beta=0d0, smu=0d0, lam=0d0, D=0d0, H=0d0, zimpedance=0d0
   ! Periodic loading
    double precision :: Tper=0d0, Aper=0d0, Omper=0d0
   ! Time solver
    double precision :: t_prev=0d0, dt_try=0d0, dt_did=0d0, &
                        dt_next=0d0, dt_max=0d0, tmax, acc
    double precision, pointer :: time => null()
    integer, pointer :: ivmax => null()
    integer :: NSTOP, itstop=-1, it=0
    double precision :: vmaxglob
   ! For outputs
    type (ot_type) :: ot
    type (ox_type) :: ox
   ! For MPI outputs
    double precision, pointer ::  tau_glob(:) => null(), dtau_dt_glob(:) => null(), &
                                  sigma_glob(:) => null(), slip_glob(:) => null(), &
                                  v_glob(:) => null(), theta_glob(:) => null(), &
                                  P_glob(:) => null(), T_glob(:) => null()
    double precision, pointer ::  tau_max_glob(:) => null(), t_rup_glob(:) => null(), &
                                  v_max_glob(:) => null(), t_vmax_glob(:) => null()
    logical :: allocated_glob = .false.
   ! QSB
    double precision :: DYN_M,DYN_th_on,DYN_th_off
    integer :: DYN_FLAG,DYN_SKIP
    ! Flag for unit testing
    logical :: test_mode = .false.

    ! SEISMIC: added structures
    type (cns_type) :: cns_params
    type (tp_type) :: tp
    type (injection_type) :: injection
    type (features_type) :: features
    type (rk45_type) :: rk45
    type (rk45_2_type) :: rk45_2
    type (test_type) :: test
  end type problem_type

end module problem_class
