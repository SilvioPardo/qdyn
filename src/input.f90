!qdyn input all parameters

module input

  implicit none
  private

  public :: read_main

contains
!=====================================================================
! Read input file qdyn.in
! If MPI, there is one qdyn*.in file per processor
!
subroutine read_main(pb)

  use problem_class
  use mesh, only : read_mesh_parameters, read_mesh_nodes, mesh_get_size
  use output, only : ot_read_stations
  use my_mpi, only : my_mpi_tag, is_MPI_parallel
  use constants, only : FAULT_TYPE, SOLVER_TYPE, FID_IN, FID_SCREEN

  type(problem_type), intent(inout) :: pb

  double precision, dimension(:), allocatable :: read_buf
  integer :: i, j, n, N_cols, N_creep

  write(FID_SCREEN, *) 'Start reading input ...'

  open(unit=FID_IN, file='qdyn'//trim(my_mpi_tag())//'.in', action='read')

  call read_mesh_parameters(FID_IN, pb%mesh)
  write(FID_SCREEN, *) '   Mesh input complete'

  if (pb%mesh%dim==1) read(FID_IN, *) pb%finite
  read(FID_IN, *) pb%itheta_law
  read(FID_IN, *) pb%i_rns_law
  read(FID_IN, *) pb%i_sigma_cpl
  ! SEISMIC: various simulation features can be turned on (1) or off (0)
  if (pb%i_rns_law == 3) then
    read(FID_IN, *) pb%cns_params%N_creep
  endif
  read(FID_IN, *) pb%features%stress_coupling, pb%features%tp, &
                  pb%features%localisation, pb%features%injection
  read(FID_IN, *) pb%ot%ntout, pb%ox%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nwout, &
                  pb%ox%nxout_dyn, pb%ox%nwout_dyn, pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(FID_IN, *) pb%beta, pb%smu, pb%lam, pb%D, pb%H, pb%ot%v_th
  read(FID_IN, *) pb%Tper, pb%Aper
  read(FID_IN, *) pb%dt_try, pb%dt_max, pb%tmax, pb%acc
  read(FID_IN, *) pb%NSTOP
  read(FID_IN, *) pb%DYN_FLAG,pb%DYN_SKIP
  read(FID_IN, *) pb%DYN_M,pb%DYN_th_on,pb%DYN_th_off
  read(FID_IN, *) FAULT_TYPE, SOLVER_TYPE
  write(FID_SCREEN, *) '  Flags input complete'

  ! MvdE: if fluid injection (analytical) is requested, read parameters
  if ((pb%features%injection == 1) .or. (pb%features%injection == 2)) then
  
    ! Check if thermal pressurisation is also active
    if (pb%features%tp == 1) then
      write(FID_SCREEN, *) &
        "Fluid injection is not compatible with thermal pressurisation"
      stop
    endif
 !   allocate(  pb%injection%p_i(n), &                                                       ! Silvio edit
  !             pb%injection%q(n), pb%injection%B(n), pb%injection%muv(n), pb%injection%k(n), &
   !            pb%injection%fi(n), pb%injection%ct(n), pb%injection%h(n), pb%injection%rw(n), &
   !            pb%injection%re(n), pb%injection%nroots(n), pb%injection%t0(n), pb%injection%Y0(n), &
   !            pb%injection%X0(n), pb%injection%an(n)  )
   ! do i=1,n
   ! read(FID_IN, *) pb%injection%p_a,pb%injection%q, pb%injection%B0, pb%injection%muv, &
   !                 pb%injection%k, pb%injection%fi, pb%injection%ct, pb%injection%h, &
   !                 pb%injection%rw, pb%injection%re, pb%injection%nroots, pb%injection%t0, &
   !                 pb%injection%Y0, pb%injection%X0, pb%injection%an, pb%injection%tf
   
   ! write(FID_SCREEN, *) "p_a, q, B0, muv, k, fi, ct, h, rw, re, nroots, t0, Y0, X0, tf ", &
   !                    pb%injection%p_a, pb%injection%q, pb%injection%B0, pb%injection%muv, pb%injection%k, & ! Silvio edits
   !                    pb%injection%fi, pb%injection%ct, pb%injection%h, pb%injection%rw, &
   !                    pb%injection%re, pb%injection%nroots, pb%injection%t0, pb%injection%Y0, &
   !                    pb%injection%X0, pb%injection%an, pb%injection%tf                    
    !end do
  !endif

    ! Source location, source time, injection rate, diffusivity ! Silvio Edit
    pb%injection%p_a = 0d0
    pb%injection%q = -0.025d0
    pb%injection%B0 = 1d0
    pb%injection%muv = 0.09d-3
    pb%injection%k = 1d-14
    pb%injection%k_min = 1d-18
    pb%injection%k_max = 1d-14            
    pb%injection%fi = 0.15d0
    pb%injection%ct = 1.45d-10
    pb%injection%h = 500
    pb%injection%rw = 0.10d0
    pb%injection%re = 1000
    pb%injection%nroots = 10
    pb%injection%t0 = 37559376000.
    pb%injection%tf = 44150400
    pb%injection%Y0 = 0d0
    pb%injection%X0 = 0d0
    pb%injection%theta = 0.8726d0
    pb%injection%an = (/0.0008, 0.0014, 0.002, 0.0027, 0.0033, 0.0039, 0.0046, 0.0052, 0.0058, 0.0064/)
    pb%injection%an = (/0.000383170608831, 0.000701558705927, 0.001017346895086, 0.001332369333348, & 
                       0.001647063218443, 0.001961586153544, 0.002276008845199, 0.002590367736046, &
                      0.002904683516422, 0.003847477786348/)   
    !pb%injection%an = (/0.000766341288526, 0.001403117645524,0.002034694279645,0.002664739504977, &
    !                  0.003294127716977, 0.003923174121984, 0.004552020133089, 0.005180738635555, &
    !                  0.005809371010047 ,  0.006437942494247/)
    !read(FID_IN, *) pb%injection%p_i, pb%injection%q, pb%injection%B0, pb%injection%muv, pb%injection%k, &
     !               pb%injection%fi, pb%injection%ct, pb%injection%h, pb%injection%rw, &
      !              pb%injection%re, pb%injection%nroots, pb%injection%t0, pb%injection%Y0, &
       !             pb%injection%X0, pb%injection%an
  endif
 ! write(FID_SCREEN, *) pb%injection%p_i, pb%injection%q, pb%injection%B0, pb%injection%muv, pb%injection%k, & ! Silvio edits
 !                      pb%injection%fi, pb%injection%ct, pb%injection%h, pb%injection%rw, &
 !                      pb%injection%re, pb%injection%nroots, pb%injection%t0, pb%injection%Y0, &
 !                      pb%injection%X0, pb%injection%tf
 !write(FID_SCREEN, *) '10 roots',pb%injection%an
 ! write(FID_SCREEN, *) '  Flags input complete'

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor
  allocate ( pb%tau(n), pb%sigma(n), pb%v(n), pb%theta(n),  &
             pb%v_star(n), pb%ot%iot(n), pb%ot%iasp(n), &
             pb%dc(n), pb%coh(n), pb%v_pl(n) )

  !JPA if MPI, read only the nodes of this processor
  ! <SEISMIC>
  ! Read input parameters for the CNS model. These parameters are (in order):
  !   sigma:            applied normal stress (effective if no TP)
  !   tau:              initial shear stress
  !   theta:            initial porosity (phi in CNS formulation)
  !   v_pl:             load-point velocity
  !   a_tilde:          coefficient of logarithmic rate dependence
  !   mu_tilde_star:    reference friction coefficient at y_gr_star
  !   y_gr_star:        reference granular fow strain rate
  !   H:                dilatancy geometric factor
  !   phi_c:            critical state porosity
  !   phi0:             lower cut-off porosity
  !
  !   (for each out of N creep mechanisms)
  !   A:                creep law kinetic constant
  !   n:                stress exponent
  !   m:                porosity exponent
  !
  !   L:                total thickness of the fault zone
  !   iot:              flag for other time series output
  !   iasp:             flag for identification purposes
  !
  ! Note that these parameters are material (gouge) properties, and are
  ! generally not spatically uniform, and hence are allocatable
  ! See friction.f90 for a description of and references to the CNS model
  ! See user manual for detailed definitions of the above parameters
  ! (this is a permanent TODO)
  ! </SEISMIC>

  ! If the CNS model is selected
  if (pb%i_rns_law == 3) then
    N_creep = pb%cns_params%N_creep
    ! Define the number of columns to read: 13 plus 3 for each creep mechanism
    N_cols = 13 + 3*N_creep
    ! Allocate read buffer
    allocate( read_buf(N_cols) )
    ! Allocate CNS parameters
    allocate( pb%cns_params%a_tilde(n), pb%cns_params%mu_tilde_star(n), &
              pb%cns_params%y_gr_star(n), pb%cns_params%H(n), &
              pb%cns_params%phi_c(n), pb%cns_params%phi0(n), &
              pb%cns_params%L(n) )
    ! Allocate CNS parameters for N creep mechanisms
    allocate( pb%cns_params%A(N_creep*n), pb%cns_params%n(N_creep*n), &
              pb%cns_params%m(N_creep*n) )
    do i=1,n

      ! Read data into buffer
      read(FID_IN, *) read_buf(:)

      ! General fault parameters
      pb%sigma(i) = read_buf(1)
      pb%tau(i) = read_buf(2)
      pb%theta(i) = read_buf(3)
      pb%v_pl(i) = read_buf(4)

      ! Granular flow parameters
      pb%cns_params%a_tilde(i) = read_buf(5)
      pb%cns_params%mu_tilde_star(i) = read_buf(6)
      pb%cns_params%y_gr_star(i) = read_buf(7)
      pb%cns_params%H(i) = read_buf(8)

      ! Porosity parameters
      pb%cns_params%phi_c(i) = read_buf(9)
      pb%cns_params%phi0(i) = read_buf(10)

      ! Creep mechanism parameters (3 per mechanism)
      do j=1,N_creep
        ! Kinetic constant
        pb%cns_params%A((i-1)*N_creep+j) = read_buf(11+(j-1)*3)
        ! Stress exponent
        pb%cns_params%n((i-1)*N_creep+j) = read_buf(12+(j-1)*3)
        ! Porosity exponent
        pb%cns_params%m((i-1)*N_creep+j) = read_buf(13+(j-1)*3)
      enddo

      ! Fault zone thickness
      pb%cns_params%L(i) = read_buf(N_cols-2)

      ! Output parameters
      pb%ot%iot(i) = int(read_buf(N_cols-1))
      pb%ot%iasp(i) = int(read_buf(N_cols))

      ! Set dummy values to unused parameters
      pb%dc(i) = 1.0
      pb%V(i) = 0.0

    end do
    ! Deallocate the read buffer
    deallocate(read_buf)
  ! Else, the RSF model is selected
  else
    allocate ( pb%a(n), pb%b(n), pb%v1(n), &
               pb%v2(n), pb%mu_star(n))
    do i=1,n
      read(FID_IN, *) pb%sigma(i), pb%v(i), pb%theta(i),  &
                      pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
                      pb%v2(i), pb%mu_star(i), pb%v_star(i), &
                      pb%ot%iot(i), pb%ot%iasp(i), pb%coh(i), pb%v_pl(i)
    end do
  endif

  ! <SEISMIC>
  ! Read input parameters for the localisation model (CNS only).
  ! These parameters and corresponding units are (in order):
  !
  if (pb%features%localisation == 1) then
    ! Raise an error if the CNS model is not selected
    if (pb%i_rns_law /= 3) then
      write(FID_SCREEN, *) &
        "Localisation of shear strain is compatible only with the CNS model (i_rns_law = 3)"
      stop
    endif

    N_creep = pb%cns_params%N_creep
    ! Define the number of columns to read: 2 plus 3 for each creep mechanism
    N_cols = 2 + 3*N_creep
    ! Allocate read buffer
    allocate( read_buf(N_cols) )
    allocate( pb%cns_params%lambda(n), pb%theta2(n) )
    allocate( pb%cns_params%A_bulk(N_creep*n), pb%cns_params%n_bulk(N_creep*n),&
              pb%cns_params%m_bulk(N_creep*n) )

    do i=1,n

      ! Read data into buffer
      read(FID_IN, *) read_buf(:)

      ! Degree of localisation (0 <= lambda <= 1)
      pb%cns_params%lambda(i) = read_buf(1)
      ! Porosity in bulk gouge
      pb%theta2(i) = read_buf(2)

      ! Creep mechanism parameters (3 per mechanism)
      do j=1,N_creep
        ! Kinetic constant
        pb%cns_params%A_bulk((i-1)*N_creep+j) = read_buf(3+(j-1)*3)
        ! Stress exponent
        pb%cns_params%n_bulk((i-1)*N_creep+j) = read_buf(4+(j-1)*3)
        ! Porosity exponent
        pb%cns_params%m_bulk((i-1)*N_creep+j) = read_buf(5+(j-1)*3)
      enddo

    end do
    ! Deallocate the read buffer
    deallocate(read_buf)
  ! Else, set lambda and theta2 to dummy values
  else
    allocate(pb%cns_params%lambda(n), pb%theta2(n))
    do i=1,n
      pb%cns_params%lambda(i) = 1d0
      pb%theta2(i) = 0d0
    end do
  endif
  ! End reading localisation model parameters
  ! </SEISMIC>

  ! <SEISMIC>
  ! Read input parameters for the thermal pressurisation (TP) model.
  ! These parameters and corresponding units are (in order):
  !
  !   rhoc:   density times specific heat capacity of host rock [J/K/m^3]
  !   beta:   bulk compressibility [1/Pa]
  !   eta:    dynamic viscocity [Pa s]
  !   w:      half-width of shear distribution profile [m]
  !   k_t:    thermal conductivity [J/s/K/m]
  !   k_p:    intrinsic hydraulic permeability [m^2]
  !   l:      nett thermal expansion coefficient [1/K]
  !   P_a:    ambient fluid pressure [Pa]
  !   T_a:    ambient temperature [K]
  !
  if (pb%features%tp == 1) then
    allocate (  pb%tp%rhoc(n), pb%tp%beta(n), pb%tp%eta(n), pb%tp%w(n), &
                pb%tp%k_t(n), pb%tp%k_p(n), pb%tp%l(n), &
                pb%tp%P_a(n), pb%tp%T_a(n), pb%tp%dilat_factor(n) )
    do i=1,n
      read(FID_IN, *) pb%tp%rhoc(i), pb%tp%beta(i), pb%tp%eta(i), pb%tp%w(i), &
                      pb%tp%k_t(i), pb%tp%k_p(i), pb%tp%l(i), &
                      pb%tp%P_a(i), pb%tp%T_a(i), pb%tp%dilat_factor(i)
    end do
  endif
  ! End reading TP model parameters
  ! </SEISMIC>

  ! Read input parameters for the fluid injection model.
  
  !   dP_dt:  constant rate of pressure increase [Pa/s]
  !   P_a:    ambient (constant) fluid pressure [Pa]
  
  if ((pb%features%injection == 1) .or. (pb%features%injection == 2)) then !Silvio
  
    if (pb%features%tp == 1) then
      write(FID_SCREEN, *) &
        "The fluid injection module is incompatible with the thermal pressurisation module"
      stop
    endif
  endif
 !   allocate (pb%injection%p_a(n), pb%injection%dP_dT(n), pb%injection%FMCk(n), pb%injection%qBmkh(n), &
 !             pb%injection%dk_dt(n))
              !pb%injection%q(n), pb%injection%B0(n), pb%injection%muv(n), &
              !pb%injection%k(n), pb%injection%fi(n), pb%injection%ct(n), pb%injection%h(n), &
              !pb%injection%rw(n), pb%injection%re(n), pb%injection%nroots(n), pb%injection%t0(n), &
              !pb%injection%Y0(n), pb%injection%X0(n), pb%injection%an(n), pb%injection%tf(n))
 !   do i=1,n
 !    read(FID_IN, *) pb%injection%p_a(i), pb%injection%dP_dT(i), pb%injection%FMCk(i), pb%injection%qBmkh(i), &
 !                    pb%injection%dk_dt(i)
                     !pb%injection%q(i), pb%injection%B0(i), pb%injection%muv(i), &
                     !pb%injection%k(i), pb%injection%fi(i), pb%injection%ct(i), pb%injection%h(i), &
                     !pb%injection%rw(i), pb%injection%re(i), pb%injection%nroots(i), pb%injection%t0(i), &
                     !pb%injection%Y0(i), pb%injection%X0(i), pb%injection%an(i), pb%injection%tf(i)
 !   end do
 !   write(FID_SCREEN, *) "p_a, dP_dt, FMCk, qBmkh, dk_dt", pb%injection%p_a, pb%injection%dP_dt, &
 !                        pb%injection%FMCk, pb%injection%qBmkh, pb%injection%dk_dt 
 ! endif
    ! Read input parameters for the fluid injection model.
  
  !   dP_dt:  constant rate of pressure increase [Pa/s] ! Silvio Edit, first attempt. Now modified reading constant values from py files
  !   P_a:    ambient (constant) fluid pressure [Pa]
!  if (pb%features%injection == 2) then 
 !   if (pb%features%tp == 1) then
  !    write(FID_SCREEN, *) &
   !     "The fluid injection module is incompatible with the thermal pressurisation module"
    !  stop
   ! endif

    !allocate (pb%injection%p_i(n), pb%injection%q(n), pb%injection%B(n), pb%injection%muv(n), pb%injection%k(n), pb%injection%fi(n), &
     !         pb%injection%ct(n), pb%injection%h(n), pb%injection%rw(n), pb%injection%re(n)) ! What about roots?
            
                   
   ! do i=1,n
   !   read(FID_IN, *) pb%injection%p_i(i), pb%injection%q(i), pb%injection%B(i), pb%injection%muv(i), pb%injection%k(i), pb%injection%fi(i), &
   !           pb%injection%ct(i), pb%injection%h(i), pb%injection%rw(i), pb%injection%re(i) ! What about roots?
   ! end do
 ! endif  
  ! End reading fluid injection model parameters

  if (pb%mesh%dim == 2 .and. is_MPI_parallel()) then
    call read_mesh_nodes(FID_IN, pb%mesh)
    ! call ot_read_stations(pb%ot)
  endif

  close(FID_IN)
  write(FID_SCREEN, *) 'Input complete'

end subroutine read_main


end module input
