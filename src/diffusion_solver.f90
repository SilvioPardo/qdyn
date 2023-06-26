!===============================================================================
!
! SEISMIC: spectral solver module for diffusion problem
!
! This module handles integration of the diffusion equations for heat and
! fluid pressure, following the method proposed by Noda & Lapusta (2010) [N&L].
! The diffusion normal to the fault plane is solved in a discretised spectral
! domain, which is numerically stable for large time steps (as opposed to
! discretisation methods in the spatial domain, such as FEM, FD, etc.). This
! procedure allows for efficient modelling of frictional heating and associated
! pressurisation of the pore fluid (thermal pressurisation: TP)
! This module is compatible with both the rate-and-state friction and the
! Chen-Niemeijer-Spiers frameworks
!
! TODO: add a coefficient that regulates the dilatancy hardening component
!
! References:
!
! [N&L]
! Noda, H., and N. Lapusta (2010): Three‐dimensional earthquake sequence
! simulations with evolving temperature and pore pressure due to shear
! heating: Effect of heterogeneous hydraulic diffusivity,
! J. Geophys. Res., 115, B12314, doi:10.1029/2010JB007780
!
!===============================================================================

module diffusion_solver

  use problem_class, only : problem_type
  use mesh, only : spectral_mesh_type

  implicit none
  private

  public :: init_tp,update_PT, update_PT_final, update_P_constant, update_P_Galis

contains

!===============================================================================
! Routine to update the fluid pressure assuming a constant pressurisation rate
!===============================================================================
subroutine update_P_constant(time, dP_dt, pb)

  type(problem_type), intent(inout) :: pb
  double precision, intent(in) :: time
  double precision, dimension(pb%mesh%nn) :: dP_dt

  ! Update the fluid pressure using a linear time model
  pb%P = time * pb%injection%dP_dt + pb%injection%p_a 
  ! Update the time derivative accordingly (constant)
  dP_dt = pb%injection%dP_dt

end subroutine update_P_constant
!===============================================================================================
! Routine to update the fluid pressure assuming a constant pressurisation rate Galis et al. 2017
!===============================================================================================
subroutine update_P_Galis(x, time,P, dP_dt, pb) ! acoeff, nroots, p_i, q, B, muv, k, fi, ct, h, rw, re)
  type(problem_type), intent(inout) :: pb
  !double precision, intent(in) :: time ! maybe intent (in) has to be used Silvio edit
  double precision :: time 
  double precision, dimension(pb%mesh%nn) :: x, y, P, dP_dt, r, r0, h, rD, rD2, pDpar, A, B, C, Ct, sump, sumpt, pD, Q, Quort
  double precision, dimension(pb%mesh%nn) :: P1,P2,P3,P4,P5,P6,PAA,ks,k, qBmkh, FMCk, tD,pippo, hr
  !double precision :: B ! double precision non gli piaceva
  double precision, parameter :: PI = 3.14159265358979d0
 ! double precision  ::  p_i,q,B,muv,k,fi,ct,h,rw,re !check how to define in input.f90
  !integer :: r_p, r_n!, nroots stanno dentro pb%injection !check how to define in input.f90
  double precision :: reD,reD2! double precision, dimension(pb%mesh%nn) non gli piaceva
  ! ask martijn because it is an array, I deleted rD ! double precision, dimension(pb%mesh%nn) non gli piaceva
  integer :: i,j
  
  !p = ps_cylindrical_noflow_v( r, t, p_i, q, B, mu, k, h, fi, ct, rw, re, a, n )
! calculated pD(rD,tD) following Eq 10 in Galis et al. 2017
!  [original formula Eq. B-69 in Pressure Transient Testing by Lee, Rolling and Spivey, 2003]

! output, p [in Pa], is pore pressure at position r [in m] at time t [in s]
! r OR t may be vector (not both at the same time!)
! pressure outside the reservoir is 0 Pa due to the noflow boundary condition
 !  p_i =  0e6 ! Pa
 ! q   = -0.025 ! m^3/s
 ! B   =  1.0
 ! muv =  0.09e-3 ! Pa.s
 ! k   =  1e-17 ! m^2
 ! fi  =  0.05
 ! ct  =  1.45e-10
 ! h   =  500  !m
 ! rw  =  0.10 !m
 ! re  =  500  !m
 ! nroots = 10


! p_i - initial pore pressure [Pa]
! q   - flow rate [m^3 / s]
! B   - formation volume factor [-]
! mu  - viscosity of the fluid [Pa.s]
! k   - permeability of the reservoir [m^2]
! h   - height of the reservoir [m]
! fi  - porosity of the reservoir [-]
! ct  - compressibility of the fluid [1/Pa]
! rw  - radius of the wellbore [m]
! re  - radius of the reservoir [m]
! a   - roots of the Eq 11 in Galis et al. 2017
!       [original equation: non-numbered equation after Eq 69 in PTT (2003)]
! n   - number of the roots (i.e., length of a array)
! note: all unis in SI

! constant P values for the discontinuity step of fluid rate
            
  if (time < pb%injection%t0) then
    P = 0
    dP_dt = 0
    !pb%injection%q = 0
  else       
      !k_min = 1d-18
      !k_max = 1d-14
      !FMCk = pb%injection%fi*pb%injection%muv*pb%injection%ct/k!pb%injection%kstar
      if (pb%mesh%dim == 1) then
        !r_p = (pb%injection%X0 + pb%injection%re)
        !r_n = (pb%injection%X0 - pb%injection%re)
        !r = x(r_n:r_p)
        !r = abs(x - pb%injection%X0)
         h = abs((x - 1000)/4) !* COS(pb%injection%theta)!(abs(x/4) - x/4)!/COS(pb%injection%theta) pb%injection%X0
         r = abs(x - pb%injection%X0) !h/COS(pb%injection%theta) !abs(x - 800)!- pb%injection%X0)
         hr = h/COS(pb%injection%theta)
         !r= max(r, 9.5)
      elseif (pb%mesh%dim == 2) then
        r = SQRT((x - pb%injection%X0)**2 + (y - pb%injection%Y0)**2)
      endif

      do i = 1, pb%mesh%nn
        if(r(i).lt.pb%injection%rw) r(i) = pb%injection%rw ! pressure inside wellbore = pressure at the wellbore-wall
      enddo
    
   !   do i = 1, pb%mesh%nn
   !     if(r(i).gt.pb%injection%re) pb%P = 0  ! pore pressure outside reservoir
   !   enddo
      r0 = SQRT(hr**2 - h**2) + 10!h/2 * SIN(pb%injection%theta)
      rD  = (r0/pb%injection%rw) 
      reD = pb%injection%re/pb%injection%rw
      rD2  = rD*rD
      reD2 = reD*reD

 ! Calcolo P alla fine del primo flow rate
      FMCk = pb%injection%fi*pb%injection%muv*pb%injection%ct/pb%injection%k
      qBmkh = (-0.002)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h) ! pb%injection%h
      tD  =  (37590912000. - pb%injection%t0) /(FMCk*pb%injection%rw*pb%injection%rw)
  
      pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

      sump = 0
      sumpt = 0
      do i= 1,int(pb%injection%nroots)
        A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
         * bessel_j0(pb%injection%an(i)*rD)
        B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
        C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*tD)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
        Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
        * EXP(-pb%injection%an(i)*pb%injection%an(i)*tD)*(bessel_j1(pb%injection%an(i)*reD))**2
        sumpt = sumpt + Ct;
        sump = sump + C;
      end do
      pD = pDpar + PI * sump
      P1 = (pb%injection%p_a - qBmkh*pD/(2*PI)) !from Galis eq 6

   ! Calcolo P alla fine del secondo flow rate

      qBmkh = (0.0000001)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)  ! pb%injection%h
      tD  =  (47304000000. - 37590912000.) /(FMCk*pb%injection%rw*pb%injection%rw)   ! Remember to change this line
      pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

      sump = 0
      sumpt = 0
      do i= 1,int(pb%injection%nroots)
        A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
         * bessel_j0(pb%injection%an(i)*rD)
        B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
        C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*tD)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
        Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
        * EXP(-pb%injection%an(i)*pb%injection%an(i)*tD)*(bessel_j1(pb%injection%an(i)*reD))**2
        sumpt = sumpt + Ct;
        sump = sump + C;
      end do
      pD = pDpar + PI * sump
      P2 = (pb%injection%p_a - qBmkh*pD/(2*PI)+P1) !from Galis eq 6

     ! do i = 3, time
     !   ks = k(i - 2) * ((-pb%v/pb%dc) * (k(i - 2) - pb%injection%k_max) &
     !          - (1/3153600.) * (k(i - 2 ) - pb%injection%k_min))
     !   k = pb%injection%k_min + (ks(i - 2) &
     !     - pb%injection%k_min)*EXP(-(pb%sigma - pb%P(i - 1))/max(pb%sigma,9.5))
     ! enddo


 ! Calcolo P alla fine del terzo flow rate

   !   qBmkh = (-0.0025)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
   !   pippo  = (189216000. - 157680000.) /(FMCk*pb%injection%rw*pb%injection%rw)
   !   pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

   !   sump = 0
   !   sumpt = 0
   !   do i= 1,int(pb%injection%nroots)
   !     A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
   !      * bessel_j0(pb%injection%an(i)*rD)
   !     B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
   !     C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
   !     Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
   !     * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
   !     sumpt = sumpt + Ct;
   !     sump = sump + C;
   !   end do
     


   !   pD = pDpar + PI * sump
   !   P3 = (pb%injection%p_i - qBmkh*pD/(2*PI)+P2) !from Galis eq 6

      !write(994,*) P1(512), P2(512),P3(512)
      
 ! Calcolo P alla fine del quarto flow rate

   !   qBmkh = (0.001)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
   !   pippo  = (252288000. - 189216000.) /(FMCk*pb%injection%rw*pb%injection%rw)
   !  pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

   !   sump = 0
   !   sumpt = 0
   !   do i= 1,int(pb%injection%nroots)
   !     A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
   !      * bessel_j0(pb%injection%an(i)*rD)
   !     B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
   !     C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
   !     Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
   !     * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
   !     sumpt = sumpt + Ct;
   !     sump = sump + C;
   !   end do
     


   !   pD = pDpar + PI * sump
   !   P4 = (pb%injection%p_i - qBmkh*pD/(2*PI)+P2) !from Galis eq 6

      !write(994,*) P1(512), P2(512),P3(512)      


 ! Calcolo P alla fine del quinto flow rate

   !   qBmkh = (-0.005)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
   !   pippo  = (283824000. - 252288000.) /(FMCk*pb%injection%rw*pb%injection%rw)
   !   pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

   !   sump = 0
   !   sumpt = 0
   !   do i= 1,int(pb%injection%nroots)
   !     A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
   !      * bessel_j0(pb%injection%an(i)*rD)
   !     B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
   !     C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
   !     Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
   !     * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
   !     sumpt = sumpt + Ct;
   !     sump = sump + C;
   !   end do
     


   !   pD = pDpar + PI * sump
   !   P5 = (pb%injection%p_i - qBmkh*pD/(2*PI)+P2) !from Galis eq 6

      !write(994,*) P1(512), P2(512),P3(512)  
      
 ! Calcolo P alla fine del quinto flow rate

   !   qBmkh = (0.001)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
   !   pippo  = (315360000. - 283824000.) /(FMCk*pb%injection%rw*pb%injection%rw)
   !   pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

   !   sump = 0
   !   sumpt = 0
   !   do i= 1,int(pb%injection%nroots)
   !     A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
   !      * bessel_j0(pb%injection%an(i)*rD)
   !     B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
   !     C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
   !     Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
   !     * EXP(-pb%injection%an(i)*pb%injection%an(i)*pippo)*(bessel_j1(pb%injection%an(i)*reD))**2
   !     sumpt = sumpt + Ct;
   !     sump = sump + C;
   !   end do
     


   !   pD = pDpar + PI * sump
   !   P6 = (pb%injection%p_i - qBmkh*pD/(2*PI)+P2) !from Galis eq 6

  !do i = 1, pb%mesh%nn  
  !  if(r(i).gt.pb%injection%re) pb%P = 0  ! pore pressure outside reservoir
  !enddo

 ! do i = 1, 630720000
 !   ks = k(i - 2) * ((-pb%v/pb%dc) * (k(i - 2) - pb%injection%k_max) &
 !              - (1/94608000.) * (k(i - 2 ) - pb%injection%k_min))
 !   k = pb%injection%k_min + (ks(i - 1) &
 !         - pb%injection%k_min)*EXP(-(pb%sigma - pb%P(i - 1))/max(pb%sigma,9.5))
 ! enddo
                      

      !write(994,*) P1(512), P2(512),P3(512)        
 ! if (time < pb%injection%t0) then
 !   P = 0
 !   dP_dt = 0
    !pb%injection%q = 0
 ! else

      !write(6,*) r

      !qBmkh = pb%injection%q*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
      !tD  =  (time - pb%injection%t0) /(FMCk*pb%injection%rw*pb%injection%rw)
      if ((time > pb%injection%t0) .and. (time <= 37590912000.)) then  ! Remember to make it time <= t1 
      qBmkh = (-0.002)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h) ! pb%injection%h   
      tD  =  (time - pb%injection%t0) /(FMCk*pb%injection%rw*pb%injection%rw)
      PAA=0.d0
   
      !pb%P = (pb%injection%p_i - qBmkh*pD/(2*PI))
      !write(994,*) time, qBmkh, tD*(FMCk*pb%injection%rw*pb%injection%rw), 1
      elseif ((time > 37590912000.) .and. (time <= 47304000000.)) then
      qBmkh = (0.0000001)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h) ! pb%injection%h
      tD  =  (time - 37590912000.) /(FMCk*pb%injection%rw*pb%injection%rw)
      PAA=P1
      
       !pb%P = (pb%injection%p_i - qBmkh*pD/(2*PI))      
      !write(994,*) time, qBmkh, tD*(FMCk*pb%injection%rw*pb%injection%rw), 2
    !  elseif ((time > 157680000.) .and. (time <= 189216000.)) then
    !  qBmkh = (-0.0025)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
    !  tD  =  (time - 157680000.) /(FMCk*pb%injection%rw*pb%injection%rw)
    !  PAA=P2
  
      !pb%P = (pb%injection%p_i - qBmkh*pD/(2*PI))      
      !write(994,*) time,qBmkh, tD*(FMCk*pb%injection%rw*pb%injection%rw), 3
    !  elseif ((time > 189216000.) .and. (time <= 252288000.)) then 
    !  qBmkh = (0.001)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
    !  tD  =  (time - 189216000.) /(FMCk*pb%injection%rw*pb%injection%rw)
    !   PAA=P3
      !pb%P = (pb%injection%p_i - qBmkh*pD/(2*PI))                             

      
    !  elseif ((time > 252288000.) .and. (time <= 283824000.)) then 
    !  qBmkh = (-0.005)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
    !  tD  =  (time - 252288000.) /(FMCk*pb%injection%rw*pb%injection%rw)
    !  PAA=P4
      !pb%P = (pb%injection%p_i - qBmkh*pD/(2*PI))                             
    
      
    !  elseif ((time > 126144000.) .and. (time <= 315360000.)) then 
    !  qBmkh = (0.001)*pb%injection%B0*pb%injection%muv/(pb%injection%k*pb%injection%h)
    !  tD  =  (time - 283824000.) /(FMCk*pb%injection%rw*pb%injection%rw)
    !  PAA=P1
      !pb%P = (pb%injection%p_i - qBmkh*pD/(2*PI))                             
      endif            

      pDpar = (2*(0.25*rD2+tD)-reD2*log(rD)-0.25*(3*reD2*reD2-4*reD2*reD2*log(reD)-2*reD2-1)/(reD2-1))/(reD2-1)

      sump = 0
      sumpt = 0
      do i= 1,int(pb%injection%nroots)
        A = bessel_j1(pb%injection%an(i)) * bessel_y0(pb%injection%an(i)*rD) - bessel_y1(pb%injection%an(i)) &
         * bessel_j0(pb%injection%an(i)*rD)
        B = pb%injection%an(i)*((bessel_j1(pb%injection%an(i)*reD))**2  - (bessel_j1(pb%injection%an(i)))**2)
        C = A/B * EXP(-pb%injection%an(i)*pb%injection%an(i)*tD)*(bessel_j1(pb%injection%an(i)*reD))**2
!    Ct = (-pb%injection%an(i)*pb%injection%an(i))*EXP(-pb%injection%an(i)*pb%injection%an(i)*tD) &
!    *(bessel_j1(pb%injection%an(i)*reD))**2*At/Bt
        Ct = A/B*(-pb%injection%an(i)*pb%injection%an(i))/(FMCk*pb%injection%rw*pb%injection%rw)  &
        * EXP(-pb%injection%an(i)*pb%injection%an(i)*tD)*(bessel_j1(pb%injection%an(i)*reD))**2
        sumpt = sumpt + Ct;
        sump = sump + C;
      end do
     
      !do i = 1, pb%mesh%nn
      !  if(r(i).gt.pb%injection%re) pb%P = 0  ! pore pressure outside reservoir
      !enddo

     ! do i = 1, pb%mesh%nn
     !   if (x(i) > r(i)) then
     !   pb%P = 0
     !   dP_dt = 0
     !   endif
     ! enddo
      

      pD = pDpar + PI * sump
  !pb%P = pb%injection%p_i - qBmkh*pD/(2*PI) Silvio edit
      P = (pb%injection%p_a - qBmkh*pD/(2*PI)+PAA) !from Galis eq 6
!  dP_dt = - qBmkh/(2*PI)*(2/(reD2-1)+PI*sumpt)
      dP_dt = - qBmkh/(2*PI)*(2/(FMCk*pb%injection%rw*pb%injection%rw*(reD2-1))+PI*sumpt)
! DA FARE PER SILVIO


     ! do i = 1, pb%mesh%nn
     ! if(r(i).gt.700) then !301.55
     !   P(i) = 0
     !   dP_dt(i) = 0
     ! end if
     ! enddo

      !do i = 1, pb%mesh%nn
      !if(r0(i).gt.pb%injection%re) then !301.55
      !  P(i) = 0
      !  dP_dt(i) = 0
      !end if
      !enddo

      do i = 1, pb%mesh%nn
      if(r(i).gt.389) then !301.55
        P(i) = 0
        dP_dt(i) = 0
      end if
      enddo     
    !  do i = 1, pb%mesh%nn
    !    if(r(i).gt.pb%injection%re) pb%P = 0  ! pore pressure outside reservoir
    !  enddo
  endif
  
!Calcolarsi dP_dt. NON può essere pb%injection%dP_dt perchè questa è costante (nel caso precedente)
!dP_dt = pb%injection%dP_dt questa formula valida solo quando P cambia in modo costante nel tempo.
end subroutine update_P_Galis


!===============================================================================
! SEISMIC: initialisation routine for the thermal pressurisation model
!===============================================================================
subroutine init_tp(pb)

  type(problem_type), intent(inout) :: pb

  call init_spectral_mesh(pb%tp%mesh)
  call init_variables(pb)
  call init_source(pb)

  ! Calculate remaining parameters (not necessarily time-constant)
  if (pb%i_rns_law == 3) then
    ! If using CNS, pass the initial porosity
    call calc_params(pb%theta, pb)
  else
    ! if using RSF, pass one
    call calc_params(pb%theta-pb%theta+1.0, pb)
  endif

end subroutine init_tp

!===============================================================================
! SEISMIC: initiate discretisation of the spectral domain for the spectral
! diffusion solver of [N&L]. Mesh parameters are hard-coded at the top of
! mesh.f90
!===============================================================================
subroutine init_spectral_mesh(sm)

  use constants, only : PI

  type(spectral_mesh_type), intent(inout) :: sm

  integer :: i
  double precision :: spi

  spi = sqrt(2.0/PI)

  allocate(sm%lw(sm%Nl), sm%F_inv(sm%Nl))

  do i=1,sm%Nl
    ! Construct logarithmic grid of dimensionless wavenumbers [N&L, Eqn. 14]
    sm%lw(i) = sm%lw_max * exp(-sm%Dlogl * (sm%Nl - i) )

    ! Construct kernel for inverse Fourier transform [N&L, Eqn. 17]
    if (i == 1) then
      sm%F_inv(i) = spi*sm%lw(i)*(1 + 0.5*sm%Dlogl)
    else if (i == sm%Nl) then
      sm%F_inv(i) = spi*sm%lw(i)*0.5*sm%Dlogl
    else
      sm%F_inv(i) = spi*sm%lw(i)*sm%Dlogl
    endif
  enddo

end subroutine init_spectral_mesh

!===============================================================================
! SEISMIC: initiate Gaussian distribution of heat source and dilatation in the
! discretised spectral domain
!===============================================================================
subroutine init_source(pb)

  use constants, only : PI

  type(problem_type), intent(inout) :: pb

  double precision :: spi = 1.0/sqrt(2*PI)
  integer :: i, j, n

  allocate (pb%tp%Omega(pb%mesh%nn*pb%tp%mesh%Nl))

  ! Loop over all fault segments
  do i=1,pb%mesh%nn
    ! Loop over all spectral elements
    do j=1,pb%tp%mesh%Nl
      n = (i-1)*pb%tp%mesh%Nl+j
      ! Fourier transform of Gaussian heat source [N&L, Eqn. 13]
      ! Note, Omega is pre-multiplied with w
      pb%tp%Omega(n) = pb%tp%w(i) * exp(-0.5*pb%tp%mesh%lw(j)**2)*spi
    enddo
  enddo

end subroutine init_source

!===============================================================================
! SEISMIC: initiate vectors of pressure and temperature (and spectral
! equivalents Pi and Theta), and vectors of previous values
!===============================================================================
subroutine init_variables(pb)

  type(problem_type), intent(inout) :: pb

  integer :: i

  ! Allocate variables
  allocate (  pb%tp%inv_w(pb%mesh%nn), &
              pb%tp%Pi(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%Theta(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%PiTheta(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%dP_dt(pb%mesh%nn), &
              pb%dtheta_dt(pb%mesh%nn), &
              pb%dtheta2_dt(pb%mesh%nn) )

  ! Allocate parameters
  allocate (  pb%tp%inv_rhoc(pb%mesh%nn), &
              pb%tp%alpha_th(pb%mesh%nn), &
              pb%tp%alpha_hy(pb%mesh%nn), &
              pb%tp%Lam(pb%mesh%nn), &
              pb%tp%Lam_prime(pb%mesh%nn), &
              pb%tp%Lam_T(pb%mesh%nn), &
              pb%tp%phi_b(pb%mesh%nn) )

  ! Allocate previous values
  allocate (  pb%tp%Theta_prev(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%PiTheta_prev(pb%mesh%nn*pb%tp%mesh%Nl) )

  do i=1,pb%mesh%nn
    pb%tp%inv_rhoc = 1.0/pb%tp%rhoc(i)
    pb%tp%inv_w(i) = 1.0/pb%tp%w(i)
    pb%P(i) = pb%tp%P_a(i)
    pb%T(i) = pb%tp%T_a(i)
    pb%tp%dP_dt(i) = 0d0
    pb%dtheta_dt(i) = 0d0
    pb%dtheta2_dt(i) = 0d0
  enddo

  do i=1,pb%mesh%nn*pb%tp%mesh%Nl
    pb%tp%Pi(i) = 0d0
    pb%tp%Theta(i) = 0d0
  enddo

end subroutine init_variables

!===============================================================================
! SEISMIC: calculate (time-variable) parameters. For the rate-and-state
! framework, phi = 1. For the CNS framework, phi is variable
!===============================================================================
subroutine calc_params(phi, pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: phi

  ! Thermal diffusivity
  pb%tp%alpha_th = pb%tp%k_t*pb%tp%inv_rhoc
  ! Hydraulic diffusivity
  pb%tp%alpha_hy = pb%tp%k_p/(pb%tp%eta*phi*pb%tp%beta)
  ! Thermal expansion moduli (by lack of a better term)
  pb%tp%Lam = pb%tp%l/pb%tp%beta
  pb%tp%Lam_prime = pb%tp%Lam * pb%tp%alpha_th/(pb%tp%alpha_hy - pb%tp%alpha_th)
  pb%tp%Lam_T = (pb%tp%Lam + pb%tp%Lam_prime)*pb%tp%inv_rhoc
  ! Specific storativity
  pb%tp%phi_b = pb%tp%dilat_factor/(phi*pb%tp%beta)

end subroutine calc_params

!===============================================================================
! SEISMIC: solve for P(t+dt) and T(t+dt) in the spatial domain
!===============================================================================
subroutine update_PT(tau_y,phi_dot,phi,v,dt,pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: tau_y, phi_dot, phi, v
  double precision :: dt

  ! Compute PiTheta and Theta for step t+dt
  call solve_spectral(tau_y, phi_dot, phi, v, dt, pb)

end subroutine update_PT

!===============================================================================
! SEISMIC: solve for P(t+dt) and T(t+dt) in the spatial domain for a full
! time step (i.e. not the intermediate RK solver steps) by a mid-point
! integration scheme
!===============================================================================
subroutine update_PT_final(dt,pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: tau_y, phi_dot, phi
  double precision :: dt

  ! Calculate mid-point values of tau_y, phi_dot, phi between t and t+dt
  tau_y = pb%tau * 0.5 * pb%v * pb%tp%inv_w
  if (pb%i_rns_law == 3) then
    ! CNS model: include porosity
    phi = pb%theta
    ! Multiply dilatancy rate with an arbitrary constant f >= 0 to control
    ! the amount of dilatancy hardening
    phi_dot = pb%dtheta_dt
  else
    ! RSF: ignore state
    phi = 1d0
    phi_dot = 0d0
  endif

  ! Compute PiTheta and Theta for step t+dt
  call solve_spectral(tau_y, phi_dot, phi, pb%v, dt, pb)

  ! Update initial values of Theta and PiTheta for next integration step
  pb%tp%Theta_prev = pb%tp%Theta
  pb%tp%PiTheta_prev = pb%tp%PiTheta

end subroutine update_PT_final

!===============================================================================
! SEISMIC: solve for P(t+dt) and T(t+dt) in the spectral domain
!===============================================================================
subroutine solve_spectral(tau_y,phi_dot,phi,v,dt,pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: tau_y, phi_dot, phi, v
  double precision :: dt, A_T, B_T, exp_T, A_P, B_P, exp_P
  double precision :: dPi, dTheta, dPiTheta
  integer :: i, j, n

  ! Update parameters in spatial domain
  call calc_params(phi, pb)

  ! Reset values of Theta and PiTheta. For each solver step, integration
  ! is done from t to t+dt, with variable dt. Hence, the initial values of
  ! Theta and PiTheta should always correspond to those at time t
  pb%tp%Theta = pb%tp%Theta_prev
  pb%tp%PiTheta = pb%tp%PiTheta_prev

  ! Loop over all fault segments
  do i=1,pb%mesh%nn
    pb%P(i) = pb%tp%P_a(i)
    pb%T(i) = pb%tp%T_a(i)
    pb%tp%dP_dt(i) = 0d0

    if (v(i) > 1e-5) then
      ! Loop over all spectral elements
      do j=1,pb%tp%mesh%Nl
        n = (i-1)*pb%tp%mesh%Nl+j

        ! Calculate F(t+dt) = B*(1-exp(-Adt))/A + F(t)*exp(-Adt)
        ! assuming constant A, B over the duration of dt [N&L, Eqn. 10]

        ! Temperature-related parameters in spectral domain
        A_T = pb%tp%alpha_th(i)*(pb%tp%mesh%lw(j)*pb%tp%inv_w(i))**2
        B_T = tau_y(i)*pb%tp%Omega(n)*pb%tp%inv_rhoc(i)
        exp_T = exp(-A_t*dt)
        ! Update Ttheta(t+dt)
        pb%tp%Theta(n) = B_T*(1.0 - exp_T)/A_T + pb%tp%Theta(n)*exp_T
        pb%T(i) = pb%T(i) + pb%tp%mesh%F_inv(j)*pb%tp%inv_w(i)*pb%tp%Theta(n)

        ! Pressure-related parameters in spectral domain
        A_P = pb%tp%alpha_hy(i)*(pb%tp%mesh%lw(j)*pb%tp%inv_w(i))**2
        B_P = ( pb%tp%Lam_T(i)*tau_y(i) - pb%tp%phi_b(i)*phi_dot(i) )*pb%tp%Omega(n)
        exp_P = exp(-A_P*dt)
        ! Update PiTtheta(t+dt)
        ! Note that PiTheta contains the spectral representation of
        ! Pi + Lambda_prime*Theta, where Pi is the Fourier transform of P
        ! and Theta is the Fourier transform of T [see N&L, Eqn. 5 and 7]
        pb%tp%PiTheta(n) = B_P*(1.0 - exp_P)/A_P + pb%tp%PiTheta(n)*exp_P
        pb%P(i) =  pb%P(i) + pb%tp%mesh%F_inv(j)*pb%tp%inv_w(i)* &
                      (pb%tp%PiTheta(n) - pb%tp%Lam_prime(i)*pb%tp%Theta(n))

        ! Update dP/dt
        dTheta = -A_T*pb%tp%Theta(n) + B_T
        dPiTheta = -A_P*pb%tp%PiTheta(n) + B_P
        dPi = dPiTheta - pb%tp%Lam_prime(i)*dTheta
        pb%tp%dP_dt(i) = pb%tp%dP_dt(i) + pb%tp%mesh%F_inv(j)*dPi*pb%tp%inv_w(i)
      enddo
    else
      do j=1,pb%tp%mesh%Nl
        n = (i-1)*pb%tp%mesh%Nl+j
        pb%tp%Theta(n) = 0d0
        pb%tp%PiTheta(n) = 0d0
      enddo
    endif
  enddo

end subroutine solve_spectral

end module diffusion_solver
