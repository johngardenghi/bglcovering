! *****************************************************************
! *****************************************************************

program algencama

  use bmalgencan, only: algencan
  use iso_c_binding, only: c_ptr, c_loc,c_f_pointer
  use modamerica
  use modcesaro
  use VorCells_Polygons

  implicit none

  ! PARAMETERS
  real(kind=8), parameter :: PI = acos( - 1.0d0 )
  ! real(kind=8), parameter :: sigma1 = 0.01d0
  ! real(kind=8), parameter :: sigma2 = 0.20d0

  ! PROBLEM DATA TYPE
  type :: pdata_type
     ! SCALARS
     integer :: nballs, npols, maxnvpols

     ! STATIC ARRAYS
     integer, dimension(5)      :: counters
     real(kind=8), dimension(2) :: dbl, dtr

     ! ALLOCATABLE ARRAYS
     integer, dimension(:), allocatable :: nvpols
     integer, dimension(:,:), allocatable :: poledges
     real(kind=8), dimension(:,:), allocatable :: xpol,ypol
  end type pdata_type

  ! LOCAL SCALARS
  character(len=10) :: problem
  integer           :: allocerr, hlnnzmax, i, ierr, istop, itrial, &
       jnnzmax, j, k, m, maxoutit, n, nextra, nincballs,           &
       nincballs_trial, nwcalls, nwtotit, outiter, p, totiter
  logical           :: corrin, exist, extallowed, rhoauto, &
       randomKappa, randomSigma, scale
  real(kind=4)      :: finish, start
  real(kind=8)      :: ang, bdsvio, bestrad, csupn, drand, epsfeas,  &
       epscompl, epsopt, f, gamma, intersection, kappa, maxdist,     &
       nlpsupn, R, Rleft, Rmin, Rright, Rtrial, rhoini, seed, ssupn, &
       theta, theta_trial, sigma1, sigma2

  ! LOCAL STRUCT
  type(pdata_type), target :: pdata

  ! LOCAL ARRAYS
  logical,      allocatable, dimension(:)   :: lind, uind
  real(kind=8),              dimension(2)   :: d, d_trial, pbl
  real(kind=8), allocatable, dimension(:)   :: c, distances, intersections, &
       lambda, lbnd, ubnd, x, xall

  ! COVERINGSECOND
  integer :: indpol, outside
  logical :: inside
  real(kind=8) :: rp(2)
  real(kind=8), allocatable, dimension(:,:) :: polpart

  write (*,*) "Problem (america | cesaro | star | stoyan | regpol):"
  read (*,*) problem

  write (*,*) 'Enter nballs (> 0): '
  read (*,*) pdata%nballs

  write (*,*) 'Enter itrial (> 0): '
  read (*,*) itrial

  ! Define the seed to random
  seed = 123456.0d0 * itrial

  write (*,*) 'Random kappa? (T|F): '
  read (*,*) randomKappa

  if (randomKappa) then
     kappa = 0.8d0 * drand(seed) + 0.1d0
  else
     write (*,*) 'Enter kappa (0..1): '
     read (*,*) kappa
  end if

  write (*,*) 'Random sigma1 and sigma2? (T|F): '
  read (*,*) randomSigma

  if (randomSigma) then
     sigma1 = 0.04d0 * drand(seed) + 0.01d0
     sigma2 = ( 0.2d0 - sigma1 ) * drand(seed) + sigma1
  else
     write (*,*) 'Enter sigma1: '
     read (*,*) sigma1

     write (*,*) 'Enter sigma2: '
     read (*,*) sigma2
  end if

  call build_polygon( problem, pdata )

  ! Number of variables (x = (c_1,...,c_nballs,r))
  n = 2 * pdata%nballs + 1

  ! Compute the initial radius and extremes of the bisection
  Rmin = sqrt( polygon_area( c_loc(pdata) ) / ( PI * pdata%nballs ) )
  Rleft  = 1.0d-01 * Rmin
  Rright = 1.0e+01 * Rmin

  ! Compute the bottom-left corner of A
  pbl(1) = huge( 1.0d0 )

  pbl(2) = huge( 1.0d0 )
  do i = 1, pdata%npols
     pbl(2) = min( pbl(2), minval( pdata%ypol(i, 1:pdata%nvpols(i)) ) )
  end do

  do i = 1, pdata%npols
     do j = 1, pdata%nvpols(i)
        if ( abs( pdata%ypol(i,j) - pbl(2) ) .le. sqrt( epsilon(1.0d0) ) ) then
           if ( pdata%xpol(i,j) .lt. pbl(1) ) then
              pbl(1) = pdata%xpol(i,j)
           end if
        end if
     end do
  end do

  write (*,*) "STARTING TO BUILD INITIAL POINT"
  write (*,*) "   problem ............: ", problem
  if ( problem .eq. "regpol" ) write (*,*) "   nvert ..............:", pdata%nvpols(1)-1
  write (*,*) "   nballs .............:", pdata%nballs
  write (*,*) "   itrial .............:", itrial
  write (*,*) "   kappa ..............:", kappa
  write (*,*) "   sigma1 .............:", sigma1
  write (*,*) "   sigma2 .............:", sigma2
  write (*,*) "   polygon area .......:", polygon_area( c_loc(pdata) )
  write (*,*) "   polygon perimeter ..:", polygon_perimeter( c_loc(pdata) )

  ! Bisection to find R
  R = 0.0d0
  do while ( Rright - Rleft .gt. 1.0d-02 * Rleft )
     Rtrial = ( Rleft + Rright ) / 2.0d0

     write (*,*)
     write (*,*) "Starting new bisection iteration"
     write (*,*) "   R = ", Rtrial

     if ( is_admissible( Rtrial, pbl, kappa, seed, c_loc(pdata), d_trial, theta_trial, nincballs_trial ) ) then
        write (*,*) "   Found admissible honeycomb"
        write (*,*) "      R     = ", Rtrial
        write (*,*) "      d     = ", d_trial
        write (*,*) "      theta = ", theta_trial

        Rleft = Rtrial

        ! Save information of a successful iteration of bisection
        if ( Rtrial .gt. R ) then
           write (*,*) "   Updating best admissible honeycomb"
           nincballs = nincballs_trial
           d(1:2) = d_trial(1:2)
           R = Rtrial
           theta = theta_trial
        end if
     else
        write (*,*) "   Not admissible honeycomb ... nincballs = ", nincballs_trial
        Rright = Rtrial
     end if
  end do

  write (*,*)
  write (*,*) "Finished radius adjustment"
  write (*,*) "   d .........: ", d(1:2)
  write (*,*) "   nincballs..: ", nincballs
  write (*,*) "   theta .....: ", theta, "rad (", 180.0d0 * theta / PI, " deg)"
  write (*,*) "   R .........: ", R

  ! Build the honeycomb
  nextra = 2*nincballs+1
  allocate( distances(nincballs), intersections(nincballs), &
       x(nextra), xall(nextra), stat=allocerr )
  if( allocerr .ne. 0 ) then
     write (*,*) "Memory allocation error"
     stop
  end if

  call build_honeycomb( R, d, theta, kappa, c_loc(pdata), x=xall )

  call drawsol( nextra, xall, "picture-solution-initial-nopert.mp", .true., c_loc(pdata), theta=theta )

  ! Compute a perturbation of the radius of the balls
  do i = 1, pdata%nballs
     distances(i) = dist_point_to_polygon_boundary( xall(2*i-1:2*i), c_loc( pdata ) )
  end do

  maxdist = maxval( distances(1:nincballs) )

  ! The perturbation is in 100 x sigma1 % (further from boundary) to
  ! 100 x sigma2 % (closer to boundary)
  do i = 1, nincballs
     gamma = sigma1 + ( sigma2 - sigma1 ) * ( 1.0d0 - distances(i) / max( 1.0d0, maxdist ) )

     xall(2*i-1) = xall(2*i-1) + &
          R * ( gamma * ( 2.0d0 * drand(seed) - 1.0d0 ) )
     xall(2*i) = xall(2*i) + &
          R * ( gamma * ( 2.0d0 * drand(seed) - 1.0d0 ) )
  end do

  ! Move extra balls to the end of x. Extra balls are those with
  ! smallest intersection with A.
  if ( nincballs .gt. pdata%nballs ) then
     do k = 1, nincballs
        intersection = ball_intersection_polygon( xall(2*k-1:2*k), R, pdata )

        do i = 1, k-1
           if ( intersection .ge. intersections(i) ) then
              do j = k, i+1, -1
                 intersections(j) = intersections(j-1)
                 x(2*j-1:2*j) = x(2*(j-1)-1:2*(j-1))
              end do
              exit
           end if
        end do

        x(2*i-1:2*i) = xall(2*k-1:2*k)
        intersections(i) = intersection
     end do
  else
     x(1:n) = xall(1:n)
  end if
  x(nextra) = R

  call drawsol( n, x, "picture-solution-initial.mp", .false., c_loc(pdata), distances=distances, nextra=nextra )

  x(n) = R
  ! -----------------------------------------------------------------------------------------
  ! COVERINGSECOND INITIAL GUESS
  ! -----------------------------------------------------------------------------------------
!   seed = 654321.0d0 * itrial
        
!   do i = 1,pdata%nballs
! 10   continue
     
!      rp(1) = pdata%dbl(1) + ( pdata%dtr(1) - pdata%dbl(1) ) * drand(seed)
!      rp(2) = pdata%dbl(2) + ( pdata%dtr(2) - pdata%dbl(2) ) * drand(seed)

!      outside = 0
!      do indpol=1,pdata%npols

!         allocate( polpart(2, pdata%nvpols(indpol)-1), stat=allocerr )

!         if ( allocerr .ne. 0 ) then
!            write(*,*) 'Allocation error.'
!            stop
!         end if

!         polpart(1,1:pdata%nvpols(indpol)-1) = pdata%xpol(indpol,1:pdata%nvpols(indpol)-1)
!         polpart(2,1:pdata%nvpols(indpol)-1) = pdata%ypol(indpol,1:pdata%nvpols(indpol)-1)
!         call polygon_contains_point_2d_convex(pdata%nvpols(indpol)-1,polpart,rp,inside)

!         deallocate( polpart, stat=allocerr )

!         if ( allocerr .ne. 0 ) then
!            write(*,*) 'Deallocation error.'
!            stop
!         end if

!         if ( inside ) then
!            exit
!         else
!            outside = outside + 1
!         end if

!      end do
!      if ( outside .eq. pdata%npols ) then
!         go to 10
!      end if
!      x(2*i-1) = rp(1)
!      x(2*i)   = rp(2)
!   end do

!   x(2*pdata%nballs+1) = ( 0.5d0 + 1.0d0 * drand(seed) ) / pdata%nballs

!   call drawsol( n, x, "picture-solution-initial.mp", .false., c_loc(pdata), distances=distances, nextra=nextra )
  ! -----------------------------------------------------------------------------------------


  ! Constraints
  m = 1
  p = 0

  allocate(lambda(m+p),c(m+p),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write (*,*) 'Allocation error.'
     stop
  end if

  ! Parameters setting
  epsfeas  = 1.0d-08
  epscompl = 1.0d-08
  epsopt   = 1.0d-08

  maxoutit = 50

  rhoauto = .true.

  if ( .not. rhoauto ) then
     rhoini = 1.0d-08
  end if

  scale = .false.
  extallowed = .true.
  corrin = .true.

  ! Set lower bounds and upper bounds
  allocate(lind(n), lbnd(n), uind(n), ubnd(n), stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write (*,*) 'Allocation error.'
     stop
  end if

  lind(1:2*pdata%nballs) = .false.
  uind(1:2*pdata%nballs+1) = .false.
  lind(n) = .true.
  lbnd(n) = 0.0d0

  ! Upper bounds on the number of sparse-matrices non-null elements
  jnnzmax = n
  hlnnzmax = huge( 1 ) / 100

  ! Initialize counters
  pdata%counters(1:5) = 0
  
  ! Initial guess and Lagrange multipliers
  lambda(1:m+p) = 0.0d0

  ! Optimize
  write (*,*)
  write (*,*) "Starting optimization"
  call cpu_time(start)
  call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax,           &
       n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
       scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
       outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))
  call cpu_time(finish)

  write (*,*) 
  write (*,*) "OPTIMIZATION SUMMARY"
  write (*,*) "R ...............: ", f
  write (*,*) "csupn ...........: ", csupn
  write (*,*) "istop ...........: ", istop
  write (*,*) "ierr ............: ", ierr
  write (*,*) "itrial ..........: ", itrial
  write (*,*) "#evalf ..........: ", pdata%counters(1)
  write (*,*) "#evalg ..........: ", pdata%counters(2)
  write (*,*) "#evalc ..........: ", pdata%counters(3)
  write (*,*) "#evalj ..........: ", pdata%counters(4)
  write (*,*) "#evalhl .........: ", pdata%counters(5)
  write (*,*) "outiter .........: ", outiter
  write (*,*) "totiter .........: ", totiter
  write (*,*) "nwcalls .........: ", nwcalls
  write (*,*) "nwtotit .........: ", nwtotit
  write (*,*) "Optimization time: ", finish - start

  if ( ierr .eq. 0 ) then
     inquire( file="bestrad.txt", exist=exist )

     if ( exist ) then
        open( 10, file="bestrad.txt" )
        read (10,*) bestrad
        close(10)
     else
        bestrad = huge( 1.0d0 )
     end if

     write (*,*)
     write (*,*) "R = ", f, "bestrad = ", bestrad
     write (*,*)

     if ( csupn .le. epsfeas .and. f .lt. bestrad * ( 1.0d0 - sqrt( epsilon(1.0d0) ) ) ) then
        open( 10, file="bestrad.txt" )
        write (10,*) f
        close(10)

        call drawsol( n, x, "picture-solution-final.mp", .false., c_loc(pdata) )

        open( 20, file="tabline.txt" )

        write (20,9000) n, kappa, sigma1, sigma2, f, csupn, ssupn,      &
             nlpsupn, bdsvio, outiter, totiter, nwcalls, nwtotit,       &
             pdata%counters(1), pdata%counters(2), pdata%counters(3),   &
             pdata%counters(4), pdata%counters(5), finish-start, istop, &
             ierr, itrial

        close(20)

        open( 20, file="output.dat" )
        write (20,*) "PROBLEM SUMMARY"
        write (20,'(A,1X,A)') "   problem .....:", problem
        if ( problem .eq. "regpol" ) write (20,'(A,1X,I0)') "   nvert .......:", pdata%nvpols(1)-1
        write (20,'(A,1X,I0)') "   nballs ......:", pdata%nballs
        write (20,'(A,1X,I0)') "   itrial ......:", itrial
        write (20,'(A,1X,0P,F4.2)') "   kappa .......:", kappa
        write (20,'(A,1X,0P,F4.2)') "   sigma1 ......:", sigma1
        write (20,'(A,1X,0P,F4.2)') "   sigma2 ......:", sigma2
        write (20,'(A,1X,1P,D22.16)') "   pol area ....:", polygon_area( c_loc(pdata) )
        write (20,'(A,1X,1P,D22.16)') "   pol perimeter:", polygon_perimeter( c_loc(pdata) )
        write (20,*)
        write (20,*) "OPTIMIZATION SUMMARY"
        write (20,'(A,1X,1P,D22.16)') "   r ......: ", f
        write (20,'(A,1X,1P,D22.16)') "   csupn ..: ", csupn
        write (20,'(A,1X,I0)') "   istop ..: ", istop
        write (20,'(A,1X,I0)') "   ierr ...: ", ierr
        write (20,'(A,1X,I0)') "   itrial .: ", itrial
        write (20,'(A,1X,I0)') "   #evalf .: ", pdata%counters(1)
        write (20,'(A,1X,I0)') "   #evalg .: ", pdata%counters(2)
        write (20,'(A,1X,I0)') "   #evalc .: ", pdata%counters(3)
        write (20,'(A,1X,I0)') "   #evalj .: ", pdata%counters(4)
        write (20,'(A,1X,I0)') "   #evalhl : ", pdata%counters(5)
        write (20,'(A,1X,I0)') "   outiter : ", outiter
        write (20,'(A,1X,I0)') "   totiter : ", totiter
        write (20,'(A,1X,I0)') "   nwcalls : ", nwcalls
        write (20,'(A,1X,I0)') "   nwtotit : ", nwtotit
        write (20,'(A,1X,0P,F4.2,1X,A)') "   opt time: ", finish - start, "seconds"
        write (20,*)
        write (20,*) "BALLS (index, center abcissa, center ordinate)"
        do i = 1, pdata%nballs
           write (20, '(I3,2(2X,1P,D23.16))') i, x(2*i-1), x(2*i)
        end do
     end if
  end if

  deallocate( c, distances, intersections, lambda, lbnd, lind,     &
       pdata%nvpols, pdata%poledges, pdata%xpol, pdata%ypol, ubnd, &
       uind, x, xall, stat=allocerr )

9000 format(1X,I6,3(0P,F7.2,1X),1P,D24.16,4(1X,1P,D7.1),4(1X,I8),5(1X,I7), &
       1X,0P,F6.2,1X,I5,1X,I4,1X,I6)

  contains

  ! *****************************************************************
  ! *****************************************************************

  logical function is_admissible(R, pbl, kappa, seed, pdataptr, d, &
       theta, nincballs)

    ! SCALAR ARGUMENTS
    integer,      intent(out) :: nincballs
    real(kind=8), intent(in)  :: kappa, R, seed
    real(kind=8), intent(out) :: theta
    type(c_ptr),  intent(in)  :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(2), intent(in)  :: pbl
    real(kind=8), dimension(2), intent(out) :: d
    
    ! PARAMETERS
    real(kind=8), parameter :: PI = acos( - 1.0d0 )

    ! LOCAL SCALARS
    integer                   :: trial
    real(kind=8)              :: drand
    type(pdata_type), pointer :: pdata

    call c_f_pointer( pdataptr, pdata )

    do trial = 1, 100
       ! Draw the angle of rotation of the honeycomb
       theta = PI * drand(seed)
       ! theta = PI * ( 0.1d0 * drand(seed) + 0.1d0 )
       ! theta = PI * 20.0d0 / 180.0d0

       ! Draw the displacement (inside the ball of radius R centered at
       ! the bottom-left of the polygon)
10     continue

       d(1) = ( 2.0d0 * drand(seed) - 1.0d0 ) * R + pbl(1)
       d(2) = ( 2.0d0 * drand(seed) - 1.0d0 ) * R + pbl(2)

       if ( sum( ( pbl - d )**2 ) .gt. R ** 2 ) go to 10

       call build_honeycomb( R, d, theta, kappa, pdataptr, nincballs )

       if ( nincballs .ge. pdata%nballs ) then
          is_admissible = .true.
          return
       end if
    end do

    is_admissible = .false.
  end function is_admissible

  ! *****************************************************************
  ! *****************************************************************

  subroutine build_honeycomb(R, d, theta, kappa, pdataptr, nincballs, x)

    ! SCALAR ARGUMENTS
    integer,      intent(out), optional :: nincballs
    real(kind=8), intent(in)  :: kappa, R
    real(kind=8), intent(in)  :: theta
    type(c_ptr),  intent(in)  :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(2), intent(in)            :: d
    real(kind=8), dimension(:), intent(out), optional :: x

    ! LOCAL SCALARS
    integer                   :: k, l, ninc
    real(kind=8)              :: a, b
    type(pdata_type), pointer :: pdata

    ! LOCAL ARRAY
    real(kind=8), dimension(2) :: bc

    call c_f_pointer( pdataptr, pdata )

    ! Build the honeycomb and check how many balls intersects A at
    ! least 100*kappa%
    ninc = 0
    do k = -50, 50
       do l = -50, 50
          ! Compute the center of a ball (honeycomb's hexagon)
          a = R / 2.0d0 * ( k + l ) * 3.0d0;
          b = R / 2.0d0 * ( k - l ) * sqrt(3.0d0);

          bc(1) = d(1) + ( cos(theta)*a - sin(theta)*b )
          bc(2) = d(2) + ( sin(theta)*a + cos(theta)*b )

          ! Include the ball (if the is space in x and if the hexagon
          ! intersects the polygon at least kappa)
          if ( ball_intersection_polygon( bc, R, pdata ) .ge. kappa ) then
             ninc = ninc + 1
             if ( present( x ) ) x(2*ninc-1:2*ninc) = bc(1:2)
          end if
       end do
    end do

    if ( present( nincballs ) ) then
       nincballs = ninc
    end if

    if ( present( x ) ) x(2*ninc+1) = R
  end subroutine build_honeycomb

  ! *****************************************************************
  ! *****************************************************************

  subroutine build_polygon(problem, pdata)
    ! SCALAR ARGUMENTS
    character(len=10), intent(in)    :: problem
    type(pdata_type),  intent(inout) :: pdata

    ! LOCAL PARAMETERS
    real(kind=8), parameter :: oct_circ_rad = 1.0d0 / (2.0d0 * dsin(PI / 8.0d0))
    real(kind=8), parameter :: dist_orig_to_bisec_segs = &
         sqrt( ( 0.5d0*(oct_circ_rad * dcos(PI / 4.0d0) + &
         oct_circ_rad * dcos(PI / 2.0d0)) )**2.0d0 + &
         ( 0.5d0*(oct_circ_rad * dsin(PI / 4.0d0) + oct_circ_rad * dsin(PI / 2.0d0)) )**2.0d0 )
    real(kind=8), parameter :: dist_orig_to_tipstar = dist_orig_to_bisec_segs + 2.0d0 * oct_circ_rad

    ! LOCAL SCALARS
    integer :: allocerr, i, iter, nvert, nsegs

    ! LOCAL ARRAYS
    real(kind=8) :: vxb(64), vyb(64)

    if ( problem .eq. "regpol" ) then
       write (*,*) "Number of vertices: "
       read (*,*) nvert

       pdata%npols = 1
       allocate( pdata%nvpols(pdata%npols), stat=allocerr )
       if ( allocerr .ne. 0 ) then
          write(*,*) 'Memory allocation error.'
          stop
       end if

       pdata%nvpols(1:pdata%npols) = (/ nvert + 1 /)
       
       allocate( pdata%poledges(1,nvert+1), pdata%xpol(1,nvert+1), &
            pdata%ypol(1,nvert+1), stat=allocerr )
       if ( allocerr .ne. 0 ) then
          write(*,*) 'Allocation error.'
          stop
       end if

       ! POLYGONS WITH SIDE 1
       ! if ( nvert .eq. 3 ) then
       !    pdata%xpol(1,1:nvert+1) = (/ 0.0d0, 1.0d0,                 0.5d0, 0.0d0 /)
       !    pdata%ypol(1,1:nvert+1) = (/ 0.0d0, 0.0d0, sqrt( 3.0d0 ) / 2.0d0, 0.0d0 /)

       ! elseif ( nvert .eq. 4 ) then
       !    pdata%xpol(1,1:nvert+1) = (/ 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0 /)
       !    pdata%ypol(1,1:nvert+1) = (/ 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0 /)

       ! POLYGONS INSIDE THE CIRCLE OF RADIUS 1
       ! else
          do i = 1,nvert + 1
             ang = ( 2.0d0 * PI * ( i - 1 ) / nvert ) + PI/2.0d0 * modulo( nvert, 2 )
             pdata%xpol(1,i) = cos( ang )
             pdata%ypol(1,i) = sin( ang )
          end do
       ! end if

       pdata%dbl(1:2) = (/ minval(pdata%xpol(1,1:pdata%nvpols(1))), minval(pdata%ypol(1,1:pdata%nvpols(1))) /)
       pdata%dtr(1:2) = (/ maxval(pdata%xpol(1,1:pdata%nvpols(1))), maxval(pdata%ypol(1,1:pdata%nvpols(1))) /)

       pdata%poledges(1,1:nvert) = 0

    elseif ( problem .eq. "america" ) then
       ! Convex polygons
       pdata%npols = 34
       pdata%maxnvpols = 12

       put(1:2, 1:np) = put(1:2, 1:np) / 20.0d0

       allocate( pdata%nvpols(pdata%npols), pdata%poledges(pdata%npols, pdata%maxnvpols), &
            pdata%xpol(pdata%npols, pdata%maxnvpols), pdata%ypol(pdata%npols, pdata%maxnvpols), stat=allocerr )

       ! Number of vertices for each polygon
       pdata%nvpols(1:pdata%npols) = (/ 7, 5, 4, 8, 5, 4, 4, 5, 6, 6, 4, 12, 4, 5, 5, 4, 6, 7, 6, 6, 5, 5, 6, 6, 5, 9, &
            5, 7, 7, 4, 7, 4, 5, 6 /)

       ! Vertices
       pdata%xpol(1, 1:pdata%nvpols(1)) = (/ put(1,80), put(1,44), put(1,78), put(1,124), put(1,45), put(1,79), put(1,80) /)
       pdata%ypol(1, 1:pdata%nvpols(1)) = (/ put(2,80), put(2,44), put(2,78), put(2,124), put(2,45), put(2,79), put(2,80) /)

       pdata%xpol(2, 1:pdata%nvpols(2)) = (/ put(1,80), put(1,79), put(1,46), put(1,43), put(1,80) /)
       pdata%ypol(2, 1:pdata%nvpols(2)) = (/ put(2,80), put(2,79), put(2,46), put(2,43), put(2,80) /)

       pdata%xpol(3, 1:pdata%nvpols(3)) = (/ put(1,83), put(1,35), put(1,64), put(1,83) /)
       pdata%ypol(3, 1:pdata%nvpols(3)) = (/ put(2,83), put(2,35), put(2,64), put(2,83) /)

       pdata%xpol(4, 1:pdata%nvpols(4)) = (/ put(1,46), put(1,83),  &
            put(1,64), put(1,102), put(1,81), put(1,82), put(1,43), &
            put(1,46) /)
       pdata%ypol(4, 1:pdata%nvpols(4)) = (/ put(2,46), put(2,83),  &
            put(2,64), put(2,102), put(2,81), put(2,82), put(2,43), &
            put(2,46) /)

       pdata%xpol(5, 1:pdata%nvpols(5)) = (/ put(1,64), put(1,63), put(1,101), put(1,102), put(1,64) /)
       pdata%ypol(5, 1:pdata%nvpols(5)) = (/ put(2,64), put(2,63), put(2,101), put(2,102), put(2,64) /)

       pdata%xpol(6, 1:pdata%nvpols(6)) = (/ put(1,102), put(1,40), put(1,81), put(1,102) /)
       pdata%ypol(6, 1:pdata%nvpols(6)) = (/ put(2,102), put(2,40), put(2,81), put(2,102) /)

       pdata%xpol(7, 1:pdata%nvpols(7)) = (/ put(1,40), put(1,127), put(1,81), put(1,40) /)
       pdata%ypol(7, 1:pdata%nvpols(7)) = (/ put(2,40), put(2,127), put(2,81), put(2,40) /)

       pdata%xpol(8, 1:pdata%nvpols(8)) = (/ put(1,40), put(1,41), put(1,42), put(1,127), put(1,40) /)
       pdata%ypol(8, 1:pdata%nvpols(8)) = (/ put(2,40), put(2,41), put(2,42), put(2,127), put(2,40) /)

       pdata%xpol(9, 1:pdata%nvpols(9)) = (/ put(1,63), put(1,97), put(1,36), put(1,100), put(1,101), put(1,63) /)
       pdata%ypol(9, 1:pdata%nvpols(9)) = (/ put(2,63), put(2,97), put(2,36), put(2,100), put(2,101), put(2,63) /)

       pdata%xpol(10, 1:pdata%nvpols(10)) = (/ put(1,100), put(1,37), put(1,38), put(1,39), put(1,101), put(1,100) /)
       pdata%ypol(10, 1:pdata%nvpols(10)) = (/ put(2,100), put(2,37), put(2,38), put(2,39), put(2,101), put(2,100) /)

       pdata%xpol(11, 1:pdata%nvpols(11)) = (/ put(1,36), put(1,99), put(1,100), put(1,36) /)
       pdata%ypol(11, 1:pdata%nvpols(11)) = (/ put(2,36), put(2,99), put(2,100), put(2,36) /)

       pdata%xpol(12, 1:pdata%nvpols(12)) = (/ put(1,35), put(1,120), put(1,56), put(1,121), put(1,33), put(1,105), put(1,106),&
            put(1,96), put(1,97), put(1,63), put(1,64), put(1,35) /)
       pdata%ypol(12, 1:pdata%nvpols(12)) = (/ put(2,35), put(2,120), put(2,56), put(2,121), put(2,33), put(2,105), put(2,106),&
            put(2,96), put(2,97), put(2,63), put(2,64), put(2,35) /)

       pdata%xpol(13, 1:pdata%nvpols(13)) = (/ put(1,106), put(1,60), put(1,96), put(1,106) /)
       pdata%ypol(13, 1:pdata%nvpols(13)) = (/ put(2,106), put(2,60), put(2,96), put(2,106) /)

       pdata%xpol(14, 1:pdata%nvpols(14)) = (/ put(1,62), put(1,98), put(1,61), put(1,60), put(1,62) /)
       pdata%ypol(14, 1:pdata%nvpols(14)) = (/ put(2,62), put(2,98), put(2,61), put(2,60), put(2,62) /)

       pdata%xpol(15, 1:pdata%nvpols(15)) = (/ put(1,106), put(1,94), put(1,62), put(1,60), put(1,106) /)
       pdata%ypol(15, 1:pdata%nvpols(15)) = (/ put(2,106), put(2,94), put(2,62), put(2,60), put(2,106) /)

       pdata%xpol(16, 1:pdata%nvpols(16)) = (/ put(1,106), put(1,34), put(1,94), put(1,106) /)
       pdata%ypol(16, 1:pdata%nvpols(16)) = (/ put(2,106), put(2,34), put(2,94), put(2,106) /)

       pdata%xpol(17, 1:pdata%nvpols(17)) = (/ put(1,106), put(1,117), put(1,107), put(1,85), put(1,34), put(1,106) /)
       pdata%ypol(17, 1:pdata%nvpols(17)) = (/ put(2,106), put(2,117), put(2,107), put(2,85), put(2,34), put(2,106) /)

       pdata%xpol(18, 1:pdata%nvpols(18)) = (/ put(1,107), put(1,86), put(1,113), put(1,114), put(1,89), put(1,85), put(1,107) /)
       pdata%ypol(18, 1:pdata%nvpols(18)) = (/ put(2,107), put(2,86), put(2,113), put(2,114), put(2,89), put(2,85), put(2,107) /)

       pdata%xpol(19, 1:pdata%nvpols(19)) = (/ put(1,114), put(1,31), put(1,32), put(1,87), put(1,89), put(1,114) /)
       pdata%ypol(19, 1:pdata%nvpols(19)) = (/ put(2,114), put(2,31), put(2,32), put(2,87), put(2,89), put(2,114) /)

       pdata%xpol(20, 1:pdata%nvpols(20)) = (/ put(1,87), put(1,115), put(1,88), put(1,103), put(1,89), put(1,87) /)
       pdata%ypol(20, 1:pdata%nvpols(20)) = (/ put(2,87), put(2,115), put(2,88), put(2,103), put(2,89), put(2,87) /)

       pdata%xpol(21, 1:pdata%nvpols(21)) = (/ put(1,33), put(1,118), put(1,119), put(1,105), put(1,33) /)
       pdata%ypol(21, 1:pdata%nvpols(21)) = (/ put(2,33), put(2,118), put(2,119), put(2,105), put(2,33) /)

       pdata%xpol(22, 1:pdata%nvpols(22)) = (/ put(1,118), put(1,104), put(1,116), put(1,119), put(1,118) /)
       pdata%ypol(22, 1:pdata%nvpols(22)) = (/ put(2,118), put(2,104), put(2,116), put(2,119), put(2,118) /)

       pdata%xpol(23, 1:pdata%nvpols(23)) = (/ put(1,31), put(1,92), put(1,91), put(1,90), put(1,32), put(1,31) /)
       pdata%ypol(23, 1:pdata%nvpols(23)) = (/ put(2,31), put(2,92), put(2,91), put(2,90), put(2,32), put(2,31) /)

       pdata%xpol(24, 1:pdata%nvpols(24)) = (/ put(1,92), put(1,93), put(1,30), put(1,29), put(1,91), put(1,92) /)
       pdata%ypol(24, 1:pdata%nvpols(24)) = (/ put(2,92), put(2,93), put(2,30), put(2,29), put(2,91), put(2,92) /)

       pdata%xpol(25, 1:pdata%nvpols(25)) = (/ put(1,49), put(1,50), put(1,51), put(1,52), put(1,49) /)
       pdata%ypol(25, 1:pdata%nvpols(25)) = (/ put(2,49), put(2,50), put(2,51), put(2,52), put(2,49) /)

       pdata%xpol(26, 1:pdata%nvpols(26)) = (/ put(1,28), put(1,27), put(1,128), put(1,59), put(1,68), put(1,76), put(1,48),&
            put(1,129), put(1,28) /)
       pdata%ypol(26, 1:pdata%nvpols(26)) = (/ put(2,28), put(2,27), put(2,128), put(2,59), put(2,68), put(2,76), put(2,48),&
            put(2,129), put(2,28) /)

       pdata%xpol(27, 1:pdata%nvpols(27)) = (/ put(1,73), put(1,77), put(1,71), put(1,72), put(1,73) /)
       pdata%ypol(27, 1:pdata%nvpols(27)) = (/ put(2,73), put(2,77), put(2,71), put(2,72), put(2,73) /)

       pdata%xpol(28, 1:pdata%nvpols(28)) = (/ put(1,129), put(1,47), put(1,125), put(1,29), put(1,30), put(1,28), put(1,129) /)
       pdata%ypol(28, 1:pdata%nvpols(28)) = (/ put(2,129), put(2,47), put(2,125), put(2,29), put(2,30), put(2,28), put(2,129) /)

       pdata%xpol(29, 1:pdata%nvpols(29)) = (/ put(1,125), put(1,17), put(1,19), put(1,123), put(1,24), put(1,29), put(1,125) /)
       pdata%ypol(29, 1:pdata%nvpols(29)) = (/ put(2,125), put(2,17), put(2,19), put(2,123), put(2,24), put(2,29), put(2,125) /)

       pdata%xpol(30, 1:pdata%nvpols(30)) = (/ put(1,19), put(1,122), put(1,123), put(1,19) /)
       pdata%ypol(30, 1:pdata%nvpols(30)) = (/ put(2,19), put(2,122), put(2,123), put(2,19) /)

       pdata%xpol(31, 1:pdata%nvpols(31)) = (/ put(1,59), put(1,130), put(1,57), put(1,58), put(1,2), put(1,69), put(1,59) /)
       pdata%ypol(31, 1:pdata%nvpols(31)) = (/ put(2,59), put(2,130), put(2,57), put(2,58), put(2,2), put(2,69), put(2,59) /)
       
       pdata%xpol(32, 1:pdata%nvpols(32)) = (/ put(1,69), put(1,67), put(1,59), put(1,69) /)
       pdata%ypol(32, 1:pdata%nvpols(32)) = (/ put(2,69), put(2,67), put(2,59), put(2,69) /)

       pdata%xpol(33, 1:pdata%nvpols(33)) = (/ put(1,67), put(1,66), put(1,4), put(1,59), put(1,67) /)
       pdata%ypol(33, 1:pdata%nvpols(33)) = (/ put(2,67), put(2,66), put(2,4), put(2,59), put(2,67) /)

       pdata%xpol(34, 1:pdata%nvpols(34)) = (/ put(1,4), put(1,3), put(1,5), put(1,68), put(1,59), put(1,4) /)
       pdata%ypol(34, 1:pdata%nvpols(34)) = (/ put(2,4), put(2,3), put(2,5), put(2,68), put(2,59), put(2,4) /)

       ! Edges ID
       pdata%poledges(1, 1:pdata%nvpols(1) - 1) = (/ 0, 0, 0, 0, 0, 1 /)
       pdata%poledges(2, 1:pdata%nvpols(2) - 1) = (/ 1, 0, 1, 0 /)
       pdata%poledges(3, 1:pdata%nvpols(3) - 1) = (/ 0, 1, 1 /)
       pdata%poledges(4, 1:pdata%nvpols(4) - 1) = (/ 0, 1, 1, 1, 0, 0, 1 /)
       pdata%poledges(5, 1:pdata%nvpols(5) - 1) = (/ 1, 1, 0, 1 /)
       pdata%poledges(6, 1:pdata%nvpols(6) - 1) = (/ 0, 1, 1 /)
       pdata%poledges(7, 1:pdata%nvpols(7) - 1) = (/ 1, 0, 1 /)
       pdata%poledges(8, 1:pdata%nvpols(8) - 1) = (/ 0, 0, 0, 1 /)
       pdata%poledges(9, 1:pdata%nvpols(9) - 1) = (/ 1, 0, 1, 1, 1 /)
       pdata%poledges(10, 1:pdata%nvpols(10) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(11, 1:pdata%nvpols(11) - 1) = (/ 0, 0, 1 /)
       pdata%poledges(12, 1:pdata%nvpols(12) - 1) = (/ 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1 /)
       pdata%poledges(13, 1:pdata%nvpols(13) - 1) = (/ 1, 0, 1 /)
       pdata%poledges(14, 1:pdata%nvpols(14) - 1) = (/ 0, 0, 0, 1 /)
       pdata%poledges(15, 1:pdata%nvpols(15) - 1) = (/ 1, 0, 1, 1 /)
       pdata%poledges(16, 1:pdata%nvpols(16) - 1) = (/ 1, 0, 1 /)
       pdata%poledges(17, 1:pdata%nvpols(17) - 1) = (/ 0, 0, 1, 0, 1 /)
       pdata%poledges(18, 1:pdata%nvpols(18) - 1) = (/ 0, 0, 0, 1, 0, 1 /)
       pdata%poledges(19, 1:pdata%nvpols(19) - 1) = (/ 0, 1, 0, 1, 1 /)
       pdata%poledges(20, 1:pdata%nvpols(20) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(21, 1:pdata%nvpols(21) - 1) = (/ 0, 1, 0, 1 /)
       pdata%poledges(22, 1:pdata%nvpols(22) - 1) = (/ 0, 0, 0, 1 /)
       pdata%poledges(23, 1:pdata%nvpols(23) - 1) = (/ 0, 1, 0, 0, 1 /)
       pdata%poledges(24, 1:pdata%nvpols(24) - 1) = (/ 0, 0, 1, 0, 1 /)
       pdata%poledges(25, 1:pdata%nvpols(25) - 1) = (/ 0, 0, 0, 0 /)
       pdata%poledges(26, 1:pdata%nvpols(26) - 1) = (/ 0, 0, 0, 1, 0, 0, 0, 1 /)
       pdata%poledges(27, 1:pdata%nvpols(27) - 1) = (/ 0, 0, 0, 0 /)
       pdata%poledges(28, 1:pdata%nvpols(28) - 1) = (/ 0, 0, 1, 1, 0, 1 /)
       pdata%poledges(29, 1:pdata%nvpols(29) - 1) = (/ 0, 0, 1, 0, 0, 1 /)
       pdata%poledges(30, 1:pdata%nvpols(30) - 1) = (/ 0, 0, 1 /)
       pdata%poledges(31, 1:pdata%nvpols(31) - 1) = (/ 0, 0, 0, 0, 0, 1 /)
       pdata%poledges(32, 1:pdata%nvpols(32) - 1) = (/ 0, 1, 1 /)
       pdata%poledges(33, 1:pdata%nvpols(33) - 1) = (/ 0, 0, 1, 1 /)
       pdata%poledges(34, 1:pdata%nvpols(34) - 1) = (/ 0, 0, 0, 1, 1 /)

       ! Rectangle containing A
       pdata%dbl(1) = minval(put(1,1:np)) - 1.0d0
       pdata%dbl(2) = minval(put(2,1:np)) - 1.0d0
       pdata%dtr(1) = maxval(put(1,1:np)) + 1.0d0
       pdata%dtr(2) = maxval(put(2,1:np)) + 1.0d0

    elseif ( problem .eq. "cesaro" ) then
       nsegs = 4
       iter = 1

       vxb(1:nsegs) = (/0.0d0, 1.0d0, 1.0d0, 0.0d0 /)
       vxb(5:64) = 0.0d0
       vyb(1:nsegs) = (/ 0.0d0, 0.0d0, 1.0d0, 1.0d0 /)
       vyb(5:64) = 0.0d0

       call includes_triang_segs(vxb,vyb,nsegs,iter)

       ! Problem data:
       ! Convex polygons
       pdata%npols = 21
       pdata%maxnvpols = 7
       allocate( pdata%nvpols(pdata%npols), pdata%poledges(pdata%npols, pdata%maxnvpols), &
            pdata%xpol(pdata%npols, pdata%maxnvpols), pdata%ypol(pdata%npols, pdata%maxnvpols), stat=allocerr )

       ! Number of vertices for each polygon (last=first)
       pdata%nvpols(1:pdata%npols) = (/ 6, 6, 6, 5, 7, 6, 6, 6, 5, 7, 6, 6, 6, 5, 7, 6, 6, 6, 5, 7, 5 /)

       ! Vertices
       pdata%xpol(1, 1:pdata%nvpols(1)) = (/ vxb(1), vxb(2), vxb(3), vxb(63), vxb(64), vxb(1) /)
       pdata%ypol(1, 1:pdata%nvpols(1)) = (/ vyb(1), vyb(2), vyb(3), vyb(63), vyb(64), vyb(1) /)

       pdata%xpol(2, 1:pdata%nvpols(2)) = (/ vxb(3), vxb(4), vxb(5), vxb(6), vxb(7), vxb(3) /)
       pdata%ypol(2, 1:pdata%nvpols(2)) = (/ vyb(3), vyb(4), vyb(5), vyb(6), vyb(7), vyb(3) /)

       pdata%xpol(3, 1:pdata%nvpols(3)) = (/ vxb(59), vxb(60), vxb(61), vxb(62), vxb(63), vxb(59) /)
       pdata%ypol(3, 1:pdata%nvpols(3)) = (/ vyb(59), vyb(60), vyb(61), vyb(62), vyb(63), vyb(59) /)

       pdata%xpol(4, 1:pdata%nvpols(4)) = (/ vxb(63), vxb(3), vxb(7), vxb(59), vxb(63) /)
       pdata%ypol(4, 1:pdata%nvpols(4)) = (/ vyb(63), vyb(3), vyb(7), vyb(59), vyb(63) /)

       pdata%xpol(5, 1:pdata%nvpols(5)) = (/ vxb(59), vxb(7), vxb(8), vxb(9), vxb(57), vxb(58), vxb(59) /)
       pdata%ypol(5, 1:pdata%nvpols(5)) = (/ vyb(59), vyb(7), vyb(8), vyb(9), vyb(57), vyb(58), vyb(59) /)

       pdata%xpol(6, 1:pdata%nvpols(6)) = (/ vxb(11), vxb(12), vxb(13), vxb(14), vxb(15), vxb(11) /)
       pdata%ypol(6, 1:pdata%nvpols(6)) = (/ vyb(11), vyb(12), vyb(13), vyb(14), vyb(15), vyb(11) /)

       pdata%xpol(7, 1:pdata%nvpols(7)) = (/ vxb(15), vxb(16), vxb(17), vxb(18), vxb(19), vxb(15) /)
       pdata%ypol(7, 1:pdata%nvpols(7)) = (/ vyb(15), vyb(16), vyb(17), vyb(18), vyb(19), vyb(15) /)

       pdata%xpol(8, 1:pdata%nvpols(8)) = (/ vxb(19), vxb(20), vxb(21), vxb(22), vxb(23), vxb(19) /)
       pdata%ypol(8, 1:pdata%nvpols(8)) = (/ vyb(19), vyb(20), vyb(21), vyb(22), vyb(23), vyb(19) /)

       pdata%xpol(9, 1:pdata%nvpols(9)) = (/ vxb(19), vxb(23), vxb(11), vxb(15), vxb(19) /)
       pdata%ypol(9, 1:pdata%nvpols(9)) = (/ vyb(19), vyb(23), vyb(11), vyb(15), vyb(19) /)

       pdata%xpol(10, 1:pdata%nvpols(10)) = (/ vxb(9), vxb(10), vxb(11), vxb(23), vxb(24), vxb(25), vxb(9) /)
       pdata%ypol(10, 1:pdata%nvpols(10)) = (/ vyb(9), vyb(10), vyb(11), vyb(23), vyb(24), vyb(25), vyb(9) /)

       pdata%xpol(11, 1:pdata%nvpols(11)) = (/ vxb(27), vxb(28), vxb(29), vxb(30), vxb(31), vxb(27) /)
       pdata%ypol(11, 1:pdata%nvpols(11)) = (/ vyb(27), vyb(28), vyb(29), vyb(30), vyb(31), vyb(27) /)

       pdata%xpol(12, 1:pdata%nvpols(12)) = (/ vxb(31), vxb(32), vxb(33), vxb(34), vxb(35), vxb(31) /)
       pdata%ypol(12, 1:pdata%nvpols(12)) = (/ vyb(31), vyb(32), vyb(33), vyb(34), vyb(35), vyb(31) /)

       pdata%xpol(13, 1:pdata%nvpols(13)) = (/ vxb(35), vxb(36), vxb(37), vxb(38), vxb(39), vxb(35) /)
       pdata%ypol(13, 1:pdata%nvpols(13)) = (/ vyb(35), vyb(36), vyb(37), vyb(38), vyb(39), vyb(35) /)

       pdata%xpol(14, 1:pdata%nvpols(14)) = (/ vxb(35), vxb(39), vxb(27), vxb(31), vxb(35) /)
       pdata%ypol(14, 1:pdata%nvpols(14)) = (/ vyb(35), vyb(39), vyb(27), vyb(31), vyb(35) /)

       pdata%xpol(15, 1:pdata%nvpols(15)) = (/ vxb(27), vxb(39), vxb(40), vxb(41), vxb(25), vxb(26), vxb(27) /)
       pdata%ypol(15, 1:pdata%nvpols(15)) = (/ vyb(27), vyb(39), vyb(40), vyb(41), vyb(25), vyb(26), vyb(27) /)

       pdata%xpol(16, 1:pdata%nvpols(16)) = (/ vxb(43), vxb(44), vxb(45), vxb(46), vxb(47), vxb(43) /)
       pdata%ypol(16, 1:pdata%nvpols(16)) = (/ vyb(43), vyb(44), vyb(45), vyb(46), vyb(47), vyb(43) /)

       pdata%xpol(17, 1:pdata%nvpols(17)) = (/ vxb(47), vxb(48), vxb(49), vxb(50), vxb(51), vxb(47) /)
       pdata%ypol(17, 1:pdata%nvpols(17)) = (/ vyb(47), vyb(48), vyb(49), vyb(50), vyb(51), vyb(47) /)

       pdata%xpol(18, 1:pdata%nvpols(18)) = (/ vxb(51), vxb(52), vxb(53), vxb(54), vxb(55), vxb(51) /)
       pdata%ypol(18, 1:pdata%nvpols(18)) = (/ vyb(51), vyb(52), vyb(53), vyb(54), vyb(55), vyb(51) /)

       pdata%xpol(19, 1:pdata%nvpols(19)) = (/ vxb(51), vxb(55), vxb(43), vxb(47), vxb(51) /)
       pdata%ypol(19, 1:pdata%nvpols(19)) = (/ vyb(51), vyb(55), vyb(43), vyb(47), vyb(51) /)

       pdata%xpol(20, 1:pdata%nvpols(20)) = (/ vxb(43), vxb(55), vxb(56), vxb(57), vxb(41), vxb(42), vxb(43) /)
       pdata%ypol(20, 1:pdata%nvpols(20)) = (/ vyb(43), vyb(55), vyb(56), vyb(57), vyb(41), vyb(42), vyb(43) /)

       pdata%xpol(21, 1:pdata%nvpols(21)) = (/ vxb(41), vxb(57), vxb(9), vxb(25), vxb(41) /)
       pdata%ypol(21, 1:pdata%nvpols(21)) = (/ vyb(41), vyb(57), vyb(9), vyb(25), vyb(41) /)

       ! Edges ID
       pdata%poledges(1, 1:pdata%nvpols(1) - 1) = (/ 0, 0, 1, 0, 0 /)
       pdata%poledges(2, 1:pdata%nvpols(2) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(3, 1:pdata%nvpols(3) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(4, 1:pdata%nvpols(4) - 1) = (/ 1, 1, 1, 1 /)
       pdata%poledges(5, 1:pdata%nvpols(5) - 1) = (/ 1, 0, 0, 1, 0, 0 /)
       pdata%poledges(6, 1:pdata%nvpols(6) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(7, 1:pdata%nvpols(7) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(8, 1:pdata%nvpols(8) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(9, 1:pdata%nvpols(9) - 1) = (/ 1, 1, 1, 1 /)
       pdata%poledges(10, 1:pdata%nvpols(10) - 1) = (/ 0, 0, 1, 0, 0, 1 /)
       pdata%poledges(11, 1:pdata%nvpols(11) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(12, 1:pdata%nvpols(12) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(13, 1:pdata%nvpols(13) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(14, 1:pdata%nvpols(14) - 1) = (/ 1, 1, 1, 1 /)
       pdata%poledges(15, 1:pdata%nvpols(15) - 1) = (/ 1, 0, 0, 1, 0, 0 /)
       pdata%poledges(16, 1:pdata%nvpols(16) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(17, 1:pdata%nvpols(17) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(18, 1:pdata%nvpols(18) - 1) = (/ 0, 0, 0, 0, 1 /)
       pdata%poledges(19, 1:pdata%nvpols(19) - 1) = (/ 1, 1, 1, 1 /)
       pdata%poledges(20, 1:pdata%nvpols(20) - 1) = (/ 1, 0, 0, 1, 0, 0 /)
       pdata%poledges(21, 1:pdata%nvpols(21) - 1) = (/ 1, 1, 1, 1 /)

       ! Rectangle containing A
       pdata%dbl(1:2) = (/ -sqrt(2.0d0) / 2.0d0, -sqrt(2.0d0) / 2.0d0 /)
       pdata%dtr(1:2) = (/ 1.0d0 + sqrt(2.0d0) / 2.0d0, 1.0d0 + sqrt(2.0d0) / 2.0d0 /)
       
    elseif ( problem .eq. "star" ) then
       ! Convex polygons
       pdata%npols = 9
       pdata%maxnvpols = 9
       allocate( pdata%nvpols(pdata%npols), pdata%poledges(pdata%npols, pdata%maxnvpols), &
            pdata%xpol(pdata%npols, pdata%maxnvpols), pdata%ypol(pdata%npols, pdata%maxnvpols), stat=allocerr )

       pdata%nvpols(1:pdata%npols) = (/ 9, 4, 4, 4, 4, 4, 4, 4, 4 /)

       pdata%xpol(1, 1:pdata%nvpols(1)) = (/ oct_circ_rad * dcos(PI / 4.0d0), &
            oct_circ_rad * dcos(PI / 2.0d0), &
            oct_circ_rad * dcos(135.0d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(PI), &
            oct_circ_rad * dcos(225.0d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(3.0d0 * PI / 2.0d0), & 
            oct_circ_rad * dcos(315.0d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(2.0d0 * PI), & 
            oct_circ_rad * dcos(PI / 4.0d0) /)
       pdata%ypol(1, 1:pdata%nvpols(1)) = (/ oct_circ_rad * dsin(PI / 4.0d0), &
            oct_circ_rad * dsin(PI / 2.0d0), &
            oct_circ_rad * dsin(135.0d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(PI), &
            oct_circ_rad * dsin(225.0d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(3.0d0 * PI / 2.0d0), & 
            oct_circ_rad * dsin(315.0d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(2.0d0 * PI), & 
            oct_circ_rad * dsin(PI / 4.0d0) /)

       pdata%xpol(2, 1:pdata%nvpols(2)) = (/ oct_circ_rad * dcos(PI / 4.0d0), & 
            oct_circ_rad * dcos(2.0d0 * PI), &
            dist_orig_to_tipstar * dcos(22.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(PI / 4.0d0) /)
       pdata%ypol(2, 1:pdata%nvpols(2)) = (/ oct_circ_rad * dsin(PI / 4.0d0), & 
            oct_circ_rad * dsin(2.0d0 * PI), &
            dist_orig_to_tipstar * dsin(22.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(PI / 4.0d0) /)

       pdata%xpol(3, 1:pdata%nvpols(3)) = (/ oct_circ_rad * dcos(PI / 2.0d0), & 
            oct_circ_rad * dcos(PI / 4.0d0), &
            dist_orig_to_tipstar * dcos(67.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(PI / 2.0d0) /)
       pdata%ypol(3, 1:pdata%nvpols(3)) = (/ oct_circ_rad * dsin(PI / 2.0d0), & 
            oct_circ_rad * dsin(PI / 4.0d0), &
            dist_orig_to_tipstar * dsin(67.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(PI / 2.0d0) /)

       pdata%xpol(4, 1:pdata%nvpols(4)) = (/ oct_circ_rad * dcos(135.0d0 * PI / 180.0d0), & 
            oct_circ_rad * dcos(PI / 2.0d0), &
            dist_orig_to_tipstar * dcos(112.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(135.0d0 * PI / 180.0d0) /)
       pdata%ypol(4, 1:pdata%nvpols(4)) = (/ oct_circ_rad * dsin(135.0d0 * PI / 180.0d0), & 
            oct_circ_rad * dsin(PI / 2.0d0), &
            dist_orig_to_tipstar * dsin(112.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(135.0d0 * PI / 180.0d0) /)

       pdata%xpol(5, 1:pdata%nvpols(5)) = (/ oct_circ_rad * dcos(PI), & 
            oct_circ_rad * dcos(135.0d0 * PI / 180.0d0), &
            dist_orig_to_tipstar * dcos(157.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(PI) /)
       pdata%ypol(5, 1:pdata%nvpols(5)) = (/ oct_circ_rad * dsin(PI), & 
            oct_circ_rad * dsin(135.0d0 * PI / 180.0d0), &
            dist_orig_to_tipstar * dsin(157.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(PI) /)

       pdata%xpol(6, 1:pdata%nvpols(6)) = (/ oct_circ_rad * dcos(225.0d0 * PI / 180.0d0), & 
            oct_circ_rad * dcos(PI), &
            dist_orig_to_tipstar * dcos(202.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(225.0d0 * PI / 180.0d0) /)
       pdata%ypol(6, 1:pdata%nvpols(6)) = (/ oct_circ_rad * dsin(225.0d0 * PI / 180.0d0), & 
            oct_circ_rad * dsin(PI), &
            dist_orig_to_tipstar * dsin(202.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(225.0d0 * PI / 180.0d0) /)

       pdata%xpol(7, 1:pdata%nvpols(7)) = (/ oct_circ_rad * dcos(3.0d0 * PI / 2.0d0), & 
            oct_circ_rad * dcos(225.0d0 * PI / 180.0d0), &
            dist_orig_to_tipstar * dcos(247.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(3.0d0 * PI / 2.0d0) /)
       pdata%ypol(7, 1:pdata%nvpols(7)) = (/ oct_circ_rad * dsin(3.0d0 * PI / 2.0d0), & 
            oct_circ_rad * dsin(225.0d0 * PI / 180.0d0), &
            dist_orig_to_tipstar * dsin(247.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(3.0d0 * PI / 2.0d0) /)

       pdata%xpol(8, 1:pdata%nvpols(8)) = (/ oct_circ_rad * dcos(315.0d0 * PI / 180.0d0), & 
            oct_circ_rad * dcos(3.0d0 * PI / 2.0d0), &
            dist_orig_to_tipstar * dcos(292.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(315.0d0 * PI / 180.0d0) /)
       pdata%ypol(8, 1:pdata%nvpols(8)) = (/ oct_circ_rad * dsin(315.0d0 * PI / 180.0d0), & 
            oct_circ_rad * dsin(3.0d0 * PI / 2.0d0), &
            dist_orig_to_tipstar * dsin(292.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(315.0d0 * PI / 180.0d0) /)

       pdata%xpol(9, 1:pdata%nvpols(9)) = (/ oct_circ_rad * dcos(2.0d0 * PI), & 
            oct_circ_rad * dcos(315.0d0 * PI / 180.0d0), &
            dist_orig_to_tipstar * dcos(337.5d0 * PI / 180.0d0), &
            oct_circ_rad * dcos(2.0d0 * PI) /)
       pdata%ypol(9, 1:pdata%nvpols(9)) = (/ oct_circ_rad * dsin(2.0d0 * PI), & 
            oct_circ_rad * dsin(315.0d0 * PI / 180.0d0), &
            dist_orig_to_tipstar * dsin(337.5d0 * PI / 180.0d0), &
            oct_circ_rad * dsin(2.0d0 * PI) /)

       pdata%poledges(1, 1:pdata%nvpols(1) - 1) = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
       pdata%poledges(2, 1:pdata%nvpols(2) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(3, 1:pdata%nvpols(3) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(4, 1:pdata%nvpols(4) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(5, 1:pdata%nvpols(5) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(6, 1:pdata%nvpols(6) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(7, 1:pdata%nvpols(7) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(8, 1:pdata%nvpols(8) - 1) = (/ 1, 0, 0 /)
       pdata%poledges(9, 1:pdata%nvpols(9) - 1) = (/ 1, 0, 0 /)

       ! Rectangle containing A
       pdata%dbl(1:2) = (/ -2.0d0 * dist_orig_to_tipstar, -2.0d0 * dist_orig_to_tipstar /)
       pdata%dtr(1:2) = (/ 2.0d0 * dist_orig_to_tipstar, 2.0d0 * dist_orig_to_tipstar /)

    elseif ( problem .eq. "stoyan" ) then
       ! Set disjoint convex polygons Aj whose union define A
       pdata%npols = 14
       allocate( pdata%nvpols(pdata%npols), stat=allocerr )
       if ( allocerr .ne. 0 ) then
          write (*,*) 'Memory allocation error.'
          stop
       end if

       pdata%nvpols(1:pdata%npols) = (/ 5, 6, 5, 5, 5, 5, 6, 5, 6, 5, 6, 7, 5, 4 /)
       pdata%maxnvpols = maxval( pdata%nvpols(1:pdata%npols) )

       allocate( pdata%poledges(pdata%npols,pdata%maxnvpols), &
            pdata%xpol(pdata%npols,pdata%maxnvpols),          &
            pdata%ypol(pdata%npols,pdata%maxnvpols),          &
            stat=allocerr )
       if ( allocerr .ne. 0 ) then
          write (*,*) 'Allocation error.'
          stop
       end if

       ! Vertices
       pdata%xpol(1, 1:pdata%nvpols(1)) = (/ 0.0d0, 0.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.0d0 /)
       pdata%ypol(1, 1:pdata%nvpols(1)) = (/ 2.0d0/3.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0 /)

       pdata%xpol(2, 1:pdata%nvpols(2)) = (/ 35.0d0/150.0d0, 2.0d0/15.0d0, 0.0d0, 0.0d0, 2.0d0/15.0d0, 35.0d0/150.0d0 /)
       pdata%ypol(2, 1:pdata%nvpols(2)) = (/ 1.0d0/3.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 0.2d0, 0.2d0, 1.0d0/3.0d0 /)

       pdata%xpol(3, 1:pdata%nvpols(3)) = (/ 0.0d0, 0.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.0d0 /)
       pdata%ypol(3, 1:pdata%nvpols(3)) = (/ 0.2d0, -4.0d0/15.0d0, -4.0d0/15.0d0, 0.2d0, 0.2d0 /)

       pdata%xpol(4, 1:pdata%nvpols(4)) = (/ 0.0d0, 0.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.0d0 /)
       pdata%ypol(4, 1:pdata%nvpols(4)) = (/ -4.0d0/15.0d0, -1.0d0/3.0d0, -1.0d0/3.0d0, -4.0d0/15.0d0, -4.0d0/15.0d0 /)

       pdata%xpol(5, 1:pdata%nvpols(5)) = (/ 2.0d0/15.0d0, 2.0d0/15.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/15.0d0 /)
       pdata%ypol(5, 1:pdata%nvpols(5)) = (/ 2.0d0/3.0d0, 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0 /)

       pdata%xpol(6, 1:pdata%nvpols(6)) = (/ 7.0d0/15.0d0, 8.0d0/15.0d0, 0.6d0, 7.0d0/15.0d0, 7.0d0/15.0d0 /)
       pdata%ypol(6, 1:pdata%nvpols(6)) = (/ 7.0d0/15.0d0, 7.0d0/15.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0, 7.0d0/15.0d0 /)

       pdata%xpol(7, 1:pdata%nvpols(7)) = (/ 8.0d0/15.0d0, 7.0d0/15.0d0, 0.3d0, 7.0d0/15.0d0, 0.6d0, 8.0d0/15.0d0 /)
       pdata%ypol(7, 1:pdata%nvpols(7)) = (/ 7.0d0/15.0d0, 7.0d0/15.0d0, 1.0d0/3.0d0, 0.2d0, 1.0d0/3.0d0, 7.0d0/15.0d0 /)

       pdata%xpol(8, 1:pdata%nvpols(8)) = (/ 7.0d0/15.0d0, 2.0d0/15.0d0, 2.0d0/15.0d0, 0.4d0, 7.0d0/15.0d0 /)
       pdata%ypol(8, 1:pdata%nvpols(8)) = (/ 0.2d0, 0.2d0, -4.0d0/15.0d0, 0.0d0, 0.2d0 /)

       pdata%xpol(9, 1:pdata%nvpols(9)) = (/ 11.0d0/15.0d0, 0.6d0, 7.0d0/15.0d0, 0.4d0, 0.6d0, 11.0d0/15.0d0 /)
       pdata%ypol(9, 1:pdata%nvpols(9)) = (/ 2.0d0/15.0d0, 1.0d0/3.0d0, 0.2d0, 0.0d0, 0.0d0, 2.0d0/15.0d0 /)

       pdata%xpol(10, 1:pdata%nvpols(10)) = (/ 13.0d0/15.0d0, 13.0d0/15.0d0, 1.0d0, 1.0d0, 13.0d0/15.0d0 /)
       pdata%ypol(10, 1:pdata%nvpols(10)) = (/ -4.0d0/15.0d0, -1.0d0/3.0d0, -1.0d0/3.0d0, -4.0d0/15.0d0, -4.0d0/15.0d0 /)

       pdata%xpol(11, 1:pdata%nvpols(11)) = (/ 13.0d0/15.0d0, 11.0d0/15.0d0, 13.0d0/15.0d0, 1.0d0, 1.0d0, 13.0d0/15.0d0 /)
       pdata%ypol(11, 1:pdata%nvpols(11)) = (/ 1.0d0/3.0d0, 2.0d0/15.0d0, -4.0d0/15.0d0, -4.0d0/15.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0 /)

       pdata%xpol(12, 1:pdata%nvpols(12)) = (/ 13.0d0/15.0d0, 0.8d0, 11.0d0/15.0d0, 13.0d0/15.0d0, 1.0d0, 1.0d0, 13.0d0/15.0d0 /)
       pdata%ypol(12, 1:pdata%nvpols(12)) = (/ 2.0d0/3.0d0, 2.0d0/3.0d0, 8.0d0/15.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 8.0d0/15.0d0, &
            2.0d0/3.0d0 /)

       pdata%xpol(13, 1:pdata%nvpols(13)) = (/ 11.0d0/15.0d0, 0.8d0, 8.0d0/15.0d0, 0.6d0, 11.0d0/15.0d0 /)
       pdata%ypol(13, 1:pdata%nvpols(13)) = (/ 8.0d0/15.0d0, 2.0d0/3.0d0, 7.0d0/15.0d0, 1.0d0/3.0d0, 8.0d0/15.0d0 /)

       pdata%xpol(14, 1:pdata%nvpols(14)) = (/ 11.0d0/15.0d0, 0.6d0, 13.0d0/15.0d0, 11.0d0/15.0d0 /)
       pdata%ypol(14, 1:pdata%nvpols(14)) = (/ 2.0d0/15.0d0, 0.0d0, -4.0d0/15.0d0, 2.0d0/15.0d0 /)

       ! Edges ID
       pdata%poledges(1, 1:pdata%nvpols(1) - 1) = (/ 0, 1, 1, 0 /)
       pdata%poledges(2, 1:pdata%nvpols(2) - 1) = (/ 0, 1, 0, 1, 0 /)
       pdata%poledges(3, 1:pdata%nvpols(3) - 1) = (/ 0, 1, 1, 1 /)
       pdata%poledges(4, 1:pdata%nvpols(4) - 1) = (/ 0, 0, 0, 1 /)
       pdata%poledges(5, 1:pdata%nvpols(5) - 1) = (/ 1, 0, 1, 0 /)
       pdata%poledges(6, 1:pdata%nvpols(6) - 1) = (/ 1, 0, 0, 1 /)
       pdata%poledges(7, 1:pdata%nvpols(7) - 1) = (/ 1, 0, 0, 1, 1 /)
       pdata%poledges(8, 1:pdata%nvpols(8) - 1) = (/ 0, 1, 0, 1 /)
       pdata%poledges(9, 1:pdata%nvpols(9) - 1) = (/ 0, 1, 1, 0, 1 /)
       pdata%poledges(10, 1:pdata%nvpols(10) - 1) = (/ 0, 0, 0, 1 /)
       pdata%poledges(11, 1:pdata%nvpols(11) - 1) = (/ 0, 1, 1, 0, 1 /)
       pdata%poledges(12, 1:pdata%nvpols(12) - 1) = (/ 0, 1, 0, 1, 0, 0 /)
       pdata%poledges(13, 1:pdata%nvpols(13) - 1) = (/ 1, 0, 1, 0 /)
       pdata%poledges(14, 1:pdata%nvpols(14) - 1) = (/ 1, 0, 1 /)

       ! Rectangle containing A
       pdata%dbl(1:2) = (/ -1.0d0/6.0d0, -0.5d0 /)
       pdata%dtr(1:2) = (/ 7.0d0/6.0d0, 5.0d0/6.0d0 /)
    end if
  end subroutine build_polygon

  ! *****************************************************************
  ! *****************************************************************

  subroutine drawsol(n, x, filename, drawhoneycomb, pdataptr, &
       distances, theta, nextra)

    ! PARAMETERS
    real(kind=8), parameter :: PI = acos( - 1.0d0 ), eps = 1.0d-04

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(in), optional :: nextra
    logical, intent(in) :: drawhoneycomb
    real(kind=8), optional, intent(in) :: theta
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    character(len=*), intent(in) :: filename
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(in), optional :: distances(n-1)

    ! LOCAL SCALARS
    logical :: inpolyfront,inthreeballs
    integer :: i,j,k,k1,k2,nballs,cont,numint
    real(kind=8) :: alpha,maxdist,dist,r
    type(pdata_type), pointer :: pdata

    ! LOCAL ARRAYS
    real(kind=8) :: p(2),q(2),v(2,2),vert1(2),vert2(2),u(2)

    call c_f_pointer(pdataptr,pdata)
    
    open(unit=10,file=trim(filename))

    ! BEGINING
    if ( present( nextra ) ) then
       write (10,10) (n-1)/2, x(nextra)
    else
       write (10,10) (n-1)/2, x(n)
    end if

    ! REGION A. 
    ! Define convex polygons as paths
    do j = 1, pdata%npols
        write (10,50) j
        do k = 1,pdata%nvpols(j)
            write (10,100) pdata%xpol(j,k),pdata%ypol(j,k)
        end do
        write (10,80)
    end do

    ! Draw the rectangle that contains the regular polygon
    ! write (10,50) 2
    ! write (10,100) pdata%dbl(1), pdata%dbl(2)
    ! write (10,100) pdata%dbl(1), pdata%dtr(2)
    ! write (10,100) pdata%dtr(1), pdata%dtr(2)
    ! write (10,100) pdata%dtr(1), pdata%dbl(2)
    ! write (10,80)
    ! write (10,130) pdata%dbl(1), pdata%dbl(2), pdata%dbl(1), pdata%dtr(2)
    ! write (10,130) pdata%dbl(1), pdata%dtr(2), pdata%dtr(1), pdata%dtr(2)
    ! write (10,130) pdata%dtr(1), pdata%dtr(2), pdata%dtr(1), pdata%dbl(2)
    ! write (10,130) pdata%dtr(1), pdata%dbl(2), pdata%dbl(1), pdata%dbl(2)
    
    ! Fill the region
    do j = 1,pdata%npols
        write (10,150) j
    end do
   
    ! Draw boundary using the partitions
    do j = 1,pdata%npols
        do k = 1,pdata%nvpols(j)
            if ( k .lt. pdata%nvpols(j) ) then
                if ( pdata%poledges(j,k) .eq. 1 ) then
                    write (10,125) pdata%xpol(j,k),pdata%ypol(j,k),pdata%xpol(j,k+1),pdata%ypol(j,k+1)
                end if
            end if
        end do
    end do

    do j = 1,pdata%npols
       do k = 1,pdata%nvpols(j)
          if ( k .lt. pdata%nvpols(j) ) then
             if ( pdata%poledges(j,k) .eq. 0 ) then
                write (10,130) pdata%xpol(j,k),pdata%ypol(j,k),pdata%xpol(j,k+1),pdata%ypol(j,k+1)
             end if
          end if
       end do
    end do

    ! BALLS DRAW
    if ( present (distances) ) then
       maxdist = maxval( distances(1:pdata%nballs) )
    end if

    if ( present( theta ) ) then
       alpha = 180.0d0 * theta / PI
    end if

    do i = 1, (n-1)/2
       write (10,200) i, x(2*i-1), i, x(2*i)

       ! if ( present( distances ) ) then
       !    write (10,220) i, 1.0d0 - ( distances(i) / maxdist ) ! fill ball
       ! end if

       write (10,230) i ! draw

       if ( drawhoneycomb ) write (10,210) alpha, alpha, i
    end do

    if ( present( nextra ) ) then
       do i = (n-1)/2 + 1, (nextra-1)/2
          write (10,200) i, x(2*i-1), i, x(2*i)
          write (10,240) i ! draw
          if ( drawhoneycomb ) write (10,250) alpha, alpha, i
       end do
    end if

    ! ==============================================================================
    ! ==============================================================================
    r = x(n)
    nballs = ( n - 1 ) / 2
    do i = 1,nballs
       p(1:2) = x(2*i-1:2*i)
       ! Check if there a polygon vertex in the frontier of the ball
       do k1 = 1,pdata%npols
          do k2 = 1,pdata%nvpols(k1) - 1
             if ( pdata%poledges(k1,k2) .eq. 0 ) then
                vert1(1:2) = (/ pdata%xpol(k1,k2), pdata%ypol(k1,k2) /)
                if ( abs( norm2( p(1:2) - vert1(1:2) ) - r ) .le. eps ) then
                   write(10,140) p(1:2),vert1(1:2)
                   write(10,141) p(1:2), p(1:2)
                   write(10,144) vert1(1:2), vert1(1:2)
                end if
             end if
          end do
       end do
       
       ! Check intersections with other balls
       do j = i + 1,nballs
          q(1:2) = x(2*j-1:2*j)
          dist = norm2( p(1:2) - q(1:2) )
          if ( dist .le. eps ) then
             write(*,*) 'Two balls centers are too close.'
             stop
          end if
          if ( abs( dist - 2.0d0 * r ) .le. eps ) then
             write(*,*) 'Two balls are nearly tangent.'
             ! stop
          end if
          if ( dist .lt. 2.0d0 * r ) then
             ! Compute the two distinct intersecting points
             call circles_intersect_points_2d(r,p,r,q,numint,v)
             if ( numint .ne. 2 ) then
                write(*,*) 'We called circles_intersect_points_2d expecting numint=2, but numint = ',numint
                stop
             end if
             do cont = 1,2
                ! Check whether v(1:2,cont) is in the frontier of any other circle
                inthreeballs = .false.
                do k = 1,nballs
                   if ( k .ne. i .and. k .ne. j ) then
                      if ( abs( norm2( v(1:2,cont) - x(2*k-1:2*k) ) - r ) .le. eps ) then
                         inthreeballs = .true.
                      end if
                   end if
                end do
                ! Check whether v(1:2,cont) is in the frontier of the
                ! polygon and whether it is a vertex of the polygon
                inpolyfront = .false.
                do k1 = 1,pdata%npols
                   do k2 = 1,pdata%nvpols(k1) - 1
                      if ( pdata%poledges(k1,k2) .eq. 0 ) then
                         vert1(1:2) = (/ pdata%xpol(k1,k2),  pdata%ypol(k1,k2)   /)
                         vert2(1:2) = (/ pdata%xpol(k1,k2+1),pdata%ypol(k1,k2+1) /)
                         call segment_contains_point_2d (vert1,vert2,v(1:2,cont),u)
                         if ( - eps .le. u(1) .and. u(1) .le. 1.0d0 + eps .and. abs( u(2) ) .le. eps ) then
                            inpolyfront = .true.
                         end if
                      end if
                   end do
                end do
                ! Check whether v(1:2,cont), that is in the intersection of two balls, is a relevant point
                ! if ( inthreeballs .or. inpolyfront ) then
                !    write(10,140) p(1:2),v(1:2,cont)
                !    write(10,141) p(1:2),p(1:2)
                !    write(10,141) v(1:2,cont),v(1:2,cont)
                !    write(10,140) q(1:2),v(1:2,cont)
                !    write(10,141) q(1:2)
                !    write(10,141) v(1:2,cont)
                ! end if

                if ( inthreeballs .or. inpolyfront ) then
                   write(10,140) p(1:2),v(1:2,cont)
                   write(10,140) q(1:2),v(1:2,cont)
                   write(10,141) p(1:2),p(1:2)
                   write(10,141) q(1:2),q(1:2)
                   if ( inthreeballs ) then
                      write(10,142) v(1:2,cont),v(1:2,cont)
                   else
                      write(10,143) v(1:2,cont),v(1:2,cont)
                   end if
                end if
             end do
          end if
       end do
    end do
    ! ==============================================================================
    ! ==============================================================================
    
    write (10,40) x(n)

    close(10)

    ! NON-EXECUTABLE STATEMENTS

10  format('prologues := 3;',/, &
         'outputformat   := "svg";',/, &
         'outputtemplate := "%j.%o";',/, &
         'input mpcolornames;',/, &
         'beginfig(',i0,');',/, &
         'u := 20cm;',/, &
         'path pols[];',/,'pair p[];',/,'R := ',f0.10,'u;',/,'D := 2*R;',/,'')
50  format('pols[',I0,'] = ')
80  format('cycle;')
100 format('(',f20.10,'u,',f20.10,'u)--')
125 format('draw (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;')
130 format('draw (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 1.0;')
140 format('draw (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10, &
         'u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;')
! Centers
! 141 format('drawdot (',f20.10,'u,',f20.10,'u) withpen pencircle scaled 5.0 withcolor Red;')
141 format('draw fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor 0.4white;',/, &
           'fill fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor 0.9red;')

! 3 circles intersection
! 142 format('drawdot (',f20.10,'u,',f20.10,'u) withpen pencircle scaled 8.0 withcolor MidnightBlue;')
142 format('draw fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor 0.4white;',/, &
           'fill fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor MidnightBlue;')

! Two circles and border
! 143 format('drawdot (',f20.10,'u,',f20.10,'u) withpen pencircle scaled 8.0 withcolor Green;')
143 format('draw fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor 0.4white;',/, &
           'fill fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor Green;')
! Circle and vertice
! 144 format('drawdot (',f20.10,'u,',f20.10,'u) withpen pencircle scaled 8.0 withcolor Yellow;')
144 format('draw fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor 0.4white;',/, &
           'fill fullcircle scaled 8.0 shifted (',f20.10,'u,',f20.10,'u) withcolor Yellow;')

150 format('fill pols[',I0,'] withcolor 0.9Dandelion;')
40  format('% r = ',f20.10,/,'endfig;',/,'end;')
    
200 format(/,'x',I0,' = ',f0.10,'u;',/,'y',I0,' = ',f0.10,'u;')
210 format(/,'k := 1;',/, &
         'for a := 0 + ',f0.10,' step 60 until 300 + ', f0.10, ':',/, &
	 '  p[k] := z',i0,' + (R*dir a);',/, &
	 '  k := k + 1;',/, &
	 'endfor',/, &
         '% Draw the hexagon with center in z1',/, &
         'for k := 1 upto 5:',/, &
	 '  draw p[k]--p[k+1] withcolor 1.0white;',/, &
	 'endfor',/, &
	 'draw p[6]--p[1] withcolor 1.0white;')
    ! '% Draw the hexagon inscribed circle with radius r',/, &
    ! 'draw fullcircle scaled D shifted z',I0,' withcolor RoyalBlue;')

220 format('fill fullcircle scaled D shifted z',i0,' withcolor ',0P, F4.3,'white;')
230 format('draw fullcircle scaled D shifted z',i0,' withcolor RoyalBlue withpen pencircle scaled 1.0;')
240 format('draw fullcircle scaled D shifted z',i0,' withcolor RoyalBlue withpen pencircle scaled 1.0 dashed evenly;')

250 format(/,'k := 1;',/, &
         'for a := 0 + ',f0.10,' step 60 until 300 + ', f0.10, ':',/, &
	 '  p[k] := z',i0,' + (R*dir a);',/, &
	 '  k := k + 1;',/, &
	 'endfor',/, &
         '% Draw the hexagon with center in z1',/, &
         'for k := 1 upto 5:',/, &
	 '  draw p[k]--p[k+1] withcolor 0.7white dashed evenly;',/, &
	 'endfor',/, &
	 'draw p[6]--p[1] withcolor 0.7white dashed evenly;')

  end subroutine drawsol

  ! *****************************************************************
  ! *****************************************************************
  
  real(kind=8) function polygon_area (pdataptr)
    ! SCALAR ARGUMENTS
    type(c_ptr), optional, intent(in) :: pdataptr

    ! LOCAL SCALARS
    integer      :: allocerr, indpol
    real(kind=8) :: convPolArea

    ! LOCAL TYPES
    type(Polygon) :: convPol
    type(pdata_type), pointer :: pdata

    call c_f_pointer( pdataptr, pdata )

    polygon_area = 0.0d0

    do indpol = 1, pdata%npols
       convPol%n = pdata%nvpols(indpol)
       convPol%deg = .false.

       allocate ( convPol%vertex(2,convPol%n), stat=allocerr )

       if ( allocerr .ne. 0 ) then
          write (*,*) 'Allocation error.'
          stop
       end if

       convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
       convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)

       call polygon_area_2d ( convPol%n - 1, convPol%vertex(1:2,1:convPol%n - 1), convPolArea )

       polygon_area = polygon_area + convPolArea

       deallocate( convPol%vertex, stat=allocerr )

       if ( allocerr .ne. 0 ) then
          write (*,*) 'Deallocation error.'
          stop
       end if
    end do
  end function polygon_area

  ! *****************************************************************
  ! *****************************************************************

  real(kind=8) function ball_intersection_polygon(bc, r, pdata)
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )

    ! SCALAR ARGUMENTS
    real(kind=8),     intent(in) :: r
    type(pdata_type), intent(in) :: pdata

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: bc(2)

    ! LOCAL SCALARS
    logical :: inside
    integer :: allocerr,indpol,v
    real(kind=8) :: area,eps,seed,theta1,theta2,totalCoveredArea

    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol

    ! LOCAL ARRAYS
    integer      :: indexes(1)
    real(kind=8) :: pD(2,4), centers(2)

    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand

    eps = 1.0d-10
    seed = 654321.0

    ! Indexes of the balls and centers
    indexes(1) = 1
    centers(1) = bc(1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
    centers(2) = bc(2) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )

    ! Points of D
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)
    
    totalCoveredArea = 0.0d0

    do indpol=1,pdata%npols

       area = 0.0d0
       convPol%n = pdata%nvpols(indpol)
       convPol%deg = .false.

       allocate ( convPol%vertex(2,convPol%n), convPol%edges(convPol%n - 1), stat=allocerr )

       if ( allocerr .ne. 0 ) then
          write (*,*) 'Allocation error.'
          stop
       end if

       convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
       convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
       convPol%edges(1:convPol%n - 1) = pdata%poledges(indpol, 1:convPol%n - 1)

       call interCellBall(centers(1:2), r, convPol, curvPol)

       if ( curvPol%deg ) then

          call polygon_contains_point_2d_convex(convPol%n - 1,convPol%vertex(1:2,1:convPol%n - 1),centers(1:2),&
               inside)

          if ( inside ) then
             area = area + (PI * r**2)
          end if
       else
          do v = 1,curvPol%n - 1

             if ( curvPol%edges(v) .eq. 0 ) then
                area = area + 0.5d0 * ( curvPol%vertex(2,v+1) - curvPol%vertex(2,v) )&
                     * ( curvPol%vertex(1,v) + curvPol%vertex(1,v+1) )
             else
                theta1 = angle_rad_2d(curvPol%vertex(1:2,v), centers(1:2), &
                     (/ centers(1) + 1.0d0, centers(2) /))

                theta2 = angle_rad_2d(curvPol%vertex(1:2,v+1), centers(1:2), &
                     (/ centers(1) + 1.0d0, centers(2) /))

                if ( theta2 .lt. theta1) then
                   theta2 = theta2 + 2.0d0 * PI
                end if

                area = area + 0.5d0 * r**2 * (  (theta2 - theta1) + &
                     (dsin(theta2) * dcos(theta2)) - (dsin(theta1) * dcos(theta1)) ) + centers(1) * r * ( &
                     dsin(theta2) - dsin(theta1) )
             end if

          end do
       end if

       totalCoveredArea = totalCoveredArea + area

       deallocate(convPol%vertex,convPol%edges, stat=allocerr)

       if ( allocerr .ne. 0 ) then
          write (*,*) 'Deallocation error.'
          stop
       end if

       if ( allocated(curCell%vertex) ) deallocate(curCell%vertex)
       if ( allocated(curCell%edges) ) deallocate(curCell%edges)
       if ( allocated(curvPol%vertex) ) deallocate(curvPol%vertex)
       if ( allocated(curvPol%vertex_id) ) deallocate(curvPol%vertex_id)
       if ( allocated(curvPol%edges) ) deallocate(curvPol%edges)

    end do

    if ( allocated(vorCell%vertex) ) deallocate(vorCell%vertex)
    if ( allocated(vorCell%edges) ) deallocate(vorCell%edges)


    ball_intersection_polygon = totalCoveredArea / ( PI * r**2 )

  end function ball_intersection_polygon

  ! *****************************************************************
  ! *****************************************************************

  real(kind=8) function dist_point_line_segment(p, a, b)

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(2), intent(in) :: a, b, p

    ! This function computes the distance between a point p and a line
    ! segment with extremes a and b

    ! LOCAL SCALARS
    real(kind=8) :: u

    ! LOCAL ARRAY
    real(kind=8), dimension(2) :: pi

    u = sum( ( p(1:2) - a(1:2) ) * ( b(1:2) - a(1:2) ) ) / sum( (b(1:2) - a(1:2))**2 )
    u = min( max( 0.0d0, u ), 1.0d0 )

    pi(1:2) = a(1:2) + u * ( b(1:2) - a(1:2) )

    dist_point_line_segment = norm2( p(1:2) - pi(1:2) )

  end function dist_point_line_segment

  ! *****************************************************************
  ! *****************************************************************

  real(kind=8) function dist_point_to_polygon_boundary(p, pdataptr)

    ! SCALAR ARGUMENTS
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(2), intent(in) :: p

    ! LOCAL SCALARS
    integer :: j, k

    ! LOCAL ARRAYS
    real(kind=8), dimension(2) :: a, b

    ! LOCAL POINTER
    type(pdata_type), pointer :: pdata

    call c_f_pointer( pdataptr, pdata )

    dist_point_to_polygon_boundary = huge( 1.0d0 )

    do j = 1,pdata%npols
       do k = 1,pdata%nvpols(j)
          if ( k .lt. pdata%nvpols(j) ) then
             if ( pdata%poledges(j,k) .eq. 0 ) then

                a(1:2) = (/ pdata%xpol(j,k), pdata%ypol(j,k) /)
                b(1:2) = (/ pdata%xpol(j,k+1), pdata%ypol(j,k+1) /)

                dist_point_to_polygon_boundary =             &
                     min( dist_point_to_polygon_boundary,    &
                          dist_point_line_segment(p, a, b) )
             end if
          end if
       end do
    end do

  end function dist_point_to_polygon_boundary

  ! *****************************************************************
  ! *****************************************************************

  real(kind=8) function polygon_perimeter(pdataptr)
    ! SCALAR ARGUMENTS
    type(c_ptr), optional, intent(in) :: pdataptr

    ! LOCAL SCALARS
    integer :: j, k

    ! LOCAL ARRAYS
    real(kind=8), dimension(2) :: a, b

    ! LOCAL POINTER
    type(pdata_type), pointer :: pdata

    call c_f_pointer( pdataptr, pdata )

    polygon_perimeter = 0.0d0
    do j = 1,pdata%npols
       do k = 1,pdata%nvpols(j)
          if ( k .lt. pdata%nvpols(j) ) then
             if ( pdata%poledges(j,k) .eq. 0 ) then
                a(1:2) = (/ pdata%xpol(j,k), pdata%ypol(j,k) /)
                b(1:2) = (/ pdata%xpol(j,k+1), pdata%ypol(j,k+1) /)

                polygon_perimeter = polygon_perimeter + &
                     norm2( a(1:2) - b(1:2) )
             end if
          end if
       end do
    end do
  end function polygon_perimeter

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf(n,x,f,inform,pdataptr)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(1) = pdata%counters(1) + 1

    f = x(n)

    inform = 0
  end subroutine evalf

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg(n,x,g,inform,pdataptr)
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr
  
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1
    
    g(1:n-1) = 0.0d0
    g(n)     = 1.0d0

    inform = 0
  end subroutine evalg

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalc(n,x,m,p,c,inform,pdataptr)
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m+p)

    ! LOCAL SCALARS
    logical :: inside
    integer :: allocerr,i,indpol,ierror,k,status,v,v_num
    real(kind=8) :: area,convPolArea,dist,eps,r,seed,theta1,theta2,totalCoveredArea,totalPolsArea
    type(pdata_type), pointer :: pdata
    
    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol
    
    ! LOCAL ARRAYS
    real(kind=8) :: pD(2,4)
  
    ! ALLOCATABLE ARRAYS
    integer, dimension(:), allocatable :: indexes
    integer, dimension(:,:), allocatable :: nod_tri, tnbr
    real(kind=8), dimension (:,:), allocatable :: centers, vor_xy
  
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(3) = pdata%counters(3) + 1

    allocate( centers(2,pdata%nballs), indexes(pdata%nballs), &
         nod_tri(3,2*pdata%nballs), tnbr(3,2*pdata%nballs),   &
         stat=allocerr )
    if ( allocerr .ne. 0 ) then
       write (*,*) "Memory allocation error in EVALC."
       stop
    end if
    
    eps = 1.0d-10
    seed = 654321.0
    
    ! Indexes of the balls and centers
    
    do i = 1,pdata%nballs
       indexes(i) = i
       centers(1,i) = x(2*i-1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
       centers(2,i) = x(2*i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
    end do
  
    ! Radius
    
    r = x(2*pdata%nballs+1)
    
    ! Points of D
    
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)
    
    ! Check whether r is positive
    if ( r .le. 0.0d0 ) then
        c(1) = 1.0d+09
        write (*,*) '*** IN EVALC SUBROUTINE: RADIUS LESS OR EQUAL THAN ZERO. ***'
        write (*,*) 'nballs = ', pdata%nballs, 'r = ', x(2*pdata%nballs+1), "n = ", n
        return
    else
        ! Check whether centers are close enough to A
        do i = 1,pdata%nballs
            call quad_point_dist_signed_2d( pD, centers(:,i), dist )
            if ( dist .gt. r ) then
                c(1) = 1.0d+09
                write (*,*) '*** IN EVALC SUBROUTINE: THERE IS A CENTER TOO FAR FROM D. ***'
                return
            end if
        end do
    end if
    
    ! Delaunay triangulation
    
    if ( pdata%nballs .ge. 3 ) then

      call dtris2(pdata%nballs, centers, indexes, v_num, nod_tri, tnbr, ierror)
      
      if ( ierror .eq. 0 ) then
        ! Obtaining the Voronoi vertices.
        allocate ( vor_xy(2,v_num), stat=allocerr )
        
        if ( allocerr .ne. 0 ) then
            write (*,*) 'Allocation error.'
            stop
        end if
        
        do k=1,v_num
            call triangle_circumcenter_2d ( centers(1:2, nod_tri(1:3, k)), vor_xy(1:2, k) )
        end do
      
      end if
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Error in dtris2'
        write ( *, '(a,i6)' ) '  IERROR = ', ierror
        return
      end if
      
    end if
    
    totalCoveredArea = 0.0d0
    
    do i = 1,pdata%nballs
    
        if ( pdata%nballs .ge. 3 .and. ierror .eq. 0 ) then
        
            call voronoi_cell(i, pdata%nballs, centers, v_num, nod_tri, tnbr, vor_xy, vorCell, status)
            
            if ( status .ne. 0 ) then
                inform = -51
                return
            end if
            
        end if
        
        do indpol=1,pdata%npols
            
            area = 0.0d0
            convPol%n = pdata%nvpols(indpol)
            convPol%deg = .false.
            
            allocate ( convPol%vertex(2,convPol%n), convPol%edges(convPol%n - 1), stat=allocerr )
            
            if ( allocerr .ne. 0 ) then
               write (*,*) 'Allocation error.'
               stop
            end if
            
            convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
            convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
            convPol%edges(1:convPol%n - 1) = pdata%poledges(indpol, 1:convPol%n - 1)
            
            if ( pdata%nballs .ge. 3 .and. ierror .eq. 0 ) then
                call VorCellInterConvPol(centers(:,i), convPol, vorCell, curCell)
            end if
            
            if ( pdata%nballs .eq. 2 .or. (pdata%nballs .ge. 3 .and. ierror .eq. 225) ) then
                call VorCellInterConvPol_Collinear(i, pdata%nballs, centers, convPol, curCell)
            end if
            
            if ( (.not. curCell%deg) .or. pdata%nballs .eq. 1 ) then
                
                if ( pdata%nballs .eq. 1 ) then
                    call interCellBall(centers(1:2,i), r, convPol, curvPol)
                else
                    call interCellBall(centers(1:2,i), r, curCell, curvPol)
                end if
                
                if ( curvPol%deg ) then
                    
                    if ( pdata%nballs .eq. 1 ) then
                        call polygon_contains_point_2d_convex(convPol%n - 1,convPol%vertex(1:2,1:convPol%n - 1),centers(1:2,i),&
                        inside)
                    else
                        call polygon_contains_point_2d_convex(curCell%n - 1,curCell%vertex(1:2,1:curCell%n - 1),centers(1:2,i),&
                        inside)
                    end if
                    
                    if ( inside ) then
                        area = area + (PI * r**2)
                    end if
                else
                    do v = 1,curvPol%n - 1
                        
                        if ( curvPol%edges(v) .eq. 0 ) then
                            area = area + 0.5d0 * ( curvPol%vertex(2,v+1) - curvPol%vertex(2,v) )&
                            * ( curvPol%vertex(1,v) + curvPol%vertex(1,v+1) )
                        else
                            theta1 = angle_rad_2d(curvPol%vertex(1:2,v), centers(1:2,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /))
                            
                            theta2 = angle_rad_2d(curvPol%vertex(1:2,v+1), centers(1:2,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /))
                            
                            if ( theta2 .lt. theta1) then
                                theta2 = theta2 + 2.0d0 * PI
                            end if
                            
                            area = area + 0.5d0 * r**2 * (  (theta2 - theta1) + &
                            (dsin(theta2) * dcos(theta2)) - (dsin(theta1) * dcos(theta1)) ) + centers(1,i) * r * ( &
                            dsin(theta2) - dsin(theta1) )
                        end if
                        
                    end do
                end if
            end if
            
            totalCoveredArea = totalCoveredArea + area
            
            deallocate(convPol%vertex,convPol%edges, stat=allocerr)
            
            if ( allocerr .ne. 0 ) then
               write (*,*) 'Deallocation error.'
               stop
            end if
            
            if ( allocated(curCell%vertex) ) deallocate(curCell%vertex)
            if ( allocated(curCell%edges) ) deallocate(curCell%edges)
            if ( allocated(curvPol%vertex) ) deallocate(curvPol%vertex)
            if ( allocated(curvPol%vertex_id) ) deallocate(curvPol%vertex_id)
            if ( allocated(curvPol%edges) ) deallocate(curvPol%edges)
            
        end do
        
        if ( allocated(vorCell%vertex) ) deallocate(vorCell%vertex)
        if ( allocated(vorCell%edges) ) deallocate(vorCell%edges)
        
    end do
    
    if ( allocated(vor_xy) ) deallocate(vor_xy)

    totalPolsArea = 0.0d0
    
    do indpol=1,pdata%npols
        
        convPol%n = pdata%nvpols(indpol)
        convPol%deg = .false.
        
        allocate ( convPol%vertex(2,convPol%n), stat=allocerr )
        
        if ( allocerr .ne. 0 ) then
            write (*,*) 'Allocation error.'
            stop
        end if
        
        convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
        convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
        
        call polygon_area_2d ( convPol%n - 1, convPol%vertex(1:2,1:convPol%n - 1), convPolArea )
        
        totalPolsArea = totalPolsArea + convPolArea
        
        deallocate(convPol%vertex, stat=allocerr)
        
        if ( allocerr .ne. 0 ) then
            write (*,*) 'Deallocation error.'
            stop
        end if

     end do

     deallocate( centers, indexes, nod_tri, tnbr, stat=allocerr )

     c(1) = totalPolsArea - totalCoveredArea
     inform = 0

  end subroutine evalc

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)
    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m+p)
    real(kind=8), intent(in) :: x(n)
    logical, intent(out) :: sorted(m+p)
    integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
    real(kind=8), intent(out) :: jval(lim)
    
    ! LOCAL SCALARS
    integer :: allocerr,i,ierror,indpol,k,status,v,v_num
    real(kind=8) :: eps,r,seed,theta1,theta2,dist
    logical :: inside
    type(pdata_type), pointer :: pdata

    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol

    ! LOCAL ARRAYS
    real(kind=8) :: pD(2,4)
    
    ! ALLOCATABLE ARRAYS
    integer, dimension(:), allocatable :: indexes
    integer, dimension(:,:), allocatable :: nod_tri, tnbr
    real ( kind = 8 ), dimension (:,:), allocatable :: centers, vor_xy
    
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(4) = pdata%counters(4) + 1

    allocate( centers(2,pdata%nballs), indexes(pdata%nballs), &
         nod_tri(3,2*pdata%nballs),tnbr(3, 2*pdata%nballs),   &
         stat=allocerr )
    if ( allocerr .ne. 0 ) then
       write (*,*) "Memory allocation error at EVALJ."
       stop
    end if
    
    ! Points of D
    
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)

    ! Only gradients of constraints j such that ind(j) = .true. need
    ! to be computed.
    
    if ( ind(1) ) then
       if ( lim .lt. n ) then
          inform = -94
          return
       end if
       
       jsta(1) = 1
       jlen(1) = n
       
       jvar(1:n) = (/ (i,i=1,n) /)

       jval(1:n) = 0.0d0
       
       eps = 1.0d-10
       seed = 654321.0
    
       ! Indexes of the balls and centers
    
       do i = 1,pdata%nballs
          indexes(i) = i
          centers(1,i) = x(2*i-1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
          centers(2,i) = x(2*i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
       end do
    
       ! Radius
       r = x(2*pdata%nballs+1)
       
       ! Check whether r is positive
       if ( r .le. 0.0d0 ) then
           write (*,*) '*** IN EVALJ SUBROUTINE: RADIUS LESS OR EQUAL THAN ZERO. ***'
           return
       else
           ! Check whether centers are close enough to A
           do i = 1,pdata%nballs
               call quad_point_dist_signed_2d( pD, centers(:,i), dist )
               if ( dist .gt. r ) then
                   write (*,*) '*** IN EVALJ SUBROUTINE: THERE IS A CENTER TOO FAR FROM D. ***'
                   return
               end if
           end do
       end if
       
       ! Delaunay triangulation
       
       if ( pdata%nballs .ge. 3 ) then
          
          call dtris2(pdata%nballs, centers, indexes, v_num, nod_tri, tnbr, ierror)
          
          if ( ierror .eq. 0 ) then
             
             ! Obtaining the Voronoi vertices.
             allocate ( vor_xy(2,v_num), stat=allocerr )
             
             if ( allocerr .ne. 0 ) then
                write (*,*) 'Allocation error.'
                stop
             end if
             
             do k=1,v_num
                call triangle_circumcenter_2d ( centers(1:2, nod_tri(1:3, k)), vor_xy(1:2, k) )
             end do
             
          end if
          
          if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'Error in dtris2'
             write ( *, '(a,i6)' ) '  IERROR = ', ierror
             return
          end if
          
       end if
       
       do i = 1,pdata%nballs
          
          if ( pdata%nballs .ge. 3 .and. ierror .eq. 0 ) then
          
             call voronoi_cell(i, pdata%nballs, centers, v_num, nod_tri, tnbr, vor_xy, vorCell, status)
             
             if ( status .ne. 0 ) then
                inform = -52
                return
             end if
             
          end if
          
          do indpol=1,pdata%npols
             
             convPol%n = pdata%nvpols(indpol)
             convPol%deg = .false.
             
             allocate ( convPol%vertex(2,convPol%n), convPol%edges(convPol%n - 1), stat=allocerr )
             
             if ( allocerr .ne. 0 ) then
                write (*,*) 'Allocation error.'
                stop
             end if
            
             convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
             convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
             convPol%edges(1:convPol%n - 1) = pdata%poledges(indpol, 1:convPol%n - 1)
             
             if ( pdata%nballs .ge. 3 .and. ierror .eq. 0 ) then
                call VorCellInterConvPol(centers(:,i), convPol, vorCell, curCell)
             end if
             
             if ( pdata%nballs .eq. 2 .or. (pdata%nballs .ge. 3 .and. ierror .eq. 225) ) then
                call VorCellInterConvPol_Collinear(i, pdata%nballs, centers, convPol, curCell)
             end if
             
             if ( .not. curCell%deg .or. pdata%nballs .eq. 1 ) then
                
                if ( pdata%nballs .eq. 1 ) then
                   call interCellBall(centers(1:2,i), r, convPol, curvPol)
                else
                   call interCellBall(centers(1:2,i), r, curCell, curvPol)
                end if
                
                if ( curvPol%deg ) then
                   
                   if ( pdata%nballs .eq. 1 ) then
                      call polygon_contains_point_2d_convex(convPol%n - 1,convPol%vertex(1:2,1:convPol%n - 1),centers(1:2,i),&
                      inside)
                   else
                      call polygon_contains_point_2d_convex(curCell%n - 1,curCell%vertex(1:2,1:curCell%n - 1),centers(1:2,i),&
                      inside)
                   end if
                   
                   if ( inside ) then
                      jval(2*pdata%nballs+1) = jval(2*pdata%nballs+1) + 2.0d0 * PI * r
                   end if
                   
                else
                   
                   do v = 1,curvPol%n - 1
                      if ( curvPol%edges(v) .eq. 1 ) then
                         
                         theta1 = angle_rad_2d(curvPol%vertex(1:2,v), centers(1:2,i), &
                         (/ centers(1,i) + 1.0d0, centers(2,i) /))
                         theta2 = angle_rad_2d(curvPol%vertex(1:2,v+1), centers(1:2,i), &
                         (/ centers(1,i) + 1.0d0, centers(2,i) /))
                         
                         if ( theta2 .lt. theta1) then
                            theta2 = theta2 + 2.0d0 * PI
                         end if
                         
                         jval(2*pdata%nballs+1) = jval(2*pdata%nballs+1) + r * (theta2 - theta1)
                         jval(2*i-1)      = jval(2*i-1) + r * (dsin(theta1) - dsin(theta2))
                         jval(2*i)        = jval(2*i)   + r * (dcos(theta2) - dcos(theta1))
                         
                      end if
                   end do
                   
                end if
                
             end if
             
             deallocate( convPol%vertex,convPol%edges, stat=allocerr )
             
             if ( allocerr .ne. 0 ) then
                write (*,*) 'Deallocation error.'
                stop
             end if
             
             if ( allocated(curCell%vertex) ) deallocate(curCell%vertex)
             if ( allocated(curCell%edges) ) deallocate(curCell%edges)
             if ( allocated(curvPol%vertex) ) deallocate(curvPol%vertex)
             if ( allocated(curvPol%vertex_id) ) deallocate(curvPol%vertex_id)
             if ( allocated(curvPol%edges) ) deallocate(curvPol%edges)
          end do
          if ( allocated(vorCell%vertex) ) deallocate(vorCell%vertex)
          if ( allocated(vorCell%edges) ) deallocate(vorCell%edges)
       end do
       jval(2*pdata%nballs+1) = -1.0d0 * jval(2*pdata%nballs+1)
       if ( allocated(vor_xy) ) deallocate(vor_xy)

       deallocate( centers, indexes, nod_tri,tnbr, stat=allocerr )

       ! Says whether the variables' indices in jvar (related to this
       ! constraint) are in increasing order. In case they are,
       ! Algencan takes advantage of this. Implement sorted gradients
       ! of constraints if you can do this in a natural (cheap)
       ! way. Under no circumnstance use a sorting algorithm. (It is
       ! better to set sorted(1) = .false. in this case.)
       
       sorted(1) = .true.
       inform = 0
    end if
    
  end subroutine evalj

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)

    ! PARAMETERS
    real(kind=8), parameter :: PI = dacos( - 1.0d0 )
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: inclf
    integer, intent(in) :: m,n,lim,p
    integer, intent(out) :: hlnnz
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: lambda(m+p),x(n)
    integer, intent(out) :: hlrow(lim),hlcol(lim)
    real(kind=8), intent(out) :: hlval(lim)
    
    ! LOCAL SCALARS
    integer :: allocerr,i,ierror,indpol,j,k,status,v,v_num
    real(kind=8) :: dist,eps,inner_prod_div,r,seed,theta1,theta2,theta_set_inter
    logical :: inside
    type(pdata_type), pointer :: pdata

    ! LOCAL TYPES
    type(VoronoiCell) :: vorCell
    type(Polygon) :: convPol,curCell
    type(CurvilinearPolygon) :: curvPol

    ! LOCAL ARRAYS
    integer :: block_diag_ind(3),block_row_vector_ind(2)
    real(kind=8) :: block_diag_val(3),block_row_vector_val(2),normal(2),pbA(2),pD(2,4)
    
    ! ALLOCATABLE ARRAYS
    integer, dimension (:,:), allocatable :: blocks_below_diag_ind, nod_tri, tnbr
    integer, dimension (:), allocatable :: indexes, intersections
    real(kind=8), dimension (:,:), allocatable :: centers, blocks_below_diag_vals,vor_xy
    
    ! FUNCTIONS
    real(kind=8) :: angle_rad_2d, drand
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(5) = pdata%counters(5) + 1

    allocate( centers(2,pdata%nballs),indexes(pdata%nballs),nod_tri(3, &
         2*pdata%nballs),tnbr(3, 2*pdata%nballs), stat=allocerr )
    if ( allocerr .ne. 0 ) then
       write (*,*) 'Allocation error in evalhl.'
       stop
    end if
    
    ! Points of D
    
    pD(:,1) = pdata%dbl(1:2)
    pD(:,2) = (/ pdata%dtr(1), pdata%dbl(2) /)
    pD(:,3) = pdata%dtr(1:2)
    pD(:,4) = (/ pdata%dbl(1), pdata%dtr(2) /)

    ! If .not. inclf then the Hessian of the objective function must not be included
    ! The Hessian of the objective function is zero.
    
    eps = 1.0d-10
    seed = 654321.0
    
    ! Indexes of the balls and centers
    
    do i = 1,pdata%nballs
       indexes(i) = i
       centers(1,i) = x(2*i-1) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
       centers(2,i) = x(2*i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 )
    end do
  
    ! Radius
    r = x(2*pdata%nballs+1)
    
    hlnnz = 0
    
    ! Check whether r is positive
    if ( r .le. 0.0d0 ) then
        write (*,*) '*** IN EVALHL SUBROUTINE: RADIUS LESS OR EQUAL THAN ZERO. ***'
        write (*,*) '*** RETURNING THE HESSIAN AS ZERO. ***'
        return
    else
        ! Check whether centers are close enough to A
        do i = 1,pdata%nballs
            call quad_point_dist_signed_2d( pD, centers(:,i), dist )
            if ( dist .gt. r ) then
                write (*,*) '*** IN EVALHL SUBROUTINE: THERE IS A CENTER TOO FAR FROM D. ***'
                write (*,*) '*** RETURNING THE HESSIAN AS ZERO. ***'
                return
            end if
        end do
    end if
    
    ! Delaunay triangulation
       
    if ( pdata%nballs .ge. 3 ) then
        
        call dtris2(pdata%nballs, centers, indexes, v_num, nod_tri, tnbr, ierror)
        
        if ( ierror .eq. 0 ) then
            
            ! Obtaining the Voronoi vertices.
            allocate ( vor_xy(2,v_num), stat=allocerr )
            
            if ( allocerr .ne. 0 ) then
               write (*,*) 'Allocation error.'
               stop
            end if
            
            do k=1,v_num
               call triangle_circumcenter_2d ( centers(1:2, nod_tri(1:3, k)), vor_xy(1:2, k) )
            end do
            
        end if
        
        if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'Error in dtris2'
            write ( *, '(a,i6)' ) '  IERROR = ', ierror
            return
        end if
        
    end if
    
    ! Note that entries of the Hessian of the Lagrangian can be
    ! repeated. If this is case, them sum of repeated entrances is
    ! considered. This feature simplifies the construction of the
    ! Hessian of the Lagrangian.
    
    if ( hlnnz + 1 .gt. lim ) then
       inform = -95
       return
    end if
    
    ! Second order derivative with respect to the radius.
    hlnnz = hlnnz + 1 
    hlrow(hlnnz) = n
    hlcol(hlnnz) = n
    hlval(hlnnz) = 0.0d0
    
    do i = 1,pdata%nballs
        
        ! These arrays are useful for computing blocks below the diagonal.
        if ( i .ne. pdata%nballs ) then
           allocate( intersections(pdata%nballs - i),blocks_below_diag_ind(4,pdata%nballs - i), &
                blocks_below_diag_vals(4,pdata%nballs - i), stat=allocerr )
            
            if ( allocerr .ne. 0 ) then
               write (*,*) 'Allocation error.'
               stop
            end if
            
            intersections(1:pdata%nballs - i) = 0
            blocks_below_diag_ind(1:4,1:pdata%nballs - i) = 0
            blocks_below_diag_vals(1:4,1:pdata%nballs - i) = 0.0d0
        end if
        
        if ( hlnnz + 5 .gt. lim ) then
           inform = -95
           return
        end if
        
        ! Diagonal and row-vector blocks.
        
        hlnnz = hlnnz + 1
        block_diag_ind(1) = hlnnz
        hlrow(hlnnz) = 2*i - 1
        hlcol(hlnnz) = 2*i - 1
        
        hlnnz = hlnnz + 1
        block_diag_ind(2) = hlnnz
        hlrow(hlnnz) = 2*i
        hlcol(hlnnz) = 2*i - 1
        
        hlnnz = hlnnz + 1
        block_diag_ind(3) = hlnnz
        hlrow(hlnnz) = 2*i
        hlcol(hlnnz) = 2*i
        
        block_diag_val(1:3) = 0.0d0
        
        hlnnz = hlnnz + 1
        block_row_vector_ind(1) = hlnnz
        hlrow(hlnnz) = n
        hlcol(hlnnz) = 2*i - 1
        
        hlnnz = hlnnz + 1
        block_row_vector_ind(2) = hlnnz
        hlrow(hlnnz) = n
        hlcol(hlnnz) = 2*i
        
        block_row_vector_val(1:2) = 0.0d0
        
        if ( pdata%nballs .ge. 3 .and. ierror .eq. 0 ) then
        
            call voronoi_cell(i, pdata%nballs, centers, v_num, nod_tri, tnbr, vor_xy, vorCell, status)
            
            if ( status .ne. 0 ) then
                inform = -53
                return
            end if
            
        end if
        
        do indpol=1,pdata%npols
            
            convPol%n = pdata%nvpols(indpol)
            convPol%deg = .false.
             
            allocate ( convPol%vertex(2,convPol%n), convPol%edges(convPol%n - 1), stat=allocerr )
             
            if ( allocerr .ne. 0 ) then
               write (*,*) 'Allocation error.'
               stop
            end if
            
            convPol%vertex(1,1:convPol%n) = pdata%xpol(indpol,1:convPol%n)
            convPol%vertex(2,1:convPol%n) = pdata%ypol(indpol,1:convPol%n)
            convPol%edges(1:convPol%n - 1) = pdata%poledges(indpol, 1:convPol%n - 1)
            
            if ( pdata%nballs .ge. 3 .and. ierror .eq. 0 ) then
                call VorCellInterConvPol(centers(:,i), convPol, vorCell, curCell)
            end if
            
            if ( pdata%nballs .eq. 2 .or. (pdata%nballs .ge. 3 .and. ierror .eq. 225) ) then
                call VorCellInterConvPol_Collinear(i, pdata%nballs, centers, convPol, curCell)
            end if
            
            if ( .not. curCell%deg .or. pdata%nballs .eq. 1 ) then
                
                if ( pdata%nballs .eq. 1 ) then
                    call interCellBall(centers(1:2,i), r, convPol, curvPol)
                else
                    call interCellBall(centers(1:2,i), r, curCell, curvPol)
                end if
                
                if ( curvPol%deg ) then
                    
                    if ( pdata%nballs .eq. 1 ) then
                        call polygon_contains_point_2d_convex(convPol%n - 1,convPol%vertex(1:2,1:convPol%n - 1),centers(1:2,i),&
                        inside)
                    else
                        call polygon_contains_point_2d_convex(curCell%n - 1,curCell%vertex(1:2,1:curCell%n - 1),centers(1:2,i),&
                        inside)
                    end if
                    
                    if ( inside ) then
                        hlval(1) = hlval(1) - 2.0d0 * PI
                    end if
                    
                else
                
                    do v = 1,curvPol%n - 1
                        
                        if ( curvPol%edges(v) .eq. 1 ) then ! Connection is an arc
                            
                            theta1 = angle_rad_2d( curvPol%vertex(:,v), centers(:,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /) )
                            
                            theta2 = angle_rad_2d( curvPol%vertex(:,v+1), centers(:,i), &
                            (/ centers(1,i) + 1.0d0, centers(2,i) /) )
                            
                            if ( theta2 .lt. theta1) then
                                theta2 = theta2 + 2.0d0 * PI
                            end if
                            
                            ! Sec. order derivative with respect to the radius (perimeter).
                            hlval(1) = hlval(1) + theta1 - theta2
                            
                            ! Integrals.
                            
                            ! Diagonal blocks:
                            block_diag_val(1) = block_diag_val(1) + dsin(theta1 - theta2)*dcos(theta1 + theta2)
                            block_diag_val(2) = block_diag_val(2) + dcos(theta2)**2 - dcos(theta1)**2
                            block_diag_val(3) = block_diag_val(3) + dsin(theta2 - theta1)*dcos(theta1 + theta2)

                            ! Two-dimensional row vector:
                            block_row_vector_val(1) = block_row_vector_val(1) + dsin(theta1) - dsin(theta2)
                            block_row_vector_val(2) = block_row_vector_val(2) + dcos(theta2) - dcos(theta1)
                            
                            ! Starting points
                            if ( curvPol%vertex_id(v) .lt. 0 ) then ! intersection with another ball
                                
                                theta_set_inter = angle_rad_2d( curvPol%vertex(:,v), centers(:,iabs(curvPol%vertex_id(v))), &
                                (/ centers(1,iabs(curvPol%vertex_id(v))) + 1.0d0, centers(2,iabs(curvPol%vertex_id(v))) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) - ((dcos(theta_set_inter - theta1) - 1.0d0) / dsin(theta_set_inter - theta1))
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) - (1.0d0 / dtan(theta_set_inter - theta1)) &
                                                    * (dcos(theta1)**2)
                                block_diag_val(2) = block_diag_val(2) - (1.0d0 / dtan(theta_set_inter - theta1)) &
                                                    * dsin(theta1) * dcos(theta1)
                                block_diag_val(3) = block_diag_val(3) - (1.0d0 / dtan(theta_set_inter - theta1)) &
                                                    * (dsin(theta1)**2)
                                
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) &
                                                            - (1.0d0 / dtan(theta_set_inter - theta1)) * dcos(theta1) &
                                                            + (dcos(theta1) / dsin(theta_set_inter - theta1))
                                block_row_vector_val(2) = block_row_vector_val(2) &
                                                            - (1.0d0 / dtan(theta_set_inter - theta1)) * dsin(theta1) &
                                                            + (dsin(theta1) / dsin(theta_set_inter - theta1))
                                
                                ! Derivative with respect to the intersecting balls centers.
                                if ( iabs(curvPol%vertex_id(v)) .gt. i ) then
                                    if ( intersections(iabs(curvPol%vertex_id(v)) - i) .eq. 0 ) then
                                    
                                        intersections(iabs(curvPol%vertex_id(v)) - i) = 1
                                        
                                        if ( hlnnz + 4 .gt. lim ) then
                                            inform = -95
                                            return
                                        end if
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(1,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v)) - 1
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(2,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v))
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(3,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v)) - 1
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(4,iabs(curvPol%vertex_id(v)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v))
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                    else
                                        
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dcos(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dcos(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                        
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v)) - i) + &
                                        ( dsin(theta1) * dsin(theta_set_inter) / dsin(theta_set_inter - theta1) )
                                    end if
                                end if
                            end if
                            
                            if ( curvPol%vertex_id(v) .eq. 0 ) then ! intersection with the boundary of A
                                
                                ! Computing the normal vector to the boundary of A.
                                if ( v .eq. 1 ) then
                                    pbA(1:2) = curvPol%vertex(1:2,curvPol%n - 1)
                                else
                                    pbA(1:2) = curvPol%vertex(1:2,v - 1)
                                end if
                                call line_exp_normal_2d ( curvPol%vertex(:,v), pbA, normal )
                                
                                inner_prod_div = dot_product( normal, (/ dcos(theta1), dsin(theta1) /) ) / &
                                                dot_product( normal, (/ - dsin(theta1), dcos(theta1) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) - inner_prod_div
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) - inner_prod_div * (dcos(theta1)**2)
                                block_diag_val(2) = block_diag_val(2) - inner_prod_div * dsin(theta1) * dcos(theta1)
                                block_diag_val(3) = block_diag_val(3) - inner_prod_div * (dsin(theta1)**2)
                            
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) - inner_prod_div * dcos(theta1)
                                block_row_vector_val(2) = block_row_vector_val(2) - inner_prod_div * dsin(theta1)
                            end if
                            
                            ! Ending points
                            if ( curvPol%vertex_id(v+1) .lt. 0 ) then ! intersection with another ball
                                
                                theta_set_inter = angle_rad_2d( curvPol%vertex(:,v+1), centers(:,iabs(curvPol%vertex_id(v+1))), &
                                (/ centers(1,iabs(curvPol%vertex_id(v+1))) + 1.0d0, centers(2,iabs(curvPol%vertex_id(v+1))) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) - ((1.0d0 - dcos(theta_set_inter - theta2)) / dsin(theta_set_inter - theta2))
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) + (1.0d0 / dtan(theta_set_inter - theta2)) &
                                                    * (dcos(theta2)**2)
                                block_diag_val(2) = block_diag_val(2) + (1.0d0 / dtan(theta_set_inter - theta2)) &
                                                    * dsin(theta2) * dcos(theta2)
                                block_diag_val(3) = block_diag_val(3) + (1.0d0 / dtan(theta_set_inter - theta2)) &
                                                    * (dsin(theta2)**2)
                                
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) &
                                                            + (1.0d0 / dtan(theta_set_inter - theta2)) * dcos(theta2) &
                                                            - (dcos(theta2)/dsin(theta_set_inter - theta2))
                                block_row_vector_val(2) = block_row_vector_val(2) &
                                                            + (1.0d0 / dtan(theta_set_inter - theta2)) * dsin(theta2) &
                                                            - (dsin(theta2)/dsin(theta_set_inter - theta2))
                                
                                ! Derivative with respect to the intersecting balls centers.
                                if ( iabs(curvPol%vertex_id(v+1)) .gt. i ) then
                                    if ( intersections(iabs(curvPol%vertex_id(v+1)) - i) .eq. 0 ) then
                                    
                                        intersections(iabs(curvPol%vertex_id(v+1)) - i) = 1
                                        
                                        if ( hlnnz + 4 .gt. lim ) then
                                            inform = -95
                                            return
                                        end if
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(1,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1)) - 1
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(2,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1))
                                        hlcol(hlnnz) = 2*i - 1
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )

                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(3,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1)) - 1
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        hlnnz = hlnnz + 1
                                        blocks_below_diag_ind(4,iabs(curvPol%vertex_id(v+1)) - i) = hlnnz
                                        hlrow(hlnnz) = 2*iabs(curvPol%vertex_id(v+1))
                                        hlcol(hlnnz) = 2*i
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                    else
                                        
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(1,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(2,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dcos(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(3,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dcos(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                        
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) = &
                                        blocks_below_diag_vals(4,iabs(curvPol%vertex_id(v+1)) - i) - &
                                        ( dsin(theta2) * dsin(theta_set_inter) / dsin(theta_set_inter - theta2) )
                                    end if
                                end if
                            end if
                            
                            if ( curvPol%vertex_id(v+1) .eq. 0 ) then ! intersection with the boundary of A
                                
                                ! Computing the normal vector to the boundary of A.
                                if ( v+1 .eq. curvPol%n ) then
                                    pbA(1:2) = curvPol%vertex(1:2,2)
                                else
                                    pbA(1:2) = curvPol%vertex(1:2,v + 2)
                                end if
                                call line_exp_normal_2d ( pbA, curvPol%vertex(:,v+1), normal )
                                
                                inner_prod_div = dot_product( normal, (/ dcos(theta2), dsin(theta2) /) ) / &
                                                dot_product( normal, (/ - dsin(theta2), dcos(theta2) /) )
                                
                                ! Sec. order derivative with respect to the radius
                                hlval(1) = hlval(1) + inner_prod_div
                                
                                ! Diagonal blocks
                                block_diag_val(1) = block_diag_val(1) + inner_prod_div * (dcos(theta2)**2)
                                block_diag_val(2) = block_diag_val(2) + inner_prod_div * dsin(theta2) * dcos(theta2)
                                block_diag_val(3) = block_diag_val(3) + inner_prod_div * (dsin(theta2)**2)
                                
                                ! Two-dimensional row vector
                                block_row_vector_val(1) = block_row_vector_val(1) + inner_prod_div * dcos(theta2)
                                block_row_vector_val(2) = block_row_vector_val(2) + inner_prod_div * dsin(theta2)
                            end if
                        end if
                    end do
                end if
            end if
            
            deallocate( convPol%vertex,convPol%edges,stat=allocerr )
            
            if ( allocerr .ne. 0 ) then
               write (*,*) 'Deallocation error.'
               stop
            end if
            
            if ( allocated(curCell%vertex) ) deallocate(curCell%vertex)
            if ( allocated(curCell%edges) ) deallocate(curCell%edges)
            if ( allocated(curvPol%vertex) ) deallocate(curvPol%vertex)
            if ( allocated(curvPol%vertex_id) ) deallocate(curvPol%vertex_id)
            if ( allocated(curvPol%edges) ) deallocate(curvPol%edges)
        end do
        if ( allocated(vorCell%vertex) ) deallocate(vorCell%vertex)
        if ( allocated(vorCell%edges) ) deallocate(vorCell%edges)
        
        hlval(block_diag_ind(1)) = block_diag_val(1)
        hlval(block_diag_ind(2)) = block_diag_val(2)
        hlval(block_diag_ind(3)) = block_diag_val(3)
        
        do j = 1,pdata%nballs - i
            if ( intersections(j) .eq. 1 ) then
                hlval(blocks_below_diag_ind(1,j)) = blocks_below_diag_vals(1,j)
                hlval(blocks_below_diag_ind(2,j)) = blocks_below_diag_vals(2,j)
                hlval(blocks_below_diag_ind(3,j)) = blocks_below_diag_vals(3,j)
                hlval(blocks_below_diag_ind(4,j)) = blocks_below_diag_vals(4,j)
            end if
        end do
        
        hlval(block_row_vector_ind(1)) = block_row_vector_val(1)
        hlval(block_row_vector_ind(2)) = block_row_vector_val(2)
        
        if ( allocated(intersections) ) deallocate(intersections)
        if ( allocated(blocks_below_diag_ind) ) deallocate(blocks_below_diag_ind)
        if ( allocated(blocks_below_diag_vals) ) deallocate(blocks_below_diag_vals)
    
    end do
    if ( allocated(vor_xy) ) deallocate(vor_xy)

    deallocate ( centers, indexes, nod_tri, tnbr, stat=allocerr )
    
    hlval(1:hlnnz) = hlval(1:hlnnz) * lambda(1)
    inform = 0
    
  end subroutine evalhl

  ! ******************************************************************
  ! ******************************************************************

end program algencama
