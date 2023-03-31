program main

  ! Get access to MPI library
  use MPI

  implicit none

  ! ================ !
  ! DEFINE VARIABLES !
  ! ================ !

  ! Fortran stuff
  integer, parameter :: WP=8! Working precision
  real(WP), parameter :: pi=3.1415926535897932385_WP

  ! MPI stuff
  integer :: iRank, iRoot, iProc, np, ndims, comm, ierror
  integer, dimension(3) :: dims, coords
  logical, dimension(3) :: isPeriodic
  logical :: reorder

  ! Grid info
  integer :: N, N_, nover
  integer :: imin, imin_, imino, imino_, imax, imax_, imaxo, imaxo_
  integer :: q, r
  real(WP), allocatable, dimension(:) :: x

  ! Local variables
  integer :: i
  real(WP), allocatable, dimension(:) :: phi, phi_deriv, ddx, phi_diff
  real(WP) :: delta_x
  real*8 :: t1, t2, t



  ! ============================== !
  ! SETUP THE PARALLEL ENVIRONMENT !
  ! ============================== !
  call MPI_Init(ierror)

  ! MPI global communicator
  comm = MPI_COMM_WORLD

  ! Get number of processors
  call MPI_Comm_Size(comm, np, ierror)

  ! Set MPI topology
  ndims = 1
  dims(1) = np
  dims(2:3) = 1
  isPeriodic = .true.
  
  !print *, dims

  reorder = .true.
  call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,isPeriodic,reorder,comm,ierror)
  call MPI_COMM_RANK(comm,irank,ierror)
  call MPI_CART_COORDS(comm,irank,ndims,coords,ierror)
  irank = irank + 1
  iproc = coords(1) + 1
  
  !print *, coords

  ! Define a root processor at coordinates (1,1,1)
  dims = 0
  call MPI_CART_RANK(comm,dims,iroot,ierror)
  iroot = iroot + 1


  
  ! =========== !
  ! HELLO WORLD !
  ! =========== !
  !print *, 'Hello, my processor rank is ', iRank


  ! ============== !
  ! SETUP THE GRID !
  ! ============== !

  ! Define number of grid points
  N=5000000

  ! Make sure number of processors < N
  if (np.gt.N) then
     print *, 'N has to be greater or equal to np'
     stop
  end if
  
  ! Number of overlapping cells
  nover=1

  ! Set the indices
  imino = 1  				! minimum global index
  imin  = imino + nover		! minimum global index, no overlap
  imax  = imin  + N - 1		! maximum global index, no overlap
  imaxo = imax  + nover		! maximum global index
  
  ! Partition the domain
  ! Set initial partitions along x
  q = N/np
  r = mod(N,np)
  
  !print *, q
  
  if (iproc<=r) then
     N_   = q+1
	 !print *, 'N_ for proc', irank, 'is: ', N_
	 !print *, 'I am in the q+1 loop for iproc',  iproc, 'My N_ is', N_
     imin_ = imin + (iproc-1)*(q+1)
  else
     N_   = q
	 !print *, 'N_ for proc', irank, 'is: ', N_
	 !print *, 'I am in the q loop for iproc',  iproc, 'My N_ is', N_
     imin_ = imin + r*(q+1) + (iproc-r-1)*q	 
  end if
  
  imax_  = imin_ + N_ - 1
  imino_ = imin_ - nover
  imaxo_ = imax_ + nover
  
  !print *, 'imin_ and imax_ is', imin_, imax_, 'for proc', iproc
  
  ! Allocate arrays
  allocate(x(imino_:imaxo_+1))
  
  ! Define the grid
  do i=imino_,imaxo_
     x(i) = 2.0_WP * pi * real(i-imin,WP)/real(N,WP)
  end do
  
  delta_x = abs(x(imino_) - x(imino_+1))
  !print*, delta_x
  !print *, x, ' for proc', iRank, 'With N_', N_
  
  ! ================================== !
  ! Compute the derivative in parallel !
  ! ================================== !

  ! Allocate arrays
  allocate(phi(imino_:imaxo_))
  allocate(ddx(imino_:imaxo_))
  allocate(phi_deriv(imino_:imaxo_))

  ! Define phi
  do i=imino_,imaxo_
     phi(i)=sin(x(i))
	 phi_deriv(i) = cos(x(i))
  end do
  
  !Allocating array
  allocate(phi_diff(imino_:imaxo_))
  
  !print *, imin_, imax_, 'For proc', iRank

  t1 = MPI_WTIME()
  ! Define 
  do i=imin_,imax_			!We go from imin_ to imax_ to not reference the i+1 place in the overall array
  	 phi_diff(i) = (phi(i+1) - phi(i))/delta_x
  end do
 
  !print *, 'Proc', iRank, ' has indices', imino_, ' to ', imaxo_
  !print *, x
  !print *, '----------------------------------------------'
  !print *, phi_diff
  !print *, '----------------------------------------------'
  !print *, phi_deriv
  
  t2 = MPI_WTIME()
  t = t2-t1
  
  print *, 'Time to compute in processor' , irank, 'is ', T * 1e6, ' ms'
 

  ! ============ !
  ! Finalize MPI !
  ! ============ !
  call MPI_FINALIZE(ierror)

end program main
