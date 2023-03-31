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
  real(WP), allocatable, dimension(:) :: phi, ddx
  
  ! Homework variables
  real(WP), allocatable, dimension(:) :: a
  real*8 :: t1, t2, t, j


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

  reorder = .true.
  call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,isPeriodic,reorder,comm,ierror)
  call MPI_COMM_RANK(comm,irank,ierror)
  call MPI_CART_COORDS(comm,irank,ndims,coords,ierror)
  irank = irank + 1
  iproc = coords(1) + 1

  ! Define a root processor at coordinates (1,1,1)
  dims = 0
  call MPI_CART_RANK(comm,dims,iroot,ierror)
  iroot = iroot + 1
  
  ! ============================================ !
  ! HW 1 Part B: Testing Program Completion Time !
  ! ============================================ !
  
  ! Define number of points
  !  For now, keep "N" a multiple of 30
  N=30000

  ! Make sure number of processors < N
  if (np.gt.N) then
     print *, 'N has to be greater or equal to np'
     stop
  end if
  
  ! Setting indices
  imino = 1
  imin  = imino + nover
  imax  = imin  + N - 1
  imaxo = imax  + nover
  
  ! Partitioning environment
  !   Note that variables with an underscore after the name pertains to 
  !    a value that is tracking the size of the partitioned portion of the entire domain
  q = N/np
  r = mod(N,np)
  if (iproc<=r) then
     N_   = q+1
     imin_ = imin + (iproc-1)*(q+1)
  else
     N_   = q
     imin_ = imin + r*(q+1) + (iproc-r-1)*q
  end if
  imax_  = imin_ + N_ - 1
  imino_ = imin_ - nover
  imaxo_ = imax_ + nover
  
  !print *, irank, N_
  
  ! Allocate arrays
  allocate(x(imino_:imaxo_))
  allocate(a(imino_:imaxo_))
  !	 There is no need for information to be communicated across indices for this portion of the homework.
  !		That being said, I chose to leave the setup code as close as possible to the given code 
  !		for my own sanity's sake. 
  !	 Long story short, I removed a "+1" from imaxo_ that was leading to some weird values in the allocated arrays
  
  !print *, 'Processor ', iRank, ' owns the array containing: ', x
  
  ! Defining the grid for later computation
  do i = imino_, imaxo_
  	x(i) = i
  end do
  !print *, x
  
  ! Actually, finally, mercifully taking a=exp(b)
  t1 = MPI_WTIME()
  do i =imino_, imaxo_
  	a(i) = exp(x(i)) 
  end do
  !print *, a
  
  t2 = MPI_WTIME()
  t = t2-t1
  
  print *, 'Time to compute in processor' , irank, 'is ', T * 1e6, ' ms'
  
    
  ! ============ !
  ! Finalize MPI !
  ! ============ !
  call MPI_FINALIZE(ierror)

end program main
