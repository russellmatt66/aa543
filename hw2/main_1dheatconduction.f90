! Matt Russell
! UW Department of Aeronautics & Astronautics
! AA543: Compressible CFD
! HW2
! 1/31/21

program ProblemOne
implicit none

  real, parameter :: xi = 0.0
  real, parameter :: xf = 1.0
  real, parameter :: dx = 0.05
  real, parameter :: dt1 = 0.0012
  real, parameter :: dt2 = 0.0013
  real, parameter :: c1 = dt1/(2.0*dx)
  real, parameter :: c2 = dt2/(2.0*dx)
  real, parameter :: pi = 4.0*atan(1.0)

  integer i, j
  integer, parameter :: N = 5
  integer, parameter :: Nx = (xf - xi)/dx + 1

  real, dimension(N) :: A_n

  real, dimension(Nx) :: x_grid, u_ic, u_fwd1_dt1, u_fwd10_dt1, u_fwd50_dt1, &
    u_fwd1_dt2, u_fwd10_dt2, u_fwd50_dt2, u_bkwd1_dt1, u_bkwd10_dt1, u_bkwd50_dt1, &
    u_bkwd1_dt2, u_bkwd10_dt2, u_bkwd50_dt2, u_anal1_dt1, u_anal10_dt1, u_anal50_dt1, &
    u_anal1_dt2, u_anal10_dt2, u_anal50_dt2

  real, dimension(Nx,Nx) :: CDstencil, M1, M2, M1inv, M2inv, Iden

  !Interface inv
  !  Function inv(A,N,M)
  !    integer,intent(in)                :: N,M
  !    real,intent(in),dimension(N,M)    :: A
  !    real,dimension(N,M)               :: inv
  !  end Function
  !end Interface

  ! print *, 'Number of grid points =', Nx
    print *, 'c1 =',c1
    print *, 'c2 =',c2

  ! initialize x_grid, initial conditions, and Central Difference stencil for
  ! implicit solve.

  CDstencil = 0.0
  Iden = 0.0
  M1 = 0.0 ! matrices whose inverse will be used for implicit update
  M2 = 0.0

  do i=1,size(x_grid)

    x_grid(i) = xi + dx*(i-1) ! dx*i should be 'real'

    if (x_grid(i) .LE. 0.5) then
      u_ic(i) = 2.0*x_grid(i)

      do j=1,N
        A_n(j) = (-1.0)**(1+j)*(4.0/(REAL(j)*pi))
      enddo

    else if (x_grid(i) .gt. 0.5) then

      u_ic(i) = 2.0 - 2.0*x_grid(i)

      do j=1,N
        A_n(j) = (4.0/(REAL(j)*pi))
      enddo
    end if

    Iden(i,i) = 1.0

    if (i == 2) then
      CDstencil(i,i) = -2.0
      CDstencil(i+1,i) = 1.0
    else if (i == Nx-1) then
      CDstencil(i,i) = -2.0
      CDstencil(i-1,i) = 1.0
    else if (i > 2 .and. i < Nx-1) then
      CDstencil(i,i) = -2.0
      CDstencil(i+1,i) = 1.0
      CDstencil(i-1,i) = 1.0
    end if
  enddo

  M1 = Iden - c1*CDstencil
  M2 = Iden - c2*CDstencil

  ! print*, 'M1 =', M1
  ! print*, 'M2 =', M2

  M1inv = inv(M1)
  M2inv = inv(M2)

  !print *, 'M1inv =', M1inv
  !print *, 'M2inv =', M2inv

  ! B.Cs: u(x=0,t) = u(x=1,t) = 0
  u_ic(1) = 0
  u_ic(size(u_ic)) = 0
  print *, 'u_ic =', u_ic
  ! print *, 'CDstencil =', CDstencil

  ! Initialize velocities w/initial conditions
  u_fwd1_dt1 = u_ic
  u_fwd10_dt1 = u_ic
  u_fwd50_dt1 = u_ic
  u_fwd1_dt2 = u_ic
  u_fwd10_dt2 = u_ic
  u_fwd50_dt2 = u_ic
  u_bkwd1_dt1 = u_ic
  u_bkwd10_dt1 = u_ic
  u_bkwd50_dt1 = u_ic
  u_bkwd1_dt2 = u_ic
  u_bkwd10_dt2 = u_ic
  u_bkwd50_dt2 = u_ic
  u_anal1_dt1 = u_ic
  u_anal10_dt1 = u_ic
  u_anal50_dt1 = u_ic
  u_anal1_dt2 = u_ic
  u_anal10_dt2 = u_ic
  u_anal50_dt2 = u_ic

  ! print *, 'velocities initialized succesfully'

  !! Single Step
  ! Explicit
  u_fwd1_dt1 = explicitstepCD(u_fwd1_dt1,c1)
  u_fwd1_dt2 = explicitstepCD(u_fwd1_dt2,c2)

  print *,'u_fwd1_dt1=', u_fwd1_dt1
  print *,'u_fwd1_dt2=', u_fwd1_dt2

  ! Implicit
  u_bkwd1_dt1 = matmul(M1inv,u_bkwd1_dt1)
  u_bkwd1_dt2 = matmul(M2inv,u_bkwd1_dt2)

  print *,'u_bkwd1_dt1 =', u_bkwd1_dt1
  print *,'u_bkwd1_dt2 =', u_bkwd1_dt2

  !! Ten Steps
  do i=1,10
    ! Explicit
    u_fwd10_dt1 = explicitstepCD(u_fwd10_dt1,c1)
    u_fwd10_dt2 = explicitstepCD(u_fwd10_dt2,c2)

    ! Implicit
    u_bkwd10_dt1 = matmul(M1inv,u_bkwd10_dt1)
    u_bkwd10_dt2 = matmul(M2inv,u_bkwd10_dt2)
  enddo

  print *,'u_fwd10_dt1=', u_fwd10_dt1
  print *,'u_fwd10_dt2=', u_fwd10_dt2

  print *,'u_bkwd10_dt1 =', u_bkwd10_dt1
  print *,'u_bkwd10_dt2 =', u_bkwd10_dt2

  !! Fifty Steps
  do i=1,50
    ! Explicit
    u_fwd50_dt1 = explicitstepCD(u_fwd50_dt1,c1)
    u_fwd50_dt2 = explicitstepCD(u_fwd50_dt2,c2)

    ! Implicit
    u_bkwd50_dt1 = matmul(M1inv,u_bkwd50_dt1)
    u_bkwd50_dt2 = matmul(M2inv,u_bkwd50_dt2)
  enddo

  print *,'u_fwd50_dt1=', u_fwd50_dt1
  print *,'u_fwd50_dt2=', u_fwd50_dt2

  print *,'u_bkwd50_dt1 =', u_bkwd50_dt1
  print *,'u_bkwd50_dt2 =', u_bkwd50_dt2

  ! Analytical Solutions:
  do j=1,N
    do i=2,Nx-1
      u_anal1_dt1(i) = u_anal1_dt1(i) + A_n(j)*exp(-(REAL(j)**2)*(pi**2)*dt1)*sin(REAL(j)*pi*x_grid(i))
      u_anal10_dt1(i) = u_anal10_dt1(i) + A_n(j)*exp(-(REAL(j)**2)*(pi**2)*10.0*dt1)*sin(REAL(j)*pi*x_grid(i))
      u_anal50_dt1(i) = u_anal50_dt1(i) + A_n(j)*exp(-(REAL(j)**2)*(pi**2)*50.0*dt1)*sin(REAL(j)*pi*x_grid(i))
      u_anal1_dt2(i) = u_anal1_dt2(i) + A_n(j)*exp(-(REAL(j)**2)*(pi**2)*dt2)*sin(REAL(j)*pi*x_grid(i))
      u_anal10_dt2(i) = u_anal10_dt2(i) + A_n(j)*exp(-(REAL(j)**2)*(pi**2)*10.0*dt2)*sin(REAL(j)*pi*x_grid(i))
      u_anal50_dt2(i) = u_anal50_dt2(i) + A_n(j)*exp(-(REAL(j)**2)*(pi**2)*50.0*dt1)*sin(REAL(j)*pi*x_grid(i))
    end do
  end do

  !! Write to output file for gnuplot
  open(7, file = 'dat/fwddata1_dt1.dat')
  open(8, file = 'dat/fwddata10_dt1.dat')
  open(9, file = 'dat/fwddata50_dt1.dat')
  open(10, file = 'dat/fwddata1_dt2.dat')
  open(11, file = 'dat/fwddata10_dt2.dat')
  open(12, file = 'dat/fwddata50_dt2.dat')

  open(13, file = 'dat/bkwddata1_dt1.dat')
  open(14, file = 'dat/bkwddata10_dt1.dat')
  open(15, file = 'dat/bkwddata50_dt1.dat')
  open(16, file = 'dat/bkwddata1_dt2.dat')
  open(17, file = 'dat/bkwddata10_dt2.dat')
  open(18, file = 'dat/bkwddata50_dt2.dat')

  open(19, file = 'dat/analdata1_dt1.dat')
  open(20, file = 'dat/analdata10_dt1.dat')
  open(21, file = 'dat/analdata50_dt1.dat')
  open(22, file = 'dat/analdata1_dt2.dat')
  open(23, file = 'dat/analdata10_dt2.dat')
  open(24, file = 'dat/analdata50_dt2.dat')

  open(25, file = 'dat/x_grid.dat')


  do i=1,Nx
    write(7,'(f14.12)') u_fwd1_dt1(i)
    write(8,'(f14.12)') u_fwd10_dt1(i)
    write(9,'(f14.12)') u_fwd50_dt1(i)
    write(10,'(f14.12)') u_fwd1_dt2(i)
    write(11,'(f14.12)') u_fwd10_dt2(i)
    write(12,'(f14.12)') u_fwd50_dt2(i)

    write(13,'(f14.12)') u_bkwd1_dt1(i)
    write(14,'(f14.12)') u_bkwd10_dt1(i)
    write(15,'(f14.12)') u_bkwd50_dt1(i)
    write(16,'(f14.12)') u_bkwd1_dt2(i)
    write(17,'(f14.12)') u_bkwd10_dt2(i)
    write(18,'(f14.12)') u_bkwd50_dt2(i)

    write(19,'(f14.12)') u_anal1_dt1(i)
    write(20,'(f14.12)') u_anal10_dt1(i)
    write(21,'(f14.12)') u_anal50_dt1(i)
    write(22,'(f14.12)') u_anal1_dt2(i)
    write(23,'(f14.12)') u_anal10_dt2(i)
    write(24,'(f14.12)') u_anal50_dt2(i)

    write(25,'(f14.12)') x_grid(i)
  end do

do j=7,25
  close(j)
end do

contains
  function explicitstepCD(u,c) result(u_out)
    implicit none
    real :: u(:), u_out(size(u))
    real :: c
    integer :: i,n

    n = size(u)
    do i=2,n-1
      u(i)=u(i) + c*(u(i+1)-2*u(i) + u(i-1))
    enddo

    u_out = u
  end function explicitstepCD

  ! https://stackoverflow.com/a/43766617/11619514
  function inv(A) result(Ainv)
    implicit none
    real,intent(in) :: A(:,:)
    real            :: Ainv(size(A,1),size(A,2))
    real            :: work(size(A,1))            ! work array for LAPACK
    integer         :: n,info,ipiv(size(A,1))     ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call SGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call SGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function inv

end program ProblemOne
