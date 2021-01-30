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

  integer i, N
  integer, parameter :: Nx = (xf - xi)/dx + 1

  real, dimension(Nx) :: x_grid, u_ic, u_fwd1_dt1, u_fwd10_dt1, u_fwd50_dt1, &
    u_fwd1_dt2, u_fwd10_dt2, u_fwd50_dt2, u_bkwd1_dt1, u_bkwd10_dt1, u_bkwd50_dt1, &
    u_bkwd1_dt2, u_bkwd10_dt2, u_bkwd50_dt2

  real, dimension(Nx,Nx) :: CDstencil, M, I

  ! print *, 'Number of grid points =', Nx

  ! initialize x_grid, initial conditions, and Central Difference stencil for
  ! implicit solve.
  CDstencil = 0.0
  I = 0.0
  M = 0.0 ! matrix whose inverse will be used for implicit update

  do i=1,size(x_grid)
    x_grid(i) = xi + dx*i ! dx*i should be 'real'
    if (x_grid(i) .LE. 0.5) then
      u_ic(i) = 2.0*x_grid(i)
    else if (x_grid(i) .gt. 0.5) then
      u_ic(i) = 2.0 - 2.0*x_grid(i)
    end if

    I(i,i) = 1

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

  ! print *, 'u_ic =', u_ic
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

  ! print *, 'velocities initialized succesfully'

  !! Single Step
  ! Explicit
  do i = 2,Nx-1
    u_fwd1_dt1(i) = u_fwd1_dt1(i) + dt1/(2.0*dx)*(u_fwd1_dt1(i+1) - 2*u_fwd1_dt1(i) &
      + u_fwd1(i-1))
    u_fwd1_dt2(i) = u_fwd1_dt2(i) + dt2/(2.0*dx)*(u_fwd1_dt2(i+1) - 2*u_fwd1_dt2(i) &
      + u_fwd2(i-1))
  enddo

  ! Implicit
  M = I - (dt1/(2.0*dx))*CDstencil


  !! Ten Steps

  !! Fifty Steps

end program ProblemOne
