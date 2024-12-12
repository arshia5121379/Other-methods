!===========================================================================
!            2D SHSLBM program for hydrodynamic problems
!                            Zhen Chen
!                  National University of Singapore
!                             2017.1
!===========================================================================

!===========================================================================
!parameter settings
!---------------------------------------------------------------------------
module param
  implicit none
  integer,parameter ::  dim = 2,  &     ! dimension
                        q = 9,    &     ! D2Q9
                        im = 251, &     ! grid number in x-direction
                        jm = 251        ! grid number in y-direction
  real(8),parameter ::  rho0 = 1.0,  &  ! reference density
                        Re=7500.0,    &  ! Reynolds number
                        eta=5.0e-8      ! convergent criteria
end module

!===========================================================================
! variables
!---------------------------------------------------------------------------
module fluid_common
  use param
  implicit none
  integer :: i, j, k, itime, i_output
  real(8) :: dt, time
  real(8) :: pi, niu, tau, cs2, err
  real(8) :: u_lid
  real(8) :: vlx, vly
  real(8) :: start, finish
  real(8), dimension(im, jm) :: x, y        ! grid position
  real(8), dimension(im, jm) :: ux, uy, uu, uxo, uyo  ! velocity in x, y direction, total velocity
  real(8), dimension(im, jm) :: rho, rhoo, str, p ! density, stream function, pressure
  real(8), dimension(im, jm) :: rho_s, ux_s, uy_s ! intermediate density and velocity
  real(8) :: ex(0:q-1), ey(0:q-1), w(0:q-1) ! lattice velocities and corresponding weights
end module


!===========================================================================
!main program
!---------------------------------------------------------------------------
Program LBMcz
  use fluid_common
  implicit none
  
  open (50, file='SHSLBMcz.log')
  
  call getdata  ! initial input

  call cpu_time(start)

  call SHSLBM  ! numerical resolution using SHSLBM

  call cpu_time(finish)

  call final_output  ! output of computational results
  
  close(50)
  
end program

!===========================================================================
!SHSLBM evolutions
!---------------------------------------------------------------------------
subroutine SHSLBM
  use fluid_common
  implicit none
  integer :: ip, jp
  real(8) :: ff, f1, f2

  time=0.0  ! initialize physical time
  itime=0  ! initialize iteration number
  err=100.0  ! set a large value of initial error to start iteration
  do while (err > eta)
    
    ! Predictor step: rho_star, u_star
    rho_s=0.0
    ux_s=0.0
    uy_s=0.0
    
    do i=2,im-1
      do j=2,jm-1
        rhoo(i,j)=rho(i,j)
        uxo(i,j)=ux(i,j)
        uyo(i,j)=uy(i,j)
        do k=0,q-1
          ip=i-nint(ex(k))
          jp=j-nint(ey(k))
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),k,ff)
          rho_s(i,j)=rho_s(i,j)+ff  ! evaluation of intermediate density
          ux_s(i,j)=ux_s(i,j)+ff*ex(k)  ! evaluation of intermediate momentum
          uy_s(i,j)=uy_s(i,j)+ff*ey(k)
        enddo
        ux_s(i,j)=ux_s(i,j)/rho_s(i,j)  ! calculate intermediate velocity
        uy_s(i,j)=uy_s(i,j)/rho_s(i,j)
      enddo
    enddo
	
	! B.C. for intermediate variables
    ! left & right walls
    do j=2,jm-1
      rho_s(1,j)=(4.0*rho_s(2,j)-rho_s(3,j))/3.0
      ux_s(1,j)=0.0
      uy_s(1,j)=0.0

      rho_s(im,j)=(4.0*rho_s(im-1,j)-rho_s(im-2,j))/3.0
      ux_s(im,j)=0.0
      uy_s(im,j)=0.0
    enddo

    ! top & bottom walls
    do i=1,im
      rho_s(i,1)=(4.0*rho_s(i,2)-rho_s(i,3))/3.0
      ux_s(i,1)=0.0
      uy_s(i,1)=0.0

      rho_s(i,jm)=(4.0*rho_s(i,jm-1)-rho_s(i,jm-2))/3.0
      ux_s(i,jm)=u_lid
      uy_s(i,jm)=0.0
    enddo

    ! Corrector step: rho_n+1, u_n+1
    do i=2,im-1
      do j=2,jm-1
        ux(i,j)=ux_s(i,j)*rho_s(i,j)
        uy(i,j)=uy_s(i,j)*rho_s(i,j)
        do k=0,q-1
          ip=i+nint(ex(k))
          jp=j+nint(ey(k))
          call equilibrium(rho_s(ip,jp),ux_s(ip,jp),uy_s(ip,jp),k,ff)  ! calculate equilibria using intermediate flow variables
          ux(i,j)=ux(i,j)+(tau-1.0)*ex(k)*ff  !momentum correction
          uy(i,j)=uy(i,j)+(tau-1.0)*ey(k)*ff
        enddo
        ux(i,j)=(ux(i,j)-(tau-1.0)*rhoo(i,j)*uxo(i,j))/rho_s(i,j) !momentum correction and calculation of velocity
        uy(i,j)=(uy(i,j)-(tau-1.0)*rhoo(i,j)*uyo(i,j))/rho_s(i,j)

        rho(i,j)=rho_s(i,j)  ! no correction on the intermediate density
        p(i,j)=cs2*rho(i,j)  ! update pressure of EoS
      enddo
    enddo

    ! Implementation of physical boundary condition
    call boundary

    ! Convergence criteria
    call convergence_test

    time=time+dt  ! update physical time
    itime=itime+1  ! update iteration number

    if(mod(itime,i_output)==0) then  !output flow field; output and record relative errors
      call output
      write(*,*) 'Relative error: ', err
      print *, ' Chen: output at t=', time, 'finished ! '
      print *, '==============================================='
      
      write(50,*) 'Relative error: ', err
      write(50,*) ' Chen: output at t=', time, 'finished ! '
      write(50,*) '==============================================='
    endif
  enddo
  
end subroutine

!===========================================================================
!boundary treatment
!---------------------------------------------------------------------------
subroutine boundary
  use fluid_common
  implicit none

  ! left & right walls
  do j=2,jm-1
    rho(1,j)=(4.0*rho(2,j)-rho(3,j))/3.0
    p(1,j)=rho(1,j)*cs2
    ux(1,j)=0.0
    uy(1,j)=0.0

    rho(im,j)=(4.0*rho(im-1,j)-rho(im-2,j))/3.0
    p(im,j)=rho(im,j)*cs2
    ux(im,j)=0.0
    uy(im,j)=0.0
  enddo

  ! top & bottom walls
  do i=1,im
    rho(i,1)=(4.0*rho(i,2)-rho(i,3))/3.0
    p(i,1)=rho(i,1)*cs2
    ux(i,1)=0.0
    uy(i,1)=0.0
    
    rho(i,jm)=(4.0*rho(i,jm-1)-rho(i,jm-2))/3.0
    p(i,jm)=rho(i,jm)*cs2
    ux(i,jm)=u_lid
    uy(i,jm)=0.0
  enddo

end subroutine

!===========================================================================
!calculating equilibrium term
!---------------------------------------------------------------------------
subroutine equilibrium(ro,u1,u2,kk,ff)
  use fluid_common
  implicit none
  real(8) :: factor, ff, ro, u1, u2
  integer :: kk

  factor=ex(kk)*u1+ey(kk)*u2
  ff=w(kk)*ro*(1.0+factor/cs2+(factor/cs2)**2/2.0-(u1**2+u2**2)/cs2/2.0)
  
end subroutine

!===========================================================================
!test of convergent criteria
!---------------------------------------------------------------------------
subroutine convergence_test
  use fluid_common
  implicit none
  real(8) :: du

  err=0.0
  do i=1,im
    do j=1,jm-1
      du=abs(uu(i,j)-sqrt(ux(i,j)**2+uy(i,j)**2))/uu(i,j)
      if(du>err) err=du
      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
    enddo
  enddo

end subroutine

!===========================================================================
!output
!---------------------------------------------------------------------------
subroutine output
  use fluid_common
  implicit none
  character(30) :: tname

  ! calculate the stream function
  do i=1,im
    str(i,1)=0.0
  enddo

  do i=1,im
    do j=2,jm
      str(i,j)=str(i,j-1)+dt*(ux(i,j)+ux(i,j-1))/2.0
    enddo
  enddo
  
  ! output data file for the flow field. The file can be opened by Tecplot.
  write(tname, '(A, F12.6, A)') 't=', time, 's.dat'
  open(7, file = tname)
  write(7, '(A, E13.6,A)') ' title="t=', time, ' s"'
  write(7,*) 'variables = x, y, u, v, p, rho, str'
22  format(7E15.6)
  write(7,'(A,A6,I,A6,I,A10)') 'zone T="BOX"', "I=",im, "J=",jm, 'F=point'
  do i=1, im
    do j=1, jm
	  write(7, 22) x(i,j), y(i,j), ux(i,j), uy(i,j), p(i,j), rho(i,j), str(i,j)
    enddo
  enddo

  close(7)
end subroutine

!===========================================================================
!final output
!---------------------------------------------------------------------------
subroutine final_output
  use fluid_common
  implicit none
  integer :: imid, jmid

  call output  ! output of converged flow field

  imid=nint((im+1)/2.0)
  jmid=nint((jm+1)/2.0)

  open(10, file = 'X-V.dat')
  open(20, file = 'Y-U.dat')

21  format(2E15.6)
  
  ! output converged velocity profile along the horizontal centerline
  do i=1, im
    write(10, 21) x(i,jmid)/vlx, uy(i,jmid)/u_lid
  enddo

  ! output converged velocity profile along the vertical centerline
  do j=1, jm
    write(20, 21) ux(imid,j)/u_lid, y(imid,j)/vly
  enddo

  ! record running time and iteration number
  write(50,*) 'Iteration steps:', itime
  write(50,*) 'Converged time t=', time, 's'
  write(50,*) 'CPU time:', finish-start, 's'

  write(*,*) 'Convergence criterion reached!'
end subroutine

!===========================================================================
!initial input
!---------------------------------------------------------------------------
subroutine getdata
  use fluid_common
  implicit none
  real(8) :: factor

  print *, '==============================================='
  print *, '   2D SHSLBM program for hydrodynamic problems'
  print *, '               version 1.0'
  print *, '           programmed by: Zhen Chen'
  print *, '       National University of Singapore'
  print *, '==============================================='

  vlx=1.0 ! length of the cavity
  vly=1.0
  u_lid=0.1 ! moving velocity of the lid
  cs2=1.0/3.0  !square of sound speed
  pi=4.0*atan(1.0)
  dt=vlx/(im-1) ! time step
  i_output=nint(500.0/dt)  ! output time interval
  niu=u_lid*vlx/Re  ! kinematic viscosity
  tau= 0.5+niu/cs2/dt  ! relaxation parameter in BGK model

  ex(0:8) = (/0,1,0,-1,0,1,-1,-1,1/)
  ey(0:8) = (/0,0,1,0,-1,1,1,-1,-1/)
  w(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)

  do i=1, im
    do j=1, jm
      x(i, j)=(i-1)*vlx/(im-1)
      y(i, j)=(j-1)*vly/(jm-1)
      if(j==jm) then
        ux(i,j)=u_lid  ! mesh point on the top lid
      else
        ux(i,j)=0.0
      endif
      uy(i,j)=0.0
      rho(i,j)=rho0
      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
    enddo
  enddo

  print *, 'Demo case: 2D lid-driven cavity flow'
  write(*,*) 'Re=', Re, 'tau=', tau
  write(50,*) 'Demo case: 2D lid-driven cavity flow'
  write(50,*) 'Re=', Re, 'tau=', tau
  
  call output

end subroutine