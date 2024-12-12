!===========================================================================
!             2D HSLBM program for hydrodynamic problems
!                            Zhen Chen
!                  National University of Singapore
!                             2020.6
!===========================================================================

!===========================================================================
!parameter settings
!---------------------------------------------------------------------------
module param
  implicit none
  integer,parameter ::  dim = 2,  &     ! dimension
                        q = 9,    &     ! D2Q9
                        im = 21, &     ! grid number in x-direction
                        jm = 21, &     ! grid number in y-direction
                        nn = 5         ! number of sample points

  real(8),parameter ::  rho0 = 1.0,  &  ! reference density
                        Re=10.0         ! Reynolds number
end module

!===========================================================================
! variables
!---------------------------------------------------------------------------
module fluid_common
  use param
  implicit none
  integer :: i, j, k, itime, i_output
  integer :: i_in(1:im,1:nn), j_in(1:jm,1:nn) ! index of sample points in Lagrangian interpolation
  real(8) :: dt, time, tmax
  real(8) :: pi, niu, tau, cs2, err
  real(8) :: u0  ! characteristic velocity
  real(8) :: vlx, vly  ! size of the computational domain
  real(8) :: start, finish
  real(8), dimension(im, jm) :: x, y        ! mesh point position
  real(8), dimension(im, jm) :: ux, uy, uu  ! velocity in x, y direction, total velocity
  real(8), dimension(im, jm) :: rho, p ! density, pressure
  real(8), dimension(im, jm) :: rho_s, ux_s, uy_s  ! intermediate properties
  real(8) :: fneq(1:im, 1:jm, 0:q-1), feq(1:im, 1:jm, 0:q-1)  ! non-equilibrium and equilibrium parts of the distribution function
  real(8) :: ex(0:q-1), ey(0:q-1), w(0:q-1) ! lattice velocities and corresponding weights
  real(8) :: ff(0:q-1)  ! interpolated properties at the streaming nodes
  real(8) :: gg(-1:1,-1:1,1:q-1)  ! properties at the sample points
  real(8) :: a(1:im, -1:1, 1:nn), b(1:jm, -1:1, 1:nn)  ! Lagrangian interpolation polynomials
end module


!===========================================================================
!main program
!---------------------------------------------------------------------------
Program HSLBMcz
  use fluid_common
  implicit none

  open (60, file='Err_history.dat')
  
  call cpu_time(start)

  call getdata  ! initial input

  call HSLBM

  call cpu_time(finish)
  
end program

!===========================================================================
!HSLBM
!---------------------------------------------------------------------------
subroutine HSLBM
  use fluid_common
  implicit none
  integer :: ip, jp, ii, jj, l, ll
  real(8) :: f1, f2, xx, yy

  time=0.0
  itime=0
  err=100.0

  ! set up interpolation parameters
  do i=1,im
    do ii=1,nn
      l=i-(nn-1)/2+ii-1
      if(l<1) then
        i_in(i,ii)=l+im-1
      elseif(l>im) then
        i_in(i,ii)=l-im+1
      else
        i_in(i,ii)=l
      endif
    enddo
  enddo

  do j=1,jm
    do jj=1,nn
      l=j-(nn-1)/2+jj-1
      if(l<1) then
        j_in(j,jj)=l+jm-1
      elseif(l>jm) then
        j_in(j,jj)=l-jm+1
      else
        j_in(j,jj)=l
      endif
    enddo
  enddo

  a=1.0d0
  i=(im+1)/2
  do k=-1, 1
    xx=x(i,2)-dt*k
    do j=1,nn
      do l=1,nn
        if(l==j) cycle
        ii=i_in(i,j)
        jj=i_in(i,l)
        f1=xx-x(jj,2)
        if(abs(f1)<0.1*dt) f1=0.0d0
        a(i,k,j)=a(i,k,j)*f1/(x(ii,2)-x(jj,2))
      enddo
    enddo
  enddo
  
  do i=1, im
    if(i==(im+1)/2) cycle
    do k=-1, 1
      do j=1,nn
        a(i,k,j)=a((im+1)/2,k,j)
      enddo
    enddo
  enddo

  b=1.0d0
  i=(jm+1)/2
  do k=-1, 1
    yy=y(2,i)-dt*k
    do j=1,nn
      do l=1,nn
        if(l==j) cycle
        ii=j_in(i,j)
        jj=j_in(i,l)
        f1=yy-y(2,jj)
        if(abs(f1)<0.1*dt) f1=0.0d0
        b(i,k,j)=b(i,k,j)*f1/(y(2,ii)-y(2,jj))
      enddo
    enddo
  enddo

  do i=1, jm
    if(i==(jm+1)/2) cycle
    do k=-1, 1
      do j=1,nn
        b(i,k,j)=b((jm+1)/2,k,j)
      enddo
    enddo
  enddo
  ! end of interpolation setup

  ! start evolution
  do while (time<tmax)
    
    ! step 1: rho_star, u_star
    rho_s=0.0
    ux_s=0.0
    uy_s=0.0

    do i=1,im
      do j=1,jm
        do k=0, q-1
          call equilibrium(rho(i,j),ux(i,j),uy(i,j),k,feq(i,j,k))
        enddo
      enddo
    enddo

    do i=1,im
      do j=1,jm
        ! interpolate equilibrium distribution functions at streaming nodes
        ff=0.0
        ff(0)=feq(i,j,0)
        do l=1,nn
          ii=i_in(i,l)
          ff(1)=ff(1)+a(i,ex(1),l)*feq(ii,j,1)
          ff(3)=ff(3)+a(i,ex(3),l)*feq(ii,j,3)
        enddo
        do l=1,nn
          jj=j_in(j,l)
          ff(2)=ff(2)+b(j,ey(2),l)*feq(i,jj,2)
          ff(4)=ff(4)+b(j,ey(4),l)*feq(i,jj,4)
        enddo

        do k=5,q-1
          do l=1,nn
            do ll=1,nn
              ii=i_in(i,l)
              jj=j_in(j,ll)
              ff(k)= ff(k)+a(i,ex(k),l)*b(j,ey(k),ll)*feq(ii,jj,k)
            enddo
          enddo
        enddo

        ! calculate intermediate density and momentum
        do k=0,q-1
          rho_s(i,j)=rho_s(i,j)+ff(k)
          ux_s(i,j)=ux_s(i,j)+ff(k)*ex(k)
          uy_s(i,j)=uy_s(i,j)+ff(k)*ey(k)
        enddo
        ! calculate intermediate velocity
        ux_s(i,j)=ux_s(i,j)/rho_s(i,j)
        uy_s(i,j)=uy_s(i,j)/rho_s(i,j)

        ! calculate non-equilibrium distribution functions
        do k=0,q-1
          call equilibrium(rho_s(i,j),ux_s(i,j),uy_s(i,j),k,f1)
          f2=ff(k)

          fneq(i,j,k)=-tau*(f1-f2)
        enddo

      enddo
    enddo

    ! step 2: rho_n+1, u_n+1
    rho=rho_s
    do i=1,im
      do j=1,jm
        ! interpolate non-equilibrium distribution functions at streaming nodes
        ff=0.0
        ff(0)=fneq(i,j,0)
        do l=1,nn
          ii=i_in(i,l)
          ff(1)=ff(1)+a(i,-ex(1),l)*fneq(ii,j,1)
          ff(3)=ff(3)+a(i,-ex(3),l)*fneq(ii,j,3)
        enddo
        do l=1,nn
          jj=j_in(j,l)
          ff(2)=ff(2)+b(j,-ey(2),l)*fneq(i,jj,2)
          ff(4)=ff(4)+b(j,-ey(4),l)*fneq(i,jj,4)
        enddo

        do k=5,q-1
          do l=1,nn
            do ll=1,nn
              ii=i_in(i,l)
              jj=j_in(j,ll)
              ff(k)= ff(k)+a(i,-ex(k),l)*b(j,-ey(k),ll)*fneq(ii,jj,k)
            enddo
          enddo
        enddo

        ! update momentum
        ux(i,j)=ux_s(i,j)*rho_s(i,j)
        uy(i,j)=uy_s(i,j)*rho_s(i,j)
        do k=0,q-1
          ux(i,j)=ux(i,j)-(1.0-1.0/tau)*ex(k)*(ff(k)-fneq(i,j,k))
          uy(i,j)=uy(i,j)-(1.0-1.0/tau)*ey(k)*(ff(k)-fneq(i,j,k))
        enddo
        ! update velocity
        ux(i,j)=ux(i,j)/rho(i,j)
        uy(i,j)=uy(i,j)/rho(i,j)
      enddo
    enddo

    ! update pressure with EoS
    do i=1,im
      do j=1,jm
        p(i,j)=cs2*rho(i,j)
      enddo
    enddo

    time=time+dt
    itime=itime+1

    if(mod(itime,10)==0) call error_analysis
    
    if(mod(itime,i_output)==0) then
      call output
    endif

  enddo
  
end subroutine


!===========================================================================
!calculating equilibrium term
!---------------------------------------------------------------------------
subroutine equilibrium(ro,u1,u2,kk,fff)
  use fluid_common
  implicit none
  real(8) :: factor, fff, ro, u1, u2
  integer :: kk

  factor=ex(kk)*u1+ey(kk)*u2
  fff=w(kk)*ro*(1.0+factor/cs2+(factor/cs2)**2/2.0-(u1**2+u2**2)/cs2/2.0)
  
end subroutine


!===========================================================================
!output
!---------------------------------------------------------------------------
subroutine output
  use fluid_common
  implicit none
  character(30) :: tname

  write(tname, '(A, F12.6, A)') 't=', time/(vlx/u0), '.dat'
  open(7, file = tname)
  write(7, '(A, E13.6, A)') ' title="t=', time/(vlx/u0), '"'
  write(7,*) 'variables = x, y, u, v, p, rho'
22  format(6E15.6)
  write(7,'(A,A6,I,A6,I,A10)') 'zone T="BOX"', "I=",im, "J=",jm, 'F=point'
  do i=1, im
    do j=1, jm
	  write(7, 22) x(i,j), y(i,j), ux(i,j), uy(i,j), p(i,j), rho(i,j)
    enddo
  enddo

  print *, ' Chen: output at t=', time/(vlx/u0), 'finished ! '
  print *, '==============================================='
  close(7)
end subroutine

!===========================================================================
!final output
!---------------------------------------------------------------------------
subroutine error_analysis
  use fluid_common
  implicit none
  integer :: imid, jmid
  real(8) :: u_pre, v_pre, p_pre, factor1, factor2

  imid=nint((im+1)/2.0)
  jmid=nint((jm+1)/2.0)

  factor1=0.0
  factor2=0.0
  do i=1, im
    do j=1, jm
      u_pre=-u0*dcos(x(i,j)*pi/vlx)*dsin(y(i,j)*pi/vlx)*dexp(-2.0d0*pi**2.0d0*u0*time/Re/vlx)
      v_pre= u0*dcos(y(i,j)*pi/vlx)*dsin(x(i,j)*pi/vlx)*dexp(-2.0d0*pi**2.0d0*u0*time/Re/vlx)

      factor1=factor1+(ux(i,j)-u_pre)**2.0d0
      factor2=factor2+(uy(i,j)-v_pre)**2.0d0
    enddo
  enddo

  factor1=factor1/im/jm
  factor2=factor2/im/jm
  
  write(60,*) time/(vlx/u0), dsqrt(factor1/u0**2), dsqrt(factor2/u0**2)

end subroutine

!===========================================================================
!initial input
!---------------------------------------------------------------------------
subroutine getdata
  use fluid_common
  implicit none

  vlx=1.0
  vly=1.0
  cs2=1.0/3.0
  pi=4.0*atan(1.0)
  dt=2.0d0*vlx/dble(im-1)
  dt=dt/5.0d0
  tau=0.75d0
  niu=(tau-0.5d0)*dt*cs2
  u0=Re*niu/vlx
  tmax=vlx/u0*4.0
  ex(0:8) = (/0,1,0,-1,0,1,-1,-1,1/)
  ey(0:8) = (/0,0,1,0,-1,1,1,-1,-1/)
  w(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)

  do i=1, im
    do j=1, jm
      x(i, j)=2.0d0*dble(i-1)*vlx/dble(im-1)-vlx
      y(i, j)=2.0d0*dble(j-1)*vly/dble(jm-1)-vly
      ux(i,j)=-u0*dcos(x(i,j)*pi/vlx)*dsin(y(i,j)*pi/vlx)
      uy(i,j)= u0*dcos(y(i,j)*pi/vlx)*dsin(x(i,j)*pi/vlx)
      p(i,j)=rho0*cs2-rho0*u0**2.0/4.0*(dcos(2.0*x(i,j)*pi/vlx)+dcos(2.0*y(i,j)*pi/vlx))
      rho(i,j)=p(i,j)/cs2
      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
    enddo
  enddo

  i_output=nint(vlx/u0/dt) ! output interval for the flow field

  call output

end subroutine