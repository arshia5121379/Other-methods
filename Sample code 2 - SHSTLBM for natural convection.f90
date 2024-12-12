!===========================================================================
!            2D SHSTLBM program for thermal flow problems
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
                        im = 251, &     ! grid number in x-direction
                        jm = 251        ! grid number in y-direction
  real(8),parameter ::  rho0 = 1.0,  &  ! reference density
                        Pr=0.71,      &  ! Prandtl number
                        Ra=1.0e6,    &  ! Rayleigh number
                        g=-9.81,     &  ! Gravitational acceleration
                        R=1.0,       &  ! Gas constant
                        eta=1.0e-8      ! convergent criteria
end module

!===========================================================================
! variables
!---------------------------------------------------------------------------
module fluid_common
  use param
  implicit none
  integer :: i, j, k, itime, i_output
  real(8) :: dt, time
  real(8) :: pi, tau, tau_c, cs2, vc, beta, err, err2, ka, niu
  real(8) :: tem_max, tem_mean, e_max
  real(8) :: vlx, vly
  real(8) :: start, finish
  real(8), dimension(im, jm) :: x, y        ! grid position
  real(8), dimension(im, jm) :: ux, uy, uu  ! velocity in x, y direction, total velocity
  real(8), dimension(im, jm) :: rho, p, tem, e ! density, pressure, temperature,  internal energy
  real(8), dimension(im, jm) :: rho_s, ux_s, uy_s, e_s, str
  real(8) :: fneq(1:im, 1:jm, 0:q-1), gneq(1:im, 1:jm, 0:q-1)
  real(8) :: ex(0:q-1), ey(0:q-1), w(0:q-1) ! lattice velocities and corresponding weights
end module


!===========================================================================
!main program
!---------------------------------------------------------------------------
Program SHSTLBMcz
  use fluid_common
  implicit none
  
  open (50, file='LBMcz.log')

  call getdata  ! initial input

  call cpu_time(start)

  call SHSTLBM

  call cpu_time(finish)

  call final_output
  
  close (50)
end program

!===========================================================================
!SHSTLBM evolution
!---------------------------------------------------------------------------
subroutine SHSTLBM
  use fluid_common
  implicit none
  integer :: ip, jp
  real(8) :: ff, gg, f1, f2, g1, g2

  time=0.0
  itime=0
  err=100.0
  do while (err > eta .or. err2 > eta)
    
    ! step 1: rho_star, u_star, e_star
    rho_s=0.0
    ux_s=0.0
    uy_s=0.0
    e_s=0.0
    do i=2,im-1
      do j=2,jm-1
        do k=0,q-1
          ip=i-nint(ex(k))
          jp=j-nint(ey(k))
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),e(ip,jp),k,ff,gg)
          rho_s(i,j)=rho_s(i,j)+ff
          ux_s(i,j)=ux_s(i,j)+ff*ex(k)
          uy_s(i,j)=uy_s(i,j)+ff*ey(k)
          e_s(i,j)=e_s(i,j)+gg
        enddo
        ux_s(i,j)=ux_s(i,j)/rho_s(i,j)
        uy_s(i,j)=uy_s(i,j)/rho_s(i,j)
        e_s(i,j)=e_s(i,j)/rho_s(i,j)

        do k=0,q-1
          call equilibrium(rho_s(i,j),ux_s(i,j),uy_s(i,j),e_s(i,j),k,f1,g1)
          
          ip=i-nint(ex(k))
          jp=j-nint(ey(k))
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),e(ip,jp),k,f2,g2)

          fneq(i,j,k)=-tau*(f1-f2)
          gneq(i,j,k)=-tau_c*(g1-g2)
        enddo

      enddo
    enddo

!    B.C. for non-equilibrium terms
    do j=2,jm-1
      do k=0,q-1
        fneq(1,j,k)=fneq(2,j,k)
        fneq(im,j,k)=fneq(im-1,j,k)
        gneq(1,j,k)=gneq(2,j,k)
        gneq(im,j,k)=gneq(im-1,j,k)
      enddo
    enddo

    ! top & bottom walls
    do i=1,im
      do k=0,q-1
        fneq(i,1,k)=fneq(i,2,k)
        fneq(i,jm,k)=fneq(i,jm-1,k)
        gneq(i,1,k)=gneq(i,2,k)
        gneq(i,jm,k)=gneq(i,jm-1,k)
      enddo
    enddo

    ! step 2: rho_n+1, u_n+1, e_n+1
    
    ! left & right walls
    do i=2,im-1
      do j=2,jm-1
        ux(i,j)=ux_s(i,j)*rho_s(i,j)
        uy(i,j)=uy_s(i,j)*rho_s(i,j)
        e(i,j)=e_s(i,j)*rho_s(i,j)
        do k=0,q-1
          ip=i+nint(ex(k))
          jp=j+nint(ey(k))
          ux(i,j)=ux(i,j)-(1-1.0/tau)*ex(k)*fneq(ip,jp,k)
          uy(i,j)=uy(i,j)-(1-1.0/tau)*ey(k)*fneq(ip,jp,k)
          e(i,j)=e(i,j)-(1-1.0/tau_c)*gneq(ip,jp,k)
        enddo
        ux(i,j)=ux(i,j)/rho_s(i,j)
        uy(i,j)=(uy(i,j)-rho(i,j)*g*beta*(tem(i,j)-tem_mean)*dt)/rho_s(i,j)
        e(i,j)=e(i,j)/rho_s(i,j)
        rho(i,j)=rho_s(i,j)
      enddo
    enddo

    ! step 3: boundary condition
    call boundary
    
    ! step 4: convergence test
    call convergence_test

    time=time+dt
    itime=itime+1

    if(mod(itime,i_output)==0) then
      print *, 'Relative errors: ', err, err2
      print *, ' Chen: output at t=', time, 'finished ! '
      print *, '==============================================='
      
      write(50,*) 'Relative errors: ', err, err2
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
  do j=1,jm
    rho(1,j)=(4.0*rho(2,j)-rho(3,j))/3.0
    ux(1,j)=0.0
    uy(1,j)=0.0
    e(1,j)=e_max

    rho(im,j)=(4.0*rho(im-1,j)-rho(im-2,j))/3.0
    ux(im,j)=0.0
    uy(im,j)=0.0
    e(im,j)=0.0
  enddo

  ! top & bottom walls
  do i=2,im-1
    rho(i,1)=rho(i,2)-dt/cs2*(-rho(i,2)*g*beta*(e(i,2)/R-tem_mean))
    ux(i,1)=0.0
    uy(i,1)=0.0
    e(i,1)=(4.0*e(i,2)-e(i,3))/3.0
    
    rho(i,jm)=rho(i,jm-1)+dt/cs2*(-rho(i,jm-1)*g*beta*(e(i,jm-1)/R-tem_mean))
    ux(i,jm)=0.0
    uy(i,jm)=0.0
    e(i,jm)=(4.0*e(i,jm-1)-e(i,jm-2))/3.0
  enddo

end subroutine

!===========================================================================
!calculating equilibrium term
!---------------------------------------------------------------------------
subroutine equilibrium(ro,u1,u2,eo,kk,ff,gg)
  use fluid_common
  implicit none
  real(8) :: factor, ff, gg, ro, u1, u2, eo
  integer :: kk

  factor=ex(kk)*u1+ey(kk)*u2
  ff=w(kk)*ro*(1.0+factor/cs2+(factor/cs2)**2/2.0-(u1**2+u2**2)/cs2/2.0)

  if(kk==0) then
    gg=-2.0/3.0*ro*eo*(u1**2+u2**2)
  elseif(kk<=4) then
    gg=ro*eo/9.0*(1.5+1.5*factor+4.5*factor**2-1.5*(u1**2+u2**2))
  else
    gg=ro*eo/36.0*(3.0+6.0*factor+4.5*factor**2-1.5*(u1**2+u2**2))
  endif
  
end subroutine

!===========================================================================
!test of convergent criteria
!---------------------------------------------------------------------------
subroutine convergence_test
  use fluid_common
  implicit none
  real(8) :: du, u_sum, dtem, tem_sum, temp

  err=0.0
  err2=0.0
  u_sum=0.0
  tem_sum=0.0

  do i=1,im
    do j=1,jm
      p(i,j)=cs2*rho(i,j)
      temp=e(i,j)/R

      dtem=abs(temp-tem(i,j))
      du=abs(uu(i,j)-sqrt(ux(i,j)**2+uy(i,j)**2))

      u_sum=u_sum+sqrt(ux(i,j)**2+uy(i,j)**2)
      tem_sum=tem_sum+abs(temp)

      err=err+du
      err2=err2+dtem

      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
      tem(i,j)=temp
    enddo
  enddo
  
  err=err/u_sum
  err2=err2/tem_sum
end subroutine

!===========================================================================
!final output
!---------------------------------------------------------------------------
subroutine final_output
  use fluid_common
  implicit none
  integer :: imid, jmid
  real(8) :: nu, factor, um, vm, xx, yy, nu0, nu_max, nu_min, yy_min, yy_max
  character(30) :: tname

  imid=nint((im+1)/2.0)
  jmid=nint((jm+1)/2.0)

  nu=0.0
  nu0=0.0
  nu_max=0.0
  nu_min=10.0
  um=0.0
  vm=0.0
  do i=1, im-1
    do j=1, jm-1
      factor=(ux(i,j)+ux(i+1,j)+ux(i,j+1)+ux(i+1,j+1))/4.0*(tem(i,j)+tem(i+1,j)+tem(i,j+1)+tem(i+1,j+1))/4.0
      factor=factor-ka*(tem(i+1,j)-tem(i,j)+tem(i+1,j+1)-tem(i,j+1))/2.0/dt
      nu=nu+factor*dt**2/ka/tem_max/vlx**2
    enddo
  enddo

  do j=1,jm
    if(ux(imid,j)>um) then
      um=ux(imid,j)
      yy=y(imid,j)
    endif
    factor=abs(4.0*tem(2,j)-3.0*tem(1,j)-tem(3,j))/2.0/dt
    
    if(factor>nu_max) then
      nu_max=factor
      yy_max=y(1,j)
    endif
    if(factor<nu_min) then
      nu_min=factor
      yy_min=y(1,j)
    endif
  enddo

  do j=1, jm-1
    factor=abs(4.0*tem(2,j)-3.0*tem(1,j)-tem(3,j))/2.0/dt+abs(4.0*tem(2,j+1)-3.0*tem(1,j+1)-tem(3,j+1))/2.0/dt
    nu0=nu0+factor/2.0*dt
  enddo

  do i=1,im
    if(uy(i,jmid)>vm) then
      vm=uy(i,jmid)
      xx=x(i,jmid)
    endif
  enddo

  write(50,*) 'Iteration steps:', itime
  write(50,*) 'Converged time t=', time, 's'
  write(50,*) 'CPU time:', finish-start, 's'
  write(50,*) 'Averaged Nu', nu
  write(50,*) 'u_max=', um/ka, 'y=', yy
  write(50,*) 'v_max=', vm/ka, 'x=', xx
  write(50,*) 'Nu_max', nu_max, 'y=', yy_max
  write(50,*) 'Nu_min', nu_min, 'y=', yy_min
  write(50,*) 'Nu0', nu0

  do j=1,jm
    str(1,j)=0.0
  enddo

  do i=2,im
    do j=1,jm
      str(i,j)=str(i-1,j)-dt*(uy(i,j)+uy(i-1,j))/2.0
    enddo
  enddo

  write(tname, '(A, F12.6, A)') 't=', time, 's.dat'
  open(7, file = tname)
  write(7, '(A, E13.6,A)') ' title="t=', time, ' s"'
  write(7,*) 'variables = x, y, u, v, p, rho, tem, e, str'
22  format(9E15.6)
  write(7,'(A,A6,I,A6,I,A10)') 'zone T="BOX"', "I=",im, "J=",jm, 'F=point'
  do i=1, im
    do j=1, jm
	  write(7, 22) x(i,j), y(i,j), ux(i,j), uy(i,j), p(i,j), rho(i,j), tem(i,j), e(i,j), str(i,j)
    enddo
  enddo

  print *, ' Chen: output at t=', time, 'finished ! '
  print *, '==============================================='

end subroutine

!===========================================================================
!initial input
!---------------------------------------------------------------------------
subroutine getdata
  use fluid_common
  implicit none
  real(8) :: factor

  vlx=1.0
  vly=1.0
  tem_max=1.0
  tem_mean=0.5
  e_max=R*tem_max
  cs2=1.0/3.0
  vc=0.1
  beta=vc**2/vlx/tem_max/abs(g)
  niu=vc*vlx*sqrt(Pr/Ra)  ! kinematic viscosity
  ka=vc*vlx/sqrt(Pr*Ra)  ! thermal diffusivity
  pi=4.0*atan(1.0)
  dt=vlx/(im-1)
  i_output=nint(50.0/dt)
  tau= 0.5+niu/cs2/dt
  tau_c=0.5+3.0*ka/2.0/dt
  ex(0:8) = (/0,1,0,-1,0,1,-1,-1,1/)
  ey(0:8) = (/0,0,1,0,-1,1,1,-1,-1/)
  w(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)

  do i=1, im
    do j=1, jm
      x(i, j)=(i-1)*vlx/(im-1)
      y(i, j)=(j-1)*vly/(jm-1)
      ux(i,j)=0.0
      uy(i,j)=0.0
      if(i==1) then
        tem(i,j)=tem_max
        e(i,j)=e_max
      else
        tem(i,j)=0.0
        e(i,j)=0.0
      endif
      rho(i,j)=rho0
      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
    enddo
  enddo

end subroutine