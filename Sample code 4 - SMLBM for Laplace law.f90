!===========================================================================
!            2D SMLBM program for multiphase flow problems
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
                        im = 301, &     ! grid number in x-direction
                        jm = 301        ! grid number in y-direction
                        
  real(8),parameter ::  eta=1.0e-8      ! convergent criteria
end module

!===========================================================================
! variables 
!---------------------------------------------------------------------------
module fluid_common
  use param
  implicit none
  integer :: i, j, k, itime, i_output
  real(8) :: dt, time, tmax
  real(8) :: pi, cs2, err
  real(8) :: rho0_d, rho0_l, den_ratio  ! reference densities of the denser and the lighter fluids, density ratio
  real(8) :: tau_d, tau_l, tau_c  ! relaxation parameters for the denser fluid, the lighter fluid and the order parameter distribution function
  real(8) :: rr, dia  ! radius and diameter of the bubble
  real(8) :: sigma, kapa, beta, epsion, M !surface tension, two constants in calculating chemical potential, interface thickness, and mobility
  real(8) :: vlx, vly  ! size of the computational domain
  real(8) :: start, finish
  real(8), dimension(im, jm) :: x, y        ! grid position
  real(8), dimension(im, jm) :: ux, uxo, uy, uyo, uu  ! velocity in x, y direction, total velocity
  real(8), dimension(im, jm) :: rho, rhoo, p, po ! density, pressure
  real(8), dimension(im, jm) :: rho_s, ux_s, uy_s, p_s, tau  ! intermediate properties and local relaxation parameter
  real(8), dimension(im, jm) :: cl, clo, miu, miuo ! order parameter, chemical potential
  real(8) :: fneq(1:im, 1:jm, 0:q-1), feq(1:im, 1:jm, 0:q-1) ! non-equilibrium and equilibrium distribution functions for the distribution function
  real(8) :: gg(1:im, 1:jm, 0:q-1) !  forcing term
  real(8) :: gneq(1:im, 1:jm, 0:q-1), geq(1:im, 1:jm, 0:q-1) ! non-equilibrium and equilibrium distribution functions for the order parameter distribution function
  real(8) :: feq_s(1:im, 1:jm, 0:q-1), geq_s(1:im, 1:jm, 0:q-1)  ! equilibria computed by intermediate properties (used in the corrector step)
  real(8) :: ex(0:q-1), ey(0:q-1), w(0:q-1) ! lattice velocities and corresponding weights
end module


!===========================================================================
!main program
!---------------------------------------------------------------------------
Program SMLBMcz
  use fluid_common
  implicit none
  
  open (50, file='LBMcz.log')
  open (55, file='pressure_gap.dat')  ! output file of the pressure difference across the bubble surface

  call getdata  ! initial input

  call cpu_time(start)

  call SMLBM  ! evolution with the SMLBM

  call cpu_time(finish)

  call final_output
  
  close(50)
  close(55)
end program

!===========================================================================
!SMLBM
!---------------------------------------------------------------------------
subroutine SMLBM
  use fluid_common
  implicit none
  integer :: ip, jp, in, jn
  integer :: imid, jmid
  real(8) :: ff, f1, f2, factor, ff1, ff2, ff3, factor1, factor2
  real(8) :: fsx, fsy
  real(8) :: rx, ry, rxx, ryy, rxy, nux, nuy, nvx, nvy, den_grad_x, den_grad_y, miu_grad_x, miu_grad_y

  time=0.0
  itime=0
  err=100.0
  imid=(im+1)/2
  jmid=(jm+1)/2
  
  do while (time<tmax)
    
    ! Predictor step
    do i=1,im
      do j=1,jm
        rho_s(i,j)=0.0
        ux_s(i,j)=0.0
        uy_s(i,j)=0.0
        p_s(i,j)=0.0
        uxo(i,j)=ux(i,j)
        uyo(i,j)=uy(i,j)
        po(i,j)=p(i,j)
        rhoo(i,j)=rho(i,j)
        miuo(i,j)=miu(i,j)
        clo(i,j)=cl(i,j)
        cl(i,j)=0.0
      enddo
    enddo

    do i=2,im-1
      do j=2,jm-1
        do k=0,q-1
          ip=i-nint(ex(k))
          jp=j-nint(ey(k))
          
          factor=ex(k)*uxo(ip,jp)+ey(k)*uyo(ip,jp)
          ff=w(k)*po(ip,jp)+w(k)*rhoo(i,j)*cs2*(factor/cs2+factor**2/cs2**2/2.0-(uxo(ip,jp)**2+uyo(ip,jp)**2)/cs2/2.0)
          p_s(i,j)=p_s(i,j)+ff
          
          feq(ip,jp,k)=w(k)*po(ip,jp)+w(k)*rhoo(ip,jp)*cs2*(factor/cs2+factor**2/cs2**2/2.0-(uxo(ip,jp)**2+uyo(ip,jp)**2)/cs2/2.0)
          ux_s(i,j)=ux_s(i,j)+feq(ip,jp,k)*ex(k)
          uy_s(i,j)=uy_s(i,j)+feq(ip,jp,k)*ey(k)

          if(k==0) then
            geq(ip,jp,k)=clo(ip,jp)-(1.0-w(k))*miuo(ip,jp)*M*3.0
          else
            geq(ip,jp,k)=w(k)*3.0*(M*miuo(ip,jp)+clo(ip,jp)*factor)
          endif
          cl(i,j)=cl(i,j)+geq(ip,jp,k)
        enddo
        if(cl(i,j)>1.0) cl(i,j)=1.0  !impose cut-off value of the order parameter
        if(cl(i,j)<0.0) cl(i,j)=0.0
        
        rho(i,j)=cl(i,j)*rho0_d+(1.0-cl(i,j))*rho0_l
        ux_s(i,j)=ux_s(i,j)/cs2/rho(i,j)
        uy_s(i,j)=uy_s(i,j)/cs2/rho(i,j)
      enddo
    enddo

    ! B.C. for intermediate variables
    do i=1,im
      cl(i,1)=cl(i,2)
      rho(i,1)=rho(i,2)
      ux_s(i,1)=(4.0*ux_s(i,2)-ux_s(i,3))/3.0
      uy_s(i,1)=(4.0*uy_s(i,2)-uy_s(i,3))/3.0
      
      cl(i,jm)=cl(i,jm-1)
      rho(i,jm)=rho(i,jm-1)
      ux_s(i,jm)=(4.0*ux_s(i,jm-1)-ux_s(i,jm-2))/3.0
      uy_s(i,jm)=(4.0*uy_s(i,jm-1)-uy_s(i,jm-2))/3.0
    enddo

    do j=1,jm
      cl(1,j)=cl(2,j)
      rho(1,j)=rho(2,j)
      ux_s(1,j)=(4.0*ux_s(2,j)-ux_s(3,j))/3.0
      uy_s(1,j)=(4.0*uy_s(2,j)-uy_s(3,j))/3.0
      
      cl(im,j)=cl(im-1,j)
      rho(im,j)=rho(im-1,j)
      ux_s(im,j)=(4.0*ux_s(im-1,j)-ux_s(im-2,j))/3.0
      uy_s(im,j)=(4.0*uy_s(im-1,j)-uy_s(im-2,j))/3.0
    enddo

  ! calculate the chemical potential 
    do i=2,im-1
      do j=2,jm-1
        ip=i-1
        in=i+1
        jp=j-1
        jn=j+1
      
        miu(i,j)=2.0*beta*cl(i,j)*(cl(i,j)-1.0)*(2.0*cl(i,j)-1.0)
        miu(i,j)=miu(i,j)-kapa*(cl(in,jn)+cl(in,jp)+cl(ip,jn)+cl(ip,jp)+4.0*cl(ip,j)+4.0*cl(in,j)+4.0*cl(i,jp)+4.0*cl(i,jn)-20.0*cl(i,j))/dt**2/6.0
      enddo
    enddo
    
    ! B.C. for the chemical potential
    do i=1,im
      miu(i,1)=miu(i,2)
      miu(i,jm)=miu(i,jm-1)
    enddo

    do j=1,jm
      miu(1,j)=miu(2,j)
      miu(im,j)=miu(im-1,j)
    enddo
    
    ! non-equilibrium terms & forcing terms
    call forcing_term
    
    do i=2,im-1
      do j=2,jm-1
        do k=0,q-1
          factor=ex(k)*ux_s(i,j)+ey(k)*uy_s(i,j)
          feq_s(i,j,k)=w(k)*p_s(i,j)+w(k)*rho(i,j)*cs2*(factor/cs2+factor**2/cs2**2/2.0-(ux_s(i,j)**2+uy_s(i,j)**2)/cs2/2.0)
          
          ip=i-nint(ex(k))
          jp=j-nint(ey(k))
        
          fneq(i,j,k)=feq_s(i,j,k)-feq(ip,jp,k)
          
          if(k==0) then
            geq_s(i,j,k)=cl(i,j)-(1.0-w(k))*miu(i,j)*M*3.0
           else
            geq_s(i,j,k)=w(k)*3.0*(M*miu(i,j)+cl(i,j)*factor)
          endif
          gneq(i,j,k)=-tau_c*(geq_s(i,j,k)-geq(ip,jp,k))
        enddo
      enddo
    enddo
    
    ! B.C. for non-equilibrium terms & forcing terms
    do i=1,im
      do k=0,q-1
        fneq(i,1,k)=fneq(i,2,k)
        gneq(i,1,k)=gneq(i,2,k)
        
        fneq(i,jm,k)=fneq(i,jm-1,k)
        gneq(i,jm,k)=gneq(i,jm-1,k)
      enddo
    enddo

    do j=1,jm
      do k=0,q-1
        fneq(1,j,k)=fneq(2,j,k)
        gneq(1,j,k)=gneq(2,j,k)
      
        fneq(im,j,k)=fneq(im-1,j,k)
        gneq(im,j,k)=gneq(im-1,j,k)
      enddo
    enddo
    
    ! Corrector step 
    do i=2,im-1
      do j=2,jm-1
        ux_s(i,j)=ux_s(i,j)*cs2*rho(i,j)
        uy_s(i,j)=uy_s(i,j)*cs2*rho(i,j)
        do k=0,q-1
          ip=i-nint(ex(k))
          jp=j-nint(ey(k))
          in=i+nint(ex(k))
          jn=j+nint(ey(k))
          if(ip<1) ip=im-1
          if(in<1) ip=im-1
          if(jp<1) jp=jm-1
          if(jn<1) jp=jm-1
          
          ux_s(i,j)=ux_s(i,j)-0.5*ex(k)*(gg(in,jn,k)-gg(ip,jp,k))*dt*tau_d+ex(k)*(tau_d-1.0)*(fneq(in,jn,k)-fneq(i,j,k))
          uy_s(i,j)=uy_s(i,j)-0.5*ey(k)*(gg(in,jn,k)-gg(ip,jp,k))*dt*tau_d+ey(k)*(tau_d-1.0)*(fneq(in,jn,k)-fneq(i,j,k))
          cl(i,j)=cl(i,j)-(1.0-1.0/tau_c)*(gneq(in,jn,k)-gneq(i,j,k))
        enddo
        
        if(cl(i,j)>1.0) cl(i,j)=1.0
        if(cl(i,j)<0.0) cl(i,j)=0.0
        rho(i,j)=cl(i,j)*rho0_d+(1.0-cl(i,j))*rho0_l
        ux_s(i,j)=ux_s(i,j)/cs2
        uy_s(i,j)=uy_s(i,j)/cs2
      enddo
    enddo

    ! B.C. for corrected order parameter and density
    do i=1,im
      cl(i,1)=cl(i,2)
      rho(i,1)=rho(i,2)
      cl(i,jm)=cl(i,jm-1)
      rho(i,jm)=rho(i,jm-1)
    enddo

    do j=1,jm
      cl(1,j)=cl(2,j)
      rho(1,j)=rho(2,j)
      cl(im,j)=cl(im-1,j)
      rho(im,j)=rho(im-1,j)
    enddo

  ! calculate the chemical potential 
    do i=2,im-1
      do j=2,jm-1
        ip=i-1
        in=i+1
        jp=j-1
        jn=j+1
      
        miu(i,j)=2.0*beta*cl(i,j)*(cl(i,j)-1.0)*(2.0*cl(i,j)-1.0)
        miu(i,j)=miu(i,j)-kapa*(cl(in,jn)+cl(in,jp)+cl(ip,jn)+cl(ip,jp)+4.0*cl(ip,j)+4.0*cl(in,j)+4.0*cl(i,jp)+4.0*cl(i,jn)-20.0*cl(i,j))/dt**2/6.0
      enddo
    enddo
    
    ! B.C. for the chemical potential
    do i=1,im
      miu(i,1)=miu(i,2)
      miu(i,jm)=miu(i,jm-1)
    enddo

    do j=1,jm
      miu(1,j)=miu(2,j)
      miu(im,j)=miu(im-1,j)
    enddo
    
    ! implement surface tension
    do i=2,im-1
      do j=2,jm-1
        ip=i-1
        in=i+1
        jp=j-1
        jn=j+1
        
        ! surface tension
        miu_grad_x=(miu(in,j)-miu(ip,j))/3.0/dt+(miu(in,jp)-miu(ip,jp))/12.0/dt+(miu(in,jn)-miu(ip,jn))/12.0/dt
        miu_grad_y=(miu(i,jn)-miu(i,jp))/3.0/dt+(miu(ip,jn)-miu(ip,jp))/12.0/dt+(miu(in,jn)-miu(in,jp))/12.0/dt
        fsx=-cl(i,j)*miu_grad_x
        fsy=-cl(i,j)*miu_grad_y
        
        ux(i,j)=(ux_s(i,j)+fsx*dt)/rho(i,j)
        uy(i,j)=(uy_s(i,j)+fsy*dt)/rho(i,j)

        p(i,j)=p_s(i,j)
      enddo
    enddo
    
    ! B.C. for the flow variables
    ! left & right walls
    do j=1, jm
      p(1,j)=(4.0*p(2,j)-p(3,j))/3.0
      rho(1,j)=rho(2,j)
      ux(1,j)=(4.0*ux(2,j)-ux(3,j))/3.0
      uy(1,j)=(4.0*uy(2,j)-uy(3,j))/3.0
    
      p(im,j)=(4.0*p(im-1,j)-p(im-2,j))/3.0
      rho(im,j)=rho(im-1,j)
      ux(im,j)=(4.0*ux(im-1,j)-ux(im-2,j))/3.0
      uy(im,j)=(4.0*uy(im-1,j)-uy(im-2,j))/3.0
    enddo

    ! top & bottom walls
    do i=1,im
      p(i,1)=(4.0*p(i,2)-p(i,3))/3.0
      rho(i,1)=rho(i,2)
      ux(i,1)=(4.0*ux(i,2)-ux(i,3))/3.0
      uy(i,1)=(4.0*uy(i,2)-uy(i,3))/3.0
    
      p(i,jm)=(4.0*p(i,jm-1)-p(i,jm-2))/3.0
      rho(i,jm)=rho(i,jm-1)
      ux(i,jm)=(4.0*ux(i,jm-1)-ux(i,jm-2))/3.0
      uy(i,jm)=(4.0*uy(i,jm-1)-uy(i,jm-2))/3.0
    enddo
    
    time=time+dt
    itime=itime+1

    if(mod(itime,i_output)==0) then
      call output
    endif

    !output pressure gap
    write(55,*) itime, p(imid,jmid)-p(10,jmid)
    
  enddo
  
end subroutine


!===========================================================================
!calculating equilibrium term
!---------------------------------------------------------------------------
subroutine equilibrium(ww,ex,ey,pp,ro,u1,u2,kk,ff)
  !use fluid_common
  implicit none
  real(8) :: factor, ff, pp, ro, u1, u2, ww, ex, ey
  integer :: kk

  factor=ex*u1+ey*u2
  ff=ww*pp+ww*ro/3.0d0*(factor*3.0+factor**2*9.0/2.0-(u1**2+u2**2)*3.0/2.0)
  
end subroutine


subroutine forcing_term
  use fluid_common
  implicit none

  real(8) :: factor, f1, f2, ff1, ff2, u1, u2
  real(8) :: rx, ry, miux, miuy
  integer :: ip, jp
  
  do i=2,im-1
    do j=2,jm-1
      ip=i-1
      jp=j-1
      
      rx=(rho(i+1,j)-rho(ip,j))/3.0/dt+(rho(i+1,j+1)-rho(ip,j+1))/12.0/dt+(rho(i+1,jp)-rho(ip,jp))/12.0/dt
      ry=(rho(i,j+1)-rho(i,jp))/3.0/dt+(rho(i+1,j+1)-rho(i+1,jp))/12.0/dt+(rho(ip,j+1)-rho(ip,jp))/12.0/dt
      miux=(miu(i+1,j)-miu(ip,j))/3.0/dt+(miu(i+1,j+1)-miu(ip,j+1))/12.0/dt+(miu(i+1,jp)-miu(ip,jp))/12.0/dt
      miuy=(miu(i,j+1)-miu(i,jp))/3.0/dt+(miu(i+1,j+1)-miu(i+1,jp))/12.0/dt+(miu(ip,j+1)-miu(ip,jp))/12.0/dt

  
      do k=0,q-1
        factor=ex(k)*ux_s(i,j)+ey(k)*uy_s(i,j)
        ff1=w(k)*(1.0+factor/cs2+factor**2/cs2**2/2.0-(ux_s(i,j)**2+uy_s(i,j)**2)/cs2/2.0)
        ff2=ff1-w(k)
        
        f1=rx*cs2*ff2-cl(i,j)*miux*ff1
        f2=ry*cs2*ff2-cl(i,j)*miuy*ff1
      
        gg(i,j,k)=(1.0-0.5/tau_d)*((ex(k)-ux_s(i,j))*f1+(ey(k)-uy_s(i,j))*f2)
      enddo
    enddo
  enddo
  
  do i=1,im
    do k=0,q-1
      gg(i,1,k)=gg(i,2,k)
      gg(i,jm,k)=gg(i,jm-1,k)
    enddo
  enddo

  do j=1,jm
    do k=0,q-1
      gg(1,j,k)=gg(2,j,k)
      gg(im,j,k)=gg(im-1,j,k)
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
  integer :: imid, jmid
  
  imid=(im+1)/2
  jmid=(jm+1)/2
  
  write(*,*) 'Reference: sigma/r=', sigma/rr
  write(*,*) 'Computed results: ', p(imid,jmid)-p(10,jmid)
  
  write(50,*) 'time=',time,'pressure difference: ', p(imid,jmid)-p(10,jmid)
  
  print *, ' Chen: output at t=', time, 'finished ! '
  print *, '==============================================='

end subroutine

!===========================================================================
!final output
!---------------------------------------------------------------------------
subroutine final_output
  use fluid_common
  implicit none
  character(30) :: tname
  integer :: imid, jmid
  
  imid=(im+1)/2
  jmid=(jm+1)/2
  
  open(70,file='pressure.dat')
  do i=1,im
    write(70,*) x(i,jmid), p(i,jmid)
  enddo
  close(70)
    
  write(tname, '(A, F10.3, A)') 't=', time, '.dat'
  open(7, file = tname)
  write(7, '(A, E13.6,A)') ' title="t=', time, ' s"'
  write(7,*) 'variables = x, y, u, v, p, rho, fraction, miu'
22  format(9E15.6)
  write(7,'(A,A6,I,A6,I,A10)') 'zone T="BOX"', "I=",im, "J=",jm, 'F=point'
  do j=1, jm
    do i=1, im
      write(7, 22) x(i,j), y(i,j), ux(i,j), uy(i,j), p(i,j), rho(i,j), cl(i,j), miu(i,j)
    enddo
  enddo
  close(7)
  
  write(50,*) 'Reference: sigma/r=', sigma/rr
  write(50,*) 'Computed results: ', p(imid,jmid)-p(10,jmid)
end subroutine

!===========================================================================
!initial input
!---------------------------------------------------------------------------
subroutine getdata
  use fluid_common
  implicit none
  real(8) :: factor, factor1, factor2, dd, h, ff
  integer :: in, ip, jp, jn

  vlx=dble(im-1)
  vly=dble(jm-1)

  rho0_d=5.0  ! density of denser fluid
  den_ratio=100.0  ! density ratio
  rho0_l=rho0_d/den_ratio
  rr=40.0  ! radius of the drop
  dia=rr*2.0
  cs2=1.0/3.0
  pi=4.0*atan(1.0)
  dt=vlx/(im-1)

  tau_d= 0.8
  tau_l= 0.8
  tau_c=1.2
  
  ex(0:8) = (/0,1,0,-1,0,1,-1,-1,1/)
  ey(0:8) = (/0,0,1,0,-1,1,1,-1,-1/)
  w(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)

  sigma=0.005  ! surface tension
  epsion=5.0*dt  ! interface thickness
  kapa=1.5*sigma*epsion
  beta=12.0*sigma/epsion
  M=0.5  ! mobility
  tmax=50000.0   !maximum iteration
  i_output=100  ! output interval

  do i=1, im
    do j=1, jm
      x(i, j)=(i-1)*vlx/(im-1)-vlx/2.0
      y(i, j)=(j-1)*vly/(jm-1)-vly/2.0
      
      cl(i,j)=0.0
      factor=sqrt(x(i,j)**2+y(i,j)**2)
      if(factor < (rr+3.0*dt)) cl(i,j)=1.0-0.5*(1.0+tanh(2.0*(factor-rr)/epsion))
      if(factor < (rr-3.0*dt)) cl(i,j)=1.0

      ux(i,j)=0.0
      uy(i,j)=0.0

      rho(i,j)=cl(i,j)*rho0_d+(1.0-cl(i,j))*rho0_l
      tau(i,j)=tau_d
      p(i,j)=1.0
      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
    enddo
  enddo

  ! calculate initial chemical potential
  do i=2,im-1
    do j=2,jm-1
      ip=i-1
      in=i+1
      jp=j-1
      jn=j+1
      
      miu(i,j)=2.0*beta*cl(i,j)*(cl(i,j)-1.0)*(2.0*cl(i,j)-1.0)
      miu(i,j)=miu(i,j)-kapa*(cl(in,jn)+cl(in,jp)+cl(ip,jn)+cl(ip,jp)+4.0*cl(ip,j)+4.0*cl(in,j)+4.0*cl(i,jp)+4.0*cl(i,jn)-20.0*cl(i,j))/dt**2/6.0
    enddo
  enddo
    
  ! B.C. for the chemical potential
  do i=1,im
    miu(i,1)=miu(i,2)
    miu(i,jm)=miu(i,jm-1)
  enddo

  do j=1,jm
    miu(1,j)=miu(2,j)
    miu(im,j)=miu(im-1,j)
  enddo

end subroutine

