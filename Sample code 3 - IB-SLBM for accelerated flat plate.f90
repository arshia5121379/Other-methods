!===========================================================================
!            2D IB-SLBM program for hydrodynamic problems
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
                        im = 701, &     ! grid number in x-direction
                        jm = 181, &     ! grid number in y-direction
                        lm = 41        ! number of lagrangian points
  real(8),parameter ::  rho0 = 1.0,  &  ! reference density
                        alpha=31900.0,    &  ! non-dimensional parameter
                        eta=1.0e-8      ! convergent criteria
end module

!===========================================================================
! variables
!---------------------------------------------------------------------------
module fluid_common
  use param
  implicit none
  integer :: i, j, k, ii, jj, itime, i_output
  integer :: i_min, i_max, j_min, j_max
  real(8) :: dt, time, tmax, dx1, dx2, dl, dth
  real(8) :: pi, tau, cs2, err, niu, a_plate
  real(8) :: u0, x_c, y_c, del1, del2
  real(8) :: h_p, dx_refine, dy_refine
  real(8) :: start, finish
  real(8), dimension(im, jm) :: x, y        ! grid position
  real(8), dimension(im, jm) :: ux, uy, uu  ! velocity in x, y direction, total velocity
  real(8), dimension(im, jm) :: rho, vor, p ! density, vorticity, pressure
  real(8), dimension(im, jm) :: rho_s, ux_s, uy_s
  real(8) :: fneq(1:im, 1:jm, 0:q-1)
  real(8) :: ex(0:q-1), ey(0:q-1), w(0:q-1) ! lattice velocities and corresponding weights
  real(8) :: gg(-1:1,-1:1,1:q-1)
  real(8) :: ff(0:q-1)
  real(8) :: aa(lm, lm), bu(lm), bv(lm), dul(lm), dvl(lm)     ! IBM matrix & velocity correction
  real(8) :: xl(lm), yl(lm), ylo(lm), uxl(lm), uyl(lm), pl(lm), rhol(lm)     ! Phisical properties of Lagrangian points
  real(8) :: ix_rho(lm),iy_rho(lm) 

end module


!===========================================================================
!main program
!---------------------------------------------------------------------------
Program IB_SLBMcz
  use fluid_common
  implicit none
  
  open (30, file='drag.dat')

  call getdata  ! initial input

  call cpu_time(start)

  call IB_SLBM

  call cpu_time(finish)
  
end program

!===========================================================================
!LBM with BGK
!---------------------------------------------------------------------------
subroutine IB_SLBM
  use fluid_common
  implicit none
  integer :: ip, jp
  real(8) :: f1, f2, a1, a2, a3, b1, b2, b3, xx, yy
  real(8) :: yp, yn
  real(8) :: drag, factor
  real(8), dimension(1:im, -1:1, 1:3) :: a
  real(8), dimension(1:jm, -1:1, 1:3) :: b

  time=0.0
  itime=0
  err=100.0
  
  ! calculate Lagrangian interpolation parameters
  do i=2, im-1
    do k=-1, 1
      xx=x(i,2)-dt*k
      a(i,k,1)=(xx-x(i,2))*(xx-x(i+1,2))/(x(i-1,2)-x(i,2))/(x(i-1,2)-x(i+1,2))
      a(i,k,2)=(xx-x(i-1,2))*(xx-x(i+1,2))/(x(i,2)-x(i-1,2))/(x(i,2)-x(i+1,2))
      a(i,k,3)=(xx-x(i-1,2))*(xx-x(i,2))/(x(i+1,2)-x(i-1,2))/(x(i+1,2)-x(i,2))
    enddo
  enddo

  do j=2, jm-1
    do k=-1, 1
      yy=y(2,j)-dt*k
      b(j,k,1)=(yy-y(2,j))*(yy-y(2,j+1))/(y(2,j-1)-y(2,j))/(y(2,j-1)-y(2,j+1))
      b(j,k,2)=(yy-y(2,j-1))*(yy-y(2,j+1))/(y(2,j)-y(2,j-1))/(y(2,j)-y(2,j+1))
      b(j,k,3)=(yy-y(2,j-1))*(yy-y(2,j))/(y(2,j+1)-y(2,j-1))/(y(2,j+1)-y(2,j))
    enddo
  enddo 

  do while (time<tmax)
    ! step 1: rho_star, u_star
    rho_s=0.0
    ux_s=0.0
    uy_s=0.0
    do i=2,im-1
      do j=2,jm-1
        
        do ii=-1,1
          do jj=-1,1
            do k=5,q-1
              ip=i+ii
              jp=j+jj
              call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),k,gg(ii,jj,k))
            enddo
          enddo
        enddo

        do ii=-1,1
          ip=i+ii
          jp=j
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),1,gg(ii,0,1))
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),3,gg(ii,0,3))
        enddo

        do jj=-1,1
          ip=i
          jp=j+jj
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),2,gg(0,jj,2))
          call equilibrium(rho(ip,jp),ux(ip,jp),uy(ip,jp),4,gg(0,jj,4))
        enddo
        
        !interpolate equilibrium distribution functions at streaming nodes
        call equilibrium(rho(i,j),ux(i,j),uy(i,j),0,ff(0))
        ff(1)=a(i,ex(1),1)*gg(-1,0,1)+a(i,ex(1),2)*gg(0,0,1)+a(i,ex(1),3)*gg(1,0,1)
        ff(2)=b(j,ey(2),1)*gg(0,-1,2)+b(j,ey(2),2)*gg(0,0,2)+b(j,ey(2),3)*gg(0,1,2)
        ff(3)=a(i,ex(3),1)*gg(-1,0,3)+a(i,ex(3),2)*gg(0,0,3)+a(i,ex(3),3)*gg(1,0,3)
        ff(4)=b(j,ey(4),1)*gg(0,-1,4)+b(j,ey(4),2)*gg(0,0,4)+b(j,ey(4),3)*gg(0,1,4)
        do k=5,q-1
          ff(k) = b(j,ey(k),1)*(a(i,ex(k),1)*gg(-1,-1,k)+a(i,ex(k),2)*gg(0,-1,k)+a(i,ex(k),3)*gg(1,-1,k))+&
                 &b(j,ey(k),2)*(a(i,ex(k),1)*gg(-1,0,k)+a(i,ex(k),2)*gg(0,0,k)+a(i,ex(k),3)*gg(1,0,k))+&
                 &b(j,ey(k),3)*(a(i,ex(k),1)*gg(-1,1,k)+a(i,ex(k),2)*gg(0,1,k)+a(i,ex(k),3)*gg(1,1,k))
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

        do k=0,q-1
          call equilibrium(rho_s(i,j),ux_s(i,j),uy_s(i,j),k,f1)
          f2=ff(k)
          fneq(i,j,k)=-tau*(f1-f2)
        enddo

      enddo
    enddo
    
    ! left & right walls
    do j=1,jm
      do k=0,q-1
        fneq(1,j,k)=fneq(2,j,k)
        fneq(im,j,k)=fneq(im-1,j,k)
      enddo
    enddo
    
    ! top & bottom walls
    do i=1,im
      do k=0,q-1
        fneq(i,1,k)=fneq(i,2,k)
        fneq(i,jm,k)=fneq(i,jm-1,k)
      enddo
    enddo

    ! step 2: rho_n+1, u_n+1
    rho=rho_s
    do i=2,im-1
      do j=2,jm-1
        
        do ii=-1,1
          do jj=-1,1
            do k=5,q-1           
              ip=i+ii
              jp=j+jj
              gg(ii,jj,k)=fneq(ip,jp,k)
            enddo
          enddo
        enddo
        
        do ii=-1,1
          ip=i+ii
          jp=j
          gg(ii,0,1)=fneq(ip,jp,1)
          gg(ii,0,3)=fneq(ip,jp,3)
        enddo

        do jj=-1,1
          ip=i
          jp=j+jj
          gg(0,jj,2)=fneq(ip,jp,2)
          gg(0,jj,4)=fneq(ip,jp,4)
        enddo

        ! interpolate non-equilibrium distribution functions at streaming nodes
        ff(0)=fneq(i,j,0)
        ff(1)=a(i,-ex(1),1)*gg(-1,0,1)+a(i,-ex(1),2)*gg(0,0,1)+a(i,-ex(1),3)*gg(1,0,1)
        ff(2)=b(j,-ey(2),1)*gg(0,-1,2)+b(j,-ey(2),2)*gg(0,0,2)+b(j,-ey(2),3)*gg(0,1,2)
        ff(3)=a(i,-ex(3),1)*gg(-1,0,3)+a(i,-ex(3),2)*gg(0,0,3)+a(i,-ex(3),3)*gg(1,0,3)
        ff(4)=b(j,-ey(4),1)*gg(0,-1,4)+b(j,-ey(4),2)*gg(0,0,4)+b(j,-ey(4),3)*gg(0,1,4)
        do k=5,q-1
          ff(k) = b(j,-ey(k),1)*(a(i,-ex(k),1)*gg(-1,-1,k)+a(i,-ex(k),2)*gg(0,-1,k)+a(i,-ex(k),3)*gg(1,-1,k))+&
                 &b(j,-ey(k),2)*(a(i,-ex(k),1)*gg(-1,0,k)+a(i,-ex(k),2)*gg(0,0,k)+a(i,-ex(k),3)*gg(1,0,k))+&
                 &b(j,-ey(k),3)*(a(i,-ex(k),1)*gg(-1,1,k)+a(i,-ex(k),2)*gg(0,1,k)+a(i,-ex(k),3)*gg(1,1,k))
        enddo

        ! update momentum
        ux(i,j)=ux_s(i,j)*rho_s(i,j)
        uy(i,j)=uy_s(i,j)*rho_s(i,j)
        do k=0,q-1
          ux(i,j)=ux(i,j)-(1-1.0/tau)*ex(k)*ff(k)
          uy(i,j)=uy(i,j)-(1-1.0/tau)*ey(k)*ff(k)
        enddo
        
        ! calculate updated velocity
        ux(i,j)=ux(i,j)/rho(i,j)
        uy(i,j)=uy(i,j)/rho(i,j)
      enddo
    enddo

    ! step 3: boundary condition
    call boundary

    ! apply EoS to update pressure
    do i=1,im
      do j=1,jm
        p(i,j)=cs2*rho(i,j)
      enddo
    enddo

    ! calculate vorticity field
    do i=2, im-1
      do j=2,jm-1
        vor(i,j)=(uy(i+1,j)-uy(i-1,j))/(x(i+1,j)-x(i-1,j))-(ux(i,j+1)-ux(i,j-1))/(y(i,j+1)-y(i,j-1))
      enddo
    enddo

    time=time+dt
    itime=itime+1
    
    ! output flow field
    if(mod(itime,i_output)==0) then
      call output
    endif
    
    ! output the time history of the drag coefficient of the plate
    if(itime>1 .and. mod(itime,20)==0) then
      drag=0.0
      do i=1, lm
        drag=drag-rho0*dul(i)/dt
      enddo
      write(30,35) a_plate*time**2/h_p, drag/(0.5*rho0*(a_plate*time)**2*h_p)
    endif
35  format(2E15.6)

  enddo
  
end subroutine

!===========================================================================
!boundary treatment
!---------------------------------------------------------------------------
subroutine boundary
  use fluid_common
  implicit none

  ! Cylinder-Boundary condition-enforced IBM
  call IBM_solver

  ! left & right walls
  do j=1,jm
    rho(1,j)=(4.0*rho(2,j)-rho(3,j))/3.0
    ux(1,j)=(4.0*ux(2,j)-ux(3,j))/3.0
    uy(1,j)=(4.0*uy(2,j)-uy(3,j))/3.0

    rho(im,j)=(4.0*rho(im-1,j)-rho(im-2,j))/3.0
    ux(im,j)=(4.0*ux(im-1,j)-ux(im-2,j))/3.0
    uy(im,j)=(4.0*uy(im-1,j)-uy(im-2,j))/3.0
  enddo
  
  ! top & bottom walls
  do i=1,im
    rho(i,1)=(4.0*rho(i,2)-rho(i,3))/3.0
    ux(i,1)=(4.0*ux(i,2)-ux(i,3))/3.0
    uy(i,1)=(4.0*uy(i,2)-uy(i,3))/3.0

    rho(i,jm)=(4.0*rho(i,jm-1)-rho(i,jm-2))/3.0
    ux(i,jm)=(4.0*ux(i,jm-1)-ux(i,jm-2))/3.0
    uy(i,jm)=(4.0*uy(i,jm-1)-uy(i,jm-2))/3.0
  enddo
end subroutine

!===========================================================================
!Boundary condition-enforced IBM
!---------------------------------------------------------------------------
subroutine IBM_solver
  use fluid_common
  implicit none
  real(8) :: factor, rmin, rr1, rr2
  integer, dimension(lm) :: nl
  integer, dimension(lm, 25) :: il, jl    
  integer :: i_start, i_stop
  
  ! identify neighboring Eulerian mesh points around the plate
  x_c=-0.5*a_plate*time**2

  i_start=i_min+nint((dx_refine-h_p*0.6+x_c)/dx1)
  i_stop =i_min+nint((dx_refine-h_p*0.4+x_c)/dx1)

  do k=1, lm
    xl(k)=x_c
    yl(k)=yl(k)
    uxl(k)=-a_plate*time
    nl(k)=0
    rmin=100.0
    do i=i_start, i_stop
      do j=j_min,j_max
        
        if(abs(xl(k)-x(i,j))<=dx1*2.0 .and. abs(yl(k)-y(i,j))<=dx1*2.0) then
            nl(k)=nl(k)+1
            il(k, nl(k))=i
            jl(k, nl(k))=j
        endif

        rr1=sqrt((xl(k)-x(i,j))**2+(yl(k)-y(i,j))**2)
        if(rr1<rmin) then
          rmin=rr1
          ix_rho(k)=i
          iy_rho(k)=j
        endif

      enddo
    enddo
  enddo

  ! Construct matrix A
  aa=0.0
  do i=1,lm
    do j=1,lm
      do ii=1, nl(i)
        call kernel(il(i,ii),jl(i,ii),i,del1)
        call kernel(il(i,ii),jl(i,ii),j,del2)
        aa(i,j)=aa(i,j)+del1*del2*dx1**2
      enddo
      aa(i,j)=aa(i,j)
    enddo
  enddo
  
  call inverse(aa, lm)
  
  ! Construct vector B
  do k=1, lm
    bu(k)=uxl(k)  
    bv(k)=uyl(k)
    do i=1, nl(k)
      ii=il(k,i)
      jj=jl(k,i)
      call kernel(ii,jj,k,factor)
      bu(k)=bu(k)-factor*ux(ii,jj)*dx1**2
      bv(k)=bv(k)-factor*uy(ii,jj)*dx1**2
    enddo
  enddo

  ! Solve velocity correction at Lagrangian mesh points
  do i=1,lm
    dul(i)=0.0
    dvl(i)=0.0
    do j=1,lm
      dul(i)=dul(i)+aa(i,j)*bu(j)
      dvl(i)=dvl(i)+aa(i,j)*bv(j)
    enddo
  enddo

  ! cast velocity corrections to Eulerian mesh points
  do k=1,lm
    do i=1,nl(k)
      ii=il(k,i)
      jj=jl(k,i)
      call kernel(ii,jj,k,factor)
      ux(ii,jj)=ux(ii,jj)+dul(k)*factor
      uy(ii,jj)=uy(ii,jj)+dvl(k)*factor
    enddo
  enddo

  do i=1,lm
    rhol(i)=rho(ix_rho(i),iy_rho(i))
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
!kernel function
!---------------------------------------------------------------------------
subroutine kernel(iii,jjj,kkk,del)
  use fluid_common
  implicit none
  real(8) :: factor, del, dd1, dd2
  integer :: iii, jjj, kkk

  dd1=abs(x(iii,jjj)-xl(kkk))/dx1
  dd2=abs(y(iii,jjj)-yl(kkk))/dx1
  
  if(dd1>2.0 .or. dd2>2.0) then
    del=0.0
  else
    del=0.25*(1.0+cos(pi*dd1/2.0))*0.25*(1.0+cos(pi*dd2/2.0))/dx1**2
  endif

end subroutine


!===========================================================================
!output of flow field
!---------------------------------------------------------------------------
subroutine output
  use fluid_common
  implicit none
  character(30) :: tname

  write(tname, '(A, F12.6, A)') 't=', time, 's.dat'
  open(7, file = tname)
  write(7, '(A, E13.6,A)') ' title="t=', time, ' s"'
  write(7,*) 'variables = x, y, u, v, p, rho, vor'
22  format(8E15.6)
  write(7,'(A,A6,I,A6,I,A10)') 'zone T="BOX"', "I=",im, "J=",jm, 'F=point'
  do j=1, jm
    do i=1, im
	  write(7, 22) x(i,j), y(i,j), ux(i,j), uy(i,j), p(i,j), rho(i,j), vor(i,j)
    enddo
  enddo

  write(7,'(A,I,A10)') 'zone I=',lm, 'F=point'
  do i=1, lm
	write(7, 22) xl(i), yl(i), uxl(i), uyl(i), pl(i), rhol(i), 0.0
  enddo

  print *, ' Chen: output at t=', time, 'finished ! '
  print *, '==============================================='

  close(7)
end subroutine

!===========================================================================
!initial input
!---------------------------------------------------------------------------
subroutine getdata
  use fluid_common
  implicit none
  integer :: ix_refine, iy_refine, ix_start, iy_start
  real(8) :: factor
  real(8) :: ratio_x, ratio_y
  real(8) :: xx, yy, dx0, dy0, rr1, rr2, rmin

  h_p=0.5  ! height of the plate
  a_plate=1.0e-3  ! acceleration of the plate
  dx_refine=11.0*h_p  ! length of refined zone in the x-direction
  dy_refine=2.0*h_p  ! height of the refined zone in the y-direction
  iy_refine=dy_refine/0.01 ! number of refined mesh points along the y-direction
  ix_refine=nint(iy_refine*dx_refine/dy_refine)  ! number of refined mesh points along the x-direction
  ix_start=(im-ix_refine-1)/2+1  ! starting index of the refined mesh point in the x-direction
  iy_start=(jm-iy_refine-1)/2+1  ! starting index of the refined mesh point in the y-direction
  i_min=ix_start
  i_max=ix_start+ix_refine
  j_min=iy_start
  j_max=iy_start+iy_refine
  ratio_x=1.052  ! mesh stretching ratio in the x-direction
  ratio_y=1.125  ! mesh stretching ratio in the y-direction
  
  u0=0.1  ! initial plate velocity
  cs2=1.0/3.0  ! square of sound speed
  pi=4.0*atan(1.0)
  ex(0:8) = (/0,1,0,-1,0,1,-1,-1,1/)
  ey(0:8) = (/0,0,1,0,-1,1,1,-1,-1/)
  w(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
  dx1=dx_refine/ix_refine  ! mesh spacing in the refined zone
  dt=dx1  ! time step is equal to the minimal mesh spacing
  tmax=101.0  ! maximum time of computation
  i_output=nint(10.0/dt)  ! output interval
  niu=dsqrt(h_p**3*a_plate/alpha)  ! kinematic viscosity
  tau= 0.5+niu/cs2/dt  ! relaxation parameter
  j_min=j_min+0.4*h_p/dx1
  j_max=j_max-0.4*h_p/dx1

  !mesh generation
  xx=0.0
  yy=0.0

  do i=1,ix_start-1
    dx0=dx1*ratio_x**i
    xx=xx+dx0
  enddo
  x(1,:)=-xx-dx_refine+h_p*0.5

  do i=2,ix_start
    dx0=dx1*ratio_x**(ix_start-i+1)
    x(i,:)=x(i-1,1)+dx0
  enddo

  do i=ix_start+1, ix_start+ix_refine
    x(i,:)=x(i-1,1)+dx1
  enddo

  do i=ix_start+ix_refine, im
    dx0=dx1*ratio_x**(i-ix_start-ix_refine)
    x(i,:)=x(i-1,1)+dx0
  enddo 

  do j=1,iy_start-1
    dx0=dx1*ratio_y**j
    yy=yy+dx0
  enddo
  y(:,1)=-yy-0.5*dy_refine

  do j=2,iy_start
    dx0=dx1*ratio_y**(iy_start-j+1)
    y(:,j)=y(1,j-1)+dx0
  enddo

  do j=iy_start+1, iy_start+iy_refine
    y(:,j)=y(1,j-1)+dx1
  enddo

  do j=iy_start+iy_refine, jm
    dx0=dx1*ratio_y**(j-iy_start-iy_refine)
    y(:,j)=y(1,j-1)+dx0
  enddo 

  ! initialize flow field
  do i=1, im
    do j=1, jm
      ux(i,j)=0.0
      uy(i,j)=0.0
      rho(i,j)=rho0
      uu(i,j)=sqrt(ux(i,j)**2+uy(i,j)**2)
    enddo
  enddo
!
  dl=h_p/(lm-1)
  do k=1, lm
    xl(k)=0.0
    yl(k)=dl*(k-1)-0.5*h_p
    rhol(k)=rho0
    pl(k)=cs2*rhol(k)
    uxl(k)=-u0
    uyl(k)=0.0
  enddo

  call output

end subroutine

!===========================================================================
!solve the inverse matrix
!---------------------------------------------------------------------------
subroutine inverse(aa,n)
  integer n,i,j,k
  real(8):: aa(n,n),b(n,n),a(n,n)
  a=aa
  do i=1,n
    b(i,i)=1.0
  enddo

  do i=1,n
    b(i,:)=b(i,:)/a(i,i)
    a(i,i:n)=a(i,i:n)/a(i,i) 
    do j=i+1,n     
      do k=1,n   
        b(j,k)=b(j,k)-b(i,k)*a(j,i)
      enddo
      a(j,i:n)=a(j,i:n)-a(i,i:n)*a(j,i)
    enddo
  enddo

  do i=n,1,-1
    do j=i-1,1,-1
      do k=1,n
        b(j,k)=b(j,k)-b(i,k)*a(j,i)
      enddo
    enddo
  enddo

  aa=b
end