      subroutine out(t,Nx,Ny,gamma,ro,flux)
      
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny),flux(Nx,Ny)




      it = nint(t/tout)

      if (Ny .gt. 1) then
      j=Ny/2
      else
      j=1
      endif

      dls=dk1/(dke0*Diffgamma)**0.5

      do i=1,Nx
        meangamma=0.d0
        do jj=1,Ny
        meangamma=meangamma+gamma(i,jj)
        enddo
        meangamma=meangamma/Ny
      write(10,*) t/dk1,i*dx/dls,gamma(i,j),ro(i,j),flux(i,j),meangamma

      enddo
      write(10,*)
      ii=1200
      flux1=depsilon*
     . (gamma(ii-1,j)+gamma(ii+1,j)-2*gamma(ii,j))/(dx**2)
     .            -(vd*(gamma(ii+1,j)-gamma(ii-1,j))/(2*dx))
!      flux1=-2.0

      f1=(1.d0+dk*gamma(ii,j))/(1.d0+gamma(ii,j))
      f2=(dL1+dk*dL2*dc*gamma(ii,j))/(1.d0+dc*gamma(ii,j))
      Y=ro(ii,j)*gamma(ii,j)/(1.d0+gamma(ii,j))
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      flux2=(s1*s2*Phi-gamma(ii,j))/depsilon

!      write(11,*) t/dk1,gamma(ii,j),gamma(ii+400,j),gamma(ii+800,j),
!     . gamma(ii+1200,j),gamma(ii+1600,j),flux1,flux2,flux1+flux2
!     . ,ro(ii,j)

      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine outFinal(t,Nx,Ny,gamma,ro)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny)



      do i=1,Nx
        do j=1,Ny
            write(42,*) t/dk1,gamma(i,j),ro(i,j)
        enddo
      enddo
      write(42,*)
      return
      end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine loadState(t,Nx,Ny,gamma,ro)

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      common /param/ gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny)


      do i=1,Nx
        do j=1,Ny
            read(7,*) t52,gamma(i,j),ro(i,j)
        enddo
      enddo
      close(7)
      return
      end

