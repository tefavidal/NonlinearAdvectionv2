      subroutine out(t,Nx,Ny,gamma,ro,vdx)
      
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny),vdx(Nx,Ny)
      double precision gLaplace(Nx,Ny), xgradeC(Nx,Ny), ygradeC(Nx,Ny)
      double precision aux
      double precision t



      dls=dk1/(dke0*Diffgamma)**0.5

       do i=1,Nx
      do j=1,Ny


       write(10,*) t/dk1,i*dx/dls, j*dy/dls,gamma(i,j),ro(i,j),vdx(i,j)
        enddo
      enddo

      write(10,*)
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

      implicit double precision (a-h, o-z)
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

