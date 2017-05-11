      subroutine out(t,Nx,Ny,gamma,ro,vdx)
      
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny),vdx(Nx,Ny)
      double precision gLaplace(Nx,Ny), xgradeC(Nx,Ny), ygradeC(Nx,Ny)
      double precision aux




      it = nint(t/tout)

      if (Ny .gt. 1) then
      j=Ny/2
      else
      j=1
      endif

      dls=dk1/(dke0*Diffgamma)**0.5

      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)


      do i=1,Nx
        meangamma=0.d0
        do jj=1,Ny
            meangamma=meangamma+gamma(i,jj)
        enddo
        meangamma=meangamma/Ny
      write(10,*) t/dk1,i*dx/dls,gamma(i,j),ro(i,j),vdx(i,j),meangamma
      enddo
      write(10,*)

!      write(6,*)'salida'
      do jj=1,Ny
        aux=depsilon*gLaplace(101,jj)-(vdx(101,jj)*xgradeC(101,jj))
        write(51,*) t/dk1,i*dx/dls,gamma(101,jj),ro(101,jj), aux

        aux=depsilon*gLaplace(105,jj)-(vdx(105,jj)*xgradeC(105,jj))
         write(55,*) t/dk1,i*dx/dls,gamma(105,jj),ro(105,jj), aux

        aux=depsilon*gLaplace(110,jj)-(vdx(110,jj)*xgradeC(110,jj))
        write(60,*) t/dk1,i*dx/dls,gamma(110,jj),ro(110,jj), aux
      enddo
      write(51,*)
      write(55,*)
      write(60,*)

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

