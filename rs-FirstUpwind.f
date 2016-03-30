      subroutine rs(t,Nx,Ny,gamma,ro,gammaprime,roprime)
      
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision gammaprime(Nx,Ny),roprime(Nx,Ny)
      double precision f1(Nx,Ny),f2(Nx,Ny)
      double precision Phi(Nx,Ny),gLaplace(Nx,Ny)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny),vdy(Nx,Ny)
      common /param/ gamma01,ro01


      call function1(t,Nx,Ny,gamma,ro,f1,f2,Phi)
      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)
      call flow(t,Nx,Ny,vdx,vdy)


      do i=1,Nx
       do j=1,Ny

       if (i .gt. 100) then
        roprime(i,j)=-f1(i,j)*ro(i,j)+f2(i,j)*(1.d0-ro(i,j))
        if (t .ge. tE) then

        gammaprime(i,j)=1.d0/depsilon*(s1*s2*Phi(i,j)-gamma(i,j))
     .                  +depsilon*gLaplace(i,j)-
     .            (vdx(i,j)*xgradeC(i,j)+vdy(i,j)*ygradeC(i,j))

        else
        gammaprime(i,j)=1.d0/depsilon*(s1*s2*Phi(i,j)-gamma(i,j))
     .                  +depsilon*gLaplace(i,j)
!
        endif
!       roprime(i,j)=-f10*(ro(i,j)-ro01)-B*ro01*(gamma(i,j)-gamma01)
!     . -f20*(ro(i,j)-ro01)+A*(1.d0-ro01)*(gamma(i,j)-gamma01)
!
!
!        gammaprime(i,j)=1/depsilon*(s1*s2*(dM*(ro(i,j)-ro01)
!     .   +dN*(gamma(i,j)-gamma01))-(gamma(i,j)-gamma01))
       else
       roprime(i,j)=0.d0
       if (t .ge. tE) then
       gammaprime(i,j)=depsilon*gLaplace(i,j)
     .            -(vdx(i,j)*xgradeC(i,j)+vdy(i,j)*ygradeC(i,j))
       else
       gammaprime(i,j)=depsilon*gLaplace(i,j)
!       gammaprime(i,j)=0.d0
       endif

       endif


       enddo
      enddo

      return

      end
!     ***********************************************************
      subroutine function1(t,Nx,Ny,gamma,ro,f1,f2,Phi)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision f1(Nx,Ny),f2(Nx,Ny),Phi(Nx,Ny)
      common /param/ gamma01,ro01

      do i=1,Nx
       do j=1,Ny
        f1(i,j)=(1.d0+dk*gamma(i,j))/(1.d0+gamma(i,j))
        f2(i,j)=(dL1+dk*dL2*dc*gamma(i,j))/(1.d0+dc*gamma(i,j))
        Y=ro(i,j)*gamma(i,j)/(1.d0+gamma(i,j))
        Phi(i,j)=(dlambda1+Y**2)/(dlambda2+Y**2)
       enddo
      enddo

      return
      end

!      ***********************************************************

      subroutine functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      double precision gamma(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      common /param/ gamma01,ro01


      do i=1,Nx
       do j=1,Ny
!       No-Flux boundary condition
       if(i .eq. 1) then
        gammaim1=gamma(i,j)
!       gammaim1=0.d0
        gammaip1=gamma(i+1,j)
       elseif(i .eq. Nx) then
        gammaip1=gamma(i,j)
!        gammaip1=gamma01
        gammaim1=gamma(i-1,j)
       else
        gammaip1=gamma(i+1,j)
        gammaim1=gamma(i-1,j)
       endif

       if(j .eq. 1) then
        gammajm1=gamma(i,j)
!        gammajm1=0.d0
        gammajp1=gamma(i,j+1)
       elseif(j .eq. Ny) then
        gammajp1=gamma(i,j)
!        gammajp1=0.d0

        gammajm1=gamma(i,j-1)
       else
        gammajp1=gamma(i,j+1)
        gammajm1=gamma(i,j-1)
       endif

        gLapX=(gammaip1+gammaim1-2*gamma(i,j))/(dx**2)
        gLapY=(gammajp1+gammajm1-2*gamma(i,j))/(dy**2)
        gLaplace(i,j)=gLapX+gLapY



!        xgradeC(i,j)=(gammaip1-gammaim1)/(2*dx)
        xgradeC(i,j)=(gamma(i,j)-gammaim1)/(dx)
        ygradeC(i,j)=(gammajp1-gammajm1)/(2*dy)
!        xgradeC(i,j)=(gammaip1-gamma(i,j))/(dx)
!        ygradeC(i,j)=(gammajp1-gamma(i,j))/(dy)
       enddo
      enddo

      return
      end
!     ******************************************************************
      subroutine flow(t,Nx,Ny,vdx,vdy)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      double precision vdx(Nx,Ny),vdy(Nx,Ny)
      common /param/ gamma01,ro01



!     %%%%%%%%%%%%%%%%parabolic flow
!      do i=1,Nx
!       do j=1,Ny
!!       k=t/0.45
!       k=1
!      vdx(i,j)=vd*
!     . (1.d0-(j-1-(Ny-1.d0)/2)**2/((Ny-1.d0)/2)**2)
!!        vdx(i,j)=(-1)**(k+1)*vd
!!         if(j .gt. Ny/2) then
!!         vdx(i,j)=vd
!!         else
!!         vdx(i,j)=vd-0.1
!!         endif
!        vdy(i,j)=0
!       enddo
!      enddo
!     %%%%%%%%%%%%%%%%Realistic Flow
      if(Ny .eq. 80)then
        do i=1,Nx
            vdx(i,1)=0
            vdx(i,2)=vd*0.7946
            vdx(i,3)=vd*0.9720
            vdx(i,4)=vd*0.9976
            vdx(i,5)=vd*0.9999

            vdx(i,76)=vd*0.9999
            vdx(i,77)=vd*0.9976
            vdx(i,78)=vd*0.9720
            vdx(i,79)=vd*0.7946
            vdx(i,80)=0

            vdy(i,1)=0
            vdy(i,2)=0
            vdy(i,3)=0
            vdy(i,4)=0
            vdy(i,5)=0

            vdy(i,76)=0
            vdy(i,77)=0
            vdy(i,78)=0
            vdy(i,79)=0
            vdy(i,80)=0

            do j=6,75
                vdx(i,j)=vd
                vdy(i,j)=0
            enddo
        enddo

      elseif(Ny .eq. 100)then
        do i=1,Nx
            vdx(i,1)=0
            vdx(i,2)=vd*0.7004
            vdx(i,3)=vd*0.9317
            vdx(i,4)=vd*0.9895
            vdx(i,5)=vd*0.9987
            vdx(i,6)=vd*0.9999

            vdx(i,95)=vd*0.9999
            vdx(i,96)=vd*0.9987
            vdx(i,97)=vd*0.9895
            vdx(i,98)=vd*0.9317
            vdx(i,99)=vd*0.7004
            vdx(i,100)=0

            vdy(i,1)=0
            vdy(i,2)=0
            vdy(i,3)=0
            vdy(i,4)=0
            vdy(i,5)=0
            vdy(i,6)=0

            vdy(i,95)=0
            vdy(i,96)=0
            vdy(i,97)=0
            vdy(i,98)=0
            vdy(i,99)=0
            vdy(i,100)=0

            do j=7,94
                vdx(i,j)=vd
                vdy(i,j)=0
            enddo
        enddo


      elseif(Ny .eq. 120)then
        do i=1,Nx
            vdx(i,1)=0
            vdx(i,2)=vd*0.6350
            vdx(i,3)=vd*0.8903
            vdx(i,4)=vd*0.9708
            vdx(i,5)=vd*0.9943
            vdx(i,6)=vd*0.9991
            vdx(i,7)=vd*0.9999

            vdx(i,114)=vd*0.9999
            vdx(i,115)=vd*0.9991
            vdx(i,116)=vd*0.9943
            vdx(i,117)=vd*0.9708
            vdx(i,118)=vd*0.8903
            vdx(i,119)=vd*0.6350
            vdx(i,120)=0

            vdy(i,1)=0
            vdy(i,2)=0
            vdy(i,3)=0
            vdy(i,4)=0
            vdy(i,5)=0
            vdy(i,6)=0
            vdy(i,7)=0

            vdy(i,114)=0
            vdy(i,115)=0
            vdy(i,116)=0
            vdy(i,117)=0
            vdy(i,118)=0
            vdy(i,119)=0
            vdy(i,120)=0

            do j=8,113
                vdx(i,j)=vd
                vdy(i,j)=0
            enddo
        enddo



      elseif(Ny .eq. 160)then
        do i=1,Nx
            vdx(i,1)=0
            vdx(i,2)=vd*0.5144
            vdx(i,3)=vd*0.7946
            vdx(i,4)=vd*0.9165
            vdx(i,5)=vd*0.9721
            vdx(i,6)=vd*0.9910
            vdx(i,7)=vd*0.9977
            vdx(i,8)=vd*0.9995
            vdx(i,9)=vd*0.9999

            vdx(i,152)=vd*0.9999
            vdx(i,153)=vd*0.9995
            vdx(i,154)=vd*0.9977
            vdx(i,155)=vd*0.9910
            vdx(i,156)=vd*0.9721
            vdx(i,157)=vd*0.9165
            vdx(i,158)=vd*0.7946
            vdx(i,159)=vd*0.5144
            vdx(i,160)=0

            vdy(i,1)=0
            vdy(i,2)=0
            vdy(i,3)=0
            vdy(i,4)=0
            vdy(i,5)=0
            vdy(i,6)=0
            vdy(i,7)=0
            vdy(i,8)=0
            vdy(i,9)=0

            vdy(i,152)=0
            vdy(i,153)=0
            vdy(i,154)=0
            vdy(i,155)=0
            vdy(i,156)=0
            vdy(i,157)=0
            vdy(i,158)=0
            vdy(i,159)=0
            vdy(i,160)=0

            do j=10,151
                vdx(i,j)=vd
                vdy(i,j)=0
            enddo
        enddo



      else
        do i=1,Nx
            do j=1,Ny
                vdx(i,j)=vd
                vdy(i,j)=0
            enddo
        enddo

      endif

      return
      end

