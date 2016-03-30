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
        endif

       else
        roprime(i,j)=0.d0
        if (t .ge. tE) then
            gammaprime(i,j)=depsilon*gLaplace(i,j)
     .            -(vdx(i,j)*xgradeC(i,j)+vdy(i,j)*ygradeC(i,j))
        else
            gammaprime(i,j)=depsilon*gLaplace(i,j)
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
        gammaim2=gamma(i+1,j)
        gammaim1=gamma(i+1,j)
        gammaip1=gamma(i+1,j)
       elseif(i .eq. 2) then
        gammaim2=-gamma(i,j)+2*gamma(i-1,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)
       elseif(i .eq. Nx) then
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i-1,j)
       else
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

       endif

       if(Ny .eq. 1) then
        gammajm1=gamma(i,j)
        gammajp1=gamma(i,j)
       elseif(j .eq. 1) then
        gammajm1=gamma(i,j+1)
        gammajp1=gamma(i,j+1)
       elseif(j .eq. Ny) then
        gammajp1=gamma(i,j-1)
        gammajm1=gamma(i,j-1)
       else
        gammajp1=gamma(i,j+1)
        gammajm1=gamma(i,j-1)
       endif

        gLapX=(gammaip1+gammaim1-2*gamma(i,j))/(dx**2)
        gLapY=(gammajp1+gammajm1-2*gamma(i,j))/(dy**2)
        gLaplace(i,j)=gLapX+gLapY



        if(gammaip1 .eq. gamma(i,j)) then
        thetai=1.d-10
        else
        thetai=(gamma(i,j)-gammaim1)/(gammaip1-gamma(i,j))+1.d-10
        endif

        if(gamma(i,j) .eq. gammaim1)then
        thetaim1=1.d-10
        else
        thetaim1=(gammaim1-gammaim2)/(gamma(i,j)-gammaim1)+1.d-10
        endif

        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))

       xgradeC(i,j)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i,j))/(dx)


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
!
!      vdx(i,j)=vd*
!     . (1.d0-(j-1-(Ny-1.d0)/2)**2/((Ny-1.d0)/2)**2)
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


      elseif(Ny .eq. 200)then
        do i=1,Nx
            vdx(i,1)=0
            vdx(i,2)=vd*0.4347
            vdx(i,3)=vd*0.7004
            vdx(i,4)=vd*0.8517
            vdx(i,5)=vd*0.9317
            vdx(i,6)=vd*0.9708
            vdx(i,7)=vd*0.9885
            vdx(i,8)=vd*0.9958
            vdx(i,9)=vd*0.9986
            vdx(i,10)=vd*0.9996
            vdx(i,11)=vd*0.9999

            vdx(i,190)=vd*0.9999
            vdx(i,191)=vd*0.9996
            vdx(i,192)=vd*0.9986
            vdx(i,193)=vd*0.9958
            vdx(i,194)=vd*0.9885
            vdx(i,195)=vd*0.9708
            vdx(i,196)=vd*0.9317
            vdx(i,197)=vd*0.8517
            vdx(i,198)=vd*0.7004
            vdx(i,199)=vd*0.4347
            vdx(i,200)=0

            vdy(i,1)=0
            vdy(i,2)=0
            vdy(i,3)=0
            vdy(i,4)=0
            vdy(i,5)=0
            vdy(i,6)=0
            vdy(i,7)=0
            vdy(i,8)=0
            vdy(i,9)=0
            vdy(i,10)=0
            vdy(i,11)=0

            vdy(i,190)=0
            vdy(i,191)=0
            vdy(i,192)=0
            vdy(i,193)=0
            vdy(i,194)=0
            vdy(i,195)=0
            vdy(i,196)=0
            vdy(i,197)=0
            vdy(i,198)=0
            vdy(i,199)=0
            vdy(i,200)=0

            do j=12,189
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


