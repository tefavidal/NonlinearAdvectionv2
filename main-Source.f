      program main
      
      implicit double precision (a-h,o-z)

!      parameter (Nx=900,Ny=20)
!      good for paper 2D
      parameter (Nx=3000,Ny=80)
!      good for paper 1D
!      parameter (Nx=6000,Ny=1)
!        parameter (Nx=1500,Ny=1)

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1
      double precision gamma(Nx,Ny),ro(Nx,Ny),gammax(1000,Nx),XGamma(Nx)
      double precision rox(1000,Nx)
      double precision gammanullc(274),ro1nullc(274),ro2nullc(274)
      double precision  ro1nullcNF(274)
      double precision XGa(1000),Xro(1000)
      double precision gLaplace(Nx,Ny),flux(Nx,Ny)
      double precision fluxinNx(1000), fluxinNx2(1000)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny), vdy(Nx,Ny)
      dimension gamma0(10),ro0(10)
      integer ier,pgbeg
      character(len=30) ct1,ct2

      ier=pgbeg(0,'OutputData2D/U15-Source-V1_5.ps/cps',1,1)

      open(10,file ='OutputData2D/U15-Source-V1_5'
     . ,status = 'unknown',form = 'formatted')
!      open(11,file ='OutputData2D/U15-Source-V1_5-iNx'
!     . ,status = 'unknown',form = 'formatted')
!      open(12,file ='OutputData2D/Dispersion-tes'
!     . ,status = 'unknown',form = 'formatted')
!      open(13,file ='OutputData2D/nullcline'
!     . ,status = 'unknown',form = 'formatted')
!      open(14,file ='OutputData2D/vdthreshold'
!     . ,status = 'unknown',form = 'formatted')
!      open(200,file ='OutputData2D/ColorMap'
!     . ,status = 'unknown',form = 'formatted')
!      open(500,file ='OutputData2D/gamma-time'
!     . ,status = 'unknown',form = 'formatted')
!      open(30,file ='OutputData2D/tes-LongOmega-k'
!     . ,status = 'unknown',form = 'formatted')
!      open(40,file ='OutputData2D/ThirdOrderEq'
!     . ,status = 'unknown',form = 'formatted')


      t=0.d0
      iPeriod=0

      call anfang(t,Nx,Ny,gamma,ro,gammanullc,ro1nullcNF,ro2nullc)
!      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)
!      do i=1,Nx
!       do j=1,Ny
!
!
!!      Y=ro(i,j)*gamma(i,j)/(1.d0+gamma(i,j))
!!      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)
!!       flux(i,j)=1.d0/depsilon*(s1*s2*Phi-gamma(i,j))-0.2
!       flux(i,j)=depsilon*gLaplace(i,j)-(vd*xgradeC(i,j))
!
!      enddo
!
!      enddo
      call out(t,Nx,Ny,gamma,ro,flux)
      call vmap(ier,t,Nx,Ny,gamma)

!     ^^^^^^^^^^^^^First data Row^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if (Ny .gt. 1) then
         jj=Ny/2
      else
         jj=1
      endif
         do i=1,Nx
          gammax(1,i)=gamma(i,jj)
          rox(1,i)=ro(i,jj)
         enddo
      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)
      call flow(t,Nx,Ny,vdx,vdy)
      fluxinNx2(1)=depsilon*gLaplace(250,jj)
     .              -(vdx(250,jj)*xgradeC(250,jj))
      fluxinNx(1)=depsilon*gLaplace(350,jj)
     .              -(vdx(350,jj)*xgradeC(350,jj))


!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 5    continue

      it1=t/tout
      it2=(t)/tpulse
      call ODE(t,Nx,Ny,gamma,ro)

      it=(t+0.000001)/tout
      it3=(t)/tpulse

      write(ct1,'(F6.2)') t
      write(ct2,'(F6.2)') t/dk1
      ct1='t= '//ct1
      ct2='real t= '//ct2
      write(6,*) ct1,ct2


!     EXTERNAL PULSES
!      if (it3 .gt. it2) then
!      write(6,*) '*********',it2,it3,t
!      do i=1,Nx
!       do j=1,Ny
!
!        if (i .lt. 1000 .and. i .gt. 970) then
!
!        gamma(i,j)=gamma(i,j)+2
!
!        else
!
!        gamma(i,j)=gamma(i,j)
!        endif
!
!       enddo
!      enddo
!      endif


      if (it .gt. it1) then
!     Flux calculation %%%%%%%%%%%%%%
!      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)
!      do i=1,Nx
!       do j=1,Ny
!
!
!!      Y=ro(i,j)*gamma(i,j)/(1.d0+gamma(i,j))
!!      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)
!!     flux(i,j)=1.d0/depsilon*(s1*s2*Phi-gamma(i,j))-0.2
!      flux(i,j)=depsilon*gLaplace(i,j)-
!     .            (vd*xgradeC(i,j))
!
!      enddo
!      enddo
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call out(t,Nx,Ny,gamma,ro,flux)
         romin = ro(1,1)
         romax = ro(1,1)
         do i = 1,Nx
          do j = 1,Ny
             romin = min(romin,ro(i,j))
             romax = max(romax,ro(i,j))
          enddo
         enddo
         write(6,*) 'ro=',romin,romax
         if (Ny .gt. 1) then
         jj=Ny/2
         else
         jj=1
         endif
         do i=1,Nx
          gammax(it+1,i)=gamma(i,jj)
          rox(it+1,i)=ro(i,jj)
         enddo


         call vmap(ier,t,Nx,Ny,gamma)

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     Flux Calculation
      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)
      call flow(t,Nx,Ny,vdx,vdy)
      fluxinNx2(it+1)=depsilon*gLaplace(250,jj)-
     .             (vdx(250,jj)*xgradeC(250,jj))
      fluxinNx(it+1)=depsilon*gLaplace(350,jj)-
     .             (vdx(350,jj)*xgradeC(350,jj))
      endif

!      if(it .eq. 98)then
!        open(42,file ='OutputData2D/WaveProfile1'
!     .   ,status = 'unknown',form = 'formatted')
!        dx1=dx*(dke0*Diffgamma)**0.5/dk1
!        do i=1,Nx
!            do j=1,Ny
!                write(42,*) t/dk1,i*dx1,j*dx1,gamma(i,j)
!            enddo
!        enddo
!        close(42)
!
!      endif

      if (t+dt .lt. tend) then
         goto 5

      endif
      call pgend
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      ier=pgbeg(0,'OutputData2D/U15-Source-V1_5-1D.ps/cps',1,1)
!      do k=1,it+1
!      do i=1,Nx
!      XGamma(i)=gammax(k,i)
!
!      enddo
!      call vmap2(ier,k,Nx,XGamma)
!      enddo
!      call pgend

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      ier=
!     . pgbeg(0,'OutputData2D/U15-Source-V1_5-traj-i1.ps/cps',1,1)
!      i=101
!
!      write(6,*) '***************************************************'
!!      dflux=depsilon*
!!     . (gamma(i-1,1)+gamma(i+1,1)-2*gamma(i,1))/(dx**2)
!!     .            -(vd*(gamma(i+1,1)-gamma(i-1,1))/(2*dx))
!!      dflux=0
!       dflux=flux(i,1);
!      call nullcline(i,dflux,gammanullc,ro1nullc,ro2nullc)
!
!      write(6,*) 'i=',i,'   passive flux=',dflux
!!
!      call SteadyState(i,dflux,nfix,gamma0,ro0)
!      write(6,*) 'nfix==',nfix
!
!      do j=1,nfix
!        write(6,*)  'gamma=',gamma0(j),'  ro=',ro0(j)
!      enddo
!      call stability(nfix,i,gamma0,ro0)
!      do k=1,it+1
!        do k1=1,k
!            Xro(k)=rox(k,i)
!            XGa(k)=gammax(k,i)
!        enddo
!        call vmap3(ier,k,XGa,Xro,gammanullc,ro1nullc,ro2nullc,
!     .  ro1nullcNF)
!      enddo
!      call pgend
!!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      ier=
!     .pgbeg(0,'OutputData2D/U15-Source-V1_5-traj-iNx2.ps/cps',1,1)
!      i=250
!
!!      dflux=depsilon*
!!     . (gamma(i-1,1)+gamma(i+1,1)-2*gamma(i,1))/(dx**2)
!!     .            -(vd*(gamma(i+1,1)-gamma(i-1,1))/(2*dx))
!!      dflux=-16.6
!!       dflux=flux(i,1);
!!      call nullcline(i,dflux,gammanullc,ro1nullc,ro2nullc)
!!      call SteadyState(i,dflux,nfix,gamma0,ro0)
!      write(6,*) '****************************************************'
!      write(6,*) 'i=',i,'   passive flux=',dflux
!      write(6,*) 'nfix==',nfix
!
!      do j=1,nfix
!      write(6,*)  'gamma=',gamma0(j),'  ro=',ro0(j)
!      enddo
!      call stability(nfix,i,gamma0,ro0)
!
!      do k=1,it+1
!        do k1=1,k
!            Xro(k)=rox(k,i)
!            XGa(k)=gammax(k,i)
!        enddo
!        dflux=fluxinNx2(k)
!        call nullcline(i,dflux,gammanullc,ro1nullc,ro2nullc)
!        call vmap3(ier,k,XGa,Xro,gammanullc,ro1nullc,ro2nullc)
!      enddo
!      call pgend
!!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      ier=pgbeg(0,'OutputData2D/U15-Source-V1_5-traj-iNx.ps/cps'
!     . ,1,1)
!      i=350
!!
!!      dflux=depsilon*
!!     . (gamma(i-1,1)+gamma(i+1,1)-2*gamma(i,1))/(dx**2)
!!     .            -(vd*(gamma(i+1,1)-gamma(i-1,1))/(2*dx))
!!      dflux=0
!!       dflux=flux(i,1);
!!      call nullcline(i,dflux,gammanullc,ro1nullc,ro2nullc)
!!      call SteadyState(i,dflux,nfix,gamma0,ro0)
!      write(6,*) '****************************************************'
!
!      write(6,*) 'i=',i,'   passive flux=',dflux
!      write(6,*) 'nfix==',nfix
!      do j=1,nfix
!      write(6,*)  'gamma=',gamma0(j),'  ro=',ro0(j)
!      enddo
!      call stability(nfix,i,gamma0,ro0)
!
!      do k=1,it+1
!        do k1=1,k
!            Xro(k)=rox(k,i)
!            XGa(k)=gammax(k,i)
!        enddo
!        dflux=fluxinNx(k)
!        call nullcline(i,dflux,gammanullc,ro1nullc,ro2nullc)
!        call vmap3(ier,k,XGa,Xro,gammanullc,ro1nullc,ro2nullc)
!      enddo
!      call pgend
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      SPACE-TIME PLOTS
      ier=pgbeg(0,'OutputData2D/U15-Source-V1_5-spacetime.ps/cps',1,1)

      call vmap5(1,Nx,it+1,gammax)

      call pgend

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
      open(42,file ='OutputData2D/Final-State'
     . ,status = 'unknown',form = 'formatted')
      call outFinal(t,Nx,Ny,gamma,ro)
      close(42)


      end
      

