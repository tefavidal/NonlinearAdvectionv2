      program main
      
      implicit double precision (a-h,o-z)

!      parameter (Nx=900,Ny=20)
!      good for paper 2D
      parameter (Nx=1000,Ny=160)
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
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny), vdy(Nx,Ny)
      dimension gamma0(10),ro0(10)
      integer ier,pgbeg
      character(len=30) ct1,ct2

      ier=pgbeg(0,'OutputData2D/U15-Source-V1_5.ps/cps',1,1)

      open(10,file ='OutputData2D/centre'
     . ,status = 'unknown',form = 'formatted')
      open(51,file ='OutputData2D/i101'
     . ,status = 'unknown',form = 'formatted')
      open(55,file ='OutputData2D/i105'
     . ,status = 'unknown',form = 'formatted')
      open(60,file ='OutputData2D/i110'
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
      call flow(t,Nx,Ny,vdx)
      call out(t,Nx,Ny,gamma,ro,vdx)
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


!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 5    continue

      it1=t/tout
      it2=(t)/tpulse
      call ODE(t,Nx,Ny,gamma,ro,vdx)

      it=(t+0.000001)/tout
      it3=(t)/tpulse

      write(ct1,'(F6.2)') t
      write(ct2,'(F6.2)') t/dk1
      ct1='t= '//ct1
      ct2='real t= '//ct2
      write(6,*) ct1,ct2




      if (it .gt. it1) then

      call out(t,Nx,Ny,gamma,ro,vdx)
         romin = ro(1,1)
         romax = ro(1,1)
         do j = 1,Ny
          do i = 1,Nx
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

      endif


      if (t+dt .lt. tend) then
         goto 5

      endif
      call pgend

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
      

