      program main
      
      implicit none


      INTEGER, PARAMETER :: Nx=2000
      INTEGER, PARAMETER :: Ny=160

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision gammanullc(300),ro1nullc(300),ro2nullc(300)
      double precision  ro1nullcNF(274)
      double precision XGa(1000),Xro(1000)
      double precision gLaplace(Nx,Ny),flux(Nx,Ny)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny), vdy(Nx,Ny)
      double precision gamma0(10),ro0(10)
      double precision t, iPeriod, it, it1
      double precision romin, romax
      integer ier,pgbeg, i, j
      integer counter
      character(len=4) ct1
      character(len=39) ct2
      character(len=43) ct3





      t=0.d0
      counter=0;
      iPeriod=0

      call anfang(t,Nx,Ny,gamma,ro,gammanullc,ro1nullcNF,ro2nullc)
      call flow(t,Nx,Ny,vdx)


          open(10,file ='/data.lfpn/evidal/2CompFlow/salida/data   0'
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,gamma,ro,vdx)
          close(10)
      ct2='/data.lfpn/evidal/2CompFlow/salida/data'


 5    continue

      call ODE(t,Nx,Ny,gamma,ro,vdx)

      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
      write(ct1,'(I4)') counter
      ct3 = ct2 // ct1

      write(6,*) ct3

!      if (mod(counter, 10) .eq. 0)then
!        write(6,*) 'Firing'
!
!         do j=75,84
!!           do j=5,14
!            if(Nx .eq. 2000)then
!               do i=105,112
!                  ro(i,j)=ro0(1)
!                  gamma(i,j)=gamma0(1)+2
!               enddo
!            elseif(Nx .eq. 1000)then
!               do i=103,106
!                  ro(i,j)=ro0(1)
!                  gamma(i,j)=gamma0(1)+2
!               enddo
!            elseif(Nx .eq. 500)then
!               do i=202,203
!                  ro(i,j)=ro0(1)
!                  gamma(i,j)=gamma0(1)+2
!               enddo
!            endif
!         enddo
!      endif


          open(10,file =ct3
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,gamma,ro,vdx)
          close(10)




      if (t+dt .lt. tend) then
         goto 5

      endif
      
      close(10)

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
      open(42,file ='OutputData2D/Final-State'
     . ,status = 'unknown',form = 'formatted')
      call outFinal(t,Nx,Ny,gamma,ro)
      close(42)


      end
      

