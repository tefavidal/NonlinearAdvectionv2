      subroutine anfang(t,Nx,Ny,gamma,ro,gammanullc,ro1nullc,ro2nullc)
      
      implicit none
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1
      double precision gamma01,ro01,Diffgamma,dke0,dk1

      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision t, depsilono, dh, dk2, dki, dKR, dkt, dlambda
      double precision dq, dsigma0, dtheta, dtprime, tendprime
      double precision tEprime, toutprime, tpulseprime, Vmax
      double precision gammanullc(300),ro1nullc(300),ro2nullc(300)
      integer Nx, Ny



      open(8,file = 'startdata',status = 'unknown', form = 'formatted')

   
      read(8,*) dc,dh,tEprime
      read(8,*) dk1,dk2,dKR
      read(8,*) dki,dke0,dkt
      read(8,*) dL1,dL2,Diffgamma
      read(8,*) dq,depsilono,dlambda
      read(8,*) dsigma0,dtheta,dalpha
      read(8,*) Vmax
      read(8,*) tendprime,toutprime,dtprime
      read(8,*) tol,tpulseprime


      close(8)

      dk=dk2/dk1
      dlambda1=dlambda*dtheta/depsilono
      dlambda2=(1+dalpha*dtheta)/(depsilono*(1+dalpha))
      s1=dq*dsigma0/(dkt+dki)*dalpha/(dalpha+1)
      s2=dkt/(dke0*dh)
      depsilon=dk1/dke0
      depsilonp=dk1/(dki+dkt)
      vd=Vmax/(dke0*Diffgamma)**0.5
      dt=dtprime*dk1
      tE=tEprime*dk1
      tend=tendprime*dk1
      tout=toutprime*dk1
      tpulse=tpulseprime*dk1


      write(6,*) 'dk=',dk
      write(6,*) 'dlambda1=',dlambda1
      write(6,*) 'dlambda2=',dlambda2
      write(6,*) 's1=',s1
      write(6,*) 's2=',s2
      write(6,*) 's=',s1*s2
      write(6,*) 'depsilon=',depsilon
      write(6,*) 'depsilonp=',depsilonp
      write(6,*) 'vd=',vd,'dt=',dt,'tend=',tend,'tE=',tE
      ! good for paper 1 dimensional
!      write(6,*) 'dimensionless xlength=',150*dk1/(dke0*Diffgamma)**0.5
!      write(6,*) 'dimensionless ylength=',4*dk1/(dke0*Diffgamma)**0.5
!      dx=150.d0/Nx*dk1/(dke0*Diffgamma)**0.5
!      dy=dx
      ! good for paper 2 dimensional
      write(6,*) 'dimensionless xlength=',75*dk1/(dke0*Diffgamma)**0.5
      write(6,*) 'dimensionless ylength=',2*dk1/(dke0*Diffgamma)**0.5
!      dx=0.025*dk1/(dke0*Diffgamma)**0.5
      dx=50.d0/Nx*dk1/(dke0*Diffgamma)**0.5
      dy=2.d0/Ny*dk1/(dke0*Diffgamma)**0.5
      write(6,*) 'Qqqqqqqqqq',(dke0*Diffgamma)**0.5/dk1
!      write(6,*) 'dimensionless xlength=',75*dk1/(dke0*Diffgamma)**0.5
!      write(6,*) 'dimensionless ylength=',2*dk1/(dke0*Diffgamma)**0.5
!      dx=75.d0/Nx*dk1/(dke0*Diffgamma)**0.5
!      dy=dx
      write(6,*) 'dx=',dx,'dy=',dy

      open(7,file = 'OutputData2D/Initial-State', err = 20,
     .       status = 'old', form = 'formatted')

      call loadState(t,Nx,Ny,gamma,ro)
      tE=0
      go to 10

 20   call ic(t,Nx,Ny,gamma,ro,gammanullc,ro1nullc,ro2nullc)

 10   return
      
      end
