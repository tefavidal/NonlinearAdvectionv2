      subroutine vmap3(ier,k1,XGa,Xro,gammanullc,ro1nullc,ro2nullc
     . ,ro1nullcNF)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob


      double precision gammanullc(274),ro1nullc(274),ro2nullc(274)
      double precision ro1nullcNF(274)
      real gammanullc1(274),ro1nullc1(274),ro2nullc1(274)
      real ro1nullcNF1(274)
      double precision XGa(1000),Xro(1000)
      real Xro1(1000),XGa1(1000)
      common /param/ gamma01,ro01,Diffgamma,dke0,dk1
      integer ier,pgbeg
      character(len=30) ct

      do i=1,274
      gammanullc1(i)=gammanullc(i)
      ro1nullc1(i)=ro1nullc(i)
      ro2nullc1(i)=ro2nullc(i)
      ro1nullcNF1(i)=ro1nullcNF(i)
      enddo
      do i=1,k1
      Xro1(i)=Xro(i)
      XGa1(i)=XGa(i)
      enddo


      t=(k1-1)*tout/dk1
      write(ct,'(F6.2)') t
      ct='t= '//ct
      ct=ct(1:12)//'min'
      call PGSLW (10)
      CALL PGSCH(2.0)
!      call pgsvp(0.2,0.8,0.4,0.95)
!      call PGSLS (1)
      CALL PGSCI(1)
      CALL PGENV(0.,1.0,0.,1.0,0,1)
      call pglab('\gg','\(0643)','')
!      CALL PGVSIZ (1.6, 8.9, 1.55, 6.7)

      call pgsci(1)
      call PGPTXT (0.4, 0.7, 0.0, 0.5, ct)
      call PGSLW (5)
      call pgsci(8)
      call PGPT (k1, XGa1, Xro1, -10)
!      call PGSLS (2)
      CALL PGSCI(5)
      call PGLINE (k1, XGa1, Xro1)
      CALL PGSCI(4)
      call PGPT (1, XGa1(k1), Xro1(k1),-5)
      call PGSLW (10)
      CALL PGSCI(2)
      call PGLINE(273,gammanullc1, ro1nullc1)
      CALL PGSCI(1)
      call PGLINE(273,gammanullc1, ro2nullc1)
!      call PGSLS (2)
      call pgsci(1)
      call PGLINE(273,gammanullc1, ro1nullcNF1)

      end

      
