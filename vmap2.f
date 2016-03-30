      subroutine vmap2(ier,k,Nx,XGamma)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision XGamma(Nx)
      real XGamma1(Nx),XPTS(Nx),gammaMax,gammaMin
      integer Nx
      integer ier,pgbeg
      character(len=30) ct

      gammaMax=XGamma(1)
      gammaMin=XGamma(1)
       do i=1,Nx
        XGamma1(i)=XGamma(i)
        XPTS(i)=i*dx
        if (XGamma1(i) .gt. gammaMax) then
        gammaMax=XGamma(i)
        endif
        if (XGamma1(i) .lt. gammaMin) then
        gammaMin=XGamma(i)
        endif
       enddo
      gammaMax=real(1.2)
      gammaMin=real(0.1)


      t=k*tout
      write(ct,'(F6.2)') t
      ct='t= '//ct
      ct=ct(1:12)//'  '
      call PGSLW (10)
      CALL PGSCH(2.0)
      call pgsvp(0.2,0.8,0.4,0.95)
      CALL PGENV(0.,real(Nx*dx),0.,9.0,0,1)
      call pglab('x','\gg','')

      call PGLINE(Nx,XPTS, XGamma1)
      call PGPTXT (2.0, 9.2, 0.0, 0.5, ct)
!      call PGPT (Nx, XPTS, XGamma1, -10)

      end

      
