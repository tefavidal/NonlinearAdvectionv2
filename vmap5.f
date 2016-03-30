      subroutine vmap5(ier,Nx,Nt,gammax)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gammax(1000,Nx)
      real gammax1(Nx,1000)
      real tr(6),gammaxmin,gammaxmax,ratio
      common /param/ gamma01,ro01,Diffgamma,dke0,dk1
      integer ier,pgbeg
      character(len=30) ct

      do i=1,Nt
       do j=1,Nx
        gammax1(j,i)=gammax(i,j)

       enddo
      enddo


         if(ier .ne. 1)stop
c        if ier is not 1 then there is an error

         ratio=0.5

         call pgpap(12.0,ratio)
!
         dx1=dx*(dke0*Diffgamma)**0.5/dk1
         dt1=tout/dk1
c        change the size of view surface pgpap(WIDTH,ASPECT)
         tr(1) = 0.
         tr(2) = real(dx1)
         tr(3) = 0.
         tr(4) = 0.
         tr(5) = 0.
         tr(6) = real(dt1)

         call pgsch(2.)
         call pgslw(3)
         call pgsvp(0.1,0.7,0.4,0.95)

         call pgswin(real(0),real(Nx*dx1),real(0),real((Nt-1)*dt1))

c        set window (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)

         call pgbox('BCTNSP',10.0,0,'BCTNSP',10.0,0)
c        draw labeled frame around viewport (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
         call pgsch(2.)
         call pglab('x(mm)','t(min)','')

         call pgsch(2.)

         call pallette(1.,0.6)
         gammaxmin = gammax1(1,1)
         gammaxmax = gammax1(1,1)

         do i = 1,Nx
          do j = 1,Nt
             gammaxmin = min(gammaxmin,gammax1(i,j))
             gammaxmax = max(gammaxmax,gammax1(i,j))
           enddo
         enddo
        write(6,*) 'gammax=',gammaxmin,gammaxmax

         
         call pgimag(gammax1,Nx,Nt,1,Nx,1,Nt,gammaxmin,gammaxmax,tr)
c        color image from a 2D data array pgimag(A, IDIM, JDIM, I1, I2, J1, J2,D12, A2, TR)
         call pgwedg('RI',0.5,3.,gammaxmin,gammaxmax,'\gg')
c        annotate an image plot with wedge PGWEDG(SIDE, DISP, WIDTH, FG, BG, LABEL)




      end


