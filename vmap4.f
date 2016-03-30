      subroutine vmap4(ier,it,Nx,gammax)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gammax(500,Nx)
      real gammax1(500,Nx)
      real tr(6),gammaxmin,gammaxmax,ratio

      integer ier,pgbeg
      character(len=30) ct

      do i=1,it
       do j=1,Nx
        gammax1(i,j)=gammax(i,j)

       enddo
      enddo




!         call pgpage

!         call pgsubp (2,2)

c        pgbeg(UNIT,file/type,NXSUB,NYSUB)

         if(ier .ne. 1)stop
c        if ier is not 1 then there is an error

         ratio=0.7

!         if (Nx .eq. Ny) then
          call pgpap(8.0,ratio)
!         else
!          call pgpap(12.5,ratio)
!         endif

c        change the size of view surface pgpap(WIDTH,ASPECT)
         tr(1) = 0.
         tr(2) = real(tout)
         tr(3) = 0.
         tr(4) = 0.
         tr(5) = 0.
         tr(6) = real(dx)
         
!         if (Nx .eq. Ny) then
          call pgsch(4.)
          call pgslw(3)
          call pgsvp(0.2,0.8,0.4,0.95)
!         else
!         call pgsch(8.)
c        set character height
!         call pgslw(3)
c        set line width
!         call pgsvp(0.1,0.7,0.4,0.95)

!         endif
         call pgswin(real(0),real(it*tout),
     .        real(0),real(Nx*dx))
c        set window (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
!         if (Nx .eq. Ny) then
         call pgbox('BCTNSP',0.2,0,'BCTNSP',5.0,0)
c        draw labeled frame around viewport (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
         call pgsch(4.)
         call pglab('t','x','')

         call pgsch(4.)
!         else
!         call pgbox('BCTNSP',4.0,0,'BCTNSP',0.5,0)
!c        draw labeled frame around viewport (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
!         call pgsch(5)
!         call pglab('x','y','')
!         call pgsch(8.)
!         endif
         call pallette(1.,0.6)
         gammaxmin = gammax1(1,1)
         gammaxmax = gammax1(1,1)

         do i = 1,it
          do j = 1,Nx
             gammaxmin = min(gammaxmin,gammax1(i,j))
             gammaxmax = max(gammaxmax,gammax1(i,j))
           enddo
         enddo

!        write(6,*) 'gamma=',gammaxmin,gammaxmax

         
         call pgimag(gammax1,it,Nx,1,it,1,Nx,gammaxmin,gammaxmax,tr)
c        color image from a 2D data array pgimag(A, IDIM, JDIM, I1, I2, J1, J2,D12, A2, TR)
         call pgwedg('RI',0.5,3.,gammaxmin,gammaxmax,'\gg')
c        annotate an image plot with wedge PGWEDG(SIDE, DISP, WIDTH, FG, BG, LABEL)
!         call pgmtxt('LV', 1.0, -0.6, 0.0, '      t=   min')
!         call pgmtxt('LV', 0.99, -0.6, 0.0, ct1)
!         call pgmtxt('LV', 1.0, -0.6, 0.0, ct)



      end


