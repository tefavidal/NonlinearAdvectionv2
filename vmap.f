      subroutine vmap(ier,t,Nx,Ny,gamma)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma(Nx,Ny)
      real gamma1(Nx,Ny)
      real tr(6),gammamin,gammamax,ratio
      common /param/ gamma01,ro01,Diffgamma,dke0,dk1
      integer ier,pgbeg
      character(len=30) ct

      do j=1,Ny
        do i=1,Nx
            gamma1(i,j)=gamma(i,j)

       enddo
      enddo
      it = nint(t/tout)


      write(ct,'(F6.2)') t/dk1

      ct='t= '//ct
      ct=ct(1:12)//'min'



!         call pgpage

!         call pgsubp (2,2)

c        pgbeg(UNIT,file/type,NXSUB,NYSUB)

         if(ier .ne. 1)stop
c        if ier is not 1 then there is an error

!         ratio=real(4*Ny)/real(Nx)
         ratio=real(320)/real(3000)

         if (Nx .eq. Ny) then
          call pgpap(5.0,ratio)
         else
          call pgpap(12.5,ratio)
         endif
         dx1=dx*(dke0*Diffgamma)**0.5/dk1
         dy1=dy*(dke0*Diffgamma)**0.5/dk1

c        change the size of view surface pgpap(WIDTH,ASPECT)
         tr(1) = 0.
         tr(2) = real(dx1)
         tr(3) = 0.
         tr(4) = 0.
         tr(5) = 0.
         tr(6) = real(dy1)
         
         if (Nx .eq. Ny) then
          call pgsch(5.)
          call pgslw(3)
          call pgsvp(0.2,0.8,0.4,0.95)
         else
         call pgsch(5.)
c        set character height
         call pgslw(3)
c        set line width
         call pgsvp(0.1,0.7,0.4,0.95)

         endif
         call pgswin(real(0),real(Nx*dx1),
     .        real(0),real(Ny*dy1))
         write(6,*) real(Ny*dy1)
c        set window (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
         if (Nx .eq. Ny) then
         call pgbox('BCTNSP',5.0,0,'BCTNSP',5.0,0)
c        draw labeled frame around viewport (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
         call pgsch(5.)
         call pglab('x(mm)','y(mm)'
     .        ,'')

         call pgsch(5.)
         else
         call pgbox('BCTNSP',5.0,0,'BCTNSP',1.0,0)
c        draw labeled frame around viewport (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
         call pgsch(5)
         call pglab('x (mm)','y (mm)','')
         call pgsch(5.)
         endif
         call pallette(1.,0.6)
         gammamin = gamma1(1,1)
         gammamax = gamma1(1,1)

        do j = 1,Ny
         do i = 1,Nx
             gammamin = min(gammamin,gamma1(i,j))
             gammamax = max(gammamax,gamma1(i,j))
           enddo
         enddo
        write(6,*) 'gamma=',gammamin,gammamax
        gammamin=0
        gammamax=6
         
         call pgimag(gamma1,Nx,Ny,1,Nx,1,Ny,gammamin,gammamax,tr)
c        color image from a 2D data array pgimag(A, IDIM, JDIM, I1, I2, J1, J2,D12, A2, TR)
         call pgwedg('RI',0.5,3.,gammamin,gammamax,'\gg')
c        annotate an image plot with wedge PGWEDG(SIDE, DISP, WIDTH, FG, BG, LABEL)
!         call pgmtxt('LV', 1.0, -0.6, 0.0, '      t=   min')
!         call pgmtxt('LV', 0.99, -0.6, 0.0, ct1)
         call pgmtxt('LV', 1.0, -0.6, 0.0, ct)



      end

      subroutine pallette(contra,bright)
      
      real contra,bright
      real l(7),r(7),g(7),b(7)
      
      data l / 0.,0.1667,0.3333,0.5,0.6667,0.8333,1. /
      data b / 1.,1.,1.,0.5,0.,0.,0. /
      data r / 0.,0.,0.,0.5,1.,1.,1. /
      data g / 0.,0.5,1.,1.,1.,0.5,0. /
      
      call pgctab(l,r,g,b,7,contra,bright)
c     install color table to be used by PGCTAB(L, R, G, B, NC, CONTRA, BRIGHT)
      
      end
