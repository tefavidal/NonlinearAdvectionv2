      subroutine ic(t,Nx,Ny,gamma,ro,gammanullc,ro1nullc,ro2nullc)
      
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      double precision gamma(Nx,Ny),ro(Nx,Ny)
      dimension gamma0(10),ro0(10)
      dimension ar(2,2),ai(2,2),wr(2),wi(2),zr(2,2),zi(2,2)
      dimension fv1(2),fv2(2),fv3(2)
      dimension dkmax(10), omega(10), wmax(10)
      common /param/ gamma01,ro01,Diffgamma,dke0,dk1
      double precision dke(Nx,Ny),dsigma(Nx,Ny)
      double precision gammanullc(274),ro1nullc(274),ro2nullc(274)
      dimension dkc0(10)

      dflux=0.d0
      ii=101
      call nullcline(ii,dflux,gammanullc,ro1nullc,ro2nullc)

      call SteadyState(ii,dflux,nfix,gamma0,ro0)
      write(6,*) 'nfix==',nfix

      do i=1,nfix
      write(6,*)  'gamma=',gamma0(i),'  ro=',ro0(i)
      enddo

      do i=1,Nx
       do j=1,Ny

!        if (i .lt. 1510 .and. i .gt. 1500) then
       if (i .lt. 110 .and. i .gt. 100) then
        ro(i,j)=ro0(1)
        if (j .lt. 50 .and. j .gt. 30) then
            gamma(i,j)=gamma0(1)+100*gamma0(1)
        else
            gamma(i,j)=gamma0(1)
        endif
!        gamma(i,j)=gamma0(1)

       elseif(i .le. 100) then
        ro(i,j)=0
        gamma(i,j)=0
       else
        ro(i,j)=ro0(1)
        gamma(i,j)=gamma0(1)

       endif
!         RANDOM CELL FIRING
!         if (i .gt. 100) then
!         ran=rand()
!      if(ran .lt. 0.0001) then
!
!        gamma(i,j)=gamma0(1)+0.1
!
!
!       else
!
!        gamma(i,j)=gamma0(1)
!       endif
!       ro(i,j)=ro0(1)
!       else
!       ro(i,j)=0
!        gamma(i,j)=0
!
!       endif
       enddo
      enddo

      do i = 1,1
         gamma01=gamma0(i)

         ro01=ro0(i)


         dkmax(i) = 0.d0
         wmax(i) = -1.d0

         do k = 1,20000
            dknum = 1.d-3*(k-1)
            call array(dknum,gamma01,ro01,ar,ai)

            call cg(2,2,ar,ai,wr,wi,1,zr,zi,fv1,fv2,fv3,ierr)


            w = max(wr(1),wr(2))

            do keig = 1,2
               if(w .eq.  wr(keig))then
                  freq = wi(keig)
               endif
            enddo
            if(k .eq. 1)then
               omega(i) = freq
            endif
            if(wmax(i) .lt. w)then

               dkmax(i) =dknum

               omega(i) = freq
               wmax(i)  = w
            endif

           write(12,*) dknum,w,freq

         enddo
      close(12)

      s=s1*s2
      Y0=ro01*gamma01/(1.d0+gamma01)
      A=(dk*dL2-dL1)*dc/(1+dc*gamma01)**2
      B=(dk-1.d0)/(1+gamma01)**2

      write(6,*) '-----------------------------------------------------'
      write(6,*) 'A=',A,'B=',B

      dM=2*ro01*gamma01**2*(dlambda2-dlambda1)
     . /((1+gamma01)**2*(dlambda2+Y0**2)**2)
      dN=2*ro01**2*gamma01*(dlambda2-dlambda1)
     . /((1+gamma01)**3*(dlambda2+Y0**2)**2)
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'M=',dM,'N=',dN,'s=',s
      f1=(1.d0+dk*gamma01)/(1.d0+gamma01)
      f2=(dL1+dk*dL2*dc*gamma01)/(1.d0+dc*gamma01)
      SS=-(f1+f2)+(dN*s-1)/depsilon
      Delta=((f1+f2)*(1.d0-dN*s)+s*dM*B*ro01
     . -s*dM*A*(1-ro01))/(depsilon)
      a11=(dN*s-1)/depsilon
      a22=-f1-f2
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'a11=(N*s-1)/epsilon=',a11,'a22=-f1-f2=',a22
      write(6,*) '-----------------------------------------------------'
      write(6,*)'Stabilty of two-variable subsystem:S<0 and Delta>0'
      write(6,*) 'Delta (determinant)=',Delta
      write(6,*) 'S (trace)=',SS

      write(6,*) '-----------------------------------------------------'
      if ((dN*s-1)/depsilon .gt. 0) then
      dkmaximum=(dN*s-1)**0.5/depsilon
      write(6,*) '********kmax=',dkmaximum

      write(6,*) '-----------------------------------------------------'
      do k = 1,20000
            dknum = 1.d-3*(k-1)
            if (dknum .lt. dkmaximum) then

       vdthreshold=((-Delta/a22+depsilon*dknum**2)*
     .(depsilon*dknum**2-SS)**2/(dknum**2*(a11-depsilon*dknum**2)))**0.5
            write(14,*) dknum,vdthreshold
            endif
      enddo
      endif
      close(14)

      enddo
      call SteadyKc(a11,a22,SS,Delta,dkc0,n)
      write(6,*) 'Numerical kc(1)=',dkc0(1),'kc(2)=',
     . dkc0(2),
     . 'number of solutions=',n
      write(6,*) 'Exact vwave(kc)/vflow='
     . ,a22/(a11+a22-depsilon*dkc0(1)**2)
     . ,a22/(a11+a22-depsilon*dkc0(2)**2)
      write(6,*) 'epsilon*kc**2=',depsilon*dkc0(1)**2,
     . depsilon*dkc0(2)**2
      do k=1,100
      dknum=k
      call functionKc(a11,a22,SS,Delta,dknum,v)
      write(30,*) dknum,v
      enddo
      call functionKc(a11,a22,SS,Delta,dkc0(1),v)

      write(6,*) '-----------------------------------------------------'
      aa=a22
      bb=a22**2-a22*a11
      cc=-Delta*(SS+a22)
      dd=SS*a11*Delta
      Del=18*aa*bb*cc*dd-4*bb**3*dd+bb**2*c**2-4*aa*cc**3-27*aa**2*dd**2
      dp1=(3*aa*cc-bb**2)/(3*aa**2)
      dq1=(2*bb**3-9*aa*bb*cc+27*aa**2*dd)/(27*aa**3)
      angle=(3*dq1)/(2*dp1)*(-3/dp1)**0.5
      dpi=3.1415
      t1=2*(-dp1/3)**0.5*cos((acos(angle))/3-2*dpi/3)
      t2=2*(-dp1/3)**0.5*cos((acos(angle))/3-4*dpi/3)
      t3=2*(-dp1/3)**0.5*cos((acos(angle))/3-6*dpi/3)
      x1=t1-bb/(3*aa)
      x2=t2-bb/(3*aa)
      x3=t3-bb/(3*aa)

      write(6,*) 'Del positive means 3 real roots: Del=',Del
      write(6,*) 't1=',t1,'t2=',t2,'t3=',t3
      write(6,*) 'x1=',x1,'x2=',x2,'x3=',x3
      write(6,*) 'Analytical kc(1)=',(x1/depsilon)**0.5,
     . 'kc(3)=',(x3/depsilon)**0.5
      vmin=((-Delta/a22+x1)*(SS-x1)**2/((x1/depsilon)*(a11-x1)))**0.5
      write(6,*) 'Exact vmin (dimentionless)=',vmin
      write(6,*) 'Exact vmin (dimentional)=',vmin*(dke0*Diffgamma)**0.5
      write(6,*) 'Analytical vwave/vflow=',a22/(-SS+x1)
      write(6,*) 'Approximation vwave/vflow=', -0.5-0.5*a22/(a11+a22)
      write(6,*) 'approximation for kc~',(-dd/(cc*depsilon))**0.5,
     . ((a11+a22)*a11/(depsilon*(2*a22+a11)))**0.5
      write(6,*) 'approximation for vmin~',
     . (depsilon*(-Delta*SS-Delta*a22+SS*a11*a22)*(SS*a11-SS**2-SS*a22)
     . **2/(SS*a11**2*a22**2*(SS+a22)))**0.5
!     .,(4*depsilon*SS*(-Delta*SS-Delta*a22+SS*a11*a22)/
!     .(a11**2*(SS+a22)))**0.5

      write(6,*) '-----------------------------------------------------'
      do i=1,40
      write(40,*) i,aa*i**3+bb*i**2+cc*i+dd,cc*i+dd
      write(40,*) -i,aa*(-i)**3+bb*(-i)**2+cc*(-i)+dd,cc*(-i)+dd
      enddo
      i=0
      write(40,*) -i,aa*(-i)**3+bb*(-i)**2+cc*(-i)+dd,cc*(-i)+dd
      return

      end

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine SteadyState(ii,dflux,n,gamma0,ro0)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01

      dimension gamma0(10),ro0(10)
!     Apliying the constrain dtro=0 looks for when dtgamma changes
!     sign, when it does, calls zero gamma
      n = 0

      fixpoint = 0.d0

      do i = 1,90
      gamma=0.001+0.0001d0*i
      call function(ii,dflux,gamma,ro,vnew)
!     returns ro such that dtro=0 for that gamma
!     vnew=dtgamma

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerogamma(ii,dflux,gamma-0.0001,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      gamma=0.01+0.001d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-0.001,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

      do i = 1,91
      gamma=0.1+0.01d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-0.01,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

      do i = 1,91
      gamma=1+0.1d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-0.1,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

       do i = 1,91
      gamma=10+1.d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-1.d0,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

       do i = 1,91
      gamma=100+10.d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-10.d0,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

       do i = 1,91
      gamma=1000+100.d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-100.d0,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo



      return
      end
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine function(ii,dflux,gamma,ro,v)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01

      x=(1.d0*ii-50.d0)/10
      s22=s2*(1+dtanh(x))/2

      f1=(1.d0+dk*gamma)/(1.d0+gamma)
      f2=(dL1+dk*dL2*dc*gamma)/(1.d0+dc*gamma)
      ro=f2/(f1+f2)

      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      v=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux


      return
      end

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine zerogamma(ii,dflux,gamma1,gamma2,gamma,ro,v)

      implicit double precision (a-h, o-z)

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01

      call function(ii,dflux,gamma1,ro1,v1)
      call function(ii,dflux,gamma2,ro2,v2)

 10   gamma=(gamma1+gamma2)/2
      call function(ii,dflux,gamma,ro,v)
      if(v1*v .lt. 0.d0) then
      gamma2=gamma

      ro2=ro
      v2=v
      else
      gamma1=gamma

      ro1=ro
      v1=v
      endif
      if (dabs(gamma2-gamma1) .lt. 1.d-12) then
        go to 25
      endif
      if(dabs(v) .gt. 1.d-6) then
         goto 10
      endif

 25   return
      end
c     *****************************************************************
      subroutine array(dknum,gamma01,ro01,ar,ai)

      implicit double precision (a-h,o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      dimension ar(2,2),ai(2,2)



      do i = 1,2
         do j = 1,2
            ar(i,j) = 0.d0
            ai(i,j) = 0.d0
         enddo
      enddo

      Y0=ro01*gamma01/(1.d0+gamma01)
      dM=2*ro01*gamma01**2*(dlambda2-dlambda1)
     . /((1+gamma01)**2*(dlambda2+Y0**2)**2)
      dN=2*ro01**2*gamma01*(dlambda2-dlambda1)
     . /((1+gamma01)**3*(dlambda2+Y0**2)**2)

      f1=(1.d0+dk*gamma01)/(1.d0+gamma01)
      f2=(dL1+dk*dL2*dc*gamma01)/(1.d0+dc*gamma01)
      A=(dk*dL2-dL1)*dc/(1+dc*gamma01)**2
      B=(dk-1)/(1+gamma01)**2

      ar(1,1) =-1.d0/depsilon-depsilon*dknum**2+dN*s1*s2/depsilon
      ai(1,1)=-vd*dknum
      ar(1,2) =s2*s1*dM/depsilon

      ar(2,1) =-B*ro01+A*(1-ro01)
      ar(2,2)=-f1-f2



      return

      end

c     *****************************************************************
      subroutine nullcline(ii,dflux,gammanullc,ro1nullc,ro2nullc)
!     calculates the two nullclines dtro=0 and dtgamma=0
!     dtro=0 is stored in ro2nullc
!     dtgamma=0 is stored in ro1nullc
!     the corresponding gamma are stored in gammanullc
!     The three vectors are writed in 13
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01
      double precision gammanullc(274),ro1nullc(274),ro2nullc(274)


      k=1

      do i = 1,91
      gamma=0.001+0.0001d0*(i-1)

      call fnullcline2(ii,dflux,gamma,ro1,ro2)
!     saves in ro2 such that dtro=0 for that gamma
!     saves in ro1 such that dtgamma=0 for that gamma
      write(13,*) gamma,ro1,ro2

      gammanullc(k)=gamma

      ro1nullc(k)=ro1
      ro2nullc(k)=ro2
      k=k+1

      enddo

      do i = 1,91
      gamma=0.01+0.001d0*(i-1)

      call fnullcline2(ii,dflux,gamma,ro1,ro2)

      write(13,*) gamma,ro1,ro2
      gammanullc(k)=gamma
      ro1nullc(k)=ro1
      ro2nullc(k)=ro2
      k=k+1

      enddo

      do i = 1,91
      gamma=0.1+0.01d0*(i-1)

      call fnullcline2(ii,dflux,gamma,ro1,ro2)

      write(13,*) gamma,ro1,ro2
      gammanullc(k)=gamma
      ro1nullc(k)=ro1
      ro2nullc(k)=ro2
      k=k+1

      enddo

      close(13)


      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine fnullcline2(ii,dflux,gamma,ro1,ro2)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01


      f1=(1.d0+dk*gamma)/(1.d0+gamma)
      f2=(dL1+dk*dL2*dc*gamma)/(1.d0+dc*gamma)
      ro2=f2/(f1+f2)

      call fnullcline1(ii,dflux,gamma,ro1)


      return
      end



!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine fnullcline1(ii,dflux,gamma,ro1)

      implicit double precision (a-h, o-z)

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      common /param/ gamma01,ro01

      x=(1.d0*ii-50.d0)/10
      s22=s2*(1+dtanh(x))/2
!     iterating in ro looks for when dtgamma changes sign
!     when it finds it calls zeroRho
      do i = 1,90
      ro=0.0001+0.00001d0*i
      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      vnew=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zeroRho(ii,dflux,gamma,ro-0.00001,ro,ro1,u)

            endif
      endif
      v=vnew
      enddo
      do i = 1,90
      ro=0.001+0.0001d0*i
      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      vnew=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zeroRho(ii,dflux,gamma,ro-0.0001,ro,ro1,u)
            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      ro=0.01+0.001d0*(i-1)
      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      vnew=(s1*s22*Phi-gamma)/depsilon+dflux
      Y0=ro01*gamma01/(1.d0+gamma01)
       dM=2*ro01*gamma01**2*(dlambda2-dlambda1)
     . /((1+gamma01)**2*(dlambda2+Y0**2)**2)
      dN=2*ro01**2*gamma01*(dlambda2-dlambda1)
     . /((1+gamma01)**3*(dlambda2+Y0**2)**2)
      vv=s1*s22*(dM*(ro-ro01)+dN*(gamma-gamma01))-(gamma-gamma01)
      if (vnew .lt. 0) then
      write(200,*) gamma,ro,-1
      else
       write(200,*) gamma,ro,1
      endif
      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
             call zeroRho(ii,dflux,gamma,ro-0.001,ro,ro1,u)
      endif
      endif
      v=vnew
      enddo

      do i = 1,91
      ro=0.1d0+0.01d0*(i-1)
      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      vnew=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux
      if (vnew .lt. 0) then
      write(200,*) gamma,ro,-1
      else
       write(200,*) gamma,ro,1
      endif
      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

             call zeroRho(ii,dflux,gamma,ro-0.01,ro,ro1,u)
            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      ro=1+0.1d0*(i-1)
      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      vnew=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux
      if (vnew .lt. 0) then
      write(200,*) gamma,ro,-1
      else
       write(200,*) gamma,ro,1
      endif
      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zeroRho(ii,dflux,gamma,ro-0.1,ro,ro1,u)
            endif
      endif
      v=vnew
      enddo
      write(200,*)

      return
      end
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zeroRho(ii,dflux,gamma,ro1,ro2,ro,v)
!     Calculates with great accuracy in ro the nullcline dtgamma=0
!     for a specific gamma.
!     the zero is located between ro1 and ro2
      implicit double precision (a-h,o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      x=(1.d0*ii-50.d0)/10
      s22=s2*(1+dtanh(x))/2
      Y=ro1*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)
      v1=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux

      Y=ro2*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)
      v2=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux

 10   ro=(ro1+ro2)/2
      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      v=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux
      if(v1*v .lt. 0.d0) then

      ro2=ro
      v2=v
      else

      ro1=ro
      v1=v
      endif
      if (dabs(ro2-ro1) .lt. 1.d-12) then
        go to 25
      endif
      if(dabs(v) .gt. 1.d-6) then
         goto 10
      endif

 25   return
      end
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine SteadyKc(a11,a22,SS,Delta,dkc0,n)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse



      dimension dkc0(10)

      n = 0

      fixpoint = 0.d0

      do i = 1,90
      dkc=0.001+0.0001d0*i
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-0.0001,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      dkc=0.01+0.001d0*(i-1)
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-0.001,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      dkc=0.1+0.01d0*(i-1)
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-0.01,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      dkc=1+0.1d0*(i-1)
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-0.1,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo

       do i = 1,91
      dkc=10+1.d0*(i-1)
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-1,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo

       do i = 1,91
      dkc=100+10.d0*(i-1)
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-10,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo

       do i = 1,91
      dkc=1000+100.d0*(i-1)
      call functionKc(a11,a22,SS,Delta,dkc,vnew)


      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerokc(a11,a22,SS,Delta,dkc-100,dkc,fixpoint,u)
            n=n+1
            dkc0(n)=fixpoint

            endif
      endif
      v=vnew
      enddo



      return
      end
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zerokc(a11,a22,SS,Delta,dkc1,dkc2,dkc,u)

      implicit double precision (a-h,o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse

      call functionKc(a11,a22,SS,Delta,dkc1,v1)

      call functionKc(a11,a22,SS,Delta,dkc2,v2)

 10   dkc=(dkc1+dkc2)/2
      call functionKc(a11,a22,SS,Delta,dkc,v)


      if(v1*v .lt. 0.d0) then

      dkc2=dkc
      else

      dkc1=dkc
      endif
      if (dabs(dkc2-dkc1) .lt. 1.d-12) then
        go to 25
      endif
      if(dabs(v) .gt. 1.d-6) then
         goto 10
      endif

 25   return
      end
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine functionKc(a11,a22,SS,Delta,dkc,v)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse


      Q=SS**2-4*Delta
      v=1/(2**0.5)*(((Q+depsilon**2*dkc**4
     . +(2*depsilon*(a22-a11)-vd**2)*dkc**2)**2+
     .(2*vd*depsilon*dkc**3-
     .2*vd*dkc*(a11-a22))**2)**0.5+Q+depsilon**2*dkc**4
     .+(2*depsilon*(a22-a11)-vd**2)*dkc**2)**0.5+(SS-depsilon*dkc**2)/2


      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine functionKc2(a11,a22,SS,Delta,dkc,v)

      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse



      Q1=a11-depsilon*dkc**2
      Q2=SS-depsilon*dkc**2
      Q3=-Delta/a22+depsilon*dkc**2
      v=depsilon*dkc**2*Q1*Q2-2*depsilon*dkc**2*Q1*Q3-Q1*Q2*Q3+
     .  depsilon*dkc**2*Q3*Q2



      return
      end
!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine stability(nfix,ii,gamma0,ro0)
      implicit double precision (a-h, o-z)
      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,tpulse
      dimension gamma0(10),ro0(10)
      dimension ar(2,2),ai(2,2),wr(2),wi(2),zr(2,2),zi(2,2)
      dimension fv1(2),fv2(2),fv3(2)
      dimension dkmax(10), omega(10), wmax(10)


      x=(1.d0*ii-50.d0)/10
      s22=s2*(1+dtanh(x))/2
       do i = 1,nfix
         gamma01=gamma0(i)

         ro01=ro0(i)


      s=s1*s22
      Y0=ro01*gamma01/(1.d0+gamma01)
      A=(dk*dL2-dL1)*dc/(1+dc*gamma01)**2
      B=(dk-1.d0)/(1+gamma01)**2

      dM=2*ro01*gamma01**2*(dlambda2-dlambda1)
     . /((1+gamma01)**2*(dlambda2+Y0**2)**2)
      dN=2*ro01**2*gamma01*(dlambda2-dlambda1)
     . /((1+gamma01)**3*(dlambda2+Y0**2)**2)

      f1=(1.d0+dk*gamma01)/(1.d0+gamma01)
      f2=(dL1+dk*dL2*dc*gamma01)/(1.d0+dc*gamma01)
      SS=-(f1+f2)+(dN*s-1)/depsilon
      Delta=((f1+f2)*(1.d0-dN*s)+s*dM*B*ro01
     . -s*dM*A*(1-ro01))/(depsilon)
      a11=(dN*s-1)/depsilon
      a22=-f1-f2

      write(6,*)'@Stabilty of two-variable subsystem:S<0 and Delta>0'
      write(6,*) 'gamma0=',gamma01
      write(6,*) 'Delta (determinant)=',Delta
      write(6,*) 'S (trace)=',SS
!      write(6,*) 'a11=',a11
!      write(6,*) 'a22=',a22
      enddo
      end
