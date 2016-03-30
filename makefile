all:SimGoldbeter
SimGoldbeter:	
	gfortran -o program5 main-Source.f anfang.f vmap.f vmap2.f vmap3-Source.f cg.f \
	vmap5.f rs-Nonlinear.f out.f ODE-Merson.f ic-Source.f \
	-L/usr/bmp/pgplot-5.2/ -lpgplot \
	-L/usr/bmp/slatec-4.1/lib -lslatec \
	-L/usr/bmp/lapack-3.4.0 \
	-lX11 \
	-lpng

	