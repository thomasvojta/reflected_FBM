      PROGRAM fbm_square
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   time-discrete reflected fractional Brownian motion on square
!
!   fbm_ring3    30 Mar 19 :    implemented soft wall
!   fbm_ring4    23 Jun 19 :    added circular obstacles 
!   fbm_ring5    23 Jun 19 :    measure 2d density averaged over several trajectories
!   fbm_square   28 Jan 20 :    square geometry 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define PARALLEL
#define VERSION 'fbm_square'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      integer(i4b), parameter     :: M=18,NT=2**M               ! number of time steps (in which mu const) 
      integer(i4b), parameter     :: NCONF=100                    ! number of obstacle configurations
      integer(i4b), parameter     :: NWALK=1                    ! number of walkers per obstacle config
      real(r8b), parameter        :: GAMMA = 0.4D0              ! FBM correlation exponent 
      
      real(r8b),parameter         :: XMAX = 100.D0             ! dimension of square x,y between -XMAX and XMAX 
      real(r8b),parameter         :: X0= 0.D0, Y0= 1.D0          ! starting point
      logical, parameter          :: RANDOMSTART=.true.          ! if true start from random position 
      
      character(4)                :: WALL='HARD'                  ! SOFT or HARD     
      real(r8b),    parameter     :: F0=0.2D0                      ! amplitude of wall force f=f0*exp(-lam*(r-RMAX))
      real(r8b),    parameter     :: LAM =1.D0                    ! decay constant of wall force   
      
      integer(i4b), parameter     :: NOBS=0                  ! number of obstacles 
      real(r8b), parameter        :: ROBS=2.D0                  ! radius of obstacles

      real(r8b), parameter        :: STEPSIG=1.0D0             ! sigma of individual step
      character(3)                :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   

      real(r8b), parameter        :: GRIDXMIN=-105.D0, GRIDXMAX=105.D0   ! grid dimensions for measuring density
      real(r8b), parameter        :: GRIDYMIN=-105.D0, GRIDYMAX=105.D0
      real(r8b), parameter        :: CELLSIZE=0.5D0                      ! size of grid cell  
      
      logical, parameter          :: WRITETRAJEC = .false.        ! write individual trajectories
      logical, parameter          :: WRITEOBSTAC = .false.        ! write obstacle file
      logical, parameter          :: WRITEDISTRIB = .true.       ! write density distribution       

      integer(i4b), parameter     :: IRINIT=1                    ! random number seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      real(r8b), parameter        :: ROBS2=ROBS*ROBS
      integer(i4b), parameter     :: NGX=NINT((GRIDXMAX-GRIDXMIN)/CELLSIZE)          ! number of grid cells in x direction
      integer(i4b), parameter     :: NGY=NINT((GRIDYMAX-GRIDYMIN)/CELLSIZE)          ! number of grid cells in y direction
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(r8b)              :: xx(0:NT), yy(0:NT)               ! walker coordinates
      real(r8b)              :: xix(1:NT), xiy(1:NT)             ! increments
      real(r8b)              :: xnew,ynew
      real(r8b)              :: rad,rad2                             ! distance of walker from origin
      real(r8b)              :: force, forcemax                   ! repulsive wall force 

      real(r8b)              :: confxx(1:NT)                     ! average of xx after each time step  
      real(r8b)              :: conf2xx(1:NT)                    
      real(r8b)              :: sumxx(1:NT),sum2xx(1:NT)         ! sums over machines in MPI version
      real(r8b)              :: auxxx(1:NT),aux2xx(1:NT) 

      integer(i4b)           :: iconf, iwalk, it                       ! configuration, walker, and time counters   
      integer(i4b)           :: totconf                         ! actual number of confs
      
      real(r8b)              :: xobs(NOBS),yobs(NOBS)           ! coordinates of obstacles
      integer(i4b)           :: iobs, jobs                             ! obstacle counter 
      real(r8b)              :: dist                             ! distance to obstacle
      logical                :: point_valid

      real(r8b)              :: confdens(NGX,NGY)                ! density array
      real(r8b)              :: normdens, lognormdens             ! normalized density 
      real(r8b)              :: cellx,celly                      ! center coordinates of each cell
      integer(i4b)           :: igx,igy                          ! grid cell counter
      real(r8b)              :: sumdens(NGX,NGY)                ! sums over machines in MPI version
      real(r8b)              :: auxdens(NGX,NGY)                ! 
                   
      external               :: kissinit 
      real(r8b),external     :: rkiss05,erfcc 

      character(18)          :: trafile = 'tra000000W0000.dat'
      character(13)          :: obsfile = 'obs000000.dat' 
      character(11)          :: disfile = 'distrib.dat'       
      
! Now the MPI stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      include 'mpif.h'
      integer(i4b)              :: ierr
      integer(i4b)              :: id,myid                  ! process index
      integer(i4b)              :: numprocs              ! total number of processes
      integer(i4b)              :: status(MPI_STATUS_SIZE)      
#endif             
  

! Start of main program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      totconf=(NCONF/numprocs)*numprocs
      
      if (myid==0) then
         print *,'Program ',VERSION,' running on', numprocs, ' processes'
         print *,'--------------------------------------------------'
      endif ! of if (myid==0)
#else
      totconf=NCONF
      print *,'Program ',VERSION,' running on single processor'
      print *,'--------------------------------------------------'
#endif 

      confxx(:)=0.D0 
      conf2xx(:)=0.D0 

      confdens(:,:)=0.D0 
                  
! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,totconf,numprocs
!      if (myid==0) print *, 'dis. conf.', iconf
#else
      disorder_loop: do iconf=1,totconf
!         print *, 'dis. conf.', iconf
#endif 

        call gkissinit(IRINIT+iconf-1)
        call initialize_obstacles                           ! place the obstacles   

        do iwalk=1,NWALK                                    ! loop over walkers for given obstacle config       
              
          call corvec(xix,NT,M,GAMMA)                         ! create x-increments (correlated Gaussian random numbers)
          call corvec(xiy,NT,M,GAMMA)                         ! create y-increments (correlated Gaussian random numbers)
        
          if (STEPDIS.eq.'BOX') then
            xix(:) = 1 - (0.5D0*erfcc(xix(:)/sqrt(2.0D0)))                       ! map onto unit inteval
            xix(:)=xix(:)-0.5D0                                                  ! center around 0
            xiy(:) = 1 - (0.5D0*erfcc(xiy(:)/sqrt(2.0D0)))                       ! map onto unit inteval
            xiy(:)=xiy(:)-0.5D0                                                  ! center around 0
          endif
          
          if (STEPDIS.eq.'BIN') then 
            xix(:)=sign(1.D0,xix(:))
            xiy(:)=sign(1.D0,xiy(:))
          endif  
         
          xix(:)=xix(:)*STEPSIG                               ! scale increments to correct sigma 
          xiy(:)=xiy(:)*STEPSIG
        
! Time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           call initialize_point  
           forcemax=0.D0
 

           if ( WALL.eq.'HARD') then           
             do it=1, NT/2
               xnew=xx(it-1)+xix(it)
               ynew=yy(it-1)+xiy(it)
               call test_point
               if (point_valid) then
                   xx(it)=xnew
                   yy(it)=ynew
               else
                   xx(it)=xx(it-1)
                   yy(it)=yy(it-1)                   
               endif 
               
               igx= 1+int((xx(it)-GRIDXMIN)/CELLSIZE)
               igy= 1+int((yy(it)-GRIDYMIN)/CELLSIZE)
               confdens(igx,igy)=confdens(igx,igy)+1.D0  
             end do
           endif  

!           if ( WALL.eq.'SOFT') then           
!             do it=1, NT/2
!               rad = sqrt( xx(it-1)*xx(it-1) + yy(it-1)*yy(it-1) )
!               force = F0*exp(LAM*(rad-RMAX)) - F0*exp(LAM*(RMIN-rad))
!               xx(it)=xx(it-1)+xix(it)-force*xx(it-1)/rad
!               yy(it)=yy(it-1)+xiy(it)-force*yy(it-1)/rad
!               print *,xx(it),yy(it),force
!               if (force.gt.forcemax) forcemax=force
!             end do
!           endif  
!           print *, 'forcemax=',forcemax

          if (WRITETRAJEC) then                          
            write(trafile(4:9),'(I6.6)') IRINIT+iconf-1
            write(trafile(11:14),'(I4.4)') iwalk
            open(2,file=trafile,status='replace')
            write(2,*) '# Program ', VERSION
            write(2,*) '# trajectory of reflected FBM on ring'
            write(2,*) '# step distribution ',STEPDIS
            write(2,*) '# STEPSIG=', STEPSIG
            write(2,*) '# GAMMMA=', GAMMA
            write(2,*) '# NT= ', NT
            write(2,*) '# XMAX=', XMAX
            write(2,*) '# NOBS=', NOBS
            write(2,*) '# ROBS=', ROBS
            write(2,*) '# IRINIT=',IRINIT
            write(2,*) '# RNG seed=', IRINIT+iconf-1
            write (2,*)'# =================================='
            write(2,*) '#    time         x         y'
            do it=0, NT/2
              Write(2,'(1X,I7,6(2X,E13.6))')  it, xx(it), yy(it)
            enddo 
            close(2) 
          endif

        enddo ! of do iwalk=1,NWALK   
        
                     
        write(obsfile(4:9),'(I6.6)') IRINIT+iconf-1
        open(2,file=obsfile,status='replace')
        write(2,*) '# Program ', VERSION
        write(2,*) '# positions of obstacles'
        write(2,*) '# step distribution ',STEPDIS
        write(2,*) '# STEPSIG=', STEPSIG
        write(2,*) '# GAMMMA=', GAMMA
        write(2,*) '# NT= ', NT
        write(2,*) '# XMAX=', XMAX
        write(2,*) '# NOBS=', NOBS
        write(2,*) '# ROBS=', ROBS
        write(2,*) '# IRINIT=',IRINIT
        write(2,*) '# RNG seed=', IRINIT+iconf-1
        write (2,*)'# =================================='
        write(2,*) '#   iobs         xobs      yobs      ROBS'
        do iobs=1,NOBS
          Write(2,'(1X,I7,6(2X,E13.6))')  iobs, xobs(iobs), yobs(iobs), ROBS
        enddo 
        close(2) 

      end do disorder_loop      ! of do incof=1,NCONF
      
#ifdef PARALLEL
     if (myid.ne.0) then                                                  ! Send data
         call MPI_SEND(confdens,NGX*NGY,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
      else
         sumdens(:,:)=confdens(:,:)
         do id=1,numprocs-1                                                   ! Receive data
            call MPI_RECV(auxdens,NGX*NGY,MPI_DOUBLE_PRECISION,id,1,MPI_COMM_WORLD,status,ierr)
            sumdens(:,:)=sumdens(:,:)+auxdens(:,:) 
         enddo
      endif        
#else          
      sumdens(:,:)=confdens(:,:)
#endif       
      
#ifdef PARALLEL
      if (myid==0) then
#endif
      if (WRITEDISTRIB) then
        open(2,file=disfile,status='replace')
#ifdef PARALLEL          
        write(2,*) '# Program ', VERSION, ' running on ', numprocs, ' CPUs'
#else          
        write(2,*) '# Program ', VERSION, ' running serially'
#endif 
        write(2,*) '# step distribution ',STEPDIS
        write(2,*) '# STEPSIG=', STEPSIG
        write(2,*) '# GAMMMA=', GAMMA
        write(2,*) '# NT= ', NT
        write(2,*) '# XMAX=', XMAX
        write(2,*) '# NOBS=', NOBS
        write(2,*) '# ROBS=', ROBS
        write(2,*) 'Number of configurations', totconf
        write(2,*) 'Number of walkers per config', NWALK
        write(2,*) '# IRINIT=',IRINIT
        write (2,*)'# =================================='
        write(2,*) '   igx  igy    x      y     distrib'
        do igx=1,NGX
        do igy=1,NGY
            cellx= GRIDXMIN + (1.D0*igx -0.5D0)*CELLSIZE
            celly= GRIDYMIN + (1.D0*igy -0.5D0)*CELLSIZE
            normdens=sumdens(igx,igy)/CELLSIZE**2 /(NWALK*totconf)/(0.5D0*NT)
            if (normdens.gt.0.D0) then
               lognormdens=log(normdens)
            else
               lognormdens=-1000
            endif   
            Write(2,'(1X,I4,1x,I4,4(1X,E13.6))')  igx, igy, cellx, celly, normdens, lognormdens
        enddo 
        enddo
        close(2)             
      endif 
#ifdef PARALLEL
      endif ! of if (myid==0)
#endif               
      
#ifdef PARALLEL
      call MPI_FINALIZE(ierr)
#endif 
      stop      
      
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine initialize_obstacles
! Places NOBS nonoverlapping circles of radius ROBS !!!!!!!!!!!!!!!!!      
      logical       ::    overlap
      iobs=0
      do while (iobs.lt.NOBS)
!        xnew=2.D0*RMAX*(rkiss05()-0.5D0)
        xnew=XMAX*rkiss05()
        ynew=2.D0*XMAX*(rkiss05()-0.5D0)
        if ( (abs(xnew).lt.XMAX) .and. (abs(ynew).lt.XMAX) ) then       ! in allowed square
           overlap=.false. 
           do jobs=1,iobs
              dist=(xnew-xobs(jobs))**2 + (ynew-yobs(jobs))**2
              if (dist.le.4.D0*ROBS2) overlap=.true.
           enddo
           if (.not.overlap) then
             iobs=iobs+1
             xobs(iobs)=xnew
             yobs(iobs)=ynew
           endif
        endif
      enddo        
      end subroutine initialize_obstacles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine initialize_point
      if (RANDOMSTART) then  
        do 
          xnew=2.D0*XMAX*(rkiss05()-0.5D0)
          ynew=2.D0*XMAX*(rkiss05()-0.5D0)
          call test_point
          if (point_valid) exit
        end do
         xx(0)=xnew
         yy(0)=ynew       
      else
         xx(0)=X0
         yy(0)=Y0  
      endif   
      end subroutine initialize_point      

                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine test_point

      point_valid=.true.
      
      if ( (abs(xnew).gt.XMAX) .or. (abs(ynew).gt.XMAX) ) then                  ! outside of allowed square
          point_valid=.false.
      else
         do iobs=1, NOBS                                  
            dist=(xnew-xobs(iobs))**2 + (ynew-yobs(iobs))**2
            if (dist.le.ROBS2) point_valid=.false.                              ! inside obstacle
         enddo
      endif
      end subroutine test_point

            
      END PROGRAM fbm_square
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine CORVEC(xr,Ns,MM)
!
! generates a 1d array of Ns Gaussian random numbers xr
! correlated according to a (translationally invariant)
! user-supplied correlation function corfunc(is,Ns)
! 
! uses Fourier filtering method
!
! history
!      v0.9         Dec  7, 2013:        first version, uses Tao Pang FFT
!      v0.91        Oct 11, 2017:        uses much faster FFT by Ooura (in fftsg.f)     
!      v0.92        Jun 23, 2019:        allocatable arrays    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE corvec(xr,Ns,MM,gam)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

      integer(i4b)           :: Ns              ! number of sites, must be power of 2 
      integer(i4b)           :: MM               ! Ns=2^M

      real(r8b)              :: xr(0:Ns-1)      ! random number array  
      real(r8b),allocatable  :: cr(:)      ! correlation function 
      integer(i4b)           :: is
      integer(i4b)           :: ipsize      
      integer(i4b),allocatable :: ip(:)   ! workspace for FFT code
      real(r8b),allocatable  :: w(:)           ! workspace for FFT code 
      
      real(r8b)              :: gam                   ! FBM exponent, pass through to correlation function   
      
      real(r8b), external    :: gkiss05,erfcc
      external               :: rdft                   ! from Ooura's FFT package
      real(r8b),external     :: fbmcorfunc 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ipsize= 2+sqrt(1.*Ns)
      allocate(cr(0:Ns-1))
      allocate(ip(0:ipsize))
      allocate(w(0:Ns/2-1))
      if (Ns.ne.2**MM) STOP 'Size indices do not match'
! Calculate correlation function 
      do is=0,Ns-1 
         cr(is)= fbmcorfunc(is,Ns,gam) 
      enddo
! Real FFT of correlation function
       ip(0)=0
       call rdft(Ns, 1, cr, ip, w)  
       
! Create array of independent Gaussian random numbers
      do is=0,Ns-1
         xr(is)= gkiss05()
      enddo
! Real FFT of input random numbers      
       call rdft(Ns, 1, xr, ip, w)  
! Filter the Fourier components of the random numbers
! as real-space correlations are symmmetric, FT of c is real
      do is=1,Ns/2-1
          xr(2*is)=xr(2*is)*sqrt(abs(cr(2*is)))*2.D0/Ns
          xr(2*is+1)=xr(2*is+1)*sqrt(abs(cr(2*is)))*2.D0/Ns
      enddo
      xr(0)=xr(0)*sqrt(abs(cr(0)))*2.D0/Ns
      xr(1)=xr(1)*sqrt(abs(cr(1)))*2.D0/Ns 
      
! FFT of filtrered random numbers (back to real space)
       call rdft(Ns, -1, xr, ip, w)  
       
! Transform from Gaussian distribution to flat distribution on (0,1)      
!      do is = 0,Ns-1
!        xr(is) = 1 -  0.5D0*erfcc(xr(is)/sqrt(2.0D0)) 
!      end do

      deallocate(cr,ip,w)      
      return
      END SUBROUTINE corvec  
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FBM correlation function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION fbmcorfunc(i,N,gam)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)              :: corfunc,fbmcorfunc
      real(r8b)              :: gam
      integer(i4b)           :: i,N
      integer(i4b)           :: dist

!      print *, 'gamma=',gam
      dist=min(i,N-i)
      if (dist.eq.0) then 
         corfunc=1.D0
      elseif (dist.lt.1000) then   
         corfunc = 0.5D0*( dble(dist+1)**(2.D0-gam) - 2.D0*(dble(dist)**(2.D0-gam)) + dble(dist-1)**(2.D0-gam) ) 
      else 
         corfunc = (2.D0-gam)*(1.D0-gam)*dble(dist)**(-gam) 
         corfunc = corfunc + (2.D0-gam)*(1.D0-gam)*(-gam)*(-1.D0-gam)*dble(dist)**(-gam-2.D0)/12.D0          
         corfunc = 0.5D0*corfunc
      endif   

      fbmcorfunc=corfunc
      
      return
      END FUNCTION fbmcorfunc 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gaussian random number generator gkiss05
!
! generates normally distributed independent random numbers
! (zero mean, variance 1) using Box-Muller method in polar form
!
! uniform random numbers provided by Marsaglia's kiss (2005 version)
!
! before using the RNG, call gkissinit(seed) to initialize
! the generator. Seed should be a positive integer.
!
!
! History:
!      v0.9     Dec  6, 2013:   first version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION gkiss05()
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)             :: gkiss05
      real(r8b), external   :: rkiss05

      real(r8b)             :: v1,v2,s,fac
      integer(i4b)          :: iset               ! switches between members of the Box-Muller pair
      real(r8b)             :: gset
      common /gausscom/gset,iset

      if (iset.ne.1) then
        do
          v1 = 2.D0 * rkiss05() - 1.D0
          v2 = 2.D0 * rkiss05() - 1.D0
          s = v1 * v1 + v2 * v2
          if ((s<1.D0) .and. (s>0.D0)) exit
        enddo
! Box-Muller transformation creates pairs of random numbers
        fac = sqrt(-2.D0 * log(s) / s)
        gset = v1 * fac
        iset = 1
        gkiss05 = v2 * fac
      else
        iset = 0
        gkiss05 = gset
      end if
      return
      END FUNCTION gkiss05


      SUBROUTINE gkissinit(iinit)
      implicit none
      integer,parameter     :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b)          :: iinit,iset
      real(r8b)             :: gset
      common /gausscom/gset,iset

      iset=0                         ! resets the switch between the members of the Box-Muller pair
      call kissinit(iinit)           ! initializes the rkiss05 RNG
      end subroutine gkissinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123
!
!
! A call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05 call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
!
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
!
!
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit
!        v0.93    Aug 13, 2012    changed integer representation test to avoid data statements
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      FUNCTION rkiss05()
      implicit none

      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers
      real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

      real(r8b)             :: rkiss05
      integer(i4b)          :: kiss
      integer(i4b)          :: x,y,z,w              ! working variables for the four generators
      common /kisscom/x,y,z,w

      x = 69069 * x + 1327217885
      y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss = ishft(x + y + ishft (z, 16) + w , -1)
      rkiss05=kiss*am
      END FUNCTION rkiss05


      SUBROUTINE kissinit(iinit)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b) idum,ia,im,iq,ir,iinit
      integer(i4b) k,x,y,z,w,c1,c2,c3,c4
      real(r8b)    rkiss05,rdum
      parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
      common /kisscom/x,y,z,w

      !!! Test integer representation !!!
      c1=-8
      c1=ishftc(c1,-3)
!     print *,c1
      if (c1.ne.536870911) then
         print *,'Nonstandard integer representation. Stoped.'
         stop
      endif

      idum=iinit
      idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
      if (idum.eq.0) idum=1
      if (idum.ge.IM) idum=IM-1

      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         x=idum+1
      else
         x=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         y=idum+1
      else
         y=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         z=idum+1
      else
         z=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         w=idum+1
      else
         w=idum
      endif

      rdum=rkiss05()

      return
      end subroutine kissinit



 FUNCTION erfcc(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates complementary error function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      real(r8b)   :: erfcc,x
      real(r8b)   :: t,z
      z=abs(x)
      t=1.D0/(1.D0+0.5D0*z)
      erfcc=t*exp(-z*z-1.26551223D0+t*(1.00002368D0+t*(.37409196D0+t*&
     &(.09678418D0+t*(-.18628806D0+t*(.27886807D0+t*(-1.13520398D0+t*&
     &(1.48851587D0+t*(-.82215223D0+t*.17087277D0)))))))))
      if (x.lt.0.D0) erfcc=2.D0-erfcc
      return
      END  