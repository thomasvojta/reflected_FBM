      PROGRAM fbm_cone2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   time-discrete reflected fractional Brownian motion
!   on 3d cone (sector of sphere with defined opening angle THETAM
!
!   fbm_sphere1 :  16 Aug 2019      first sphere version, builds on fbm_disk6 
!   fbm_cone1   :  13 Jan 2020      first cone version   
!   fbm_cone2   :  13 Jan 2020      add "narrow" density using THETAM/NARRFAC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define PARALLEL
#define VERSION 'fbm_cone2'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real(r8b), parameter        :: PI=3.14159265358979323D0

      integer(i4b), parameter     :: M=20, NT=2**M               ! number of time steps (in which mu const) 
      integer(i4b), parameter     :: NCONF=1e5                    ! number of walkers
      real(r8b), parameter        :: GAMMA = 0.4D0              ! FBM correlation exponent 
      
      real(r8b),parameter         :: RMAX = 1.D4                 ! radius
      real(r8b),parameter         :: THETAM=Pi/4                 ! opening angle of cone, allowed angles between 0 and THETAM
      real(r8b),parameter         :: X0= 0.D0, Y0=0.D0, Z0=0.01D0         ! starting point

      real(r8b), parameter        :: STEPSIG=1.D0                 ! sigma of individual step
      character(3)                :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   

      character(4)                :: WALL='HARD'                  ! SOFT or HARD     
      real(r8b),    parameter     :: F0=0.25D0*STEPSIG                      ! amplitude of wall force f=f0*exp(-lam*(r-RMAX))
      real(r8b),    parameter     :: LAM = 0.25D0/STEPSIG            ! decay constant of wall force            
            
      integer(i4b), parameter     :: NOUT = 500                    ! max number of output times, used as dimension of sumarrays  
      logical,parameter           :: WRITEDISTRIB = .TRUE.        ! write radial distribution    
      integer(i4b), parameter     :: NBIN =  50000                   ! number of bins for radial distribution
      real(r8b),parameter         :: NARRFAC = 5.D0                 ! narrow distribution between theta=0 and THETAM/NARRFAC
	  
      logical,parameter           :: WRITEPHIDISTRIB = .TRUE.        ! write angular distribution    
      integer(i4b), parameter     :: NPHIBIN =  200               ! number of bins for angular distribution
      logical,parameter           :: WRITETHETADISTRIB = .TRUE.        ! write angular distribution    
      integer(i4b), parameter     :: NTHETABIN =  200               ! number of bins for angular distribution
      integer(i4b), parameter     :: NTSTART=0.5*2**M             ! begin and end of measuring distribution
      integer(i4b), parameter     :: NTEND=2**M 

      integer(i4b), parameter     :: IRINIT=1                    ! random number seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	  real(r8b), parameter        :: PI2=3.14159265358979323D0/2.D0
      real(r8b), parameter        :: RMAX2=RMAX**2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(r8b)              :: xx(0:NT),yy(0:NT),zz(0:NT)                   ! walker coordinates
      real(r8b)              :: xix(1:2*NT),xiy(1:2*NT),xiz(1:2*NT)            ! increments

      integer(i4b)           :: outtime(1:NOUT)                    ! it indices of output times       
      real(r8b)              :: confrr(1:NOUT)                     ! average of rr at each output time
      real(r8b)              :: conf2rr(1:NOUT)                    
      real(r8b)              :: sumrr(1:NOUT),sum2rr(1:NOUT)         ! sums over machines in MPI version
      real(r8b)              :: auxrr(1:NOUT),aux2rr(1:NOUT) 

      integer(i4b)           :: iconf, it, ibin                    ! configuration, and time counters   
      integer(i4b)           :: itout, itoutmax
      integer(i4b)           :: totconf                         ! actual number of confs

      real(r8b)              :: rrdis(0:NBIN)                    ! density histogram 
      real(r8b)              :: sumdis(0:NBIN)
      real(r8b)              :: rrdisnarr(0:NBIN)                ! narrow density histogram 
      real(r8b)              :: sumdisnarr(0:NBIN)
      real(r8b)              :: PP,PPvol,PPvolerr                   ! P(r) 
      real(r8b)              :: rad,rad2,volume                     ! distance from origin
      real(r8b)              :: forcer 
                              
      real(r8b)              :: phidis(-NPHIBIN:NPHIBIN)                    ! angular histogram in xy plane
      real(r8b)              :: sumphidis(-NPHIBIN:NPHIBIN)
      real(r8b)              :: phi,PPphierr                      

      real(r8b)              :: thetadis(0:NTHETABIN)                    ! angular histogram in xy plane
      real(r8b)              :: sumthetadis(0:NTHETABIN)
      real(r8b)              :: theta,PPthetaerr                      


      external               :: kissinit 
      real(r8b),external     :: rkiss05,erfcc 

      character(15)          :: avrfile = 'avr00000000.dat'
      character(15)          :: disfile = 'dis00000000.dat' 
      character(15)          :: narfile = 'nar00000000.dat' 
      character(15)          :: phifile = 'phi00000000.dat' 
      character(15)          :: thetafile = 'the00000000.dat' 
            
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
         print *,'Program ',VERSION,' runing on', numprocs, ' processes'
         print *,'--------------------------------------------------'
      endif ! of if (myid==0)
#else
      totconf=NCONF
      print *,'Program ',VERSION,' runing on single processor'
      print *,'--------------------------------------------------'
#endif 

! Precalculate output times
      it=1
      itout=1
      do while (it.lt.NT)
         outtime(itout)=it
         it = max(it+1,int(it*1.05D0))
         itout=itout+1
      enddo      
      itoutmax=itout
      outtime(itoutmax)=NT 

      confrr(:)=0.D0 
      conf2rr(:)=0.D0
      rrdis(:)=0.D0 
      rrdisnarr(:)=0.D0	  
      phidis(:)=0.D0
	  thetadis(:)=0.D0

! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,totconf,numprocs
!      if (myid==0) print *, 'dis. conf.', iconf
#else
      disorder_loop: do iconf=1,totconf
      if (100*(iconf/100).eq.iconf)   print *, 'dis. conf.', iconf
#endif 

        call gkissinit(IRINIT+iconf-1)
        call corvec(xix,2*NT,M+1,GAMMA)                         ! create x-increments (correlated Gaussian random numbers)
        call corvec(xiy,2*NT,M+1,GAMMA)                         ! create y-increments (correlated Gaussian random numbers)
        call corvec(xiz,2*NT,M+1,GAMMA)                         ! create z-increments (correlated Gaussian random numbers)
        
        if (STEPDIS.eq.'BOX') then
          xix(:) = 1 - (0.5D0*erfcc(xix(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xix(:)=xix(:)-0.5D0                                                  ! center around 0
          xiy(:) = 1 - (0.5D0*erfcc(xiy(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xiy(:)=xiy(:)-0.5D0                                                  ! center around 0
          xiz(:) = 1 - (0.5D0*erfcc(xiz(:)/sqrt(2.0D0)))                       ! map onto unit inteval
          xiz(:)=xiz(:)-0.5D0                                                  ! center around 0
        endif
        
        if (STEPDIS.eq.'BIN') then 
          xix(:)=sign(1.D0,xix(:))
          xiy(:)=sign(1.D0,xiy(:))
          xiz(:)=sign(1.D0,xiz(:))
        endif  
        
        xix(:)=xix(:)*STEPSIG                               ! scale increments to correct sigma 
        xiy(:)=xiy(:)*STEPSIG                               ! scale increments to correct sigma 
        xiz(:)=xiz(:)*STEPSIG                               ! scale increments to correct sigma 
        
! Time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        xx(0)=X0
        yy(0)=Y0
        zz(0)=Z0
        if ( WALL.eq.'HARD') then                                                ! hard reflecting boundaries
           do it=1, NT
              xx(it)=xx(it-1)+xix(it)
              yy(it)=yy(it-1)+xiy(it)
              zz(it)=zz(it-1)+xiz(it)
              rad2= xx(it)**2 + yy(it)**2 +zz(it)**2
              phi= atan2( yy(it), xx(it) )
			  theta= atan2( sqrt(xx(it)**2 + yy(it)**2), zz(it) )
              if ( (rad2.gt.RMAX2).or.(theta.gt.THETAM) ) then
                 xx(it)=xx(it-1)
                 yy(it)=yy(it-1)
                 zz(it)=zz(it-1)
              end if
           end do  
        endif    
        if ( WALL.eq.'SOFT') then
           do it=1, NT
              rad = sqrt( xx(it-1)**2 + yy(it-1)**2 + zz(it-1)**2 )
              forcer = F0*exp(LAM*(rad-RMAX))
              if (forcer.gt.1.D-6*STEPSIG) then 
                xx(it)=xx(it-1) +xix(it) -forcer*xx(it-1)/rad
                yy(it)=yy(it-1) +xiy(it) -forcer*yy(it-1)/rad 
                zz(it)=zz(it-1) +xiz(it) -forcer*zz(it-1)/rad 
              else
                xx(it)=xx(it-1) +xix(it)
                yy(it)=yy(it-1) +xiy(it)
                zz(it)=zz(it-1) +xiz(it)
              endif  
           end do  
        endif     
        
           do itout=1,itoutmax
              rad2= xx(outtime(itout))**2 + yy(outtime(itout))**2 + zz(outtime(itout))**2
              confrr(itout)=confrr(itout) + sqrt(rad2)
              conf2rr(itout)=conf2rr(itout) + rad2
           enddo
           
           if (WRITEDISTRIB) then
              do it = NTSTART,NTEND 
                 rad=sqrt( xx(it)**2 + yy(it)**2 + zz(it)**2)
     			 theta= atan2( sqrt(xx(it)**2 + yy(it)**2), zz(it) )
                 ibin=nint( rad*NBIN/RMAX )
                 if ( (ibin.ge.0) .and. (ibin.le.NBIN)) then
                   rrdis(ibin)=rrdis(ibin)+1.D0 
                   if( theta.lt.THETAM/NARRFAC )  then
                     rrdisnarr(ibin)=rrdisnarr(ibin)+1.D0
                   endif  
                 endif    
              enddo  
           endif 

           if (WRITEPHIDISTRIB) then
              do it = NTSTART,NTEND 
                 phi=atan2( yy(it), xx(it) )                                         ! angle in zy plane
                 ibin=nint( phi*NPHIBIN/Pi )
                 if ( (ibin.ge.-NPHIBIN) .and. (ibin.le.NPHIBIN)) then
                   phidis(ibin)=phidis(ibin)+1.D0 
                 endif    
              enddo  
           endif 

           if (WRITETHETADISTRIB) then
              do it = NTSTART,NTEND 
     			 theta= atan2( sqrt(xx(it)**2 + yy(it)**2), zz(it) )
                 ibin=nint( theta*NTHETABIN/PI )
                 if ( (ibin.ge.0) .and. (ibin.le.NTHETABIN)) then
                   thetadis(ibin)=thetadis(ibin)+1.D0 
                 endif    
              enddo  
           endif 
                      
      end do disorder_loop      ! of do incof=1,NCONF

! Now collect and analyze data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef PARALLEL
      call MPI_REDUCE(confrr,sumrr,NOUT,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)   
      call MPI_REDUCE(conf2rr,sum2rr,NOUT,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)   
      call MPI_REDUCE(rrdis,sumdis,NBIN+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)   
      call MPI_REDUCE(rrdisnarr,sumdisnarr,NBIN+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)   
      call MPI_REDUCE(phidis,sumphidis,2*NPHIBIN+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)   
      call MPI_REDUCE(thetadis,sumthetadis,NTHETABIN+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)   
#else          
      sumrr(:)=confrr(:)
      sum2rr(:)=conf2rr(:)
      sumdis(:)=rrdis(:)
      sumdisnarr(:)=rrdisnarr(:)
      sumphidis(:)=phidis(:)
      sumthetadis(:)=thetadis(:)
#endif          

#ifdef PARALLEL
      if (myid==0) then
#endif       
        write(avrfile(4:11),'(I8.8)') NT       
        open(2,file=avrfile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'average displacement of reflected FBM'
        write(2,*) 'step distribution ',STEPDIS
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'WALL=', WALL
        write(2,*) 'F0=',F0,' LAM=',LAM 
        write(2,*) 'RMAX=', RMAX
        write(2,*) 'THETAM=', THETAM
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=', totconf
        write (2,*)'=================================='
        write(2,*) '   time         <r>         <r^2>'
        do itout=1,itoutmax
          Write(2,'(1X,I8,6(2X,E13.6))')  outtime(itout), sumrr(itout)/totconf, sum2rr(itout)/totconf 
        enddo 
        close(2) 

      if (WRITEDISTRIB) then
        write(disfile(4:11),'(I8.8)') NT 

        open(2,file=disfile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'Radial distribution between NTSTART and NTEND'
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'WALL=', WALL
        write(2,*) 'F0=',F0,' LAM=',LAM 
        write(2,*) 'NTSTART= ', NTSTART
        write(2,*) 'NTEND= ', NTEND
        write(2,*) 'RMAX=', RMAX
        write(2,*) 'THETAM=', THETAM
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=',totconf
        write (2,*)'=================================='
        write(2,*) 'ibin  r=(RMAX*ibin)/NBIN  r/RMAX  2*pi*r**2*(1-cos(THETAM))*P(r)  P(r)  P(r)*RMAX**3  RMAX-r  P(r)_err'
        do ibin=1,NBIN 
          rad= (RMAX*ibin)/NBIN
          PP= (sumdis(ibin)*NBIN)/(RMAX*totconf*(NTEND-NTSTART+1))
          volume=(RMAX*(ibin+0.5D0)/NBIN)**3  -  (RMAX*(ibin-0.5D0)/NBIN)**3  
          volume=(2.D0/3.D0)*Pi*(1.D0-cos(THETAM))*volume
          PPvol= sumdis(ibin)/(volume*totconf)/(NTEND-NTSTART+1)
          PPvolerr= sqrt(sumdis(ibin))/(volume*totconf)/(NTEND-NTSTART+1)
          Write(2,'(1X,I7,9(2X,E13.6))') ibin, rad, rad/RMAX, PP, PPvol, PPvol*RMAX**3, RMAX-rad,PPvolerr, sumdis(ibin)
        enddo 
        close(2) 

        write(narfile(4:11),'(I8.8)') NT 

        open(2,file=narfile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'Narrow Radial distribution between NTSTART and NTEND'
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'WALL=', WALL
        write(2,*) 'F0=',F0,' LAM=',LAM 
        write(2,*) 'NTSTART= ', NTSTART
        write(2,*) 'NTEND= ', NTEND
        write(2,*) 'RMAX=', RMAX
		write(2,*) 'NARRFAC=', NARRFAC
        write(2,*) 'THETAM=', THETAM
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=',totconf
        write (2,*)'=================================='
        write(2,*) 'ibin  r=(RMAX*ibin)/NBIN  r/RMAX  2*pi*r**2*(1-cos(THETAM))*P(r)  P(r)  P(r)*RMAX**3  RMAX-r  P(r)_err'
        do ibin=1,NBIN 
          rad= (RMAX*ibin)/NBIN
          PP= (sumdisnarr(ibin)*NBIN)/(RMAX*totconf*(NTEND-NTSTART+1))
          volume=(RMAX*(ibin+0.5D0)/NBIN)**3  -  (RMAX*(ibin-0.5D0)/NBIN)**3  
          volume=(2.D0/3.D0)*Pi*(1.D0-cos(THETAM/NARRFAC))*volume
          PPvol= sumdisnarr(ibin)/(volume*totconf)/(NTEND-NTSTART+1)
          PPvolerr= sqrt(sumdisnarr(ibin))/(volume*totconf)/(NTEND-NTSTART+1)
          Write(2,'(1X,I7,9(2X,E13.6))') ibin, rad, rad/RMAX, PP, PPvol, PPvol*RMAX**3, RMAX-rad,PPvolerr, sumdisnarr(ibin)
        enddo 
        close(2) 

      endif

      if (WRITEPHIDISTRIB) then
        write(phifile(4:11),'(I8.8)') NT 

        open(2,file=phifile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'Phi distribution between NTSTART and NTEND'
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'WALL=', WALL
        write(2,*) 'F0=',F0,' LAM=',LAM 
        write(2,*) 'NTSTART= ', NTSTART
        write(2,*) 'NTEND= ', NTEND
        write(2,*) 'RMAX=', RMAX
        write(2,*) 'THETAM=', THETAM
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=',totconf
        write (2,*)'=================================='
        write(2,*) '   ibin    phi=(Pi*ibin)/NPHIBIN  P(phi)  P(phi)_err'
        do ibin=-NPHIBIN,NPHIBIN 
          phi= (Pi*ibin)/NPHIBIN
          PP= (sumphidis(ibin)*NPHIBIN)/(Pi*totconf)/(NTEND-NTSTART+1)
          PPphierr = (sqrt(sumphidis(ibin))*NPHIBIN)/(Pi*totconf)/(NTEND-NTSTART+1)
          Write(2,'(1X,I7,8(2X,E13.6))') ibin, phi, PP, PPphierr
        enddo 
        close(2) 
      endif

      if (WRITETHETADISTRIB) then
        write(thetafile(4:11),'(I8.8)') NT 

        open(2,file=thetafile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'Theta distribution between NTSTART and NTEND'
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'WALL=', WALL
        write(2,*) 'F0=',F0,' LAM=',LAM 
        write(2,*) 'NTSTART= ', NTSTART
        write(2,*) 'NTEND= ', NTEND
        write(2,*) 'RMAX=', RMAX
        write(2,*) 'THETAM=', THETAM
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=',totconf
        write (2,*)'=================================='
        write(2,*) '   ibin    theta=(Pi/2*ibin)/NTHETABIN  P(theta)  P(theta)_err'
        do ibin=0,NTHETABIN 
          theta= (PI*ibin)/NTHETABIN
          PP= (sumthetadis(ibin)*NTHETABIN)/(PI*totconf)/(NTEND-NTSTART+1)/sin(theta)
          PPthetaerr = (sqrt(sumthetadis(ibin))*NTHETABIN)/(PI*totconf)/(NTEND-NTSTART+1)/sin(theta)
          Write(2,'(1X,I7,8(2X,E13.6))') ibin, theta, PP, PPthetaerr
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
      END PROGRAM fbm_cone2

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine CORVEC(xr,Ns,MM)
!
! allocatable array version (prevents seg faults)
! generates a 1d array of Ns Gaussian random numbers xr
! correlated according to a (translationally invariant)
! user-supplied correlation function corfunc(is,Ns)
! 
! uses Fourier filtering method
!
! history
!      v0.9         Dec  7, 2013:        first version, uses Tao Pang FFT
!      v0.91        Oct 11, 2017:        uses much faster FFT by Ooura (in fftsg.f)        
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