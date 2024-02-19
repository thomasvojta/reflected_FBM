      PROGRAM fbm5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   time-discrete reflected fractional Brownian motion
!
!   fbm1   10 Sep 2017   first fbm version 
!   fbm2   10 Sep 2017   combined parallel vs serial version using preprocessor 
!   fbm3   10 Sep 2017   also write scaled distrib 
!   fbm4   10 Sep 2017   longer filename to allow longer times 
!   fbm5   11 Oct 2017   uses new corvec with Ooura's FFT (much faster!) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #define PARALLEL
#define VERSION 'fbm5'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      integer(i4b), parameter     :: M=14,NT=2**M               ! number of time steps (in which mu const) 
      integer(i4b), parameter     :: NCONF=1000                 ! number of disorder configurations
      real(r8b), parameter        :: V0=0.D0                       ! drift velocity

      character(3)                :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   

      logical,parameter           :: WRITEDISTRIB = .TRUE.        ! write final xx distribution    
      integer(i4b), parameter     :: NBIN =  1000                   ! number of bins for distribution
      real(r8b), parameter        :: MAXXX = 5.D3              ! largest - xx for histogram  

      real(r8b), parameter        :: outtimefac=2**0.25D0          ! factor for consecutive output times 

      integer(i4b), parameter     :: IRINIT=1                    ! random number seed
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(r8b)              :: y(1:NT)                           ! correlated random numbers
      real(r8b)              :: xx                               ! position
      real(r8b)              :: confxx(1:NT)                     ! average of xx after each time step  
      real(r8b)              :: conf2xx(1:NT)                    
      real(r8b)              :: sumxx(1:NT),sum2xx(1:NT)         ! sums over machines in MPI version
      real(r8b)              :: auxxx(1:NT),aux2xx(1:NT) 

      integer(i4b)           :: iconf, it, ibin                 ! configuration, and time counters   
      integer(i4b)           :: totconf                         ! actual number of confs
      
      real(r8b)              :: xxdis(-NBIN:NBIN)                    ! density histogram 
      real(r8b)              :: sumdis(-NBIN:NBIN)
      real(r8b)              :: auxdis(-NBIN:NBIN)
      real(r8b)              :: PP, scalefac                     ! P(x)

      external               :: kissinit 
      real(r8b),external     :: rkiss05,erfcc 

      character(13)          :: rhofile = 'avx000000.dat'
      character(13)          :: disfile = 'dis000000.dat'
      
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

      confxx(:)=0.D0 
      conf2xx(:)=0.D0 
      xxdis(:)=0.D0

! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,totconf,numprocs
!      if (myid==0) print *, 'dis. conf.', iconf
#else
      disorder_loop: do iconf=1,totconf
!         print *, 'dis. conf.', iconf
#endif 

        call gkissinit(IRINIT+iconf-1)
        call corvec(y,NT,M)                         ! creates correlated Gaussian random numbers
        
        if (STEPDIS.eq.'BOX') then
!          y(:) = 1 - (0.5D0*erfcc(y(:)/sqrt(2.0D0)))                       ! map onto unit inteval
!          y(:)=y(:)-0.5D0                                                  ! center around 0
        endif
        
        if (STEPDIS.eq.'BIN') then 
          y(:)=sign(1.D0,y(:))
        endif  
        
! Time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           xx=0.D0  
           do it=1, NT/2
              xx=xx+y(it)
              if (xx.lt.0) xx=-xx                            ! reflecting wall   
              confxx(it)=confxx(it) + xx   
              conf2xx(it)=conf2xx(it) + xx*xx   
           end do  

           if (WRITEDISTRIB) then
              ibin=nint( xx*NBIN/MAXXX )
              if ( (ibin.ge.-NBIN) .and. (ibin.le.NBIN)) then
                 xxdis(ibin)=xxdis(ibin)+1.D0 
              endif    
!              print *, xx,ibin
           endif

      end do disorder_loop      ! of do incof=1,NCONF


! Now collect and analyze data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
     if (myid.ne.0) then                                                  ! Send data
         call MPI_SEND(confxx,NT,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
         call MPI_SEND(conf2xx,NT,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,ierr)
         call MPI_SEND(xxdis,2*NBIN+1,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
      else
         sumxx(:)=confxx(:)
         sum2xx(:)=conf2xx(:)
         sumdis(:)=xxdis(:)
         do id=1,numprocs-1                                                   ! Receive data
            call MPI_RECV(auxxx,NT,MPI_DOUBLE_PRECISION,id,1,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(aux2xx,NT,MPI_DOUBLE_PRECISION,id,2,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(auxdis,2*NBIN+1,MPI_DOUBLE_PRECISION,id,3,MPI_COMM_WORLD,status,ierr)
            sumxx(:)=sumxx(:)+auxxx(:) 
            sum2xx(:)=sum2xx(:)+aux2xx(:) 
            sumdis(:)=sumdis(:)+auxdis(:)
         enddo
      endif        
#else          
      sumxx(:)=confxx(:)
      sum2xx(:)=conf2xx(:)
      sumdis(:)=xxdis(:)
#endif
         
#ifdef PARALLEL
      if (myid==0) then
#endif       
        write(rhofile(4:9),'(I6.6)') NT/2       
        open(2,file=rhofile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'average displacement of reflected FBM'
        write(2,*) 'step distribution ',STEPDIS
        write(2,*) 'NT= ', NT
        write(2,*) 'V0=', V0
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=', totconf
        write (2,*)'=================================='
        write(2,*) '   time         <x>         <x^2>'
        it=1
        do while(it.le.NT/2)  
          Write(2,'(1X,I7,6(2X,E13.6))')  it, sumxx(it)/totconf, sum2xx(it)/totconf 
          it=max(it+1,nint(outtimefac*it))
        enddo 
        close(2) 

      if (WRITEDISTRIB) then
        write(disfile(4:9),'(I6.6)') NT/2 

        open(2,file=disfile,status='replace')
        write(2,*) 'Program ', VERSION
        write(2,*) 'Density distribution at time', NT/2 
        write(2,*) 'step distribution ',STEPDIS        
        write(2,*) 'NT= ', NT
        write(2,*) 'V0=', V0
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=',totconf
        write (2,*)'=================================='
        write(2,*) '   ibin    x=(MAXX*ibin)/NBIN   P(x)   x/sigma   P*sigma'
        do ibin=0,NBIN 
          xx= (MAXXX*ibin)/NBIN
          PP= (sumdis(ibin)*NBIN)/(MAXXX*totconf)
          scalefac=sqrt(sum2xx(NT/2)/totconf)
          Write(2,'(1X,I7,5(2X,E13.6))') ibin, xx, PP, xx/scalefac, PP*scalefac   
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
      END PROGRAM fbm5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine CORVEC(xr,Ns,M)
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE corvec(xr,Ns,M)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

      integer(i4b)           :: Ns              ! number of sites, must be power of 2 
      integer(i4b)           :: M               ! Ns=2^M

      real(r8b)              :: xr(0:Ns-1)      ! random number array  
      real(r8b)              :: cr(0:Ns-1)      ! correlation function 
      integer(i4b)           :: is
      
      integer(i4b)           :: ip(0:2+int(sqrt(1.*Ns)))   ! workspace for FFT code
      real(r8b)              :: w(0:Ns/2-1)           ! workspace for FFT code 
      
      real(r8b), external    :: gkiss05,erfcc
      external               :: rdft                   ! from Ooura's FFT package
      real(r8b),external     :: corfunc 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (Ns.ne.2**M) STOP 'Size indices do not match'
! Calculate correlation function 
      do is=0,Ns-1 
         cr(is)= corfunc(is,Ns) 
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
      
      return
      END SUBROUTINE corvec 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Desired correlation function of the random numbers in real space
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION corfunc(i,N)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)              :: corfunc
      real(r8b)              :: gam
      integer(i4b)           :: i,N
      integer(i4b)           :: dist

      gam = 0.8D0    ! Gamma controlling FBM

      dist=min(i,N-i)
      if (dist.eq.0) then 
         corfunc=1.D0
      else   
         corfunc = 0.5D0*( (dist+1)**(2.D0-gam) - 2.D0*(dist**(2.D0-gam)) + (dist-1)**(2.D0-gam) ) 
      endif   

      return
      END FUNCTION corfunc


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