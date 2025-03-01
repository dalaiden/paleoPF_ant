subroutine cost_euclide(fcost,area,weight,undef,MeanDiffModObs,methode,cov,Ci,CiUndef,sigma,NbStdDev,ln_geoweight,Tmask)
  IMPLICIT NONE

  interface
    subroutine inv_matrix(a)
      REAL, DIMENSION(:,:), ALLOCATABLE :: a !(n,n) 
    end subroutine inv_matrix
  end interface

  INTEGER nlon,nlat,npoint,indx ! grid size (lon, lat, valid point, point larger than bar)
  INTEGER :: ii, jj, i          ! loop index
  INTEGER methode               ! type of methode FAP BEST
  INTEGER, DIMENSION(:), ALLOCATABLE :: s, indy !(nlon*nlat)

  REAL :: fcost, undef, sigma, CiUndef, NbStdDev
  REAL :: weight              ! variable weight
  DOUBLE PRECISION :: radian,sqrt2pi,d2pi,bar
  DOUBLE PRECISION :: accu1,accu2         ! accumulator for fcost
  DOUBLE PRECISION :: dmax,dmin           ! max/min dist to obs 

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rgeoweight
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: area !(nlon,nlat)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rgeoweight_indy
  REAL, DIMENSION(:)    , ALLOCATABLE :: Ci !(0:nlat)
  REAL, DIMENSION(:)    , ALLOCATABLE :: state_array !(nlon*nlat)
  REAL, DIMENSION(:)    , ALLOCATABLE :: dummy, rms !(nlon*nlat)
  REAL, DIMENSION(:,:)  , ALLOCATABLE :: covT !(nlon*nlat,nlon*nlat)
  REAL, DIMENSION(:,:)  , ALLOCATABLE :: cov_inv !(nlon*nlat,nlon*nlat)
  REAL, DIMENSION(:,:)  , ALLOCATABLE :: MeanDiffModObs !(nlon,nlat)
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: cov !(nlon*nlat,nlon*nlat)

  LOGICAL :: ln_geoweight                      ! flag for geoweight (1 if false) 

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tmask
  INTEGER :: nmask, nTotmask! (nTotmask = number of grid cell which isn't mask)

  nlon=size(MeanDiffModObs,1)
  nlat=size(MeanDiffModObs,2)
  radian=2*acos(0.0)/180.0
  sqrt2pi=sqrt(4*acos(0.0))
  d2pi=4*acos(0.0)


  !START VALUE FOR ACCUMULATOR
  SELECT CASE(methode)
  CASE(1)
    write(*,*)"    |Methode Likelihood"
    accu1=0
    accu2=0

  CASE(2)
    write(*,*)"    |Methode old euclide"
    accu1=0
    accu2=0

  CASE default
    write(*,*)"    |Unknow Methode ",methode
    STOP 0
  END SELECT
  ! --- geographic weight ---
  allocate(rgeoweight(size(area,1),size(area,2)))
  allocate(rgeoweight_indy(size(area,1)*size(area,2)))
  IF ( .NOT. ln_geoweight ) then
    rgeoweight=1.
  else
    rgeoweight=area
  endif
  ! ---
  SELECT CASE(methode)
  CASE(1)
    ! Reduce the full cov. matrix to a box 
    allocate(state_array(nlon*nlat))    !npoint (nombre de grille avec au - 1 proxy)
    allocate(indy(nlon*nlat))           !idem
    allocate(s(nlon*nlat))              !idem
    allocate(covT(nlon*nlat,nlon*nlat)) !(sum(Tmask)**2
    allocate(rms(nlon*nlat))            !npoint

    !Test if size cov = sum(Tmask)
    nTotmask=sum(Tmask)


    if(size(cov)==0) then  ! FK
      write(*,*) "No covariance matrice loaded in cost_euclide.f90"
    else 
      write(*,*) "Covariance matrix loaded in cost_euclide.f90"
      IF( (size(cov,1).ne.size(cov,2)) .or. (size(cov,1).ne.nTotmask) ) THEN
        write(*,*) "    |W.A.R.N.I.N.G  -  Wrong size of cov matrix", size(cov,1), "x", size(cov,2)
        write(*,*) "    |It should be ((nlat*nlon)-nTotmask)", nTotmask, "x", nTotmask
        stop
      ELSE
        write(*,*) "    |size cov = nlat*nlon - nTotmask   [OK]"
      ENDIF
    endif

    rms(:)=0.
    i=1
    nmask=0
    DO jj=1,nlat  
      DO ii=1,nlon
        ! A.M. Count number of grid that is continent (nmask) in the case of ocean.
        ! If it is not the ocean case, nmask = 0 because mask are fill with 1.
        ! If Tmask(ii,jj) = 0 then continent. We see if < 0.5 because of fortran issue.
        if(Tmask(ii,jj)<0.5) nmask=nmask+1;
        if((Tmask(ii,jj)<0.5).AND.(MeanDiffModObs(ii,jj).ne.undef)) then
          write(*,*) "    |W.A.R.N.I.N.G - Proxy ", "lat", jj, "lon", ii, "is located within the mask",Tmask(ii,jj),nmask
        endif
        if((Tmask(ii,jj)>0.5).AND.(MeanDiffModObs(ii,jj).ne.undef).AND.(weight*rgeoweight(ii,jj)).ne.0.0) then
          state_array(i)=MeanDiffModObs(ii,jj)
          indy(i)=jj
          rgeoweight_indy(i)=rgeoweight(ii,jj)
          s(i)=(jj-1)*nlon+ii-nmask

          if(size(cov) .NE. 0) covT(i,1:size(cov,2))=cov(s(i),:,1)

          rms(i)=Ci(s(i))
          if (rms(i).eq.undef) rms(i)=CiUndef
          PRINT *, ii, jj, state_array(i), rms(i)
          i=i+1
        endif
      ENDDO      
    ENDDO

    ! --- number of valid point (not undef) ---
    npoint=i-1
    write(*,*) npoint, "proxy available"

    if(size(cov)==0) then  ! FK No cov

      allocate(dummy(npoint))    
      write(*,*) size(dummy,1)       
      DO ii=1,npoint
        write(*,*) ii
        dummy(ii) = state_array(ii) / sqrt(rms(ii))
      ENDDO

    else ! FK Cov
      ! ---
      allocate(cov_inv(npoint,npoint))
      DO ii=1,npoint
        cov_inv(:,ii)=covT(1:npoint,s(ii))*sigma
        !write(*,*) sigma, rms(ii)
      ENDDO
      deallocate(s); deallocate(covT);  

      ! Add measurement error
      DO ii=1,npoint
        cov_inv(ii,ii)=cov_inv(ii,ii)+rms(ii)
      ENDDO  

      write(*,*)'    |Calculating inverse of the cov matrix, size',npoint


      write(*,*) size(cov_inv,1)

      call inv_matrix(cov_inv)

      ! Multiply inverse covariance matrix by a state vector
      allocate(dummy(npoint))    
      DO ii=1,npoint
        dummy(ii)=0.
        DO jj=1,npoint
          dummy(ii)=dummy(ii)+cov_inv(ii,jj)*state_array(jj)
        ENDDO
        !PRINT *, ii, dummy(ii)
        !write(*,*) "dummy=", state_array(ii)
      ENDDO      
    endif

    ! --- Calculate max and min distance between model and observations ---
    dmax=abs(state_array(1))
    dmin=abs(state_array(1))

    DO ii=2,npoint
      if (dmax.lt.abs(state_array(ii))) dmax=abs(state_array(ii))
      if (dmin.gt.abs(state_array(ii))) dmin=abs(state_array(ii))
    ENDDO
    write(*,*) '    ||mod-obs|_max  = ', dmax 
    write(*,*) '    ||mod-obs|_min  = ', dmin 
    ! ---
    ! --- compute term (d-H)TC-1(d-H) (eq. 4, Dubinkina 2010) ---
    indx=0
    bar=NbStdDev*sqrt(maxval(rms))
    DO ii=1,npoint
      !ACCUMULATION
      dmax=min(abs(dummy(ii)),bar)
      accu1=accu1+weight*rgeoweight_indy(ii)*dmax**2
      accu2=accu2+rgeoweight_indy(ii)
      if (dmax.ne.bar) indx=indx+1
    ENDDO
    write(*,*) '    |number of dmax smaller than barrier', indx
    write(*,*) '    |barrier = ',bar
    deallocate(state_array, dummy, rms, indy)
    ! ---  
  CASE(2)
    DO ii=1,nlon
      DO jj=1,nlat
        IF((MeanDiffModObs(ii,jj).ne.undef).AND.((weight*rgeoweight(ii,jj)).ne.0.0)) THEN
          !ACCUMULATION
          accu1=accu1+weight*rgeoweight(ii,jj)*MeanDiffModObs(ii,jj)**2
          accu2=accu2+weight*rgeoweight(ii,jj)
        ENDIF
      ENDDO
    ENDDO
  END SELECT


  !RETURN VALUE
  SELECT CASE(methode)
  CASE(1)
    ! --- compute fcost (eq. 4, Dubinkina 2010) ---
    fcost=exp(-accu1/2.)
    write(*,*) "    |accu",accu1
  CASE(2)
    ! --- compute mean rms ---
    if (accu2.ne.0.0) then
      fcost=sqrt(accu1/accu2)
    else
      fcost=10.0
    endif
  END SELECT

  if(size(cov)==0) then  ! FK No cov

    write(*,*) "NO COV"

  else ! FK Cov

    write(*,*) "COV"

  endif

END SUBROUTINE cost_euclide



!These subroutines for finding inverse
SUBROUTINE inv_matrix(a)
  IMPLICIT NONE

  REAL, DIMENSION(:,:), ALLOCATABLE :: a !(m,n)

  REAL, DIMENSION(:,:), ALLOCATABLE :: v !(n,n)
  REAL, DIMENSION(:), ALLOCATABLE :: w !(n)
  INTEGER n,i,j
  REAL TINY1,maxy,minl
  PARAMETER (TINY1=1.0e-8) !Threshhold for smallest eigenvalue

  INTEGER LWORK,INFO,LWMAX, idxFirstBig, idxLastZero
  REAL, DIMENSION(:), ALLOCATABLE :: WORK !(MAX(1,LWORK))

  n=size(a,1)
  allocate(w(n))
  allocate(v(n,n))
  v=a

  LWMAX=n
  allocate(WORK(LWMAX))

  LWORK = -1
  CALL SSYEV( 'V', 'U', n, v, n, w, WORK, LWORK, INFO )
  LWORK = INT( WORK( 1 ) )
  deallocate(WORK)
  allocate(WORK(LWORK))
  CALL SSYEV( 'V', 'U', n, v, n, w, WORK, LWORK, INFO )

  if(INFO.eq.0) write(*,*)'    |LAPACK Eigen value decomposition was done successfully'
  if(INFO.lt.0) write(*,*)'    |In LAPACK Eigen value decomposition the argument',INFO,'had an illegal value'
  if(INFO.gt.0) write(*,*)'    |LAPACK Eigen value decomposition did not converge'

  write(*,*) '    |8 smallest eigenvalues after shift'
  write(*,*) '    |', w(1:4)
  write(*,*) '    |', w(5:8)
  write(*,*) '    |8 largest eigenvalues'
  write(*,*) '    |', w(size(w,1)-7:size(w,1)-4)
  write(*,*) '    |', w(size(w,1)-3:size(w,1))

  idxLastZero=-1; idxFirstBig=-1
  !DIR$ novector
  do i=1,n
    if(w(i).lt.0.) then
      w(i)=0.
      idxLastZero=i
    endif
  enddo

  if(idxLastZero.eq.-1) then
    idxFirstBig=1
  else
    if(idxLastZero+1.le.size(w,1)) then
      idxFirstBig=idxLastZero+1
    else
      write(*,*) "    |ERROR w matrix is completely null !!"
      stop
    endif
  endif

  if(idxLastZero.gt.0) then
    write(*,*) '    |Zero eigenvalue',w(i)
    write(*,*) '    |Smallest nonzero eigenvalue',w(i+1)
    minl=w(idxFirstBig)-TINY1
    w(:)=w(:)+minl
    write(*,*) '    |8 smallest eigenvalues after shift'
    write(*,*) '    |', w(1:4)
    write(*,*) '    |', w(5:8)
    write(*,*) '    |8 largest eigenvalues'
    write(*,*) '    |', w(size(w,1)-7:size(w,1)-4)
    write(*,*) '    |', w(size(w,1)-3:size(w,1))
  endif

  maxy=0.
  DO i=1,n
    DO j=1,n
      a(i,j)=v(j,i)/sqrt(w(i))
      if(abs(a(i,j)).gt.maxy) maxy=abs(a(i,j))
    ENDDO
  ENDDO
  write(*,*)'    |Max element in square root of the inv.cov.matrix =',maxy

  deallocate(v);  deallocate(w); deallocate(WORK);


END SUBROUTINE inv_matrix


