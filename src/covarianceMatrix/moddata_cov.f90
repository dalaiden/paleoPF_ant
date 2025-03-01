PROGRAM moddata_cov
  use netcdf
  IMPLICIT NONE
  !This program reads model and data files, process them, compute a cost function 

  !Interface for pointer manipulation in these to function
  interface
    subroutine readvarcdf(filename, varName, pvalue, pdimName, pdimSize, dimtolimit, limit)
      REAL, DIMENSION(:), POINTER :: pvalue
      CHARACTER(len=*) :: varName, filename
      CHARACTER(len=256), DIMENSION(:), POINTER :: pdimName
      INTEGER, DIMENSION(:), POINTER :: pdimSize
      INTEGER, DIMENSION(:), OPTIONAL :: limit
      CHARACTER(len=*), OPTIONAL :: dimtolimit
    end subroutine readvarcdf
    subroutine average_mod(undef, inputtab, outputtab, timedimidx, first, last)
      REAL :: undef
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: inputtab
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: outputtab
      INTEGER :: timedimidx,first,last
    end subroutine average_mod
    subroutine average_ref(undef, inputtab, outputtab, timedimidx, monthF, monthL)
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: inputtab
      REAL, DIMENSION(:,:), ALLOCATABLE :: outputtab
      INTEGER :: timedimidx,monthF,monthL
      REAL :: undef
    end subroutine average_ref
    subroutine covariance(cov,state_array)
      REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: cov !(nlon*nlat,nlon*nlat) 
      REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: state_array !(nlon*nlat,ntmaxs)
    end subroutine covariance
    subroutine write_cov(filename,data_out,data_name)
      CHARACTER(len=256) :: filename
      REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: data_out !(NX,NY)
      CHARACTER(len=*) :: data_name
    end subroutine write_cov
  end interface

  !Definition of variable to use these previous function
  REAL, DIMENSION(:), POINTER :: pvalue
  CHARACTER(len=256), DIMENSION(:), POINTER :: pdimName
  INTEGER, DIMENSION(:), POINTER :: pdimSize


  !DEFINITION OF VARIABLES
  INTEGER :: nlat !=32 number of latitudes
  INTEGER :: nlon !=64 number of longitudes
  INTEGER :: ntmax !=2100  ! number of times
  INTEGER, PARAMETER :: ntmaxs=150  ! number of times
  INTEGER :: itime, ilat, ilon
  REAL(kind=4) :: undef, maxy

  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: modsmooth, moddsigF !(nlon,nlat,ntmaxs)
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: ref3D, mod3D !(nlon,nlat,ntmaxs)
  REAL(kind=4), DIMENSION(:,:,:,:), ALLOCATABLE :: ref4D, mod4D !(nlon,nlat,P_U3,ntmaxs)
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE ::   meanmod !(nlon,nlat)

  !   REAL(kind=4) ::                                       wTs,         wPp,         wGeopg,         wSss,         wSst
  CHARACTER(len=256) ::                modvar!,          modTs,       modPp,       modGeopg,       modSss,       modSst
  CHARACTER(len=256) ::                refmodvar!,       refmodTs,    refmodPp,    refmodGeopg,    refmodSss,    refmodSst
  CHARACTER(len=256) ::                covvar!,          covTs,       covPp,       covGeopg,       covSss,       covSst

  CHARACTER(len=256) :: varName

  INTEGER :: i,j,ii,jj,k,l,gridtype,netcdfDimension,ncid

  INTEGER :: monthF, monthL
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: cov !(nlon*nlat,nlon*nlat)
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: state_array !(nlon*nlat,ntmaxs)
  CHARACTER(len=10) :: methode


  integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in

  ! GET ARGUMENT FOR SELECT METHODE
  call getarg(1,methode)   

  ! READ MODDATA_CO.PARAM
  open(61,file='moddata_cov.param')
  !remove header
  read(61,*)
  read(61,*)
  !Variable weight
  read(61,*)
  read(61,*) varName
  write(*,*) 'Variable name = ', trim(varName)

  ! First and last month
  read(61,*)
  read(61,*) monthF
  read(61,*) monthL
  write(*,*) 'Average over a period; the first and the last months of the period are :'
  write(*,*) monthF, monthL


  ! Name of the data file
  read(61,*)
  read(61,'(A256)') modvar
  read(61,'(A256)') refmodvar
  read(61,'(A256)') covvar
  write(*,*) 'modTs    = ', trim(modvar)
  write(*,*) 'refmodTs    = ', trim(refmodvar)
  write(*,*) 'covTs    = ', trim(covvar)

  if ((monthF.lt.1).or.(monthF.gt.12).or.(monthL.lt.1).or.(monthL.gt.12).or.(monthF.gt.monthL)) then
    write(*,*)'Error: check the months over which obs will be assimilated!'
    stop
  endif

  ! CHECK THE VARIABLE CONSISTENCY
  ! Check if the variable is in the ocean or in the atmosphere
  ! gridtype=0 -> ATMOSPHERE    gridtype=1 -> OCEAN
  gridtype=10000;
  if((varName.eq.'sst').or.(varName.eq.'sss').or.(varName.eq.'sim').or.(varName.eq.'albq')) then
    gridtype=1
  else if((varName.eq.'ts').or.(varName.eq.'pp').or.(varName.eq.'geopg')) then
    gridtype=0
  endif

  if(gridtype.eq.0) then
    write(*,*) 'Variable ',trim(varName),' is on continental grid'
  else if(gridtype.eq.1) then
    write(*,*) 'Variable ',trim(varName),' is on oceanic grid'
  else if((gridtype.gt.1).or.(gridtype.lt.0)) then
    write(*,*) 'Variable ',trim(varName),' is unknown by this software'
    stop
  endif

  !DEFINITION OF CONSTANTS AND OF GRID CHARACTERISTICS
  if(gridtype.eq.0) then
    undef=-1E+32
    nlat=32
    nlon=64
  else
    undef=-1E+32
    nlat=65
    nlon=120
  endif
  write(*,*) "DEFINITION OF CONSTANTS AND OF GRID CHARACTERISTICS"
  write(*,*) "nlat =", nlat, "nlon =", nlon, "undef =", undef

  ! Get the dimension of the netcdf to know if 3D or 4D variables
  call check( nf90_open(trim(refmodvar), NF90_NOWRITE, ncid) )
  call check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )
  netcdfDimension=ndims_in
  write(*,*) 'netcdf Dimension = ', netcdfDimension,'D'
  call check( nf90_close(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!     READ OBS REF       !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) " "
  write(*,*) 'Reading ref file for the mod '
  write(*,*) '- - - - - - - - - - - - - - -'
  if(len(trim(refmodvar))==0) then 
    write(*,*) "No ref for ts"
    allocate(meanmod(nlon,nlat))
    meanmod=0
  else
    write(*,*) trim(refmodvar)
    call readvarcdf(refmodvar, trim(varName), pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/1,12/))

    if(netcdfDimension.eq.3) then! FOR 3D VARIABLE
      allocate(ref3D(pdimSize(1),pdimSize(2),pdimSize(3)))
      ref3D=reshape(pvalue, (/pdimSize(1),pdimSize(2),pdimSize(3)/))
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
    endif
    if(varName.eq.'geopg') then! FOR 4D VARIABLE
      allocate(ref4D(pdimSize(1),pdimSize(2),pdimSize(3),pdimSize(4)))
      allocate(ref3D(pdimSize(1),pdimSize(2),pdimSize(4)))
      ref4D=reshape(pvalue, (/pdimSize(1),pdimSize(2),pdimSize(3),pdimSize(4)/))
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      ref3D(:,:,:)=ref4D(:,:,3,:)
    endif
    write(*,*)'Average from', monthF, 'month to', monthL, 'month'
    allocate(meanmod(size(ref3D,1),size(ref3D,2)))
    call average_ref(undef,ref3D,meanmod,3,monthF,monthL)
  endif
  write(*,*) "meanmod INFORMATION"
  write(*,*) "meanmod size", size(meanmod,1), size(meanmod,2)
  write(*,*) "meanmod(1,1)=",meanmod(1,1)




  if(allocated(ref3D)) deallocate(ref3D)
  if(allocated(ref4D)) deallocate(ref4D)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!     READ OBS       !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) " "
  !for 3d variable (lon,lat,t)
  write(*,*) 'reading of the mod file ', trim(varName)
  write(*,*) '- - - - - - - - - - - - - - -'
  write(*,*) trim(modvar)
  call readvarcdf(modvar, trim(varName), pvalue, pdimName, pdimSize) !limit = (start, count)
  if(netcdfDimension.eq.3) then! FOR 3D VARIABLE
    allocate(mod3D(pdimSize(1),pdimSize(2),pdimSize(3)))
    allocate(modsmooth(pdimSize(1),pdimSize(2),pdimSize(3)/12))
    do j=0, pdimSize(1)-1
      do k=0, pdimSize(2)-1
        do l=0, pdimSize(3)-1
          mod3D(j+1,k+1,l+1)=pvalue(j+pdimSize(1)*k+pdimSize(1)*pdimSize(2)*l+1)
        enddo
      enddo
    enddo
  endif
  if(varName.eq.'geopg') then! FOR 4D VARIABLE
    allocate(mod4D(pdimSize(1),pdimSize(2),pdimSize(3),pdimSize(4)))
    allocate(mod3D(pdimSize(1),pdimSize(2),pdimSize(4)))
    allocate(modsmooth(pdimSize(1),pdimSize(2),pdimSize(4)/12))
    do j=0, pdimSize(1)-1
      do k=0, pdimSize(2)-1
        do l=0, pdimSize(3)-1
          do ii=0, pdimSize(4)-1
            mod4D(j+1,k+1,l+1,ii+1)=pvalue(j+pdimSize(1)*k+pdimSize(1)*pdimSize(2)*l+pdimSize(1)*pdimSize(2)*pdimSize(3)*ii+1)
          enddo
        enddo
      enddo
    enddo
    mod3D(:,:,:)=mod4D(:,:,3,:)
  endif
  write(*,*)'Size of mod3D',pdimSize(:)


  write(*,*)'Average from', monthF, 'month to', monthL, 'month'
  call average_mod(undef,mod3D,modsmooth,3,monthF,monthL) !3 is the time index in pdimSize

  deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize); deallocate(mod3D)
  if(size(mod4D) .ne. 0) then ;deallocate(mod4D); endif



    write(*,*) "modsmooth INFORMATION"
    write(*,*) "modsmooth size", size(meanmod,1), size(meanmod,2)
    write(*,*) "modsmooth(1,1,1)=",modsmooth(1,1,1)



    !compute the ensemble means and normalized anomalies
    allocate(moddsigF(size(modsmooth,1), size(modsmooth,2), size(modsmooth,3)))
    DO itime=1,size(modsmooth,3) !Choise_cov=size(modsmooth,3)=size(moddsig,3)
      if(varName.eq.'geopg') then
        moddsigF(:,:,itime)=(modsmooth(:,:,itime)-meanmod(:,:))/9.81
      else 
        do ilon=1, size(modsmooth,1)
          do ilat=1, size(modsmooth,2)
            if (modsmooth(ilon,ilat,itime).ne.undef .AND. meanmod(ilon,ilat).ne.undef) then
              moddsigF(ilon,ilat,itime)=modsmooth(ilon,ilat,itime)-meanmod(ilon,ilat)
            else
              moddsigF(ilon,ilat,itime)=undef
            endif
          enddo
        enddo
      endif
    ENDDO
    ntmax=size(moddsigF,3)

    write(*,*) " "
    if(gridtype.eq.0) then
      write(*,*) "              (1,32)         (57,24)"
      write(*,*) 'modsmooth ',modsmooth(1,32,1),'   ',modsmooth(57,24,1)
      write(*,*) 'moddsigF  ',moddsigF(1,32,1), '   ',moddsigF(57,24,1)
    else if(gridtype.eq.1) then
      write(*,*) "              (1,32)         (57,24)       (120,52)"
      write(*,*) 'modsmooth ',modsmooth(1,32,1),'   ',modsmooth(57,24,1),modsmooth(120,52,1)
      write(*,*) 'moddsigF  ',moddsigF(1,32,1), '   ',moddsigF(57,24,1),moddsigF(120,52,1)
    endif
    write(*,*) " "

    !COMPUTATION OF THE COVARAINCE MATRIX 
    write(*,*) " "
    write(*,*)'COMPUTATION OF THE COVARAINCE MATRIX '
    write(*,*)'- - - - - - - - - - - -- - - - - - - '
    !Define the dim of cov matrix
    write(*,*) 'moddsigF(1,1,1)=',moddsigF(1,1,1)
    write(*,*) 'undef=',undef
    i=1
    DO jj=1,nlat
      DO ii=1,nlon
        if(moddsigF(ii,jj,1).ne.undef) i=i+1
      ENDDO
    ENDDO
    allocate(state_array(i-1,size(moddsigF,3)))
    allocate(cov(i-1,i-1))

    write(*,*)'Size of the cov matrix',size(cov,1)
    write(*,*)'Dim in time',size(state_array,2)
    write(*,*)'Dim in time should be equal or greater than size of the cov matrix for nonsingularity'



    if(gridtype.eq.0) then
      !Fill in the state array
      state_array(:,:)=0.
      cov(:,:)=0.
      do j=1,size(state_array,2)
        i=1
        DO jj=1,nlat
          DO ii=1,nlon
            if(moddsigF(ii,jj,j).ne.undef) then 
              state_array(i,j)=moddsigF(ii,jj,j)
              i=i+1
            endif
          ENDDO
        ENDDO
      enddo
    else if(gridtype.eq.1) then
      state_array(:,:)=0.
      cov(:,:)=0.
      undef=-undef
      maxy=0.
      do j=1,size(state_array,2)
        i=1
        DO jj=1,nlat
          DO ii=1,nlon
            if(abs(moddsigF(ii,jj,j)).ne.undef) then 
              state_array(i,j)=moddsigF(ii,jj,j)
              if(abs(state_array(i,j)).gt.maxy) maxy=state_array(i,j)
              i=i+1
            endif
          ENDDO
        ENDDO
      enddo
      write(*,*)'Max in st. array',maxy
    endif

    CALL covariance(cov,state_array)
    deallocate(state_array)  

    write(*,*) " "
    write(*,*) 'cov ',cov(4,4),'   ',cov(100,10)
    write(*,*) 'cov ',cov(10,100),'   ',cov(2000,2000)
    write(*,*) " "

    write(*,*)'The covariance matrix was written to'
    write(*,*) covvar
    CALL write_cov(covvar,cov,"cov")

    deallocate(cov)
    deallocate(modsmooth)
    deallocate(moddsigF)

  END PROGRAM moddata_cov

