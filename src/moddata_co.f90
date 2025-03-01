!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  THIS PROGRAM READS MODEL AND DATA FILES TO COMPUTRE A COST FUNCTION  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM moddata_co
  USE NETCDF
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!     DEFINE INTERFACE FOR POINTER MANIPULATION IN FUNCTION     !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine readvarcdf(filename, varName, pvalue, pdimName, pdimSize, missval, dimtolimit, limit)
      REAL, DIMENSION(:), POINTER :: pvalue
      CHARACTER(len=*) :: varName, filename
      CHARACTER(len=256), DIMENSION(:), POINTER :: pdimName
      INTEGER, DIMENSION(:), POINTER :: pdimSize
      INTEGER, DIMENSION(:), OPTIONAL :: limit
      CHARACTER(len=*), OPTIONAL :: dimtolimit
      REAL(kind=4), OPTIONAL :: missval
    end subroutine readvarcdf
    subroutine average(undef, inputtab, outputtab, timedimidx)
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: inputtab
      REAL, DIMENSION(:,:), ALLOCATABLE :: outputtab
      INTEGER :: timedimidx
      REAL :: undef
    end subroutine average
    subroutine averageNew(undef, obs, mod, outputtab, timedimidx)
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: mod
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: obs
      REAL, DIMENSION(:,:), ALLOCATABLE :: outputtab
      INTEGER :: timedimidx
      REAL :: undef
    end subroutine averageNew
    subroutine filter_scale(var4t,undef,y_coord,lambda,lscale,nst)
      REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: var4t !(nlon,nlat)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: y_coord !table of lat in clio format (nlon,nlat)
      INTEGER :: nst
      REAL(kind=4) :: undef, lambda, lscale
    end subroutine filter_scale
    subroutine cost_euclide(fcost,y_coord,weight,undef,MeanDiffModObs,methode,cov,Ci,CiUndef,sigma,NbStdDev,ln_geoweight,Tmask)
      REAL(kind=4) :: fcost, undef, sigma, CiUndef, NbStdDev
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: y_coord !(nlon,nlat)
      REAL(kind=4), DIMENSION(:), ALLOCATABLE :: Ci !(0:nlat*nlon)
      REAL(kind=4) :: weight
      REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: MeanDiffModObs !(nlon,nlat)
      REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: cov !(nlon*nlat,nlon*nlat)
      INTEGER :: methode
      LOGICAL :: ln_geoweight
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tmask
    end subroutine cost_euclide
    subroutine write_cov(filename,data_out,data_name)
      CHARACTER(len=256) :: filename
      REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: data_out !(NX,NY)
      CHARACTER(len=*) :: data_name
    end subroutine write_cov
    subroutine read_cov(filename,data_in,data_name, t)
      CHARACTER(len=256) :: filename
      REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: data_in !(NX,NY)
      CHARACTER(len=*) :: data_name
      INTEGER :: t
    end subroutine read_cov
    subroutine ExtractUsefullDimension(DimAsk, DimInSize, DataIn, DataOut)
      character(len=1), dimension(4) :: DimAsk
      INTEGER, DIMENSION(:), POINTER :: DimInSize
      REAL, DIMENSION(:), POINTER :: DataIn
      REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: DataOut
    end subroutine ExtractUsefullDimension
    subroutine MakeBox(tlon, tlat, lon_min, lon_max, lat_min, lat_max, Mask)
      double precision ::  lon_min, lon_max, lat_min, lat_max
      integer, DIMENSION(:,:), ALLOCATABLE :: Mask
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tlon, tlat
    end subroutine MakeBox
    subroutine AverageBoxAllocateMod(data, Tlat, latproxy, lonproxy, latmin, latmax, lonmin, lonmax)
      REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: data !(nlon,nlat,ntmaxs)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tlat
      integer, DIMENSION(:), ALLOCATABLE :: latproxy, lonproxy, latmin, latmax, lonmin, lonmax
    end subroutine AverageBoxAllocateMod
    subroutine AverageBoxMask(data, Tlat, latproxy, lonproxy, IndMask)
      REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: data !(nlon,nlat,ntmaxs)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tlat
      integer, DIMENSION(:), ALLOCATABLE :: latproxy, lonproxy
      Integer, DIMENSION(:,:), ALLOCATABLE :: IndMask  ! (nlon,nlat)
    end subroutine AverageBoxMask
    subroutine nanmean(vector,iilat,iilon,undef,nanmeanvalue)
      REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: vector
      integer, dimension(1) :: iilat, iilon
      REAL(kind=4) :: undef,nanmeanvalue
    end subroutine nanmean
  end interface

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!     VARIABLE DEFINITION       !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TYPE variable
    character(nf90_max_name) :: name
    character(len=10) :: type
    character(len=6) :: box_type    
    character(len=3) :: ErrObs_dim    
    character(len=3) :: Data_weights_dim    
    character(len=3) :: Model_weights_dim    
    character(len=1), dimension(4) :: dim
    real :: weight
    real :: Ci
    real :: sigma
    integer :: nb_box ! FK
    integer :: startDateY ! FK
    integer :: startDateD ! FK    
    logical :: refmodfixed
    integer :: startDateYobs, startDateDobs
    integer :: startDateYrefmod, startDateDrefmod
    character(len=256) :: obs
    character(len=256) :: refobs
    character(len=256) :: errobs
    character(len=256) :: refmod
    character(len=256) :: cov
    character(len=256) :: modelweightspath
    character(len=256) :: dataweightspath
    integer, DIMENSION(:), ALLOCATABLE :: latmin
    integer, DIMENSION(:), ALLOCATABLE :: latmax
    integer, DIMENSION(:), ALLOCATABLE :: lonmin
    integer, DIMENSION(:), ALLOCATABLE :: lonmax
    integer, DIMENSION(:), ALLOCATABLE :: latproxy
    integer, DIMENSION(:), ALLOCATABLE :: lonproxy
  END TYPE variable

  type(variable), dimension(10) :: var

  character(nf90_max_name), DIMENSION(:,:), ALLOCATABLE :: invar
  character(nf90_max_name), DIMENSION(:), ALLOCATABLE :: boxvar

  ! LOOP INDEX
  INTEGER :: ii, jj, i, j
  ! DEBUG INDEX
  !INTEGER :: i1debug, i2debug, j1debug, j2debug, buffer

  !Definition of variable to use these previous function
  REAL, DIMENSION(:), POINTER :: pvalue
  CHARACTER(len=256), DIMENSION(:), POINTER :: pdimName
  INTEGER, DIMENSION(:), POINTER :: pdimSize

  INTEGER :: nlat, nlon, ninvar
  REAL(kind=4) :: undef, radian, lambda, NbStdDev, pundef
  REAL(kind=4) :: lscale, fcost, fcostTotal

  DOUBLE PRECISION :: lat_min,lat_max,lon_min,lon_max
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tlon, Tlat, TlonDeg, TlatDeg
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tmask

  INTEGER(kind=4), DIMENSION(:,:), ALLOCATABLE :: Mask, OceanArea      !(nlon,nlat) FK on le tappe en real..
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE ::  DiffModObs  !(nlon,nlat,ntmaxs)
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE ::    MeanDiffModObs, MeanDiffModObsF !(nlon,nlat)
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE ::  obsAnomaly, modAnomaly, refobs3D, refmod3D, obs3D, mod3D, OceanAreaFK !(nlon,nlat,ntmaxs)

  CHARACTER(len=256) :: newdata_atmos, newdata_atmos2, newdata_atmos3, newdata_atmos4, newdata_atmos5, newdata_atmos6, newdata_ocean, newdata_evolu, OceanAreaPath, OceanAreaVar

  INTEGER :: nbvariable, idx, count
  INTEGER :: t,v,nst,nttime,idxLat,idxLon, t_ref
  INTEGER :: t1_y, t1_d, totalt_ref
  INTEGER :: totalt, debut_mois, fin_mois, totalt_month, frequence, frequence_ref ! FK
  INTEGER :: fin_mois_model, debut_mois_model ! FK

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: rms !(nlon,nlat)
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: dataweights !(nlon,nlat)
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: modelweights !(nlon,nlat)

  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: rms_3D !(time,nlon,nlat)
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: dataweights_3D !(time,nlon,nlat)
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: modelweights_3D !(time,nlon,nlat)

  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: cov !(nlon*nlat,nlon*nlat)


  REAL(kind=4) :: CiUndef
  REAL(kind=4), DIMENSION(:), ALLOCATABLE :: Ci !(nlon*nlat,nlon*nlat)
  CHARACTER(len=10) :: methode, param, modestr
  INTEGER :: periode
  LOGICAL :: mode_online=.true.

  LOGICAL :: ln_geow=.true.  
  LOGICAL :: lexist_atmos, lexist_atmos2, lexist_atmos3, lexist_atmos4, lexist_atmos5, lexist_atmos6, lexist_ocean    ! flag for file presence

  CHARACTER(len=256) :: fcostEPFfile, box_file
  REAL(kind=4) :: fcostEPF
  INTEGER :: iflag

  ! Variablie pour allocate cov dim
  integer :: nlonCov, nlatCov

  ! Variable pour les log
  INTEGER(kind=4), DIMENSION(:,:), ALLOCATABLE :: year_BP_modref, month_BP_modref      !(ntmaxs,1)
  INTEGER(kind=4), DIMENSION(:,:), ALLOCATABLE :: year_BP_obs, month_BP_obs            !(ntmaxs,1)

  ! Sea-ice Month (FK)
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: FK_calc !(nlon,nlat)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!     DEF NAMELIST MODDATA_CO     !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NAMELIST/namglob/ lambda, nst, lscale, NbStdDev, ln_geow, debut_mois, fin_mois, debut_mois_model, fin_mois_model,  &
    frequence, fcostEPFfile, newdata_atmos, newdata_atmos2, newdata_atmos3, newdata_atmos4, newdata_atmos5, newdata_atmos6, newdata_ocean, &
    newdata_evolu, ninvar, t1_y, t1_d, OceanAreaPath, OceanAreaVar, box_file
  NAMELIST/namdom/ lon_min, lon_max, lat_min, lat_max
  NAMELIST/namvar/ invar
  NAMELIST/nambox/ boxvar
  !!----------------------------------------------------------------------

  ! Mask for boxes
  Integer, DIMENSION(:,:), ALLOCATABLE :: IndMask  ! (nlon,nlat)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!     GET ARGUMENT TO SELECT ASSIM METHODE       !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,methode)
  call getarg(2,param)
  call getarg(3,modestr)
  read(param,*) periode
  write(*,*) "##########################"
  write(*,*) "#### Start Moddata_co ####"
  write(*,'(a1,a14,a4,a8)') ' ', "#### methode: ", trim(methode), "    ####"
  write(*,'(a1,a14,i4,a8)') ' ', "#### periode: ", periode, "    ####"
  write(*,'(a1,a14,a4,a8)') ' ', "#### mode:    ", trim(modestr), "    ####"
  write(*,*) "##########################"
  if(modestr=='off') mode_online=.false.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!     READ MODDATA_CO.PARAM       !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(161,file='moddata_co.namelist',status='old',form='formatted')
  read(161,NML=namglob)
  !.Laplacian filter
  write(*,*) 'Filter parameters : ',lambda,nst,lscale
  !.Barrier
  write(*,*) 'Barrier',NbStdDev
  !.Geo w
  write(*,*) 'Geo. weight = ',ln_geow
  !.Domain
  read(161,NML=namdom)
  write(*,'(1X,A,4F10.2)') 'Domain = ',lon_min, lon_max, lat_min, lat_max

  if((lat_min.eq.-1).or.(lat_max.eq.-1)) then
    lat_min=-1; lat_max=-1
    write(*,*) '  No limitation for latitude'
  else
    write(*,*) '  Latitude limitation is ', lat_min, lat_max
  endif
  if((lon_min.eq.-1).or.(lon_max.eq.-1)) then
    lon_min=-1; lon_max=-1
    write(*,*) '  No limitation for longitude'
  else
    write(*,*) '  Longitude limitation is ', lon_min, lon_max
  endif

  !Read variable parameter
  allocate(invar(28,ninvar))
  read(161,NML=namvar)
  nbvariable=0

  write(*,*) "Variable parameters : "
  do v=1,ninvar
    var(v)%name=invar(1,v)
    var(v)%type=invar(2,v)
    var(v)%dim(1)=invar(3,v)
    var(v)%dim(2)=invar(4,v)
    var(v)%dim(3)=invar(5,v)
    var(v)%dim(4)=invar(6,v)
    read(invar(7,v),*) var(v)%weight
    read(invar(8,v),*) var(v)%Ci
    var(v)%Ci=var(v)%Ci*var(v)%Ci
    read(invar(9,v),*) var(v)%sigma
    read(invar(10,v),*) var(v)%nb_box
    allocate(var(v)%latproxy(var(v)%nb_box)); allocate(var(v)%lonproxy(var(v)%nb_box)); allocate(var(v)%latmin(var(v)%nb_box))
    allocate(var(v)%latmax(var(v)%nb_box)); allocate(var(v)%lonmin(var(v)%nb_box)); allocate(var(v)%lonmax(var(v)%nb_box))
    nbvariable=max(nbvariable,var(v)%nb_box)
    read(invar(11,v),*) var(v)%box_type
    read(invar(12,v),*) var(v)%startDateY
    read(invar(13,v),*) var(v)%startDateD
    var(v)%refmodfixed=(invar(14,v)=='1')
    read(invar(15,v),*) var(v)%startDateYobs
    read(invar(16,v),*) var(v)%startDateDobs
    read(invar(17,v),*) var(v)%startDateYrefmod
    read(invar(18,v),*) var(v)%startDateDrefmod
    var(v)%obs=invar(19,v)
    var(v)%refobs=invar(20,v)
    var(v)%errobs=invar(21,v)
    read(invar(22,v),*) var(v)%ErrObs_dim
    var(v)%refmod=invar(23,v)
    var(v)%cov=invar(24,v)
    var(v)%dataweightspath=invar(25,v)
    read(invar(26,v),*) var(v)%Data_weights_dim
    var(v)%modelweightspath=invar(27,v) 
    read(invar(28,v),*) var(v)%Model_weights_dim

    write(*,*) "  ",trim(var(v)%name), "(",var(v)%type,")"
    write(*,*) "  ------------------"
    write(*,*) "  Weight=",var(v)%weight, "Ci=",var(v)%Ci, "Sigma=",var(v)%sigma
    call newwrite2char("  obs:    ", trim(var(v)%obs))
    call newwrite2char("  refobs: ", trim(var(v)%refobs))
    call newwrite2char("  errobs: ", trim(var(v)%errobs))
    call newwrite2char("  cov:    ", trim(var(v)%cov))
    call newwrite2char("  refmod: ", trim(var(v)%refmod))
    call newwrite2char("  Data weights: ", trim(var(v)%dataweightspath))
    call newwrite2char("  Model weights: ", trim(var(v)%modelweightspath))   
    if(.not.mode_online) write(*,*) "  ",var(v)%startDateY, var(v)%startDateD ! FK chgmts des dates
    if(mode_online) write(*,*) "  ",var(v)%startDateYobs, var(v)%startDateDobs
  enddo

  write(*,*) "Model data :"
  if(mode_online) then
    call newwritechar2int("  mod: "//trim(newdata_atmos),t1_y,t1_d)
    call newwritechar2int("  mod: "//trim(newdata_atmos2),t1_y,t1_d)
    call newwritechar2int("  mod: "//trim(newdata_atmos3),t1_y,t1_d)
    call newwritechar2int("  mod: "//trim(newdata_atmos4),t1_y,t1_d) 
    call newwritechar2int("  mod: "//trim(newdata_atmos5),t1_y,t1_d)
    call newwritechar2int("  mod: "//trim(newdata_atmos6),t1_y,t1_d)
    call newwritechar2int("  mod: "//trim(newdata_ocean),t1_y,t1_d)
    call newwritechar2int("  mod: "//trim(newdata_evolu),t1_y,t1_d)
  else
    ! FK on change de nom de t1_y et t1_d vers debut_mois et fin_mois
    call newwritechar2int("  mod: "//trim(newdata_atmos),debut_mois,fin_mois)
    call newwritechar2int("  mod: "//trim(newdata_atmos2),debut_mois,fin_mois)
    call newwritechar2int("  mod: "//trim(newdata_atmos3),debut_mois,fin_mois)
    call newwritechar2int("  mod: "//trim(newdata_atmos4),debut_mois,fin_mois)
    call newwritechar2int("  mod: "//trim(newdata_atmos5),debut_mois,fin_mois)
    call newwritechar2int("  mod: "//trim(newdata_atmos6),debut_mois,fin_mois)  
    call newwritechar2int("  mod: "//trim(newdata_ocean),debut_mois,fin_mois)
    call newwritechar2int("  mod: "//trim(newdata_evolu),debut_mois,fin_mois)
  endif
  write(*,*) "Clio Area Path : ", trim(OceanAreaPath)
  write(*,*) "Clio Area Var : ", trim(OceanAreaVar)

  if(.not.mode_online) then
    write(*,*) "**********************"
    write(*,*) "Month model"  
    write(*,*) debut_mois_model,fin_mois_model
    write(*,*) "**********************"

    write(*,*) "**********************"
    write(*,*) "Month data"  
    write(*,*) debut_mois,fin_mois
    write(*,*) "**********************"  
  endif

  ! Assigne nombre de boxes
  if(nbvariable.ne.0) then
    allocate(boxvar(ninvar*(nbvariable*6+1)))
    read(161,NML=nambox)
    write(*,*) "Box parameters : "
    idx=1
    ! Regarde si il y a des box pour les variables et alloue les dimensions des box pour chaque variable
    do v=1,ninvar
      if((var(v)%nb_box.ne.0).and.(var(v)%name==boxvar(idx))) then
        count=idx+6*var(v)%nb_box
        i=idx+1
        nbvariable=1
        do while (i <= count)
          !write(*,*) idx, i, nbvariable, count
          read(boxvar(i),*) var(v)%latproxy(nbvariable)
          read(boxvar(i+1),*) var(v)%lonproxy(nbvariable)
          read(boxvar(i+2),*) var(v)%latmin(nbvariable)
          read(boxvar(i+3),*) var(v)%latmax(nbvariable)
          read(boxvar(i+4),*) var(v)%lonmin(nbvariable)
          read(boxvar(i+5),*) var(v)%lonmax(nbvariable)
          i=i+6
          nbvariable=nbvariable+1
        enddo
      endif
      idx=1+idx+6*var(v)%nb_box
    enddo
    if(idx>1) then
      write(*,*) 'CHECKBOX : average are computed over area before compute diffModObs'
    else
      write(*,*) 'CHECKBOX : no average are computed over area before compute diffModObs'
    endif
    close(161)
  endif

  !number of months between 2 assimilations  
  totalt_month=t1_y*12+t1_d/30
  write(*,*) "Assimilation frequence (mois): ", frequence

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   READ CONSTANT & FCOST INITIAL VALUE   !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  radian=2*acos(0.0)/180.0
  if((methode.eq."FAP").or.(methode.eq."NeTF")) then
    fcostTotal=1.
    write(*,*) "/!\ W.A.R.N.I.N.G methode is FAP, then fcostTotal is initialized to ",fcostTotal
  else
    fcostTotal=0.
    write(*,*) "/!\ W.A.R.N.I.N.G methode is not FAP, then fcostTotal is initialized to ",fcostTotal
  endif

  ! Test if atmos and ocean files from simulation have
  ! been saved and then exists ! If not : write in fcostopt.dat
  ! and STOP moddata_co
  call newwritechar("  "//trim(newdata_atmos))
  INQUIRE(FILE=newdata_atmos, EXIST=lexist_atmos)

  call newwritechar("  "//trim(newdata_atmos2))
  INQUIRE(FILE=newdata_atmos2, EXIST=lexist_atmos2)

  call newwritechar("  "//trim(newdata_atmos3))
  INQUIRE(FILE=newdata_atmos3, EXIST=lexist_atmos3)

  call newwritechar("  "//trim(newdata_atmos4))
  INQUIRE(FILE=newdata_atmos4, EXIST=lexist_atmos4)

  call newwritechar("  "//trim(newdata_atmos5))
  INQUIRE(FILE=newdata_atmos5, EXIST=lexist_atmos5)

  call newwritechar("  "//trim(newdata_atmos6))
  INQUIRE(FILE=newdata_atmos6, EXIST=lexist_atmos6)

  call newwritechar("  "//trim(newdata_ocean))
  INQUIRE(FILE=newdata_ocean, EXIST=lexist_ocean)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! LOOP FOR MULTIVARIABLE ASSIMILATION !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do v=1, ninvar

    IF ((var(v)%type=="ocean" .AND. .NOT. lexist_ocean) .OR. (var(v)%type=="atmos" .AND. .NOT. lexist_atmos)) THEN
      !CREATE OUTPUT FCOSTOPT.DAT
      open (37, FILE='fcostopt.dat')
      write(37,*) 'Cost function of the simulation'
      write(37,*) 0.
      close(37)
      WRITE(*,*) "W A R N I N G : error in getfile. fcost=0."
      STOP
    END IF

    ! IF PROXY EMPTY, DON'T DO ANYTHING !!!
    if ((var(v)%obs(len(trim(var(v)%obs))-13:len(trim(var(v)%obs)))) == 'ProxysEmpty.nc') then
      write(*,*) 'Because your file is named ProxysEmpty.nc we assumed that there is no data in your assimilation box' 
      !CREATE OUTPUT FCOSTOPT.DAT
      open (37, FILE='fcostopt.dat')
      write(37,*) 'Cost function of the simulation'
      write(37,*) 1.
      close(37)
      WRITE(*,*) "W A R N I N G : nodata in your assimilation box. fcost=1"
      STOP
    endif

    write(*,*) "########### ",trim(var(v)%name)," ###########"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!     READ MOD       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Mod :', var(v)%type
    if(var(v)%type=="atmos") then
      call newwritechar("  "//trim(newdata_atmos))
      if(mode_online) call readvarcdf(newdata_atmos, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_atmos, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    elseif(var(v)%type=="atmos2") then
      call newwritechar("  "//trim(newdata_atmos2))
      if(mode_online) call readvarcdf(newdata_atmos2, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_atmos2, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    elseif(var(v)%type=="atmos3") then
      call newwritechar("  "//trim(newdata_atmos3))
      if(mode_online) call readvarcdf(newdata_atmos3, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_atmos3, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    elseif(var(v)%type=="atmos4") then
      call newwritechar("  "//trim(newdata_atmos4))
      if(mode_online) call readvarcdf(newdata_atmos4, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_atmos4, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    elseif(var(v)%type=="atmos5") then
      call newwritechar("  "//trim(newdata_atmos5))
      if(mode_online) call readvarcdf(newdata_atmos5, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_atmos5, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    elseif(var(v)%type=="atmos6") then
      call newwritechar("  "//trim(newdata_atmos6))
      if(mode_online) call readvarcdf(newdata_atmos6, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_atmos6, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    elseif(var(v)%type=="ocean") then
      if(var(v)%name=="sim") then
        call newwritechar("  "//trim(newdata_ocean))
        if(mode_online) call readvarcdf(newdata_ocean, "albq", pvalue, pdimName, pdimSize)
        if(.not.mode_online) call readvarcdf(newdata_ocean, "albq", pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
      elseif (var(v)%name/="sim") then
        call newwritechar("  "//trim(newdata_ocean))
        if(mode_online) call readvarcdf(newdata_ocean, var(v)%name, pvalue, pdimName, pdimSize)
        if(.not.mode_online) call readvarcdf(newdata_ocean, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
      endif
    elseif(var(v)%type=="evolu") then
      call newwritechar("  "//trim(newdata_evolu))
      if(mode_online) call readvarcdf(newdata_evolu, var(v)%name, pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_evolu, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/debut_mois_model,frequence/))  !limit = (start, count)
    endif
    call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, mod3D)
    deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
    totalt=size(mod3D,3)
    write(*,*) 'Mod size=',size(mod3D,1),size(mod3D,2),size(mod3D,3)
    write(*,*) 'totalt1',totalt

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !!!!!     READ OBS       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    pundef=0.
    nttime=(t1_y-var(v)%startDateYobs)*12+(t1_d-var(v)%startDateDobs)/30+1
    if(mode_online) then
      write(*,*) 'Obs : select from nttime',nttime,'and count',totalt_ref
      write(*,*) var(v)%startDateYobs,var(v)%startDateDobs
    else
      write(*,*) 'Obs : select from debut_mois',debut_mois,'and count',totalt_ref
      write(*,*) var(v)%startDateY,var(v)%startDateD
    endif
    call newwritechar("  "//trim(var(v)%obs))
    if(mode_online) call readvarcdf(var(v)%obs, var(v)%name, pvalue, pdimName, pdimSize, pundef, dimtolimit = "time", limit = (/nttime,totalt/))
    if(.not.mode_online) call readvarcdf(var(v)%obs, var(v)%name, pvalue, pdimName, pdimSize, pundef, dimtolimit = "time", limit = (/debut_mois,frequence/))  !limit = (start, count)
    call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, obs3D)
    deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
    write(*,*) 'Obs size',size(obs3D,1),size(obs3D,2),size(obs3D,3)

    if ( .not. var(v)%refmodfixed ) THEN

      allocate(year_BP_obs(totalt,1))
      call readvarcdf(var(v)%obs, "year_BP", pvalue, pdimName, pdimSize, pundef, dimtolimit = "year_BP", limit = (/nttime,totalt/))
      year_BP_obs(1:totalt,1:1)=reshape(pvalue, (/totalt,1/))
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);

      allocate(month_BP_obs(totalt,1))
      call readvarcdf(var(v)%obs, "month_BP", pvalue, pdimName, pdimSize, pundef, dimtolimit = "month_BP", limit = (/nttime,totalt/))
      month_BP_obs(1:totalt,1:1)=reshape(pvalue, (/totalt,1/))
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);

    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! CHOOSE GRID/UNDEF PARAMETER !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(var(v)%type=="ocean") then
      ! ATTENTION GRILLE REGULIERE OU PAS ???
      if(mode_online) call readvarcdf(newdata_ocean, "tlon", pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_ocean, "lon", pvalue, pdimName, pdimSize)
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      Tlon=reshape(pvalue, (/size(mod3D,1),size(mod3D,2)/))
      if(mode_online) call readvarcdf(newdata_ocean, "tlat", pvalue, pdimName, pdimSize)
      if(.not.mode_online) call readvarcdf(newdata_ocean, "lat", pvalue, pdimName, pdimSize)
      Tlat=reshape(pvalue, (/size(mod3D,1),size(mod3D,2)/))
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat 
      allocate(OceanArea(nlon,nlat))
      call readvarcdf(OceanAreaPath, OceanAreaVar, pvalue, pdimName, pdimSize) ! ici que ça foire       
      OceanArea=reshape(pvalue, (/nlon,nlat/))
      Tlon(:,:)=OceanArea! OceanArea is already normalize between
      Tlat(:,:)=OceanArea! 0 and 1
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize); deallocate(OceanArea)
    elseif(var(v)%type=="atmos") then
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      call readvarcdf(newdata_atmos, "lon", pvalue, pdimName, pdimSize)
      Tlon(1:size(mod3D,1),1:1)=reshape(pvalue, (/size(mod3D,1),1/))
      do i=2,size(mod3D,2) 
        Tlon(:,i)=Tlon(:,1) 
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      call readvarcdf(newdata_atmos, "lat", pvalue, pdimName, pdimSize)
      Tlat(1:1,1:size(mod3D,2))=reshape(pvalue, (/1,size(mod3D,2)/))
      do i=2,size(mod3D,1)
        Tlat(i,:)=Tlat(1,:)
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat
      Tlon(:,:)=cos(TlonDeg(:,:)*radian)
      Tlat(:,:)=cos(TlatDeg(:,:)*radian)
    elseif(var(v)%type=="atmos2") then
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      call readvarcdf(newdata_atmos2, "lon", pvalue, pdimName, pdimSize)
      Tlon(1:size(mod3D,1),1:1)=reshape(pvalue, (/size(mod3D,1),1/))
      do i=2,size(mod3D,2) 
        Tlon(:,i)=Tlon(:,1) 
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      call readvarcdf(newdata_atmos2, "lat", pvalue, pdimName, pdimSize)
      Tlat(1:1,1:size(mod3D,2))=reshape(pvalue, (/1,size(mod3D,2)/))
      do i=2,size(mod3D,1)
        Tlat(i,:)=Tlat(1,:)
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat
      Tlon(:,:)=cos(TlonDeg(:,:)*radian)
      Tlat(:,:)=cos(TlatDeg(:,:)*radian)
    elseif(var(v)%type=="atmos3") then
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      call readvarcdf(newdata_atmos3, "lon", pvalue, pdimName, pdimSize)
      Tlon(1:size(mod3D,1),1:1)=reshape(pvalue, (/size(mod3D,1),1/))
      do i=2,size(mod3D,2) 
        Tlon(:,i)=Tlon(:,1) 
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      call readvarcdf(newdata_atmos3, "lat", pvalue, pdimName, pdimSize)
      Tlat(1:1,1:size(mod3D,2))=reshape(pvalue, (/1,size(mod3D,2)/))
      do i=2,size(mod3D,1)
        Tlat(i,:)=Tlat(1,:)
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat
      Tlon(:,:)=cos(TlonDeg(:,:)*radian)
      Tlat(:,:)=cos(TlatDeg(:,:)*radian)
    elseif(var(v)%type=="atmos4") then
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      call readvarcdf(newdata_atmos4, "lon", pvalue, pdimName, pdimSize)
      Tlon(1:size(mod3D,1),1:1)=reshape(pvalue, (/size(mod3D,1),1/))
      do i=2,size(mod3D,2) 
        Tlon(:,i)=Tlon(:,1) 
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      call readvarcdf(newdata_atmos4, "lat", pvalue, pdimName, pdimSize)
      Tlat(1:1,1:size(mod3D,2))=reshape(pvalue, (/1,size(mod3D,2)/))
      do i=2,size(mod3D,1)
        Tlat(i,:)=Tlat(1,:)
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat
      Tlon(:,:)=cos(TlonDeg(:,:)*radian)
      Tlat(:,:)=cos(TlatDeg(:,:)*radian)
    elseif(var(v)%type=="atmos5") then
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      call readvarcdf(newdata_atmos5, "lon", pvalue, pdimName, pdimSize)
      Tlon(1:size(mod3D,1),1:1)=reshape(pvalue, (/size(mod3D,1),1/))
      do i=2,size(mod3D,2) 
        Tlon(:,i)=Tlon(:,1) 
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      call readvarcdf(newdata_atmos5, "lat", pvalue, pdimName, pdimSize)
      Tlat(1:1,1:size(mod3D,2))=reshape(pvalue, (/1,size(mod3D,2)/))
      do i=2,size(mod3D,1)
        Tlat(i,:)=Tlat(1,:)
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat
      Tlon(:,:)=cos(TlonDeg(:,:)*radian)
      Tlat(:,:)=cos(TlatDeg(:,:)*radian)
    elseif(var(v)%type=="atmos6") then
      allocate(Tlon(size(mod3D,1),size(mod3D,2)))
      allocate(Tlat(size(mod3D,1),size(mod3D,2)))
      Tlon=-1000; Tlat=-1000
      call readvarcdf(newdata_atmos6, "lon", pvalue, pdimName, pdimSize)
      Tlon(1:size(mod3D,1),1:1)=reshape(pvalue, (/size(mod3D,1),1/))
      do i=2,size(mod3D,2) 
        Tlon(:,i)=Tlon(:,1) 
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      call readvarcdf(newdata_atmos6, "lat", pvalue, pdimName, pdimSize)
      Tlat(1:1,1:size(mod3D,2))=reshape(pvalue, (/1,size(mod3D,2)/))
      do i=2,size(mod3D,1)
        Tlat(i,:)=Tlat(1,:)
      enddo
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      nlon=size(mod3D,1); nlat=size(mod3D,2); undef=-99.99
      if(pundef/=0.) undef=pundef
      write(*,*) "Missing Value =",undef
      allocate(TlonDeg(size(mod3D,1),size(mod3D,2)))
      allocate(TlatDeg(size(mod3D,1),size(mod3D,2)))
      TlonDeg=Tlon
      TlatDeg=Tlat
      Tlon(:,:)=cos(TlonDeg(:,:)*radian)
      Tlat(:,:)=cos(TlatDeg(:,:)*radian)
    elseif(var(v)%type=="evolu") then
      nlon=1; nlat=1
      allocate(Tlon(1,1))
      allocate(Tlat(1,1))
      Tlon=-1000; Tlat=-1000
    else
      write(*,*) "Unknow type ",var(v)%type, " for variable ",trim(var(v)%name)
      stop
    endif

    write(*,*) "Grid : nlat =", nlat, "nlon =", nlon
    write(*,*) "Tlon =", size(Tlon,1),size(Tlon,2)
    write(*,*) "Tlat =", size(Tlat,1),size(Tlat,2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!     READ MOD REF       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( var(v)%refmodfixed ) THEN           !Si la reference du modele est fixee
      write(*,*) 'Mod ref :    Model reference is fixed'  
      nttime=(t1_d-1)/30+1
      totalt_ref=totalt
      IF (totalt .GT. 12) THEN
        IF (MOD(totalt,12) .NE. 0) THEN
          WRITE(*,*) "totalt has to be < 12 or multiple of 12"
          STOP
        ELSE
          WRITE(*,*) "totalt is < 12 or multiple of 12"
          totalt_ref=12
        END IF
      END IF
      ! FK on enlève kes reference bougeante, et je rajoute une conditions: si pas de fichier de refmod, on prend en compte zero.
      ! FK si la fréquence d'assimilation est plus grande qu'un an, la taille de la ref est égale à 12. On doit mettre cette option 
      ! pcq si frequence = mensuelle, taille de la ref=1.
      frequence_ref = frequence
      IF (frequence_ref .GT. 12) THEN ! FK 
        frequence_ref=12
      END IF ! FK 
      if(len(trim(var(v)%refmod))==0) then  ! FK
        call newwritechar("  No ref") ! FK
        allocate(refmod3D(nlon,nlat,frequence_ref)) ! FK
        refmod3D=0 ! FK
      else ! FK
        call newwritechar("  "//trim(var(v)%refmod))
        if(mode_online) call readvarcdf(var(v)%refmod, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/nttime,totalt_ref/)) ! Online
        if(.not.mode_online) call readvarcdf(var(v)%refmod, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/1,frequence_ref/)) ! FK
        call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, refmod3D)
        deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize)
      endif
      write(*,*) 'RefMod size',size(refmod3D,1),size(refmod3D,2),size(refmod3D,3)
    ELSEIF ( .not. var(v)%refmodfixed ) THEN !Si la reference du modele evolue
      write(*,*) 'Mod ref :    Model reference evolve trough time'
      nttime=((t1_y-var(v)%startDateYrefmod)*12+(t1_d-var(v)%startDateDrefmod)/30+1)

      totalt_ref=totalt
      write(*,*) 't1_y',t1_y
      write(*,*) 'var(v)%startDateYrefmod)',var(v)%startDateYrefmod
      write(*,*) 't1_d',t1_d
      write(*,*) 'var(v)%startDateDrefmod)',var(v)%startDateDrefmod
      write(*,*) 'nttime', nttime
      write(*,*) 'totalt_ref', totalt_ref
      call newwritechar("  "//trim(var(v)%refmod))
      call readvarcdf(var(v)%refmod, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/nttime,totalt_ref/))
      call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, refmod3D)
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      allocate(year_BP_modref(totalt_ref,1))
      call readvarcdf(var(v)%refmod, "year_BP", pvalue, pdimName, pdimSize, dimtolimit = "year_BP", limit = (/nttime,totalt_ref/))
      year_BP_modref(1:totalt_ref,1:1)=reshape(pvalue, (/totalt_ref,1/))
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      allocate(month_BP_modref(totalt_ref,1))
      call readvarcdf(var(v)%refmod, "month_BP", pvalue, pdimName, pdimSize, dimtolimit = "month_BP", limit = (/nttime,totalt_ref/))
      month_BP_modref(1:totalt_ref,1:1)=reshape(pvalue, (/totalt_ref,1/))
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      write(*,*) 'RefMod size',size(refmod3D,1),size(refmod3D,2),size(refmod3D,3)
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!  transfer albq to sim  !!!!! pour assimiler sea-ice month
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     The following calculates the mean number of months per year with a sea-ice
    !     concentration greater than 50% in case of the loaded proxy variable is 'sim'.    
    IF (trim(var(v)%name) .EQ. 'sim') then
      IF (totalt .NE. 12) THEN ! FK
        WRITE(*,*) "The assimilation has to be annual to assimilate the variable sea-ice month." ! FK
        STOP ! FK
      ELSE    
        ! Pour modèle
        allocate(FK_calc(nlon,nlat))
        where (mod3D(:,:,:).lt. 0 .or. mod3D(:,:,:).gt. 1)
          mod3D(:,:,:) = undef
        elsewhere (mod3D(:,:,:).ge. 0.5 .and. mod3D(:,:,:).ne. undef)
          mod3D(:,:,:) = 1
        elsewhere (mod3D(:,:,:).lt. 0.5 .and. mod3D(:,:,:).ne. undef)
          mod3D(:,:,:) = 0
        endwhere
        FK_calc(:,:) = sum(mod3D,3) / (totalt/12)
        mod3D(:,:,:) = spread(FK_calc,3,totalt) 
        where (mod3D(:,:,:).lt. -12)
          mod3D(:,:,:) = undef
        endwhere
        deallocate(FK_calc)

        ! Pour ref modèle
        allocate(FK_calc(nlon,nlat))
        where (refmod3D(:,:,:).lt. 0 .or. refmod3D(:,:,:).gt. 1)
          refmod3D(:,:,:) = undef
        elsewhere (refmod3D(:,:,:).ge. 0.5 .and. refmod3D(:,:,:).ne. undef)
          refmod3D(:,:,:) = 1
        elsewhere (refmod3D(:,:,:).lt. 0.5 .and. refmod3D(:,:,:).ne. undef)
          refmod3D(:,:,:) = 0
        endwhere
        FK_calc(:,:) = sum(refmod3D,3) / (totalt/12)
        refmod3D(:,:,:) = spread(FK_calc,3,totalt) 
        where (refmod3D(:,:,:).lt. -12)
          refmod3D(:,:,:) = undef
        endwhere
        deallocate(FK_calc)
      END IF ! FK       
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!     READ OBS REF       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
    nttime=(t1_d-1)/30+1
    if(len(trim(var(v)%refobs))==0) then 
      call newwritechar("Obs ref : No ref for "//trim(var(v)%name))
      if(mode_online) allocate(refobs3D(nlon,nlat,totalt_ref))
      if(.not.mode_online) allocate(refobs3D(nlon,nlat,frequence_ref))
      refobs3D=0
    else
      call newwritechar("  "//trim(var(v)%refobs))
      if(mode_online) call readvarcdf(var(v)%refobs, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/nttime,totalt_ref/))
      if(.not.mode_online) call readvarcdf(var(v)%refobs, var(v)%name, pvalue, pdimName, pdimSize, dimtolimit = "time", limit = (/1,frequence_ref/)) ! FK
      call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, refobs3D)
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
    endif

    if ( .not. var(v)%refmodfixed ) THEN
      write(*,*) '-------     Time correspondance      -------'
      write(*,'(A15,A13,I6,A6,I3,A11,I6,A6,I3,A16,I4)') '| MODEL REF: ','FROM yearBP', year_BP_modref(1,1), 'month', month_BP_modref(1,1), 'TO yearBP', year_BP_modref(size(year_BP_modref,1),1), 'month', month_BP_modref(size(year_BP_modref,1),1),' -> TOTAL MONTH =',totalt_ref
      write(*,'(A15,A13,I6,A6,I3,A11,I6,A6,I3,A16,I4)') '| PROXY    : ','FROM yearBP', year_BP_obs(1,1), 'month', month_BP_obs(1,1), 'TO yearBP', year_BP_obs(size(year_BP_obs,1),1), 'month', month_BP_obs(size(year_BP_obs,1),1),' -> TOTAL MONTH =',totalt
      write(*,*) '--------------------------------------------'
      deallocate(year_BP_modref); deallocate(month_BP_modref); deallocate(year_BP_obs); deallocate(month_BP_obs)
    endif

    IF (size(refobs3D,3) .GT. 12) THEN ! FK
      WRITE(*,*) "size month refmod has to be max 12" ! FK
      STOP ! FK
    END IF   ! FK

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! COMPUTE ANOMALY AND MEAN !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "Compute anomalies"
    !write(*,*) 'totalt2',totalt
    allocate(obsAnomaly(nlon, nlat, totalt))
    allocate(modAnomaly(nlon, nlat, totalt))
    allocate(DiffModObs(nlon, nlat, totalt))

    do t=1, totalt
      t_ref=MOD(t,frequence); IF (t_ref==0) t_ref=frequence   !assume ref have only 12 time slide ! FK je comprends pas
      if (t_ref.gt.12) then 
        t_ref = t_ref-12
      endif
      do idxLon=1,nlon
        do idxLat=1,nlat
          modAnomaly(idxLon,idxLat,t)=(mod3D(idxLon,idxLat,t)-refmod3D(idxLon,idxLat,t_ref))*var(v)%weight
          if((obs3D(idxLon,idxLat,t).eq.undef).or.(mod3D(idxLon,idxLat,t).eq.undef).or.(refobs3D(idxLon,idxLat,t_ref).eq.undef).or.(refmod3D(idxLon,idxLat,t_ref).eq.undef)) then
            obsAnomaly(idxLon,idxLat,t)=undef
          else
            obsAnomaly(idxLon,idxLat,t)=(obs3D(idxLon,idxLat,t)-refobs3D(idxLon,idxLat,t_ref))*var(v)%weight
          endif
        enddo
      enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!     READ MOD and data weights       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read weight to multiply with data and model results FK
    ! LOAD MODELS WEIGHTS (either 2D or 3D)
    if(len(trim(var(v)%modelweightspath))==0) then  ! FK
      write(*,*) "No model weight loaded"
    else
      write(*,*)'Load model weights ' ! FK

      CiUndef=-99.99
      ALLOCATE(modelweights(nlon,nlat))

      IF (trim(var(v)%Model_weights_dim) .EQ. '2D') then

        call readvarcdf(var(v)%modelweightspath, "model_weights", pvalue, pdimName, pdimSize)
        modelweights=reshape(pvalue, (/nlon,nlat/))

        ! test file dimension
        IF (SIZE(pdimSize) .NE. 2) THEN
          PRINT *, "error in ",trim(var(v)%modelweightspath)," dimension :"
          PRINT *, SIZE(pdimSize)," are available for rms while it is defined as 2D"
          STOP
        END IF  

      ELSEIF (trim(var(v)%Model_weights_dim) .EQ. '3D') then

        call readvarcdf(var(v)%modelweightspath, "model_weights", pvalue, pdimName, pdimSize, pundef, dimtolimit = "time", limit = (/debut_mois,frequence/))
        call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, modelweights_3D)        
        call average(CiUndef,modelweights_3D, modelweights, 3)

      END IF
      DEALLOCATE(pvalue,pdimName, pdimSize)

      do t=1, totalt
        do idxLon=1,nlon
          do idxLat=1,nlat
            if(modAnomaly(idxLon,idxLat,t).ne.undef) then
              modAnomaly(idxLon,idxLat,t)=modAnomaly(idxLon,idxLat,t)*modelweights(idxLon,idxLat)
            endif
          enddo
        enddo
      enddo

    endif ! FK

    if(len(trim(var(v)%dataweightspath))==0) then  ! FK
      write(*,*) "No data weight loaded"
    else
      write(*,*)'Load data weights ' ! FK

      ! DATA WEIGHTS
      ALLOCATE(dataweights(nlon,nlat))

      IF (trim(var(v)%Data_weights_dim) .EQ. '2D') then

        call readvarcdf(var(v)%dataweightspath, "data_weights", pvalue, pdimName, pdimSize)
        dataweights=reshape(pvalue, (/nlon,nlat/))

        ! test file dimension
        IF (SIZE(pdimSize) .NE. 2) THEN
          PRINT *, "error in ",trim(var(v)%dataweightspath)," dimension :"
          PRINT *, SIZE(pdimSize)," are available for rms while it is defined as 2D"
          STOP
        END IF  

      ELSEIF (trim(var(v)%Data_weights_dim) .EQ. '3D') then

        call readvarcdf(var(v)%dataweightspath, "data_weights", pvalue, pdimName, pdimSize, pundef, dimtolimit = "time", limit = (/debut_mois,frequence/))
        call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, dataweights_3D)        
        call average(CiUndef,dataweights_3D, dataweights, 3)

      END IF
      DEALLOCATE(pvalue,pdimName, pdimSize)

      do t=1, totalt
        do idxLon=1,nlon
          do idxLat=1,nlat
            if(obsAnomaly(idxLon,idxLat,t).ne.undef) then
              obsAnomaly(idxLon,idxLat,t)=obsAnomaly(idxLon,idxLat,t)*dataweights(idxLon,idxLat)
            endif
          enddo
        enddo
      enddo

    endif ! FK

    write(*,*) 'box allocate proxy = ',var(v)%nb_box
    IF ( var(v)%nb_box.ne.0 ) then

      if ( var(v)%box_type == 'square' ) then

        write(*,*) "square box method"
        call AverageBoxAllocateMod(modAnomaly, Tlat, var(v)%latproxy, var(v)%lonproxy, var(v)%latmin, var(v)%latmax, var(v)%lonmin, var(v)%lonmax)

      elseif ( var(v)%box_type == 'file' ) then

        write(*,*) "file box method"
        !!! NEW WAY FOR BOXES 
        !! INCLUDE A Test based on a variable read in .PARAM
        !! READ FILE DEFINING BOXES: FILL In IndMask
        allocate(IndMask(nlon,nlat))
        call readvarcdf(box_file, "ireg", pvalue, pdimName, pdimSize) 
        IndMask=reshape(pvalue, (/nlon,nlat/))
        !write(*,*) 'IndMask',IndMask(180,152),IndMask(258,142),IndMAsk(nlon,nlat)
        deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize)

        !! CALL AVERAGE OVER BOXES
        !write(*,*) 'before box',modAnomaly(319,159,1)

        call AverageBoxMask(modAnomaly, Tlat, var(v)%latproxy, var(v)%lonproxy, IndMask)
        !write(*,*) 'after box',modAnomaly(319,159,1)
        !! END TEST CHOICE BOXES

      endif
    ENDIF

    ! Calcul du diffModObs
    write(*,*) "Compute diffModObs"
    do t=1, totalt
      do idxLon=1,nlon
        do idxLat=1,nlat
          if (obsAnomaly(idxLon,idxLat,t).ne.undef) then
            DiffModObs(idxLon,idxLat,t)=modAnomaly(idxLon,idxLat,t)-obsAnomaly(idxLon,idxLat,t)
          else
            DiffModObs(idxLon,idxLat,t)=undef
          endif
        enddo
      enddo
    enddo 

    IF ( var(v)%nb_box.ne.0 ) THEN
      if (allocated(var(v)%latproxy)) deallocate(var(v)%latproxy)
      if (allocated(var(v)%lonproxy)) deallocate(var(v)%lonproxy)
      if (allocated(var(v)%latmin)) deallocate(var(v)%latmin)
      if (allocated(var(v)%latmax)) deallocate(var(v)%latmax)
      if (allocated(var(v)%lonmin)) deallocate(var(v)%lonmin)
      if (allocated(var(v)%lonmax)) deallocate(var(v)%lonmax)
    ENDIF

    write(*,*) "Compute mean"
    allocate(MeanDiffModObs(nlon, nlat))
    call average(undef,DiffModObs, MeanDiffModObs, 3) !3 is the time index in pdimSize
    !call averageNew(undef,obsAnomaly, modAnomaly, MeanDiffModObs, 3) !3 is the time index in pdimSize

    write(*,*) "Compute laplacian filter"
    allocate(MeanDiffModObsF(nlon, nlat))
    MeanDiffModObsF=MeanDiffModObs
    if(var(v)%type/="evolu") CALL filter_scale(MeanDiffModObsF,undef,Tlat,lambda,lscale,nst) ! Diff OK

    deallocate(obs3D); deallocate(mod3D); deallocate(refobs3D); deallocate(refmod3D)
    deallocate(obsAnomaly); deallocate(modAnomaly); deallocate(DiffModObs); deallocate(MeanDiffModObs)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !!!!!    DEFINE DOMAIN      !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! we use extend* for the case of the box is for example (60 10 30 2)
    write(*,*) 'Apply domain'
    allocate(Mask(nlon, nlat))
    Mask=0
    write(*,'(4F7.2)') lon_min, lon_max, lat_min, lat_max
    call MakeBox(TlonDeg, TlatDeg, lon_min, lon_max, lat_min, lat_max, Mask)
    MeanDiffModObsF(:,:)=Mask(:,:)*MeanDiffModObsF(:,:)+(Mask(:,:)-1)*(-undef)
    deallocate(Mask)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!    LOAD NOISE/COV    !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! A.M. Load mask if ocean else, fill mask on nlon*nlat with value 1
    write(*,*) "type=", var(v)%type
    write(*,*) "periode=", periode

    ! FK Ici on charge d'abord un masque. Avant, on allait chopper la variable tmask ds clio mais pas portable.. Ici, on regarde plutôt si un fichier d'airte est chargé.
    ! Si c'est le cas, on charge le nc et on fait *0+1 => on obtient notre masque. Si pas de fichier d'aires mentionné, tmask=1 et à la même taille que les lon et lat du modèle chargé.
    allocate(Tmask(nlon,nlat));
    Tmask(:,:)=1;
    if((len(trim(OceanAreaPath))/=0).and.(var(v)%type=="ocean")) then ! FK
      call readvarcdf(OceanAreaPath, OceanAreaVar, pvalue, pdimName, pdimSize) ! FK
      call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, OceanAreaFK) ! FK
      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);  ! FK
      where (OceanAreaFK(:,:,1)==0)
        Tmask(:,:) = 0
      elsewhere
        Tmask(:,:) = 1
      endwhere
    endif

    ! FK on rajoute la possibilité qu'il n'y ait pas de matrice de cov
    if(len(trim(var(v)%cov))==0) then  ! FK
      write(*,*) "No matrix de covariance taken into account"
    else    
      write(*,*)'Load Covariance Matrix ' ! FK
      ! Find dim of cov matrix and load it from NetCDF-files 
      call readvarcdf(var(v)%cov, "cov", pvalue, pdimName, pdimSize)
      nlonCov=pdimSize(1); nlatCov=pdimSize(2);

      if ( size(pdimSize) .lt. 3) then
        periode=1
        write(*,*)'  W.A.R.N.I.N.G : Cov matrix have not dim time, so periode=', periode
      else
        if ( pdimSize(3) .lt. periode ) then
          periode=1
          write(*,*)'  W.A.R.N.I.N.G : Cov matrix dim time is too small, so periode=', periode
        endif
      endif

      deallocate(pvalue); deallocate(pdimName); deallocate(pdimSize);
      allocate(cov(nlonCov,nlatCov,1))
      allocate(Ci(nlonCov))
      CALL read_cov(var(v)%cov, cov,"cov", periode)
    endif ! FK

    PRINT *, 'ALLOCATE ',nlon,nlat
    ALLOCATE(rms(nlon,nlat))
    rms(:,:)=0.
    write(*,*) "type=", var(v)%type, "taille", size(rms)

    IF (len(trim(var(v)%errobs))==0) THEN
      CiUndef=var(v)%Ci
      write(*,*) '  No files for measurement error, therefore it will be Ci',var(v)%Ci
      Ci=var(v)%Ci
      rms(:,:)=var(v)%Ci
    ELSE

      ! FK Load error, either 2D or 3D
      write(*,*) '  Measurement error file is', trim(var(v)%errobs)
      CiUndef=-99.99
      IF (trim(var(v)%ErrObs_dim) .EQ. '2D') then

        CALL readvarcdf(trim(var(v)%errobs), "rms", pvalue, pdimName, pdimSize)
        rms=reshape(pvalue, (/nlon,nlat/))

        ! test file dimension
        IF (SIZE(pdimSize) .NE. 2) THEN
          PRINT *, "error in ",trim(var(v)%errobs)," dimension :"
          PRINT *, SIZE(pdimSize)," are available for rms while it is defined as 2D"
          STOP
        END IF  

      ELSEIF (trim(var(v)%ErrObs_dim) .EQ. '3D') then

        call readvarcdf(trim(var(v)%errobs), "rms", pvalue, pdimName, pdimSize, pundef, dimtolimit = "time", limit = (/debut_mois,frequence/))
        call ExtractUsefullDimension(var(v)%Dim, pdimSize, pvalue, rms_3d)        
        call average(CiUndef,rms_3d, rms, 3)

      END IF
      DEALLOCATE(pvalue,pdimName, pdimSize)

      ! FK on doit rajouter ça pour avoir la bonne taille quand pas de cov
      if(len(trim(var(v)%cov))==0) then ! FK
        nlonCov=nlon*nlat  ! FK
        allocate(Ci(nlonCov))  ! FK
      endif ! FK

      ! Fill the Ci matrix (only the diagonal)
      Ci(:)=CiUndef
      j=0
      DO jj=1,nlat
        DO ii=1,nlon
          if  (Tmask(ii,jj).gt.0.5) then! It is needeed for the ocean (AM 05/2012)
            j=j+1
            IF (ABS(rms(ii,jj)).NE.ABS(undef)) Ci(j)=rms(ii,jj)*rms(ii,jj)
          endif
        ENDDO
      ENDDO

      !! END P.M.
    END IF

    DEALLOCATE(rms)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !!!!!   EFFICIENT PF    !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) " "
    write(*,*) 'Get fcostEPF for the <<efficient particle filter>> part'
    fcostEPF=1.0
    if(len(trim(fcostEPFfile)).ne.0) then
      open(10,file=fcostEPFfile,IOSTAT=iflag,status='old',form='formatted')
      if (iflag.eq.0) then
        read(10,*) fcostEPF
      else
        write(*,*) 'The fcostEPF file does not exist!!! Therefore fcostEPF takes a value 1'
      endif
      close(10)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!    COMPUTE FCOST    !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( var(v)%type=="evolu" ) Tlat=0
    !Select cost function
    SELECT CASE(methode)
    CASE("FAP") 
      ! CALL cost_euclide(fcost,Tlat,var(v)%weight,undef,MeanDiffModObsF,2,cov,Ci,CiUndef,var(v)%sigma,NbStdDev,ln_geow,Tmask) !old
      CALL cost_euclide(fcost,Tlat,var(v)%weight,undef,MeanDiffModObsF,1,cov,Ci,CiUndef,var(v)%sigma,NbStdDev,ln_geow,Tmask) !likelihood
      write(*,*)'  Likelihood=',fcost
      write(*,*)'  fcostEPF=',fcostEPF
      fcostTotal=fcostTotal*var(v)%weight*fcost*fcostEPF
    CASE("NeTF")
      ! CALL cost_euclide(fcost,Tlat,var(v)%weight,undef,MeanDiffModObsF,2,cov,Ci,CiUndef,var(v)%sigma,NbStdDev,ln_geow,Tmask) !old
      CALL cost_euclide(fcost,Tlat,var(v)%weight,undef,MeanDiffModObsF,1,cov,Ci,CiUndef,var(v)%sigma,NbStdDev,ln_geow,Tmask) !likelihood
      write(*,*)'  Likelihood=',fcost
      write(*,*)'  fcostEPF=',fcostEPF
      fcostTotal=fcostTotal*var(v)%weight*fcost*fcostEPF
    CASE("Best")
      ! CALL cost_euclide(fcost,Tlat,var(v)%weight,undef,MeanDiffModObsF,1,cov,Ci,CiUndef,var(v)%sigma,NbStdDev,ln_geow,Tmask) !old
      CALL cost_euclide(fcost,Tlat,var(v)%weight,undef,MeanDiffModObsF,2,cov,Ci,CiUndef,var(v)%sigma,NbStdDev,ln_geow,Tmask) !likelihood
      write(*,*)'fcost=',fcost
      fcostTotal=fcostTotal+var(v)%weight*fcost

    CASE default
      write(*,*)"Unknow methode : ", methode
      STOP 0
    END SELECT

    if (allocated(cov))     deallocate(cov); 
    if (allocated(Ci))     deallocate(Ci)
    if (allocated(MeanDiffModObsF))     deallocate(MeanDiffModObsF)
    if (allocated(Tlon))     deallocate(Tlon); 
    if (allocated(Tlat))     deallocate(Tlat)
    if (allocated(TlonDeg))     deallocate(TlonDeg); 
    if (allocated(TlatDeg))     deallocate(TlatDeg)
    if (allocated(Tmask))     deallocate(Tmask)

  enddo

  write(*,*) '###############################'  
  write(*,*) "fcostTotal=",fcostTotal
  write(*,*) '###############################'  

  !CREATE OUTPUT FCOSTOPT.DAT
  open (37, FILE='fcostopt.dat')
  write(37,*) 'Cost function of the simulation'
  write(37,*) fcostTotal

  if(.not.mode_online) then
    ! FK on sauve un fichier texte reprenant l'année prise en compte, qui sera notée ds output_fcost via Partifilter.f90
    ! !   annee=fin_mois/12 
    ! !   annee=debut_mois+frequence
    ! !   
    ! !   annee_debut=debut_mois
    ! !   annee_fin=fin_mois
    open (38, FILE='mois_start.dat')
    write(38,*) debut_mois_model
    open (39, FILE='mois_end.dat')
    write(39,*) fin_mois_model
  endif

END PROGRAM moddata_co

SUBROUTINE ExtractUsefullDimension(DimAsk, DimInSize, DataIn, DataOut)
  IMPLICIT NONE
  character(len=1), dimension(4) :: DimAsk
  INTEGER, DIMENSION(:), POINTER :: DimInSize
  REAL, DIMENSION(:), POINTER :: DataIn
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: DataOut

  REAL(kind=4), DIMENSION(:,:,:,:), ALLOCATABLE :: DataTmp
  INTEGER, DIMENSION(4) :: DimOutSize, DimOutStart, DimOutEnd
  INTEGER :: SizeX, SizeY, SizeT, i

  DimOutStart=1; DimOutSize=1; SizeX=1; SizeY=1; SizeT=1
  do i=1, size(DimInSize)
    if((DimAsk(i)=='X').or.(DimAsk(i)=='Y').or.(DimAsk(i)=='T')) then
      DimOutSize(i)=DimInSize(i)
      DimOutStart(i)=1
      if(DimAsk(i)=='X') SizeX=DimInSize(i)
      if(DimAsk(i)=='Y') SizeY=DimInSize(i)
      if(DimAsk(i)=='T') SizeT=DimInSize(i)
    else
      DimOutSize(i)=1
      read(DimAsk(i),*) DimOutStart(i)
    endif
  enddo

  !write(*,'(A,I2,A,I2,A,I2,A,I2,A)') "DataTmp=reshape(DataIn, (/",DimInSize(1),",",DimInSize(2),",",DimInSize(3),",",DimInSize(4),"/))"
  if(size(DimInSize).eq.1) then
    allocate(DataTmp(DimInSize(1), 1, 1, 1))
    DataTmp=reshape(DataIn, (/DimInSize(1),1,1,1/))
  elseif(size(DimInSize).eq.2) then
    allocate(DataTmp(DimInSize(1), DimInSize(2), 1, 1))
    DataTmp=reshape(DataIn, (/DimInSize(1),DimInSize(2),1,1/))
  elseif(size(DimInSize).eq.3) then
    allocate(DataTmp(DimInSize(1), DimInSize(2), DimInSize(3), 1))
    DataTmp=reshape(DataIn, (/DimInSize(1),DimInSize(2),DimInSize(3),1/))
  elseif(size(DimInSize).eq.4) then
    allocate(DataTmp(DimInSize(1), DimInSize(2), DimInSize(3), DimInSize(4)))
    DataTmp=reshape(DataIn, (/DimInSize(1),DimInSize(2),DimInSize(3),DimInSize(4)/))
  endif

  if(allocated(DataOut)) deallocate(DataOut)
  allocate(DataOut(SizeX,SizeY,SizeT))
  DimOutEnd=DimOutStart+DimOutSize-1
  if(DimOutSize(4)==1) then
    !write(*,*) "4 size 1"
    !write(*,*) "DataOut(1:",SizeX,",1:",SizeY,",1:",SizeT,")=DataTmp(",DimOutStart(1),":",DimOutEnd(1),",",DimOutStart(2),":",DimOutEnd(2),",",DimOutStart(3),":",DimOutEnd(3),",",DimOutStart(4),")"
    DataOut(1:SizeX,1:SizeY,1:SizeT)=DataTmp(DimOutStart(1):DimOutEnd(1), DimOutStart(2):DimOutEnd(2), DimOutStart(3):DimOutEnd(3), DimOutStart(4))
  elseif(DimOutSize(3)==1) then
    !write(*,*) "3 size 1"
    !write(*,*) "DataOut(1:",SizeX,",1:",SizeY,",1:",SizeT,")=DataTmp(",DimOutStart(1),":",DimOutEnd(1),",",DimOutStart(2),":",DimOutEnd(2),",",DimOutStart(3),",",DimOutStart(4),":",DimOutEnd(4),")"
    DataOut(1:SizeX,1:SizeY,1:SizeT)=DataTmp(DimOutStart(1):DimOutEnd(1), DimOutStart(2):DimOutEnd(2), DimOutStart(3), DimOutStart(4):DimOutEnd(4))
  elseif(DimOutSize(2)==1) then
    !write(*,*) "2 size 1"
    DataOut(1:SizeX,1:SizeY,1:SizeT)=DataTmp(DimOutStart(1):DimOutEnd(1), DimOutStart(2), DimOutStart(3):DimOutEnd(3), DimOutStart(4):DimOutEnd(4))
  elseif(DimOutSize(1)==1) then
    !write(*,*) "1 size 1"
    DataOut(1:SizeX,1:SizeY,1:SizeT)=DataTmp(DimOutStart(1), DimOutStart(2):DimOutEnd(2), DimOutStart(3):DimOutEnd(3), DimOutStart(4):DimOutEnd(4))
  else
    write(*,*) "Problem with ExtractUsefullDimension"
    stop
  endif

  deallocate(DataTmp)
END SUBROUTINE ExtractUsefullDimension



SUBROUTINE MakeBox(tlon, tlat, lon_min, lon_max, lat_min, lat_max, Mask)! Re-writen by AM (05/2012)
  IMPLICIT NONE
  double precision ::  lon_min, lon_max, lat_min, lat_max
  integer, DIMENSION(:,:), ALLOCATABLE :: Mask
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tlon, tlat
  WHERE (tlon .GT. 360.)  ! Due to CLIO grid which doesn't varies
    tlon=tlon-360.       ! between 0 to 360 degrees of LONGITUDE
  END WHERE
  IF(lon_max.le.lon_min) THEN
    WHERE((tlat(:,:).ge.lat_min).AND.(tlat(:,:).le.lat_max).AND.(tlon(:,:).ge.lon_min))
      Mask(:,:)=1
    end where
    WHERE((tlat(:,:).ge.lat_min).AND.(tlat(:,:).le.lat_max).AND.(tlon(:,:).le.lon_max))
      Mask(:,:)=1
    end where
  ELSE
    WHERE((tlat(:,:).ge.lat_min).AND.(tlat(:,:).le.lat_max).AND.(tlon(:,:).ge.lon_min).AND.(tlon(:,:).le.lon_max))
      Mask(:,:)=1
    end where
  ENDIF
END SUBROUTINE MakeBox





SUBROUTINE AverageBoxAllocateMod(data, Tlat, latproxy, lonproxy, latmin, latmax, lonmin, lonmax)
  IMPLICIT NONE

  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: data ! (nlon,nlat,ntmaxs)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tlat  ! (nlon,nlat)
  REAL(kind=4) :: somme, moyenne
  INTEGER :: t, i, j, nbbox
  INTEGER, DIMENSION(:), ALLOCATABLE :: latproxy, lonproxy, latmin, latmax, lonmin, lonmax

  write(*,*) "Compute Average for VIAU and DAVIS box. Allocate the average to the proxy grid"

  DO nbbox=1,size(latproxy),1
    write(*,*) latproxy(nbbox), lonproxy(nbbox), latmin(nbbox), latmax(nbbox), lonmin(nbbox), lonmax(nbbox)

    do t=1, size(data,3)
      moyenne=0.
      somme=0.
      if(lonmin(nbbox).gt.lonmax(nbbox)) then
        do i=lonmin(nbbox), size(data,1)
          do j=latmin(nbbox), latmax(nbbox)
            moyenne=moyenne+(data(i,j,t)*Tlat(i,j))
            somme=somme+Tlat(i,j)
          enddo
        enddo
        do i=1, lonmax(nbbox)
          do j=latmin(nbbox), latmax(nbbox)
            moyenne=moyenne+(data(i,j,t)*Tlat(i,j))
            somme=somme+Tlat(i,j)
          enddo
        enddo
      else
        do i=lonmin(nbbox), lonmax(nbbox)
          do j=latmin(nbbox), latmax(nbbox)
            moyenne=moyenne+(data(i,j,t)*Tlat(i,j))
            somme=somme+Tlat(i,j)
          enddo
        enddo
      endif

      moyenne=moyenne/somme

      ! Re-assignation de la valeur moyenne a la lat lon du proxy
      data(lonproxy(nbbox),latproxy(nbbox),t)=moyenne

    enddo
  ENDDO

END SUBROUTINE AverageBoxAllocateMod

SUBROUTINE nanmean(vector,iilat,iilon,undef,nanmeanvalue)
  IMPLICIT NONE
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: vector
  REAL(kind=4), DIMENSION(:), ALLOCATABLE :: smallvector,futuremean
  integer, dimension(1) :: iilat, iilon
  REAL(kind=4) :: undef,nansomme,nanmeanvalue

  allocate(smallvector(size(vector,3)))
  allocate(futuremean(size(vector,3)))
  smallvector=RESHAPE(vector(iilon,iilat,:),(/size(vector,3)/))
  where(smallvector.eq.undef)
    futuremean = 0
  elsewhere
    futuremean = 1
  end where

  nansomme=SUM(smallvector, MASK=smallvector.ne.undef)
  nanmeanvalue=nansomme/sum(futuremean)

  deallocate(smallvector)
  deallocate(futuremean)

END SUBROUTINE nanmean

SUBROUTINE AverageBoxMask(data, Tlat, latproxy, lonproxy, IndMask)
  IMPLICIT NONE

  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: data ! (nlon,nlat,ntmaxs)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tlat  ! (nlon,nlat)
  DOUBLE PRECISION :: Deg_rad
  Integer, DIMENSION(:,:), ALLOCATABLE :: IndMask  ! (nlon,nlat)
  REAL(kind=4) :: somme(0:50), moyenne(0:50)
  INTEGER :: t, i, j, nbbox
  INTEGER, DIMENSION(:), ALLOCATABLE :: latproxy, lonproxy



  write(*,*) "Compute Average on box and allocate to proxy grid" 
  ! write(*,*) 'IndMask',IndMask(180,152),IndMask(258,142)
  ! DO nbbox=1,size(latproxy),1
  !   write(*,*) latproxy(nbbox), lonproxy(nbbox)
  ! ENDDO
  !  Deg_rad=3.141592654/180.0
  Deg_rad=1.0

  do t=1, size(data,3)
    moyenne=0.
    somme=0.
    do i=1, size(data,1)
      do j=1, size(data,2)
        if ((IndMask(i,j).gt.0).and.(data(i,j,t).gt.-1E10)) then
          moyenne(IndMask(i,j))=moyenne(IndMask(i,j))+data(i,j,t)*cos(Tlat(i,j)*Deg_rad)
          somme(IndMask(i,j))=somme(IndMask(i,j))+cos(Tlat(i,j)*Deg_rad)
          !             write(*,*) IndMask(i,j),data(i,j,t),Tlat(i,j)
        endif
      enddo
    enddo

    ! Re-assignation de la valeur moyenne a la lat lon du proxy
    DO nbbox=1,size(latproxy),1
      moyenne(nbbox)=moyenne(nbbox)/somme(nbbox)
      write(*,*) nbbox,moyenne(nbbox),Tlat(1,nbbox)
      data(lonproxy(nbbox),latproxy(nbbox),t)=moyenne(nbbox)
    ENDDO
  enddo
END SUBROUTINE AverageBoxMask
