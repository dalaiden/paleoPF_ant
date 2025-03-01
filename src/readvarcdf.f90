subroutine readvarcdf(filename, varName, pvalue, pdimName, pdimSize, missval, dimtolimit, limit)
  USE NETCDF
  IMPLICIT NONE
  
  REAL, DIMENSION(:), POINTER :: pvalue
  CHARACTER(len=*) :: varName, filename
  CHARACTER(len=256), DIMENSION(:), POINTER :: pdimName
  INTEGER, DIMENSION(:), POINTER :: pdimSize
  
  INTEGER, DIMENSION(:), OPTIONAL :: limit
  CHARACTER(len=*), OPTIONAL :: dimtolimit
  REAL, OPTIONAL :: missval
  
  INTEGER :: ncid, status
  INTEGER :: nbDim, unlimDimID, varID, nbvarDim, dimtolimitID
  INTEGER, DIMENSION(nf90_max_var_dims) :: varDimID

  INTEGER :: i
  INTEGER, DIMENSION(:), allocatable :: start
  
  status = nf90_open(filename, nf90_nowrite, ncid)
  if(status /= nf90_NoErr) then
    write(*,*) "error 1", NF90_STRERROR(status)
    stop
  endif
    
  status = nf90_inquire(ncid, nbDim, unlimitedDimId = unlimDimID)
  if(status /= nf90_NoErr) then
    write(*,*) "error 2", NF90_STRERROR(status)
    stop
  endif
  !write(*,*) "File have ", nbDim, " dimension(s)"
  !write(*,*) "Unlimited dimension id is ", unlimDimID

     
  status = nf90_inq_varid(ncid, varName, varID)
  if(status /= nf90_NoErr) then
    write(*,*) "  error 4", NF90_STRERROR(status)
    stop
  endif
  !write(*,*)trim(varName), " id is ", varID
  
  status = nf90_inquire_variable(ncid, varID, ndims = nbvarDim, dimids = varDimID)
  if(status /= nf90_NoErr) then
    write(*,*) "  error 5", NF90_STRERROR(status)
    stop
  endif
  
  if(present(missval)) status = nf90_get_att(ncid, varID, "missing_value", missval)

  allocate(pdimName(nbvarDim))
  allocate(pdimSize(nbvarDim))
  allocate(start(nbvarDim))
  start=1
  dimtolimitID=-1
  do i=1,nbvarDim
    status = nf90_inquire_dimension(ncid, varDimID(i), name = pdimName(i), len = pdimSize(i))
    if(status /= nf90_NoErr) then
      write(*,*) "   error 3", NF90_STRERROR(status)
      stop
    endif
    
    if(pdimName(i) == dimtolimit) dimtolimitID=i
  enddo


  if(present(limit)) then
    if(dimtolimitID == -1) then
      write(*,*) "   Can't find dimension to limit: ", "'", dimtolimit, "'"
      stop
    else
      write(*,'(A,A,A,I1,A,I5,A,I5)') "   limit dimension = ", trim(dimtolimit), "(", dimtolimitID, ") between ", limit(1), " and ", limit(1)+limit(2)-1
      pdimSize(dimtolimitID)=limit(2)
      start(dimtolimitID)=limit(1)
    endif
  endif
    
  allocate(pvalue(product(pdimSize)))
  status = nf90_get_var(ncid, varID, pvalue, start = start, count = pdimSize)
  if(status /= nf90_NoErr) then
    write(*,*) "   error 6", NF90_STRERROR(status)
    write(*,*) "   start=", start
    write(*,*) "   pdimSize=", pdimSize
    stop
  endif
  
  status = nf90_close(ncid);
  if(status /= nf90_NoErr) then
    write(*,*) "   error 7", NF90_STRERROR(status)
    stop
  endif

END SUBROUTINE  readvarcdf
