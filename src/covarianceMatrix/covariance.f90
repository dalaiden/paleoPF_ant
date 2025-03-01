! These subroutines are used for constructing the covariance matrix
SUBROUTINE covariance(cov,state_array)
  IMPLICIT NONE
    
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: cov !(grid,grid) 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: state_array !(grid,ntmaxs)

  INTEGER ntmaxs,i,j,k,n
  REAL(kind=4) :: mean,maxy
    
  n=size(state_array,1)
  ntmaxs=size(state_array,2)
     
  write(*,*)'    |Subtracting mean from state array for cov matrix'
  ! Do I actually need to subtract the mean?
  DO i=1,n
    !Calculate the mean in time
    mean=0.
    DO j=1,ntmaxs
      mean=mean+state_array(i,j)
    ENDDO
    !Subtract the mean
    DO j=1,ntmaxs
      state_array(i,j)=state_array(i,j)-mean/ntmaxs
    ENDDO
  ENDDO
 
  write(*,*)'    |Calculating the cov matrix'
  maxy=0.     
  DO i=1,n  
    DO k=i,n
      cov(i,k)=0.
      DO j=1,ntmaxs   
	cov(i,k)=cov(i,k)+state_array(i,j)*state_array(k,j)
      ENDDO      
      if (ntmaxs.ne.1) cov(i,k)=cov(i,k)/(ntmaxs-1)
      if (abs(cov(i,k)).gt.maxy) maxy=abs(cov(i,k))
    ENDDO      
    if (i.lt.n) then
      cov(i+1:n,i)=cov(i,i+1:n) 
    endif
  ENDDO
  write(*,*)'    |Max in the cov=',maxy
  
  DO i=1,n  
    if (abs(cov(i,i)).lt.maxy) maxy=abs(cov(i,i))
  ENDDO
  write(*,*)'    |Min diag element of the cov=',maxy

END SUBROUTINE covariance

SUBROUTINE write_cov(filename,data_out,data_name)
    
    USE NETCDF
    IMPLICIT NONE

    interface
      subroutine check(status)
	INTEGER, INTENT (IN) :: status
      end subroutine check
    end interface

    CHARACTER(len=256) :: filename
    REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: data_out !(NX,NY)
    CHARACTER(len=*) :: data_name
    
    INTEGER :: NX, NY
    INTEGER :: ncid, varid, dimids(2)
    INTEGER :: x_dimid, y_dimid

    NX=size(data_out,1)
    NY=size(data_out,2)

    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    CALL check( nf90_def_dim(ncid, "X", NX, x_dimid) )
    CALL check( nf90_def_dim(ncid, "Y", NY, y_dimid) )
    dimids =  (/x_dimid,y_dimid/)
    CALL check( nf90_def_var(ncid, data_name, NF90_FLOAT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, data_out) )
    CALL check( nf90_close(ncid) )
    write(*,*) "    |Writing the cov.matrix into the file was done successfully"
END SUBROUTINe write_cov

SUBROUTINE read_cov(filename,data_in,data_name, t)
    USE NETCDF
    IMPLICIT NONE

    interface
      subroutine check(status)
	INTEGER, INTENT (IN) :: status
      end subroutine check
    end interface

    CHARACTER(len=256) :: filename
    REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: data_in !(NX,NY)
    CHARACTER(len=*) :: data_name
    integer :: t, ncid, varid
    integer, dimension(nf90_max_var_dims) :: varDimID
    integer :: nbvarDim

    if (len(trim(filename))==0) then
      data_in=0.
      write(*,*) "    |No cov.matrix file : cov=0"
    else
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, data_name, varid) )
      call check(nf90_inquire_variable(ncid, varid, ndims=nbvarDim, dimids=varDimID))
      call check( nf90_get_var(ncid, varid, data_in, (/1,1,t/), (/size(data_in,1), size(data_in,2), size(data_in,3)/) ))
      call check( nf90_close(ncid) )
      write(*,*) "    |Reading the cov.matrix from the file was done successfully"
    endif
END SUBROUTINE read_cov

SUBROUTINE check(status)
    USE NETCDF
    IMPLICIT NONE
    
    INTEGER, INTENT (IN) :: status
    if(status /= nf90_noerr) then 
      write(*,*)"    |Error in writing the cov.matrix to the file, ", trim(nf90_strerror(status))
      stop 
    end if
END SUBROUTINE check  


