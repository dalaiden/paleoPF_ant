subroutine average_mod(undef, inputtab, outputtab, timedimidx, first, last)
  IMPLICIT NONE
  
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: inputtab
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: outputtab
  INTEGER :: timedimidx, first, last
  REAL :: undef
  
  INTEGER :: i,j,t,k,maxi,maxj
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nbpoint
  
  if(size(inputtab,timedimidx) < 12) then
    write(*,*) "time length is too short ", size(inputtab,timedimidx), " < 12"
    stop
  endif
  
  maxi=-1
  do i=1,3
    if(i /= timedimidx .and. maxi == -1) maxi=size(inputtab,i)
    if(i /= timedimidx .and. maxi /= -1) maxj=size(inputtab,i)
  enddo

  allocate(nbpoint(maxi,maxj))
  outputtab=0.0
  do t=1,(size(inputtab,timedimidx)/12)
    nbpoint=0
    do i=1,maxi
      do j=1,maxj
        do k=first,last
          if(inputtab(i,j,k+(t-1)*12).ne.undef) then
            outputtab(i,j,t)=outputtab(i,j,t)+inputtab(i,j,k+(t-1)*12)
	    nbpoint(i,j)=nbpoint(i,j)+1
          endif
	enddo
        if(nbpoint(i,j).ne.0) then
         outputtab(i,j,t)=outputtab(i,j,t)/real(nbpoint(i,j))
        else
         outputtab(i,j,t)=undef
        endif
      enddo
    enddo
  enddo
END SUBROUTINE  average_mod

subroutine average_ref(undef, inputtab, outputtab, timedimidx, monthF, monthL)
  IMPLICIT NONE
  
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: inputtab
  REAL, DIMENSION(:,:), ALLOCATABLE :: outputtab
  INTEGER :: timedimidx,monthF,monthL
  REAL :: undef

  INTEGER :: i,j,t,maxi,maxj
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nbpoint
    
  maxi=-1
  do i=1,3
    if(i /= timedimidx .and. maxi == -1) maxi=size(inputtab,i)
    if(i /= timedimidx .and. maxi /= -1) maxj=size(inputtab,i)
  enddo
   
  
  allocate(nbpoint(maxi,maxj))
  outputtab=0.0
  nbpoint=0
  do i=1,maxi
    do j=1,maxj
      do t=monthF,monthL
        if(inputtab(i,j,t).ne.undef) then
          outputtab(i,j)=outputtab(i,j)+inputtab(i,j,t)
          nbpoint(i,j)=nbpoint(i,j)+1
        endif  
      enddo
      if(nbpoint(i,j).ne.0) then
        outputtab(i,j)=outputtab(i,j)/real(nbpoint(i,j))
      else
        outputtab(i,j)=undef
      endif
    enddo
  enddo
END SUBROUTINE  average_ref

