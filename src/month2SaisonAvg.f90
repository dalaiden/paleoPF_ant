subroutine average(undef, inputtab, outputtab, timedimidx)
  IMPLICIT NONE
  
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: inputtab
  REAL, DIMENSION(:,:), ALLOCATABLE :: outputtab
  INTEGER :: timedimidx
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
      do t=1,(size(inputtab,timedimidx))
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
END SUBROUTINE  average


subroutine averageNew(undef, obs, mod, outputtab, timedimidx)
  IMPLICIT NONE
  
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: mod
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: obs
  REAL, DIMENSION(:,:), ALLOCATABLE :: outputtab
  INTEGER :: timedimidx, nbpoint, nbpointtotalutilise
  REAL :: undef, currentAccu, oldvalue

  INTEGER :: i,j,t,maxi,maxj
 
    
  maxi=-1
  do i=1,3
    if(i /= timedimidx .and. maxi == -1) maxi=size(mod,i)
    if(i /= timedimidx .and. maxi /= -1) maxj=size(mod,i)
  enddo

  do i=1,maxi
    do j=1,maxj
      outputtab(i,j)=0.0
      oldvalue=obs(i,j,1)
      nbpoint=0
      nbpointtotalutilise=0
      currentAccu=0.0
      do t=1,size(obs,timedimidx)
        if((obs(i,j,t).ne.oldvalue).or.(t.eq.size(obs,timedimidx))) then
          if(nbpoint.ne.0) then
            outputtab(i,j)=outputtab(i,j)+currentAccu/real(nbpoint)
            nbpointtotalutilise=nbpointtotalutilise+1
          endif
          nbpoint=0
          oldvalue=obs(i,j,t)
          currentAccu=0.0
        endif

        if(obs(i,j,t).ne.undef) then
          nbpoint=nbpoint+1
          currentAccu=currentAccu+abs(obs(i,j,t)-mod(i,j,t))
        endif

        if((t.eq.size(obs,timedimidx)).and.(nbpoint.ne.0)) then
          outputtab(i,j)=outputtab(i,j)+currentAccu/real(nbpoint)
          nbpointtotalutilise=nbpointtotalutilise+1
        endif
      enddo
      if(nbpointtotalutilise.eq.0) then
        outputtab(i,j)=undef
      else
        outputtab(i,j)=outputtab(i,j)/real(nbpointtotalutilise)
      endif

    enddo
  enddo
END SUBROUTINE  averageNew

