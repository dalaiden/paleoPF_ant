subroutine filter_scale(var4t,undef,area,lambda,lscale,nst)
  IMPLICIT NONE
       
  REAL, DIMENSION(:,:), ALLOCATABLE :: var4t !(nlon,nlat)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: area
  INTEGER :: nst
  REAL :: undef, lambda, lscale
      
  INTEGER :: dy,n,ii,jj
  REAL, DIMENSION(:,:), ALLOCATABLE :: var4to !(nlon,nlat)
  REAL, DIMENSION(:,:), ALLOCATABLE :: temp !(0:nlon+1,0:nlat+1)
  REAL :: lrel,lrel2,dx,fx1,fx2,fy1,fy2,totalinc

  INTEGER :: nlat  ! number of latitudes
  INTEGER :: nlon  ! number of longitudes

  nlon=size(var4t,1)
  nlat=size(var4t,2)
!   write(*,*) "YOOOOOOOOOOOOOOOOOOO size=", nlon,nlat
  allocate(var4to(nlon,nlat))
  allocate(temp(0:nlon+1,0:nlat+1))
       

! note lambda mean eq diffus*pseudo-time
  dy=1.
!  radian=3.141592654/180.0
  lrel=lscale/621600.0
  lrel2=lrel*lrel
  var4to=var4t
  DO n=1,nst
      do ii=1,nlon
        do jj=1,nlat
          temp(ii,jj)=var4t(ii,jj)
        enddo
      enddo
      do jj=1,nlat
        temp(0,jj)=var4t(nlon,jj)
        temp(nlon+1,jj)=var4t(1,jj)
      enddo
      do ii=1,nlon/2
        temp(ii,0)=temp(nlon/2+ii,1)
        temp(ii,nlat+1)=temp(nlon/2+ii,nlat)
      enddo
      do ii=1+nlon/2,nlon
        temp(ii,0)=temp(nlon+1-ii,1)
        temp(ii,nlat+1)=temp(nlon+1-ii,nlat)
      enddo

      do ii=1,nlon
        do jj=1,nlat
          if (temp(ii,jj).ne.undef) then
            dx=max(1.0*area(ii,jj),0.2)
            if (temp(ii+1,jj).ne.undef) then 
              fx1=lambda/dx*(temp(ii+1,jj)-temp(ii,jj))
            else
              fx1=0
            endif
            if (temp(ii-1,jj).ne.undef) then
              fx2=lambda/dx*(temp(ii-1,jj)-temp(ii,jj))
            else
              fx2=0
            endif
            if (temp(ii,jj+1).ne.undef) then
              fy1=lambda/dy*(temp(ii,jj+1)-temp(ii,jj))
            else
              fy1=0
            endif
            if (temp(ii,jj-1).ne.undef) then
              fy2=lambda/dy*(temp(ii,jj-1)-temp(ii,jj))
            else
              fy2=0
            endif
            totalinc=(fx1+fx2)/dx+(fy1+fy2)/dy
            totalinc=totalinc-lambda*(temp(ii,jj)-var4to(ii,jj))/lrel2
            var4t(ii,jj)=temp(ii,jj)+totalinc
!           temp(ii,jj)=temp(ii,jj)+totalinc
          else
            var4t(ii,jj)=undef
          endif
        enddo
      enddo
  ENDDO ! end of do loop over n diffusion pseudo time steps

  return

 end
