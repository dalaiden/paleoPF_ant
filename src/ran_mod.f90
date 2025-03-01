module ran_mod
  integer, parameter:: b8 = selected_real_kind(14)
  real(b8), parameter :: pi = 3.141592653589793239_b8 
! module contains three functions
! ran1 returns a uniform random number between 0-1
! spread returns random number between min - max
! normal returns a normal distribution

contains
    subroutine init_random_seed(seedtxt)
        CHARACTER(len=255) :: seedtxt
        integer, dimension(:), allocatable :: seed
        integer :: i, nbgraine

        nbgraine=1
        do i=1,len(seedtxt)
          if(seedtxt(i:i)==",") nbgraine=nbgraine+1
        enddo
        allocate(seed(nbgraine))
        read(seedtxt,*) seed
        if(seed(1)==-1) then
          deallocate(seed)
          call random_seed
          call random_seed(SIZE = i)
          allocate(seed(i))
          call random_seed(GET = seed)
        else
          call random_seed(PUT = seed)
        endif
        write(*,*) "Seed : ",seed
    end subroutine init_random_seed

    function ran1()  !returns random number between 0 - 1
        implicit none
        real(b8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
    end function ran1

    function rmultinom(N, prob, sizeProb)
        !return a vector of length K of integers in 0:N
        !N specifying the total number of objects that are put into K boxes in the typical multinomial experiment
        !prob is a numeric non-negative vector of length K, specifying the probability for the K classes; is internally normalized to sum 1.
        !sizeProb is the size of prob vector
        INTEGER :: N, sizeProb, i, j
        REAL(b8) :: alea
	REAL, dimension(:) :: prob
        REAL, dimension(sizeProb) :: probSum, rmultinom
        
        
        do j=1, sizeProb
          if(j>1) then
            probSum(j) = probSum(j-1)+(prob(j)/sum(prob))
          else
            probSum(j) = (prob(j)/sum(prob))
          endif
        enddo
        
        rmultinom = 0
        do i=1,N
          alea = ran1()
          do j=1,sizeProb
            if(j>1) then
	      if(alea<probSum(j) .AND. alea>probSum(j-1)) then
	        rmultinom(j) = rmultinom(j)+1
	        exit
	      endif
	    else
	      if(alea<probSum(j)) then
	        rmultinom(j) = rmultinom(j)+1
	        exit
	      endif
	    endif   
          enddo
        enddo
    end function rmultinom

    function spread(min,max)  !returns random number between min - max
        implicit none
        real(b8) spread
        real(b8) min,max
        spread=(max - min) * ran1() + min
    end function spread

    function normal(mean,sigma) !returns a normal distribution
        implicit none
        real(b8) normal,tmp
        real(b8) mean,sigma
        integer flag
        real(b8) fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0_b8
            do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
                r1=2.0_b8*ran1()-1.0_b8
                r2=2.0_b8*ran1()-1.0_b8
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0_b8*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
    end function normal

end module ran_mod 
