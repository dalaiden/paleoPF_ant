  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Use exemple :
  ! ./PartFilter /home/astr/sallaz/Repository/trunk/LOVECLIM/wkdir/co_?/ ic001850_000 3 FAP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   function intYO(X)
!     REAL :: X, espylon
!     INTEGER :: intYO, intX
!     parameter(espylon=1E-7)
!     intX=int(X)
!     if(int((X+espylon)).gt.intX) then
!       intYO=intX+1
!     else
!       intYO=intX
!     endif
!     return
!   end function intYO

  PROGRAM Filtre_a_particules
    USE ran_mod
    IMPLICIT NONE

    
    !Variable definition
    REAL, DIMENSION(:), allocatable :: Likelyhood, prob, niq
    INTEGER, DIMENSION(:), allocatable :: rsam, Weight, newWeight, TargWeight
    INTEGER :: nbEnsemble, numExp, i, j
    CHARACTER(len=10), DIMENSION(:), allocatable :: NameExp, NewNameExp
    CHARACTER(len=10) :: parm, methode, mois_start, mois_end
    CHARACTER(len=20) :: numExpTxt
    CHARACTER(len=255) :: expdir, seedtxt
    CHARACTER(len=300) :: fcostFileName, mois_start_FileName, mois_end_FileName
!     CHARACTER(len=300) :: weightFileName
    LOGICAL::fichierexiste
    REAL :: sumprob
    
    !Get exe arguments
    call getarg(1,expdir)
        
    call getarg(3,parm)
    read(parm,*)nbEnsemble
    
    call getarg(4,methode)

    call getarg(5,seedtxt)
    
    allocate(TargWeight(nbEnsemble))
    allocate(newWeight(nbEnsemble))
    allocate(Weight(nbEnsemble))
    allocate(Likelyhood(nbEnsemble))
    allocate(NameExp(nbEnsemble))
    allocate(NewNameExp(nbEnsemble))

    !Read cost and weight
    i=1
    do numExp = 1, 1+nbEnsemble-1

      write(numExpTxt,*) numExp
      write(fcostFileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"/"//"fcostopt.dat"
      !write(weightFileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"P/"//trim(icnum)//"/weight.dat"
      inquire( file=fcostFileName, exist=fichierexiste)

      if(fichierexiste) then
        open(61,file=fcostFileName)
        read(61,*)
        read(61,*) Likelyhood(i)
        close(61)
        !open(62,file=weightFileName)
        !read(62,*) Weight(i)
        Weight(i) = 1
        !close(62)
      else
        call newwrite2char("Impossible to find : ", fcostFileName)
        Select case(methode)
          case("FAP")
            Likelyhood(i)=0
          case("Best")
            Likelyhood(i)=666
        endselect
        Weight(i)=0
      endif
      write(NameExp(i),*) trim(adjustl(numExpTxt))//""
      i=i+1
    enddo
    
    !Methode selection
    Select case(methode)
      case("FAP")
         !write(*,'(A27)',advance='no') "Methode Filtre a Particule "
         allocate(prob(nbEnsemble))
         allocate(niq(nbEnsemble))
         allocate(rsam(nbEnsemble))
     
         !call init_random_seed(seedtxt)
    
         !residual methode
         prob = Weight * Likelyhood
         sumprob=sum(prob)
         if(sumprob == 0) then
           prob = 1./nbEnsemble
         else
           prob = prob / sumprob
         endif
         niq = nbEnsemble * prob
         if(minval(Likelyhood)==maxval(Likelyhood)) then
           rsam=1
         else
           rsam = int(niq)
         endif
         !write(*,*)(niq)
         !write(*,*)(niq - rsam)
         if(any(niq - rsam > 0)) then
            rsam = rsam + rmultinom(nbEnsemble-sum(rsam), niq-rsam, size(niq-rsam))
         endif
         TargWeight=rsam   
         NewNameExp=NameExp
         
         newWeight=TargWeight         
         do i = 1, nbEnsemble
            if(newWeight(i) == 0) then
               do j = 1, nbEnsemble
                  if(newWeight(j) > 1) then
                     newWeight(j)=newWeight(j)-1
                     newWeight(i)=newWeight(i)+1
                     NewNameExp(i)=NewNameExp(j)
					 exit
                  endif
               enddo
            endif
         enddo
		 
      case("Best")
        write(*,'(A20)') "Methode Select Best "
        call sort(NameExp, Weight, Likelyhood, nbEnsemble, 'min2max')
        TargWeight=0
        TargWeight(1)=nbEnsemble
        newWeight=1
        NewNameExp=NameExp(1)
        
        
      case default
        write(*,*)"Unknow methode"
        STOP 0
    End Select

    i=1
    write(*,'(A99)')'NAME    START    END          FCOST        CURRENT WEIGHT  TARGET WEIGHT RESTART COPY    NEW WEIGHT'
    do numExp = 1, 1+nbEnsemble-1
      write(numExpTxt,*) numExp
      !write(weightFileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"P/"//trim(icnum)//"/weight.dat"
!       write(nameFileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"P/"//"name.dat" Ca sert à rien ce truc
      !open(63,file=weightFileName)
      !write(63,*) newWeight(i)
      !close(63)
!       open(64,file=nameFileName)
!       write(64,*) NewNameExp(i)
!       close(64)

      write(mois_start_FileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"/"//"mois_start.dat" 
      open(63,file=mois_start_FileName)
      read(63,*) mois_start
      close(63)
      write(mois_end_FileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"/"//"mois_end.dat" 
      open(64,file=mois_end_FileName)
      read(64,*) mois_end
      close(64)
      
      
        
!         write(fcostFileName,'(A)') trim(expdir)//trim(adjustl(numExpTxt))//"P/"//"fcostopt.dat"
! 
!         open(61,file=fcostFileName)
!         read(61,*)
!         read(61,*) Likelyhood(i)
!         close(61)


      write(*,'(A10,A7,A7,E,A5,I5,A12,I5,A10,A10,A6,I5)')adjustl(NameExp(i)),mois_start,mois_end,Likelyhood(i),' ',Weight(i),' ',TargWeight(i),' ',adjustl(NewNameExp(i)),' ',NewWeight(i)
      i=i+1
    enddo

  END PROGRAM Filtre_a_particules

  
  SUBROUTINE sort(NameExp, Weight, Likelyhood, nbEnsemble, sens)
    IMPLICIT NONE
    INTEGER :: nbEnsemble
    REAL, DIMENSION(nbEnsemble) :: Likelyhood
    INTEGER, DIMENSION(nbEnsemble) :: Weight
    CHARACTER(len=10), DIMENSION(nbEnsemble) :: NameExp
    INTEGER :: pos, tmpint
    INTEGER, DIMENSION(1):: posmax
    REAL :: tmpreal
    CHARACTER(len=10) :: tmpchar
    CHARACTER(len=7) :: sens
    LOGICAL, DIMENSION(nbEnsemble) :: mask
    mask=1
    
    do pos = 1, nbEnsemble-1
        Select case(sens)
          case('min2max')
             posmax=minloc(Likelyhood, mask)
          case('max2min')
             posmax=maxloc(Likelyhood, mask)
          case default
            write(*,*)"Unknow option ", sens
            STOP 0
        End Select
    	
    	mask(pos)=0
    	
    	tmpreal=Likelyhood(pos)
    	Likelyhood(pos)=Likelyhood(posmax(1))
    	Likelyhood(posmax(1))=tmpreal
    	
    	tmpint=Weight(pos)
    	Weight(pos)=Weight(posmax(1))
    	Weight(posmax(1))=tmpint
     	
    	tmpchar=NameExp(pos)
    	NameExp(pos)=NameExp(posmax(1))
    	NameExp(posmax(1))=tmpchar
    enddo
    
  END SUBROUTINE sort
