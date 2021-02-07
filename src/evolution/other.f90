subroutine eliminateLoci(nanim,genotypes,nloci,nblock,istore,eliminated,oldnumber,nremained,newnblock)
! eliminate loci from the genotype array.  The loci to be eliminate are given in a list
! after it returns the old position of the new reordered loci.  This is beacuse the
!  sampling recombined gamete need to consider that loci are not equidistance anymore

  implicit none

  INTEGER                            , INTENT(IN)     :: istore  ! how SNP are stored

  INTEGER                            , INTENT(IN)     :: nanim
  INTEGER                            , INTENT(IN)     :: nloci,nblock     !nblock and nloci are the same if storage is as integer
  INTEGER          , DIMENSION(:,:,:), INTENT(INOUT)  :: genotypes
  INTEGER          , DIMENSION(:)    , INTENT(IN)     :: eliminated
  INTEGER          , DIMENSION(:)    , INTENT(OUT)    :: oldnumber
  INTEGER                            , INTENT(OUT)    :: nremained
  INTEGER                            , INTENT(OUT)    :: newnblock

  INTEGER :: i,j,id,iloci
  integer :: iblck1,iblck2,ibit1,ibit2,nbits
  DOUBLE PRECISION :: itot

  nbits=32 ! 32 bits  in an integer variable

  oldnumber(1:nloci)=0

  IF(istore == 1) then  ! todo: read how elim. works for binary storage
    iblck1 =0  !current block being used to pass loci
    iblck2 =0  !current block being chek
    ibit1  =0  !last position used to pass a loci
    ibit2  =0  !last SNP position
    iloci=0
    do i=1,nloci
      ibit2=ibit2-1                  ! position of the new loci to be check
      if(ibit2 <0) then
        iblck2 = iblck2 + 1          ! a new block
        ibit2  = nbits  - 1          ! start from the first bit
      endif
      if(eliminated(i) == 0) then    ! the loci will be kept.
        iloci=iloci+1                ! counting
        ibit1=ibit1-1                ! new position of the loci to be kept
        if(ibit1 <0) then
          iblck1 = iblck1 + 1        ! a new block
          ibit1  = nbits  - 1        ! start from the first bit
        endif

        do id=1,nanim
            if(     Btest(genotypes(id,1,iblck2),ibit2) ) genotypes(id,1,iblck1)= ibset(genotypes(id,1,iblck1),ibit1)
            if(.not.Btest(genotypes(id,1,iblck2),ibit2) ) genotypes(id,1,iblck1)= ibclr(genotypes(id,1,iblck1),ibit1)

            if(     Btest(genotypes(id,2,iblck2),ibit2) ) genotypes(id,2,iblck1)= ibset(genotypes(id,2,iblck1),ibit1)
            if(.not.Btest(genotypes(id,2,iblck2),ibit2) ) genotypes(id,2,iblck1)= ibclr(genotypes(id,2,iblck1),ibit1)
        enddo
        oldnumber(iloci)        = i            !this SNP in new position <iloci> used to be in position <i> (before elimination)
      endif
    enddo
    nremained=iloci
    newnblock=iblck1

  ! the remainig bit which were left unused in the last block
    ibit1 =ibit1-1  !the first unused bit
    IF(ibit1 >= 0 ) then
      do WHILE (ibit1 >=0)
        do id=1,nanim
          genotypes(id,1,iblck1)= IBCLR(genotypes(id,1,iblck1),ibit1)   !make them homozygote 00
          genotypes(id,2,iblck1)= IBCLR(genotypes(id,2,iblck1),ibit1)   !make them homozygote 00
        enddo
        ibit1 =ibit1-1
      ENDDO
    endif

  endif

!===========================================================================================
!if storage is done as integer
!===========================================================================================
  IF(istore==0) then
    iloci=0
    do i=1,nloci
      if(eliminated(i) == 0) then    ! the loci will be kept.
        iloci=iloci+1                ! counting
        do id=1,nanim
            genotypes(id,1,iloci)= genotypes(id,1,i)
            genotypes(id,2,iloci)= genotypes(id,2,i)
        enddo
        oldnumber(iloci)        = i            !this SNP in new position <iloci> used to be in position <i> (before elimination)
      endif
    enddo
    nremained=iloci
    newnblock=iloci
    IF(iloci < nloci) then    !  zeroeing the reming unused blocks
      do i=(iloci+1),nloci
        do id=1,nanim
          genotypes(id,1,i)= 0
          genotypes(id,2,i)= 0
        enddo

      end do

    endif
  endif



  return
END subroutine

!=============================================================================
subroutine calcPedNRM(nanim,pedigree,ParAmat,Amat)
!---------------------------------------------------------------------------------------------
! it calculates the relationship between all individuals in one generation 
! it does not calculate the relationship between the parents(ancestors) and offspring generation
! the subroutine is useful when assuming a discrete generation and relationship of
! individual from previous generations are not needed (therefore note saved)
!
! it assumes that the parents NRM (input) and the offspring one (output) are 
! half stored (lower) diagonal matrix
!----------------------------------------------------------------------------------------------

  USE trsm_module

  implicit none


  INTEGER                          , INTENT(IN)  :: nanim
  INTEGER          , DIMENSION(:,:), INTENT(IN)  :: Pedigree
  DOUBLE PRECISION , DIMENSION(:  ), INTENT(IN)  :: ParAmat
  DOUBLE PRECISION , DIMENSION(:  ), INTENT(OUT) :: Amat

  INTEGER :: id1,id2,ipar1,ipar2,ipos1,ipos2,ipos,i,j,k

  DOUBLE PRECISION:: val1,val2,val3

  Amat(:)=0
  do id1=1,nanim
    ipos=lowerPos(id1,1)

     !looping over the two parents of id1
    do i=1,2
      ipar1=pedigree(id1,i) ! the parent to be considered (if it is not missing i.e. base animal)

      IF(ipar1>0) then
        !looping over the two parents of id2
        do j=1,2
          ipos1=ipos
          do id2=1,(id1-1),1
            ipar2=pedigree(id2,j)
            IF(ipar2>0) then
              ipos2=lowerPos(ipar1,ipar2)
              Amat(ipos1)= Amat(ipos1)+ParAmat(ipos2)/4.d0  !the relationship due to the ipar side
            endif
            ipos1 = ipos1+1
          enddo
        end do
      endif

    end do

    !the diagonal of id1
    ipos1=(id1+1)*id1/2
    Amat(ipos1)= 1.d0
    IF(pedigree(id1,1) > 0 .AND. pedigree(id1,2) > 0) then
      ipos2       = lowerPos(pedigree(id1,1), pedigree(id1,2))  !position storing the relationship between both parents
      Amat(ipos1) = 1.d0 + ParAmat(ipos2)/2.d0      !diag(id1,id1)=1+f
    endif

  end do

  !open ( 20, file = "pedA.txt", status = 'unknown' )
  !open ( 21, file = "pedB.txt", status = 'unknown' )
  !val1=Amat(1)
  !ipos=1
  !write ( 20, * ) 1, 1, aMAT( 1 )
  !write ( 21, * ) 1, 1, ParaMAT( 1 )
  !do id1 = 2, nanim
  !    do id2 = 1, (id1-1),1
  !        ipos=ipos+1
  !        val1=val1+Amat(ipos)*2
  !        write ( 20, * ) ID1, ID2, aMAT( IPOS )
  !        write ( 21, * ) ID1, ID2, ParaMAT( IPOS )
  !    enddo
  !    ipos = ipos + 1
  !    val1 = val1 + Amat( ipos ) 
  !    write ( 20, * ) ID1, ID1, aMAT( IPOS )
  !    write ( 21, * ) ID1, ID1, ParaMAT( IPOS )
  !enddo

  !call calcMeanCoancestry(nanim, Amat,val2,val3,i)
  !write ( 20, * ) "# ", val1 / dble( nanim * nanim ), val2 / 2.D0, val3 - 1.D0

!open(20,file="pedF.txt", status='unknown', position ='append') 
!write(20,*)val1/dble(nanim*nanim),val2/2.d0, val3-1.d0
!close(20)
!  val1=0.d0
!  val2=0.d0
!  k=0
!  ipos=0
!do id1=1,nanim
!do id2=1,(id1-1),1
!ipos=ipos+1
!val2=val2+Amat(ipos)*2
!k=k+1
!
!IF(Amat(ipos) > 2.d0 .OR. Amat(ipos) .lt.0.d0)then
!WRITE(*,*)' wrong off diagonal ',id1,id2,ipos,k,Amat(ipos)
!WRITE(*,*)'pedigree',id1, pedigree(id1,1),pedigree(id1,2)
!ipos1=lowerPos(pedigree(id1,1),id2)
!ipos2=lowerPos(pedigree(id1,2),id2)
!WRITE(*,*)Amat(ipos),ParAmat(ipos1),ParAmat(ipos2)
!endif
!
!
!
!
!
!end do
!ipos=ipos+1
!k=k+1
!val2=val2+Amat(ipos)
!val1=val1+Amat(ipos)
!
!IF(Amat(ipos) > 2.d0 .OR. Amat(ipos) .lt.1.d0)then
!WRITE(*,*)' wrong diagonal ',id1,ipos,Amat(ipos)
!WRITE(*,*)'pedigree',id1, pedigree(id1,1),pedigree(id1,2),ipos
!ipos1=lowerPos(pedigree(id1,1),pedigree(id1,2))
!WRITE(*,*)Amat(ipos),ParAmat(ipos1)
!endif
!enddo
!
!  WRITE(*,'(a30,i,2f,2i)') ' mean diagonal NRM', 0,val1/DBLE(nanim),val2/DBLE(nanim*nanim),k,nanim
!
!  call calcMeanCoancestry(nanim, Amat,val2,val1,i)
!  WRITE(*,'(a30,i,2f,2i)') ' mean diagonal NRM', 1,val1,val2,k,nanim
!
  return

END subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcNRM(nanim,genotypes,nloci,nblock,istore,Amat)
!---------------------------------------------------------------------------------------------
! it calculates the relationship between all individuals using genotypes
!
! genotype
! it assumes that there is no missnig genotypes



  implicit none

  INTEGER                            , INTENT(IN)     :: istore  ! how SNP are stored

  INTEGER                            , INTENT(IN)     :: nanim
  INTEGER                            , INTENT(IN)     :: nloci,nblock     !nblock and nloci are the same if storage is as integer
  INTEGER          , DIMENSION(:,:,:), INTENT(INOUT)  :: genotypes

  DOUBLE PRECISION , DIMENSION(:)    ,INTENT(OUT) :: Amat

  INTEGER                             :: i,j,id,iloci,a1,a2,ihap
  integer                             :: iblck1,iblck2,ibit1,ibit2,nbits
  DOUBLE PRECISION                    :: itot,p1,p2,q1,q2,x11
  CHARACTER (LEN=2 ), dimension (2,2) :: genocode
  CHARACTER (LEN=30)                  :: formatoid,formatosnp

  DOUBLE PRECISION , DIMENSION(3,3) :: IBS
  INTEGER, DIMENSION(:), ALLOCATABLE :: igenotype

ALLOCATE(igenotype(nanim))


  IBS(1,1)=2.d0
  IBS(1,2)=1.d0
  IBS(1,3)=0.d0

  IBS(2,1)=1.d0
  IBS(2,2)=1.d0
  IBS(2,3)=1.d0

  IBS(3,1)=0.d0
  IBS(3,2)=1.d0
  IBS(3,3)=2.d0
  i=(nanim)*(nanim+1)/2
  Amat(1:i)=0.d0


  nbits=32 ! 32 bits  in an integer variable

100 FORMAT(a3)



  iblck1 =0        !block for loci i
  ibit1  =0        !bit for loci i
  do iloci=1,nloci
    IF(istore == 1 ) then
      ibit1=ibit1-1                  ! position of the new loci i to be check
      if(ibit1 <0 ) then
       iblck1 = iblck1 + 1          ! a new block
       ibit1  = nbits  - 1          ! start from the first bit
      endif

      do id=1,nanim
         a1=1 ! first allele is clear bit
         a2=1 ! second allele is clear bit
         if(Btest(genotypes(id,1,iblck1),ibit1) ) a1=2    ! the
         if(Btest(genotypes(id,2,iblck1),ibit1) ) a2=2
         igenotype(id) = a1+a2-1
      enddo
    ELSEIF(istore==0) then
     do id=1,nanim
       a1=genotypes(id,1,iloci)
       a2=genotypes(id,2,iloci)
       igenotype(id) = a1+a2-1
     end do

    endif

    j=0
    do id=1,nanim
      do i=1,id
        j=j+1
        Amat(j)=Amat(j)+IBS(igenotype(id),igenotype(i))
      end do
    end do
  enddo
  i=(nanim)*(nanim+1)/2
  Amat(1:i)=Amat(1:i)/DBLE(nloci)

  DEALLOCATE(igenotype)
  return
END subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcNRMselected(nanim,genotypes,nloci,nblock,istore,Amat,usedloci)
!---------------------------------------------------------------------------------------------
! it calculates the relationship between all individuals using genotypes
!
! genotype
! it assumes that there is no missnig genotypes



  implicit none

  INTEGER                            , INTENT(IN)     :: istore  ! how SNP are stored

  INTEGER                            , INTENT(IN)     :: nanim
  INTEGER                            , INTENT(IN)     :: nloci,nblock     !nblock and nloci are the same if storage is as integer
  INTEGER          , DIMENSION(:,:,:), INTENT(INOUT)  :: genotypes

  DOUBLE PRECISION , DIMENSION(:)    ,INTENT(OUT) :: Amat
  INTEGER       , DIMENSION(:), INTENT(IN)  :: usedloci
  INTEGER                             :: i,j,id,iloci,a1,a2,ihap
  integer                             :: iblck1,iblck2,ibit1,ibit2,nbits
  DOUBLE PRECISION                    :: itot,p1,p2,q1,q2,x11
  CHARACTER (LEN=2 ), dimension (2,2) :: genocode
  CHARACTER (LEN=30)                  :: formatoid,formatosnp

  DOUBLE PRECISION , DIMENSION(3,3) :: IBS
  INTEGER, DIMENSION(:), ALLOCATABLE :: igenotype
  integer :: n
ALLOCATE(igenotype(nanim))


  IBS(1,1)=2.d0
  IBS(1,2)=1.d0
  IBS(1,3)=0.d0

  IBS(2,1)=1.d0
  IBS(2,2)=1.d0
  IBS(2,3)=1.d0

  IBS(3,1)=0.d0
  IBS(3,2)=1.d0
  IBS(3,3)=2.d0
  i=(nanim)*(nanim+1)/2
  Amat(1:i)=0.d0


  nbits=32 ! 32 bits  in an integer variable

100 FORMAT(a3)

  n=0
  iblck1 =0        !block for loci i
  ibit1  =0        !bit for loci i
  do iloci=1,nloci
      ibit1=ibit1-1                  ! position of the new loci i to be check
      if(ibit1 <0 ) then
       iblck1 = iblck1 + 1          ! a new block
       ibit1  = nbits  - 1          ! start from the first bit
      endif

    if(usedLoci(iloci) == 1) then
    n=n+1
    IF(istore == 1 ) then
      do id=1,nanim
         a1=1 ! first allele is clear bit
         a2=1 ! second allele is clear bit
         if(Btest(genotypes(id,1,iblck1),ibit1) ) a1=2    ! the
         if(Btest(genotypes(id,2,iblck1),ibit1) ) a2=2
         igenotype(id) = a1+a2-1
      enddo
    ELSEIF(istore==0) then
     do id=1,nanim
       a1=genotypes(id,1,iloci)
       a2=genotypes(id,2,iloci)
       igenotype(id) = a1+a2-1
     end do

    endif

    j=0
    do id=1,nanim
      do i=1,id
        j=j+1
        Amat(j)=Amat(j)+IBS(igenotype(id),igenotype(i))
      end do
    end do
    endif
  enddo
  i=(nanim)*(nanim+1)/2
!  Amat(1:i)=Amat(1:i)/DBLE(nloci)
  Amat(1:i)=Amat(1:i)/DBLE(n)
!WRITE(*,*)' nloci used for IBS calculation',n


  DEALLOCATE(igenotype)
  return
END subroutine


!==============================================================================
subroutine calcMeanCoancestry(nanim, amat,meanNRM,meanDiag,ifail)

  implicit none


  INTEGER                            ,INTENT(IN)  :: nanim
  DOUBLE PRECISION ,DIMENSION(:)     ,INTENT(IN ) :: Amat
  DOUBLE PRECISION                   ,INTENT(OUT) :: meanNRM
  DOUBLE PRECISION                   ,INTENT(OUT) :: meanDiag
  INTEGER                            ,INTENT(OUT) :: ifail

  DOUBLE PRECISION :: val1
  INTEGER          :: ipos,maxpos,i


  ifail=1
  maxpos =nanim*(nanim+1)/2
  val1=0.d0
  do ipos=1,maxpos
    val1=val1+Amat(ipos)
  end do
  meanNRM=val1*2.d0  ! all values time 2 because matrix half stored (but diag was added twice)

  ! but diagonal should be only time 1
  val1=0.d0
  ipos=0
  do i=1,nanim
    ipos=ipos+i
    val1=val1 +Amat(ipos)
  end do
  meanNRM  = (meanNRM-val1)/DBLE(nanim*nanim)
  meanDiag = val1/DBLE(nanim)
  ifail    = 0

  return
END subroutine

!=========================================================================
subroutine createSDP(nanim, Ainv, sex,iun,ifail)
implicit none

  INTEGER                        , INTENT(IN)  :: nanim
  INTEGER                        , INTENT(IN)  :: iun
  DOUBLE PRECISION , DIMENSION(:), INTENT(IN)  :: Ainv
  INTEGER          , DIMENSION(:), INTENT(IN)  :: sex
  INTEGER                        , INTENT(OUT)  :: ifail

  CHARACTER (LEN=100) :: formato1,formato2
  INTEGER :: n, m,ipos,i,j,iblock

  m=6

  WRITE(iun,*)nanim+1, "= nDIM"   ! variable to optimise contributions and auxliary variable
  WRITE(iun,*)m,       "= nBLOCK"
  WRITE(iun,'(7i)') nanim+1,1,1,1,1,nanim
  WRITE(iun,'(a1)',ADVANCE='no') "{"

  do i=1,nanim
  WRITE(iun,'(i1,a2)',ADVANCE='no') 0,","
  end do
  WRITE(iun,'(i1,a2)',ADVANCE='no') 1,"}"
  WRITE(iun,*)

  WRITE(iun,*)

! writing -F0
100 format(i,2x,i,2x,2i,f,a)
101 format(i,2x,i,2x,2i,i,a)

  iblock=1
  ipos=0
  do i=1,nanim
    do j = 1,i
      ipos=ipos+1
      WRITE(iun,100)0,iblock, i,j, -Ainv(ipos)  ! remmeber negative F0
    end do
  end do

  WRITE(iun,100)0,2,  1,1, -(-0.5d0)
  WRITE(iun,100)0,3,  1,1, -( 0.5d0)
  WRITE(iun,100)0,4,  1,1, -(-0.5d0)
  WRITE(iun,100)0,5,  1,1, -( 0.5d0)

  WRITE(iun,*)

! writing Fi (i=1,nanim)

  iblock=1
  i= nanim+1
  do j=1,nanim
    WRITE(iun,101)j,1,   j,i, 1
    IF(sex(j) == 1) then  ! sex 1 = male
      WRITE(iun,101)j,2,   1,1,  1 ! flag for male
      WRITE(iun,101)j,3,   1,1, -1 ! neg flag for male
    else
      WRITE(iun,101)j,4,   1,1,  1 ! flag for female
      WRITE(iun,101)j,5,   1,1, -1 ! neg flag for female
    endif
    WRITE(iun,101) j,6,   j,j,  1  ! minimum value
    WRITE(iun,*)
  end do

  i=nanim+1
  WRITE(iun,101) i,1, i,i, 2  ! minimum value
  WRITE(iun,*)


  ifail=0
  return
end subroutine


!===========================================================================
!=========================================================================
subroutine createSDPres ( nanim, Ainv, sex, iun, ifail,nres,resAinv,resF )
    implicit none

    INTEGER, INTENT ( IN ) :: nanim
    INTEGER, INTENT ( IN ) :: iun
    DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN ) :: Ainv
    INTEGER, DIMENSION ( : ), INTENT ( IN ) :: sex
    INTEGER, INTENT ( OUT ) :: ifail

    integer, intent(inout) ::nres
    DOUBLE PRECISION, DIMENSION ( :,: ), INTENT ( IN ) :: resAinv
    DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN ) :: resF

    CHARACTER ( LEN = 100 ) :: formato1, formato2
    INTEGER :: n, m, ipos, i, j,k, iblock

   if(nres <=0) nres=0
    m = 6+nres
write(*,*)' parmeters'
    WRITE ( iun, * ) nanim + 1, "= nDIM"   ! variable to optimise contributions and auxliary variable
    WRITE ( iun, * ) m, "= nBLOCK"
    WRITE ( iun, '(6i)',advance="no" ) nanim + 1, 1, 1, 1, 1, nanim
    do i=1,nres
        WRITE ( iun, '(i)', advance = "no" ) nanim+1
    enddo
    write(iun,*)
    WRITE ( iun, '(a1)', ADVANCE = 'no' ) "{"

    do i = 1, nanim
        WRITE ( iun, '(i1,a2)', ADVANCE = 'no' ) 0, ","
    end do
    WRITE ( iun, '(i1,a2)', ADVANCE = 'no' ) 1, "}"
    WRITE ( iun, * )

    WRITE ( iun, * )

! writing -F0
100 format ( i, 2x, i, 2x, 2i, f, a )
101 format ( i, 2x, i, 2x, 2i, i, a )
    write ( *, * ) ' block 1-6'

    iblock = 1
    ipos = 0
    do i = 1, nanim
        do j = 1, i
            ipos = ipos + 1
            WRITE ( iun, 100 ) 0, iblock, i, j, - Ainv( ipos )  ! remmeber negative F0
        end do
    end do

    WRITE ( iun, 100 ) 0, 2, 1, 1, - ( - 0.5D0 )
    WRITE ( iun, 100 ) 0, 3, 1, 1, - ( 0.5D0 )
    WRITE ( iun, 100 ) 0, 4, 1, 1, - ( - 0.5D0 )
    WRITE ( iun, 100 ) 0, 5, 1, 1, - ( 0.5D0 )

    WRITE ( iun, * )
!----------------
    write ( *, * ) ' block for restriction',nres
    do k = 1, nres
        iblock = 6+k
        ipos = 0
        do i = 1, nanim
            do j = 1, i
                ipos = ipos + 1
                WRITE ( iun, 100 ) 0, iblock, i, j, - resAinv( ipos,k )  ! remmeber negative F0
            end do
        end do
        WRITE ( iun, 100 ) 0, iblock, (nanim+1), (nanim+1), - 2.d0*resF( k )  ! remmeber negative F0

        WRITE ( iun, * )

    end do







! writing Fi (i=1,nanim)

    iblock = 1
    i = nanim + 1
    do j = 1, nanim
        WRITE ( iun, 101 ) j, 1, j, i, 1
        IF ( sex ( j ) == 1 ) then  ! sex 1 = male
            WRITE ( iun, 101 ) j, 2, 1, 1, 1 ! flag for male
            WRITE ( iun, 101 ) j, 3, 1, 1, - 1 ! neg flag for male
        else
            WRITE ( iun, 101 ) j, 4, 1, 1, 1 ! flag for female
            WRITE ( iun, 101 ) j, 5, 1, 1, - 1 ! neg flag for female
        end if
        WRITE ( iun, 101 ) j, 6, j, j, 1  ! minimum value
        WRITE ( iun, * )

        do k = 1, nres
            iblock=6+k
            WRITE ( iun, 101 ) j, iblock, j, i, 1
        enddo
        
        
    end do

    i = nanim + 1
    WRITE ( iun, 101 ) i, 1, i, i, 2  ! minimum value
    WRITE ( iun, * )


    ifail = 0
    return
end subroutine





