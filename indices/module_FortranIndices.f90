MODULE module_FortranIndices
! Module with the structure for the computation of 1D indices using Fortran
!   Using developments from PyNCplot (https://git.cima.fcen.uba.ar/lluis.fita/pyncplot/-/wikis/home)
!

!!!!!!!! Contents

! Structure
! index1D_from3D_d3: Subroutine to compute 1D indices using 3rd dimension from 3D arrays producing 2D 
!   arrays
! stats1D_from3D_d3: Subroutine to compute 1D statistics using 3rd dimension from 3D arrays producing 
!   2D arrays

  USE module_definitions
  USE module_basic
  USE module_fortranindices1d

  CONTAINS
  
!! STRUCTURE -- structure 

  SUBROUTINE index1D_from3D_d3(d1, d2, d3, indexn, matv, maskv, missv, opmatv2Da, opmatv2Db,          &
    opmatv2Dc, opmatv3Da, opmatv3Db, opmatv3Dc, indexv, op2DIindexa, op2DIindexb, op2DIindexc,        &
    op2DRindexa, op2DRindexb, op3DIindexa, op3DIindexb, op3DIindexc, op3DRindexa, op3DRindexb)
  ! Subroutine to compute 1D indices using 3rd dimension from 3D arrays producing 2D arrays
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: d1, d2, d3
    CHARACTER(len=*), INTENT(in)                         :: indexn
    REAL, INTENT(in)                                     :: missv
    REAL, DIMENSION(d1,d2), OPTIONAL, INTENT(in)         :: opmatv2Da, opmatv2Db, opmatv2Dc
    REAL, DIMENSION(d1,d2,d3), INTENT(in)                :: matv
    REAL, DIMENSION(d1,d2,d3), OPTIONAL, INTENT(in)      :: opmatv3Da, opmatv3Db, opmatv3Dc
    LOGICAL, DIMENSION(d1,d2), INTENT(in)                :: maskv
    REAL, DIMENSION(d1,d2), INTENT(out)                  :: indexv
    INTEGER, DIMENSION(d1,d2), OPTIONAL, INTENT(out)     :: op2DIindexa, op2DIindexb, op2DIindexc
    REAL, DIMENSION(d1,d2), OPTIONAL, INTENT(out)        :: op2DRindexa, op2DRindexb
    INTEGER, DIMENSION(d1,d2,d3), OPTIONAL, INTENT(out)  :: op3DIindexa, op3DIindexb, op3DIindexc
    REAL, DIMENSION(d1,d2,d3), OPTIONAL, INTENT(out)     :: op3DRindexa, op3DRindexb

    ! Local
    INTEGER                                              :: i, j, k
    INTEGER                                              :: Navailindexn
    INTEGER                                              :: Ia
    CHARACTER(len=Sm), DIMENSION(:), ALLOCATABLE         :: availindexn
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: I1Da, I1Db
    INTEGER, DIMENSION(:,:), ALLOCATABLE                 :: I2Da, I2Db
    REAL, DIMENSION(:,:), ALLOCATABLE                    :: R2Da, R2Db
    REAL, DIMENSION(:,:,:), ALLOCATABLE                  :: R3Da, R3Db
    CHARACTER(len=Sm)                                    :: fname
    
!!!!!!! Variables
! d1, d2, d3:  shape of the 3D array
! indexn: name of the index to compute
! matv: 3D array with the values
! missv: missing value of matv
! maskv: 2D mask whether to compute or not in the grid point (.FALSE.: compute)
! optmatv2D[a-c]: optional 2D variable to be used to compute the index
! optmatv3D[a-c]: optional 3D variable to be used to compute the index

  fname = 'index1D_from3D_d3'
  
  Navailindexn = 1
  IF (ALLOCATED(availindexn)) DEALLOCATE(availindexn)
  ALLOCATE(availindexn(Navailindexn))
  
  availindexn(1) = 'mean'
  
  indexv = fillValueR
  
  SELECT CASE (TRIM(indexn))
  
    CASE ('hcwi_hot')
      ! Heat wave
      !   Necessary variables: tasn, tasx, tasn11drunq90, tasx11drunq90
      !   to be assumed to be: maskv, optmatv3Da, optmatv2Da, optmatv2Db

      IF (.NOT.PRESENT(opmatv3Da)) THEN
        msg = "To compute '" // TRIM(indexn) // "' is necessary at least 'daily tas max', tasmax"
        CALL StopRun(msg, fname, -1)
      END IF
    
      IF (ALLOCATED(R2Da)) DEALLOCATE(R2Da)
      ALLOCATE(R2Da(d1,d2))
      IF (ALLOCATED(R2Db)) DEALLOCATE(R2Db)
      ALLOCATE(R2Db(d1,d2))
      ! Computing tasmin90 [30-year baseline period (1981-2010)]
      ! 90th percentile ('Q90') of the 330 respective temperature values in an 11-day window centred 
      !   on that day, for all years in the baseline period
      IF (.NOT.PRESENT(opmatv2Da)) THEN
        IF (ALLOCATED(R3Da)) DEALLOCATE(R3Da)
        ALLOCATE(R3Da(d1,d2,Npercents))
        
        CALL stats1D_from3D_d3(d1, d2, d3, 'percentiles', matv, maskv, missv, R2Da, R3Da)
        R2Da = R3Da(:,:,19)
      ELSE
        R2Da = opmatv2Da
      END IF
    
      ! Computing tasmax90
      IF (.NOT.PRESENT(opmatv2Da)) THEN
        IF (ALLOCATED(R3Da)) DEALLOCATE(R3Da)
        ALLOCATE(R3Da(d1,d2,Npercents))
        
        CALL stats1D_from3D_d3(d1, d2, d3, 'percentiles', matv, maskv, missv, R2Da, R3Da)
        R2Db = R3Da(:,:,19)
      ELSE
        R2Db = opmatv2Db
      END IF

      indexv = 0.
      ! Amount of values
      ALLOCATE(I2Da(d3,2))
      op2DIindexa = 0
    
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            CALL hcwi_hot(d3, matv(i,j,:), opmatv3Da(i,j,:), R2Da(i,j), R2Db(i,j), missv,             &
              op2DIindexa(i,j), I2Da)
            op3DIindexa(i,j,1:op2DIindexa(i,j)) = I2Da(1:op2DIindexa(i,j),1)
            op3DIindexb(i,j,1:op2DIindexa(i,j)) = I2Da(1:op2DIindexa(i,j),2)
          END IF
        END DO
      END DO

    CASE DEFAULT
    
      msg = "1D index '" // TRIM(indexn) // "' not ready !!"

      CALL StopRunAvail(msg, fname, Navailindexn, availindexn)
    
    END SELECT

    DEALLOCATE(availindexn)
    IF (ALLOCATED(I1Da)) DEALLOCATE(I1Da)
    IF (ALLOCATED(I1Db)) DEALLOCATE(I1Db)
    IF (ALLOCATED(I2Da)) DEALLOCATE(I2Da)
    IF (ALLOCATED(I2Db)) DEALLOCATE(I2Db)
    IF (ALLOCATED(R2Da)) DEALLOCATE(R2Da)
    IF (ALLOCATED(R2Db)) DEALLOCATE(R2Db)
    IF (ALLOCATED(R3Da)) DEALLOCATE(R3Da)
    IF (ALLOCATED(R3Db)) DEALLOCATE(R3Db)

    RETURN

  END SUBROUTINE index1D_from3D_d3

  SUBROUTINE stats1D_from3D_d3(d1, d2, d3, statn, matv, maskv, missv, statv, percents)
  ! Subroutine to compute 1D statistics using 3rd dimension from 3D arrays producing 2D arrays
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: d1, d2, d3
    CHARACTER(len=*), INTENT(in)                         :: statn
    REAL, INTENT(in)                                     :: missv
    REAL, DIMENSION(d1,d2,d3), INTENT(in)                :: matv
    LOGICAL, DIMENSION(d1,d2), INTENT(in)                :: maskv
    REAL, DIMENSION(d1,d2), INTENT(out)                  :: statv
    REAL,DIMENSION(d1,d2,Npercents), INTENT(out), OPTIONAL :: percents

    ! Local
    INTEGER                                              :: i, j, k
    INTEGER                                              :: Navailstatn
    INTEGER                                              :: Ia
    REAL                                                 :: Ra, Rb
    REAL, DIMENSION(d1,d2)                               :: R2Da, R2Db 
    CHARACTER(len=Sm), DIMENSION(:), ALLOCATABLE         :: availstatn
    CHARACTER(len=Sm)                                    :: fname
    
!!!!!!! Variables
! d1, d2, d3:  shape of the 3D array
! statn: name of the statistics to compute
! matv: 3D array with the values
! missv: missing value of matv
! maskv: 2D mask whether to compute or not in the grid point (.FALSE.: compute)

  fname = 'stats1D_from3D_d3'
  
  Navailstatn = 6
  IF (ALLOCATED(availstatn)) DEALLOCATE(availstatn)
  ALLOCATE(availstatn(Navailstatn))
  
  availstatn(1) = 'min'
  availstatn(2) = 'max'
  availstatn(3) = 'mean'
  availstatn(4) = 'mean2'
  availstatn(5) = 'stddev'
  availstatn(6) = 'percentiles'
  
  statv = fillValueR
  
  SELECT CASE (TRIM(statn))
  
    CASE ('min')
    
      statv = fillValueR
         
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            Ia = 0
            DO k=1, d3
              IF (matv(i,j,k) /= missv) THEN
                IF (matv(i,j,k) < statv(i,j)) statv(i,j) = matv(i,j,k)
              END IF
            END DO
          END IF
        END DO
      END DO
  
    CASE ('max')
    
      statv = -fillValueR
         
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            Ia = 0
            DO k=1, d3
              IF (matv(i,j,k) /= missv) THEN
                IF (matv(i,j,k) > statv(i,j)) statv(i,j) = matv(i,j,k)
              END IF
            END DO
          END IF
        END DO
      END DO
  
    CASE ('mean')
    
      statv = 0.
      ! Amount of values
      Ia = 0
    
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            Ia = 0
            DO k=1, d3
              IF (matv(i,j,k) /= missv) THEN
                statv(i,j) = statv(i,j) + matv(i,j,k)
                Ia = Ia + 1
              END IF
            END DO
            statv(i,j) = statv(i,j) / Ia
          END IF
        END DO
      END DO
  
    CASE ('mean2')
    
      statv = 0.
      ! Amount of values
      Ia = 0
    
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            Ia = 0
            DO k=1, d3
              IF (matv(i,j,k) /= missv) THEN
                statv(i,j) = statv(i,j) + matv(i,j,k)*matv(i,j,k)
                Ia = Ia + 1
              END IF
            END DO
            statv(i,j) = statv(i,j) / Ia
          END IF
        END DO
      END DO

    CASE ('stddev')
    
      statv = 0.
      R2Da = 0.
      R2Db = 0.
      
      ! Amount of values
      Ia = 0
    
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            Ia = 0
            DO k=1, d3
              IF (matv(i,j,k) /= missv) THEN
                Ra = Ra + matv(i,j,k)
                Rb = Rb + matv(i,j,k)*matv(i,j,k)
                Ia = Ia + 1
              END IF
            END DO
            Ra = Ra / Ia
            Rb = Rb / Ia
            statv(i,j) = SQRT(Rb-Ra*Ra)
          END IF
        END DO
      END DO

    CASE ('percentiles')
    
      IF (PRESENT(percents)) THEN
    
        percents = 0.
         
        DO i=1, d1
          DO j=1, d2
            IF (.NOT.maskv(i,j)) THEN
              CALL percentilesR(d3, matv(i,j,:), missv, Npercents, percents(i,j,:))
            END IF
          END DO
        END DO
        
      END IF

    CASE DEFAULT
    
      msg = "1D index '" // TRIM(statn) // "' not ready !!"

      CALL StopRunAvail(msg, fname, Navailstatn, availstatn)
    
    END SELECT

    DEALLOCATE(availstatn)

    RETURN

  END SUBROUTINE stats1D_from3D_d3

END MODULE module_FortranIndices

