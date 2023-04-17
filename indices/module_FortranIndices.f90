MODULE module_FortranIndices
! Module with the instructions for the computation of 1D indices using Fortran
!   Using developments from PyNCplot (https://git.cima.fcen.uba.ar/lluis.fita/pyncplot/-/wikis/home)
!

!!!!!!!! Contents
! Indices
! hcwi_hot: Heat Wave Index (from HCWI), which is implemented in the Copernicus European Drought 
!   Observatory (EDO)
! percentilesR: Subroutine to provide the percentiles of a given set of values of type real

! Structure
! index1D_from3D_d3: Subroutine to compute 1D indices using 3rd dimension from 3D arrays producing 2D 
!   arrays
! stats1D_from3D_d3: Subroutine to compute 1D statistics using 3rd dimension from 3D arrays producing 
!   2D arrays

  USE module_definitions
  USE module_basic

  CONTAINS

!! Indices
  SUBROUTINE hcwi_hot(dt, tasn, tasx, tasnq90, tasxq90, missv, Nheatwaves, heatwaves)
  !  Heat Wave Index (from HCWI), which is implemented in the Copernicus European Drought 
  !    Observatory (EDO)
  ! FROM: https://edo.jrc.ec.europa.eu/documents/factsheets/factsheet_heatColdWaveIndex.pdf
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: dt
    REAl, INTENT(in)                                     :: missv
    REAL, INTENT(in)                                     :: tasnq90, tasxq90
    REAL, DIMENSION(dt), INTENT(in)                      :: tasn, tasx
    INTEGER, INTENT(out)                                 :: Nheatwaves
    INTEGER, DIMENSION(dt,2)                             :: heatwaves
    
    ! Local
    INTEGER                                              :: it
    INTEGER                                              :: Nheatwaves0
    INTEGER, DIMENSION(dt,2)                             :: heatwaves0
    LOGICAL, DIMENSION(dt)                               :: hotday
    CHARACTER(len=Sm)                                    :: fname
    
!!!!!!! Variables
! dt: amound time-steps to use
! tasn: daily minimum temperature
! tasx: daily maximum temperature
! tasnq90: 90 percentile of daily minimum temperature
! tasxq90: 90 percentile of daily maximum temperature
! missv: missing value

    fname = 'hcwi_hot'
    
    ! Hot wave as: there are at least three consecutive days with both Tmin and Tmax above (for 
    !   heatwaves) or below (for cold waves) their daily threshold values (defined as described 
    !   previously). When two successive heat or cold waves are separated in time by one day, these 
    !   are considered to be mutually dependent events, and so are merged (“pooled”) as a single event.
    hotday = .FALSE.
    DO it=1, dt
      IF (tasn(it) > tasnq90 .AND. tasx(it) > tasxq90) THEN
        hotday(it) = .TRUE.
      END IF
    END DO
    
    Nheatwaves0 = 0
    heatwaves0 = -1
    ! index beginning heat-wave
    ih = -1
    ! index ending heat-wave
    eh = -1
    ! Look for heatwaves
    DO it=1, dt
      IF (hotday(it)) THEN
        IF (ih == -1) THEN
          ih = it
        ELSE
          eh = it
        END IF
      ELSE
        ! Heat wave if consecutive days >= 3
        IF (ih /= -1 .AND. eh /= -1 .AND. eh - ih >= 2) THEN
          Nheatwaves0 = Nheatwaves0 + 1
          heatwaves0(Nheatwaves,1) = ih
          heatwaves0(Nheatwaves,2) = eh
          ih = -1
          eh = -1
        END IF
      END IF
    END DO        

    ! Pooling Heat waves
    Nheatwaves = 0
    heatwaves = -1
    DO it=1, Nheatwaves-1
      ! Should be checked if the it+2, ..., it+N also should be pooled?
      IF (heatwaves(it+1,1) - heatwaves(it,2) == 1) THEN
        Nheatwaves = Nheatwaves + 1
        heatwaves(Nheatwaves,1) = heatwaves(it,1)
        heatwaves(Nheatwaves,2) = heatwaves(it+1,2)        
      ELSE
        Nheatwaves = Nheatwaves + 1
        heatwaves(Nheatwaves,1) = heatwaves(it,1)
        heatwaves(Nheatwaves,2) = heatwaves(it,2)      
      END IF
    
    END DO

    RETURN
  
  END FUNCTION hcwi_hot
  
  SUBROUTINE percentilesR(Nvals, vals, missv, Npercents, percents, optypepercent)
  ! Subroutine to provide the percentiles of a given set of values of type real

  IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: Nvals, Npercents
    REAL, INTENT(in)                                     :: missv
    REAL, DIMENSION(Nvals), INTENT(in)                   :: vals
    CHARACTER(len=Sm), INTENT(in), OPTIONAL              :: optypepercent
    REAL, DIMENSION(Npercents), INTENT(out)              :: percents

    ! Local
    INTEGER                                              :: i, Nmiss, percenv
    INTEGER                                              :: floorv, modv
    CHARACTER(len=Sm)                                    :: typepercent
    REAL                                                 :: dpercent
    REAL, DIMENSION(Nvals)                               :: sortedvals
    CHARACTER(len=Sm)                                    :: fname

!!!!!!! Variables
! Nvals: number of values
! vals: values
! Npercents: number of percents
! optypepercent: type of percentiles 
!  FROM: https://en.wikipedia.org/wiki/Percentile
!  rank: Using the nearest-rank method (default)
!  C0: Using the interpolation

    fname = 'percentilesR'
    
    percents = fillValueR
    
    IF (PRESENT(optypepercent)) THEN
      typepercent = optypepercent
    ELSE
      typepercent = 'rank'
    END IF

    sortedvals = vals
    ! Using from: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
    CALL SortR(sortedvals, Nvals, missv)
    
    ! Missed values will be at the end of the sortedvals
    Nmiss = 0
    DO i=1, Nvals
      IF (sortedvals(i) == missv) Nmiss = Nmiss + 1
    END DO
    
    ! FROM: https://en.wikipedia.org/wiki/Percentile
    dpercent = (Nvals-Nmiss)/(Npercents-1)
    percents(1) = sortedvals(1)

    IF (TRIM(typepercent) == 'rank') THEN
      ! Using the nearest-rank method
      DO i=2, Npercents-1
        percents(i) = sortedvals(INT((i-1)*dpercent))
      END DO

    ELSE
      ! Using the interpolation C = 1/2
      DO i=2, Npercents-1
        percenv = (i-1)*5
        floorv = INT(FLOOR(percenv*(Nvals-Nmiss)/100.))
        modv = INT(MOD(percenv*(Nvals-Nmiss)*1.,100.))
        IF (modv /= 0) THEN
          percents(i) = sortedvals(floorv+1)
        ELSE
          percents(i) = (sortedvals(floorv) + sortedvals(floorv+1))/2.
        END IF
      END DO
    END IF
    
    percents(Npercents) = sortedvals(Nvals-Nmiss)

    RETURN

  END SUBROUTINE percentilesR
  
!! STRUCTURE -- structure 

  SUBROUTINE index1D_from3D_d3(d1, d2, d3, indexn, matv, maskv, missv, opmatv2Da, opmatv2Db,          &
    opmatv2Dc, opmatv3Da, opmatv3Db, opmatv3Dc, indexv)
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

    ! Local
    INTEGER                                              :: i, j, k
    INTEGER                                              :: Navailindexn
    INTEGER                                              :: Ia
    CHARACTER(len=Sm), DIMENSION(:), ALLOCATABLE         :: availindexn
    REAL, DIMENSION(:,:), ALLOCATABLE                    :: R2Da, R2Db
    REAL, DIMENSION(:,:,:), ALLOCATABLE                  :: R2Da, R2Db
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

      IF (.NOT.PRESENT(optmatv3Da)) THEN
        msg = "To compute '" // TRIM(indexn) // "' is necessary at least 'daily tas max', tasmax"
        CALL StopRun(msg, fname)
      END IF
    
      IF (ALLOCATED(R2Da)) DEALLOCATE(R2Da)
      ALLOCATE(R2Da(d1,d2))
      IF (ALLOCATED(R2Db)) DEALLOCATE(R2Db)
      ALLOCATE(R2Db(d1,d2))
      ! Computing tasmin90 [30-year baseline period (1981-2010)]
      ! 90th percentile ('Q90') of the 330 respective temperature values in an 11-day window centred 
      !   on that day, for all years in the baseline period
      IF (.NOT.PRESENT(optmatv2Da)) THEN
        IF (ALLOCATED(R3Da)) DEALLOCATE(R3Da)
        ALLOCATE(R3Da(d1,d2,Npercents))
        
        CALL stats1D_from3D_d3(d1, d2, d3, 'percentiles', matv, maskv, missv, R2Da, R3Da)
        R2Da = R3Da(:,:,19)
      ELSE
        R2Da = optmatv2Da
      END IF
    
      ! Computing tasmax90
      IF (.NOT.PRESENT(optmatv2Da)) THEN
        IF (ALLOCATED(R3Da)) DEALLOCATE(R3Da)
        ALLOCATE(R3Da(d1,d2,Npercents))
        
        CALL stats1D_from3D_d3(d1, d2, d3, 'percentiles', matv, maskv, missv, R2Da, R3Da)
        R2Db = R3Da(:,:,19)
      ELSE
        R2Db = optmatv2Db
      END IF

      indexv = 0.
      ! Amount of values
      Ia = 0
    
      DO i=1, d1
        DO j=1, d2
          IF (.NOT.maskv(i,j)) THEN
            indexv(i,j) = hcwi_hot(d3, matv, optmatv3Da, R2Da, R2Db, missv)
          END IF
        END DO
      END DO

    CASE DEFAULT
    
      msg = "1D index '" // TRIM(indexn) // "' not ready !!"

      CALL StopRunAvail(msg, fname, Navailindexn, availindexn)
    
    END SELECT

    DEALLOCATE(availindexn)

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
