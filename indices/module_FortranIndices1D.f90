MODULE module_FortranIndices1D
! Module with the instructions for the computation of 1D indices using Fortran
!   Using developments from PyNCplot (https://git.cima.fcen.uba.ar/lluis.fita/pyncplot/-/wikis/home)
!

!!!!!!!! Contents
! Indices
! hcwi_hot: Heat Wave Index (from HCWI), which is implemented in the Copernicus European Drought 
!   Observatory (EDO)
! percentilesR: Subroutine to provide the percentiles of a given set of values of type real

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
    INTEGER                                              :: it, ih, eh
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
  
  END SUBROUTINE hcwi_hot
  
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

END MODULE module_FortranIndices1D

