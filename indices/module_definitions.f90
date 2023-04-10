MODULE module_definitions
! Module with the definitions for the computation of 1D indices using Fortran
!   Using developments from PyNCplot (https://git.cima.fcen.uba.ar/lluis.fita/pyncplot/-/wikis/home)
!

  IMPLICIT NONE
  
  ! small huge
  INTEGER, PARAMETER                                     :: Sh = 500
  ! small large
  INTEGER, PARAMETER                                     :: Sl = 100
  ! small character
  INTEGER, PARAMETER                                     :: Sm = 50
  ! small tiny
  INTEGER, PARAMETER                                     :: St = 30
  ! small micro
  INTEGER, PARAMETER                                     :: Su = 10

  ! Error title
  CHARACTER(len=32), PARAMETER                           :: errormsg = 'ERROR -- error -- ERROR -- error'
  ! message
  CHARACTER(len=Sh)                                      :: msg
  
  ! Real FillValue
  REAL, PARAMETER                                        :: fillValueR = 1.e20
  
  ! Amount of percentiles
  ! For minval, 5 to 95%, maxval
  INTEGER, PARAMETER                                     :: Npercents = 21

END MODULE module_definitions 
