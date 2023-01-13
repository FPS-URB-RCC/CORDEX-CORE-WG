MODULE module_definitions

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

END MODULE module_definitions 
