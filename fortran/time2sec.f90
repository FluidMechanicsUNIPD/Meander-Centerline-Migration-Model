      SUBROUTINE TIME2SEC(VALUES, TIME_SEC)
!
!------------------------------------------------------------------------
!  PURPOUSE: convert time from the form `hhmmss.ss' (hours, minutes,
!            seconds and milliseconds) into seconds
!------------------------------------------------------------------------
!
      IMPLICIT NONE

      INTEGER YY, MO, DD, HH, MM, SS, MS, VALUES(8)
      CHARACTER(8) DATE
      CHARACTER(10) TIME
      CHARACTER(5) ZONE
      REAL*8 TIME_SEC
!
!------------------------------------------------------------------------
      CALL DATE_AND_TIME(DATE, TIME, ZONE, VALUES)
!
      YY = REAL(VALUES(1))
      MO = REAL(VALUES(2))
      DD = REAL(VALUES(3))
      HH = REAL(VALUES(5))
      MM = REAL(VALUES(6))
      SS = REAL(VALUES(7))
      MS = REAL(VALUES(8))

      TIME_SEC = (((((YY * 12.d0 + MO) * 30.d0 + DD) * 24.d0 + HH) * &
     &                   60.d0 + MM) * 60.d0 + SS) * 1000.d0 + MS
!
!------------------------------------------------------------------------
      RETURN
      END
