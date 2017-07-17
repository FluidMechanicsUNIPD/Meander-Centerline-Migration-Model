!-----------------------------------------------------------------------
! PRINT CPU TIME AND DATE
!-----------------------------------------------------------------------
      subroutine printtime(tempo, jfile, jind)
      implicit none

      integer tempo(8), jfile, jind, DD, MO, YY, HH, MM, SS

!-----------------------------------------------------------------------
! print on file
      if (jind.eq.0) then
        write(jfile,85,advance='no') 'Simulation start : '
        write(6,85,advance='no') 'Simulation start : '
       else
        write(jfile,85,advance='no') 'Simulation end   : '
        write(6,85,advance='no') 'Simulation end   : '
      end if

      DD = tempo(3)
      MO = tempo(2)
      YY = tempo(1)
      HH = tempo(5)
      MM = tempo(6)
      SS = tempo(7)

! stampa l'istante nel formato DD/MM/YYYY HH/MM/SS
!      write(jfile,20) DD, '/', MO, '/', SS, HH, ':', MM, ':', SS

! giorno
      if (DD.le.9) then
        write(jfile,30,advance='no') '0',DD,'/'
        write(6,30,advance='no') '0',DD,'/'
       else
        write(jfile,40,advance='no') DD,'/'
        write(6,40,advance='no') DD,'/'
      end if

! mese
      if (MO.le.9) then
        write(jfile,30,advance='no') '0',MO,'/'
        write(6,30,advance='no') '0',MO,'/'
       else
        write(jfile,40,advance='no') MO,'/'
        write(6,40,advance='no') MO,'/'
      end if

! anno
      write(jfile,31,advance='no') YY
      write(6,31,advance='no') YY

! ore
      if (HH.le.9) then
        write(jfile,30,advance='no') '0',HH,':'
        write(6,30,advance='no') '0',HH,':'
       else
        write(jfile,40,advance='no') HH,':'
        write(6,40,advance='no') HH,':'
      end if

! minuti
      if (MM.le.9) then
        write(jfile,30,advance='no') '0',MM,':'
        write(6,30,advance='no') '0',MM,':'
       else
        write(jfile,40,advance='no') MM,':'
        write(6,40,advance='no') MM,':'
      end if

! secondi
      if (jind.eq.0) then
        if (SS.le.9) then
          write(jfile,33) '0',SS
          write(6,33) '0',SS
         else
          write(jfile,43) SS
          write(6,43) SS
        end if
       else
        if (SS.le.9) then
          write(jfile,32) '0',SS
          write(6,32) '0',SS
         else
          write(jfile,42) SS
          write(6,42) SS
        end if
      end if

!-----------------------------------------------------------------------  
85    format(a22,2x)
20    format(2(i2,a1),i4,5x,2(i2,a1),i2)
30    format(a1,i1,a1)
31    format(i4,3x)
32    format(a1,i1,/)
33    format(a1,i1)
40    format(i2,a1)
42    format(i2,/)
43    format(i2)
!-----------------------------------------------------------------------
      return
      end
