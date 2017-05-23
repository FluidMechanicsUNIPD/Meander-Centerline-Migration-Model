c--------------------------------------------------------------
c       subroutines calcolo integrali
c--------------------------------------------------------------
        subroutine CavSimpds(dzz,ff,sum,NN)
c--------------------------------------------------------------
c       Integrazione con Cavalieri Simpson di una funzione reale
c       passo dzz costante
c--------------------------------------------------------------
        implicit none
        integer NN,j
        real*8 dzz,ff,sum 
        dimension ff(NN)
c--------------------------------------------------------------
        sum=0.D0
        do j=2,NN-1
                sum=sum+ff(j)
        end do
        sum=dzz*(sum+(ff(1)+ff(NN))/2.0d0)

        return
        end

c--------------------------------------------------------------
        subroutine CavSimds2(dzz,ff1,ff2,sum1,sum2,NN)
c--------------------------------------------------------------
c       Integrazione con Cavalieri Simpson delle parti reale ed immaginaria di una funzione
c       passo dzz costante
c--------------------------------------------------------------
      implicit none
        integer NN,j
        real*8 dzz,ff1,ff2,sum1,sum2  
        dimension ff1(NN),ff2(NN)
c--------------------------------------------------------------
        sum1=0.D0
        sum2=0.D0
        do j=2,NN-1
                sum1=sum1+ff1(j)
                sum2=sum2+ff2(j)
        end do
        sum1=dzz*(sum1+(ff1(1)+ff1(NN))/2.0d0)
        sum2=dzz*(sum2+(ff2(1)+ff2(NN))/2.0d0)

        return
        end

c--------------------------------------------------------------
        subroutine CavSimp(zz,ff,sum,NN)
c--------------------------------------------------------------
c       Integrazione con Cavalieri Simpson di una funzione reale
c       passo dz variabile
c--------------------------------------------------------------
      implicit none
        integer NN,j
        real*8 zz,ff,sum,dzz  
        dimension zz(NN),ff(NN)
c--------------------------------------------------------------
        sum=0.D0
        do j=1,NN-1
                dzz=zz(j+1)-zz(j)
                sum=sum + dzz*(ff(j) +  ff(j+1))/2.D0
        end do
        return
        end

c--------------------------------------------------------------
        subroutine CavSimp2(zz,ff1,ff2,sum1,sum2,NN)
c--------------------------------------------------------------
c       Integrazione con Cavalieri Simpson di due funzioni reali
c       passo dz variabile
c--------------------------------------------------------------
      implicit none
        integer NN,j
        real*8 zz,ff1,ff2,sum1,sum2,dzz  
        dimension zz(NN),ff1(NN),ff2(NN)
c--------------------------------------------------------------
        sum1=0.D0
        sum2=0.D0
        do j=1,NN-1
                dzz=zz(j+1)-zz(j)
                sum1=sum1 + dzz*(ff1(j) + ff1(j+1))/2.D0
                sum2=sum2 + dzz*(ff2(j) + ff2(j+1))/2.D0
        end do
        return
        end

c--------------------------------------------------------------
        subroutine CavSimp3(zz,ff1,ff2,ff3,sum1,sum2,sum3,NN)
c--------------------------------------------------------------
c       Integrazione con Cavalieri Simpson di 3 funzioni reali
c       passo dz variabile
c--------------------------------------------------------------
      implicit none
        integer NN,j
        real*8 zz,ff1,ff2,ff3,sum1,sum2,sum3,dzz  
        dimension zz(NN),ff1(NN),ff2(NN),ff3(NN)
c--------------------------------------------------------------
        sum1=0.D0
        sum2=0.D0
        sum3=0.D0
        do j=1,NN-1
                dzz=zz(j+1)-zz(j)
                sum1=sum1 + dzz*(ff1(j) + ff1(j+1))/2.D0
                sum2=sum2 + dzz*(ff2(j) + ff2(j+1))/2.D0
                sum3=sum3 + dzz*(ff3(j) + ff3(j+1))/2.D0
        end do
        return
        end
