c#######################################################################
c#                                MARKS                                #
c#######################################################################
c#    c - call crossAmp in transElt                                    #
c#    a - crossAmp subroutine                                          #
c#    g - grabData subroutine                                          #
c#######################################################################
      include 'inputs.f'
c
      program transElt
      use m_input
c
      real(8) :: res, deltaE, tmp
      double precision :: f3j, f6j
      real(8) :: mTot, j, jp, l
      real(8), dimension(:,:,:), allocatable :: amp
      integer :: jStep, jpstep, lStep, eStep, epStep, nEig, nEigp, nR
      integer :: param(4,7), paramp(4,7), i, mTStep,iin
      real(8), allocatable :: eig(:), eigp(:)
      integer, allocatable :: ke(:,:), kep(:,:) 
      character(14) :: name11, name12, name21, name22
      call reader
c
      open(40)
      param = 0.d0
      paramp = 0.d0
      name11 = '_kei1_kep1.dat'
      name12 = '_kei1_kep2.dat'
      name21 = '_kei2_kep1.dat'
      name22 = '_kei2_kep2.dat'
      open(16,file='transAmp.dat')
      open(17,file='transAmp2.dat')
      iin = 3
      if (jtot .ne. 0) then
         open(11,file=fileName // name11)
         open(12,file=fileName // name12)
         iin = 1
      end if
      open(13,file=fileName // name21)
      open(14,file=fileName // name22)
      do i =iin,4
         read(10 + i,*)
         read(10 + i,*) param(i,:)
      enddo
      nEig = param(1,7)+param(2,7)+param(3,7)+param(4,7)
      nR = param(3,4)
      do i =iin,4
         close(10+i)
      enddo

      open(21,file=fileNamep // name11)
      open(22,file=fileNamep // name12)
      open(23,file=fileNamep // name21)
      open(24,file=fileNamep // name22)
      do i =1,4
         read(20 + i,*)
         read(20 + i,*) paramp(i,:)
      enddo
      nEigp = paramp(1,7)+paramp(2,7)+paramp(3,7)+paramp(4,7)
      jpTot = dfloat(paramp(1,1))
      do i =1,4
         close(20+i)
      enddo

      allocate(eig(nEig))
      allocate(eigp(nEigp))
      allocate(ke(nEig,2))
      allocate(kep(nEigp,2))
      
      allocate(amp(0:25,0:25+jTot,0:25))

      write(16,*) '             E',"             E'",'        Delta E',
     &'    <psi|D|psi>'
      write(17,*) '             E',"             E'",'        Delta E',
     &'    <psi|D|psi>', '    m'
      do eStep = 1,nEig
      write(16,*) eStep
      write(17,*) eStep
      do epStep = 1, nEigp
      print*, estep,epstep
         res = 0.d0
         call crossAmp(dfloat(jTot),dfloat(jpTot),eStep,epStep,nR,nEig,
     &   nEigp,amp,deltaE,eig, eigp, ke, kep,fileName,fileNamep)
         do mTStep = -1*min(jTot,jpTot), min(jTot,jpTot)
         mTot = dfloat(mTStep)
         tmp = 0.d0
            do jStep = 0, 24
               j = dfloat(jStep)
               do lStep = abs(jTot-jStep), jTot + jStep
                  l = dfloat(lStep)
                  do jpStep = abs(jStep - 1), jStep + 1,2
                     jp = dfloat(jpStep)
                     tmp = tmp + amp(jStep, lStep,
     &     jpStep)*canalElt(dfloat(jTot), dfloat(jpTot), mTot, j, jp, l)
                  enddo
               enddo
            enddo
            tmp = tmp*dsqrt((2.d0*jpTot+1)*(2.d0*jTot+1))*(-1.d0)**mTot*
     &     f3j(dfloat(jpTot),1.d0, dfloat(jTot), -1.d0*mTot, 0.d0, mTot)
      if (dabs(tmp) .gt. 1.d-8) then
      write(17,'(4f15.8,i5)') eig(eStep), eigp(epStep),deltaE,tmp,mTStep
      end if
            res = res + tmp**2
         end do ! mTStep
         if (dabs(res) .gt. 1.d-8) then
            write(16,'(4f15.8)') eig(eStep), eigp(epStep),deltaE, res
         end if
      enddo ! epStep
      write(16,*)
      write(17,*)
      enddo ! eStep
      close(16)
      close(17)

      deallocate(amp,eig,eigp,ke,kep)
      close(40)
      end program
c#######################################################################

      real(8) function canalElt(jTot,jpTot,mTot,j,jp,l) result(res)
      real(8), intent(in) :: jTot, jpTot, mTot, j, jp, l
      real(8) :: f3j, f6j

      res = dsqrt((2*jp+1.d0)*(2*j+1.d0))*(-1.d0)**(l-jp)*f3j(jp,1.d0,j,
     & 0.d0, 0.d0, 0.d0)*f6j(jp,1.d0,j,jTot,l,jpTot)
  
      end function

      subroutine faclog
c#######################################################################
c#    initialisation of logarithms of factorials array                 #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /logfac/ fct(nfctmx)
      data ntimes /0/
c
      ntimes = ntimes+1
      if (ntimes .gt. 1) return
      fct(1) = 0.d0
      do 10 i = 1,nfctmx-1
         ai = i
         fct(i+1) = fct(i)+dlog(ai)
 10   continue
c
      return
      end


      double precision function f3j (fj1,fj2,fj3, fm1,fm2,fm3)
c#######################################################################
c#    calculates 3j coefficients from racah formula                    #
c#    (messiah: t2, p 910; formula 21) .                               #
c#    clebsch-gordan coefficients are given by (p. 908, formula 12) :  #
c#                         j -j +m                |j    j     j|       #
c#    <j j m m |j m> = (-1) 1  2   (2*j+1)**(0.5) | 1    2     |       #
c#      1 2 1 2                                   |m    m    -m|       #
c#                                                | 1    2     |       #
c#---------------------------------------------------------------------#
c#    has been tested for j up to 200.                                 #
c#    logfac should contain the logarithms of factorials               #
c#    fj-fm integer not checked                                        #
c#    j.m.l. (1975)                                                    #
c#######################################################################
      implicit double precision (a-h,o-z)
      integer t,tmin,tmax
      parameter (nfctmx=1001)
      data tiny,zero,one /0.01d0,0.d0,1.d0/ ,ntimes /1/
      common /logfac/ fct(nfctmx)
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      cc = zero
      if (fj3 . gt. (fj1+fj2+tiny))      go to 100
      if (dabs(fj1-fj2) .gt. (fj3+tiny)) go to 100
      if (dabs(fm1+fm2+fm3) .gt. tiny)   go to 100
      if (dabs(fm1) .gt. (fj1+tiny))     go to 100
      if (dabs(fm2) .gt. (fj2+tiny))     go to 100
      if (dabs(fm3) .gt. (fj3+tiny))     go to 100
      fk1 = fj3-fj2+fm1
      fk2 = fj3-fj1-fm2
      fk3 = fj1-fm1
      fk4 = fj2+fm2
      fk5 = fj1+fj2-fj3
      fk1m = fk1-tiny
      fk2m = fk2-tiny
      fk1p = fk1+tiny
      fk2p = fk2+tiny
      if (fk1m .lt. zero) k1 = fk1m
      if (fk1p .gt. zero) k1 = fk1p
      if (fk2m .lt. zero) k2 = fk2m
      if (fk2p .gt. zero) k2 = fk2p
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      tmin = 0
      if (k1+tmin .lt. 0) tmin = -k1
      if (k2+tmin .lt. 0) tmin = -k2
      tmax = k3
      if (k4-tmax .lt. 0) tmax = k4
      if (k5-tmax .lt. 0) tmax = k5
      n1 = fj1+fj2-fj3+one+tiny
      n2 = fj2+fj3-fj1+one+tiny
      n3 = fj3+fj1-fj2+one+tiny
      n4 = fj1+fm1+one+tiny
      n5 = fj2+fm2+one+tiny
      n6 = fj3+fm3+one+tiny
      n7 = fj1-fm1+one+tiny
      n8 = fj2-fm2+one+tiny
      n9 = fj3-fm3+one+tiny
      n10 = fj1+fj2+fj3+2.d0+tiny
      x = fct(n1)+fct(n2)+fct(n3)+fct(n4)+fct(n5)+fct(n6)
     &   +fct(n7)+fct(n8)+fct(n9)-fct(n10)
      x = 0.5d0*x
      do 10  t = tmin,tmax
         phase = one
         if (mod(t,2) .ne. 0) phase = -one
         cc = cc+phase*dexp(-fct(t+1)   -fct(k1+t+1)-fct(k2+t+1)
     &                      -fct(k3-t+1)-fct(k4-t+1)-fct(k5-t+1)+x)
 10   continue
      fsp = dabs(fj1-fj2-fm3)+tiny
      ns = fsp
      if (mod(ns,2) .gt. 0) cc = -cc
 100  f3j = cc
      return
      end


      double precision function f6j (fj1,fj2,fj3,fl1,fl2,fl3)
c#######################################################################
c#    calculation of 6j-coefficients                                   #
c#---------------------------------------------------------------------#
c#    checked by j.m.l. (1980), modified 11/1981, not checked          #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /logfac/ fct(nfctmx)
      data tiny /.01/ ,ntimes /1/
c
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      d = fdelta (fj1,fj2,fj3)
      d = d*fdelta (fj1,fl2,fl3)
      d = d*fdelta (fl1,fj2,fl3)
      d = d*fdelta (fl1,fl2,fj3)
      f6j = 0.d0
      if (dabs(d) .eq. 0.d0) return
c
      fk1 = fj1+fj2+fj3
      fk2 = fj1+fl2+fl3
      fk3 = fl1+fj2+fl3
      fk4 = fl1+fl2+fj3
      fk5 = fj1+fj2+fl1+fl2
      fk6 = fj2+fj3+fl2+fl3
      fk7 = fj3+fj1+fl3+fl1
      fmin = dmin1 (fk5,fk6,fk7)
      fmax = dmax1 (fk1,fk2,fk3,fk4)
      min = fmin+tiny
      max = fmax+tiny
      k1 = fk1+tiny
      k2 = fk2+tiny
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      k6 = fk6+tiny
      k7 = fk7+tiny
      if (min-max) 1000,3,3
 3    if (max) 1000,4,4
 4    if (min) 1000,5,90
 5    k1 = -k1
      k2 = -k2
      k3 = -k3
      k4 = -k4
      bot = fct(k1+1)+fct(k2+1)+fct(k3+1)+fct(k4+1)+fct(k5+1)+fct(k6+1)
     &     +fct(k7+1)
      bot = dexp(bot)
      f6j = d/bot
      return
c
 90   f6j = 0.
      do 100 i = max,min
         boite = ((-1.)**i)
         iz = i+1
         m1 = i-k1
         m2 = i-k2
         m3 = i-k3
         m4 = i-k4
         m5 = k5-i
         m6 = k6-i
         m7 = k7-i
         dot = fct(iz+1)
         bot = fct(m1+1)+fct(m2+1)+fct(m3+1)+fct(m4+1)+fct(m5+1)
     &        +fct(m6+1)+fct(m7+1)
         b1 = dot-bot
         boite = boite*dexp(b1)
         f6j = f6j+boite
 100  continue
      f6j = f6j*d
 1000 return
c
      end
      double precision function fdelta (fl1,fl2,fl3)
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /logfac/ fct(nfctmx)
      data   eps /.01/
c
      ia=fl1+fl2+fl3 +eps
      a=2.*(fl1+fl2+fl3)+1.
      ib=a +eps
      ib=ib/2
      if(ib-ia)1,6,1
    6 continue
      ik1=fl1+fl2-fl3+eps
      ik2=fl2+fl3-fl1+eps
      ik3=fl3+fl1-fl2+eps
      kk=fl1+fl2+fl3+1+eps
      if(ik1)1,2,2
    2 if(ik2)1,3,3
    3 if(ik3)1,4,4
    4 d1=fct(kk+1)
      d2=fct(ik1+1)+fct(ik2+1)+fct(ik3+1)
      d3 = (d2 - d1) / 2.d0
      fdelta = dexp (d3)
      go to 5
    1 fdelta=0.
    5 return
c
      end

      subroutine crossAmp(jTot, jpTot, e, ep, nR, nEig, nEigp, amp,
     & deltaE, eig, eigp, ke, kep,fileName,fileNamep)
c#######################################################################
c#  Computes the crossed terms F_{J,M,j,l}(R)*F_{J',M',j',l'}(R) and   #
c#    sums them over R using the weights w(R) and the 1/R² factor.     #
c#    Uses psi_{E,J,M} = 1/R*SUM_{j,l}Y_{J,M,j,l}*F_{E,J,M,j,l}(R)     #
c#######################################################################
      real(8), dimension(:), allocatable :: rr, ww
      real(8), intent(inout) :: jTot, jpTot 
      real(8), intent(out) :: amp(0:25,0:25+int(jTot),0:25), deltaE
      real(8), allocatable :: uuNew(:,:,:,:), uuNewp(:,:,:,:)
      integer, intent(inout) :: e, ep, nEig, nEigp
      integer :: rStep, jStep, lStep, jpStep, nR
      real(8), intent(inout) :: eig(nEig), eigp(nEigp)
      integer ::  ke(nEig,2), kep(nEigp,2) 
      character(20), intent(in) :: fileName,fileNamep
c
      allocate(rr(nR))
      allocate(ww(nR))
      allocate(uuNew(nEig,nR,0:25,0:25+int(jTot)))
      allocate(uuNewp(nEigp,nR,0:25,0:25+int(jpTot)))
c
      call grabData(fileName,rr, ww, eig, ke, uuNew, nR, nEig,jtot)
      call grabData(fileNamep,rr, ww, eigp, kep, uuNewp,nR,nEigp,jptot)
c
      amp = 0.d0
c
      write(30,*), eig(e), eigp(ep)
      do jStep = 0, 24
         do lStep = abs(jStep-int(jTot)), jStep+int(jTot)
            do jpStep = abs(jStep-1), jStep+1,2
               do rStep = 1,nR
                  amp(jStep, lStep,jpStep) = amp(jStep, lStep,jpStep) +
     & uuNew(e,rStep,jStep,lStep)*uuNewp(ep,rStep,jpStep,lStep)
c     & /rr(rStep)**2
               end do
            write(30,*), jStep, jpStep, lStep, amp(jStep,lStep,jpStep)
            end do
         end do
      end do
      write(30,*)
      deltaE = eigp(ep)-eig(e)
c
      deallocate(rr, ww,uunew,uunewp)
      end subroutine

      subroutine grabData(fileName,rr,ww,eig,ke,uuNew,nR,nEig,jTot)
c#######################################################################
c#      Reads data from the 4 files with kei = 1,2 and kep = 1,2.      #
c#######################################################################
      integer, intent(in) :: nEig, nR
      real(8) :: jTot
      real(8), intent(inout) :: uuNew(nEig,nR,0:25,0:25+int(jTot))
      integer :: jStep, lStep, rStep, eStep, i, j, chStep, junk2, c
      integer, dimension(4,7) :: param
      integer, intent(inout) :: ke(nEig,2)
      character(len=4) :: junk1
      real(8), intent(inout) :: rr(nR), ww(nR), eig(nEig)
      integer :: index(nEig), ketmp(nEig,2), iin
      real(8) :: eigtmp(nEig), uunewtmp(nEig,nR,0:25,0:25+int(jTot))
      character(20) :: fileName
      character(14) :: name11, name12, name21, name22
      name11 = '_kei1_kep1.dat'
      name12 = '_kei1_kep2.dat'
      name21 = '_kei2_kep1.dat'
      name22 = '_kei2_kep2.dat'
c
      iin=3
      if (int(jTot) .ne. 0) then
         iin = 1
         open(11,file=trim(fileName) // trim(name11))
         open(12,file=fileName // name12)
      end if
         open(13,file=fileName // name21)
         open(14,file=fileName // name22)
      do i =iin,4
         read(10 + i,*)
         read(10 + i,*) param(i,:)
      enddo
c
      uuNew = 0.d0
      c = 1
      do i = iin, 4
         read(10+i,*)
         read(10+i,*)
         do rStep = 1,nR
            read(10+i,*) rr(rStep), ww(rStep)
         enddo
         read(10+i,*)
         read(10+i,*)
         do chStep = 1,param(i,5)
            read(10+i,*)
         enddo 
         do eStep = 1, param(i,7)
            read(10+i,*)
            read(10+i,*)junk1, junk2, eig(c)
            ke(c,1) = param(i,2)
            ke(c,2) = param(i,3)
            do jStep = 2-param(i,3), 24, 2
               do lStep = abs(jStep-int(jTot)),jStep+int(jTot)
                  if (mod(param(i,2)+param(i,3),2) .eq. mod(lStep,2))
     & then
                     do rStep = 1,nR
                        read(10+i,*) uuNew(c,rStep,jStep,lStep)
                     enddo
                  endif
               enddo
            enddo
            c = c+1
         enddo
      enddo
c
      call sortrx(nEig,eig,index)
      do eStep = 1,nEig
         eigtmp(estep)         = eig(index(estep))
         ketmp(estep,:)        = ke(index(estep),:)
         uunewtmp(estep,:,:,:) = uunew(index(estep),:,:,:)
      end do
      do eStep = 1,nEig
         eig(estep)         = eigtmp(estep)        
         ke(estep,:)        = ketmp(estep,:)
         uunew(estep,:,:,:) = uunewtmp(estep,:,:,:) 
      end do
      do i = iin,4
         close(10+i)
      end do
      call maxValue(nEig,nR,int(jTot),uuNew)
      
      end subroutine
c
c A Viel 31 March 2008 - double precision
C From Leonard J. Moss of SLAC:

C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.  

      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL*8
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      REAL*8    DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL*8    DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
C===================================================================
C
C     All done
 
      END
      subroutine reader
      use   m_input
      read  (5,input)
c      write (6,input)
      return
      end

      subroutine maxValue(nEig,nR,jTot,uuNew)
      integer, intent(in) :: nEig, nR,jTot
      real(8), intent(inout) :: uuNew(nEig,nR,0:25,0:25+jTot)
      integer :: e, r, j, l, jm, lm, rm
      real(8) :: mVal
      do e = 1,nEig
         mVal=0.d0
         do j = 0,25
            do l = abs(j-jTot), j+jTot
               do r=1,nR
                  if (dabs(uuNew(e,r,j,l)) .gt. mVal) then
                     jm = j
                     lm = l
                     rm = r
                     mVal = dabs(uuNew(e,r,j,l))
                  end if
               end do
            end do
         end do
         uuNew(e,:,:,:) = 0.d0
         uuNew(e,1,jm,lm) = 1.d0
         write(40,*)e,jm,lm
      end do
      write(40,*)
      end subroutine
