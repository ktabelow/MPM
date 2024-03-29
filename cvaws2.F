C
C    Copyright (C) 2015 Weierstrass-Institut fuer
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines are part of the MPM toolbox and contain
C  FORTRAN 77 code needed in functions smoothESTATICS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cvaws2(y,nv,nvd,n1,n2,n3,nmask,iii,jjj,hakt,lambda,
     1                theta,si2,bi,bin,thnew,ncores,lwght,wght,swjy)
C
C   y        observed values of regression function
C   nv       number of vector components
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   kritical value
C   theta    estimates from last step   (input)
C   si2      inverse covariance matrices of vectors in y
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   ncores   number of cores
C
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)

      integer nv,nvd,n1,n2,n3,ncores,nmask
      integer*4 iii(nmask),jjj(*)
      logical aws
      real*8 y(nv,nmask),theta(nv,nmask),bi(nmask),thnew(nv,nmask),
     1  si2(nvd,nmask),lambda,wght(2),hakt,lwght(*),swjy(nv,ncores),
     2  bin(nmask)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr,imask,jmask
      real*8 bii,biinv,sij,swj,z1,z2,z3,wj,hakt2,hmax2,w1,w2,spmb,spf
      external lkern,KLdistsi
      real*8 lkern,KLdistsi
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=4.d0/3.d0
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
C      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,bin,nv,nvd,n1,n2,n3,hakt2,hmax2,theta,
C$OMP& ih3,lwght,wght,y,swjy,nmask,iii,jjj,si2)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(imask,jmask,i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO imask=1,nmask
         iind = iii(imask)
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
C    nothing to do, final estimate is already fixed by control
         bii=bi(imask)/lambda
         biinv=1.d0/bii
         spmb=0.25d0/bii
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jmask = jjj(jind)
                  if(jmask.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=KLdistsi(theta(1,imask),theta(1,jmask),
     1                            si2(1,imask),nv)
                     IF (sij.ge.biinv) CYCLE
                     IF (sij.gt.spmb) wj=wj*(1.d0-spf*(bii*sij-0.25d0))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jmask)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,imask)=swjy(k,thrednr)/swj
         END DO
         bin(imask)=swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bin)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=exp(-xsq*8.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN
      END

      real*8 function KLdistsi(thi,thj,si2,nv)
      implicit logical (a-z)
      integer nv
      real*8 thi(nv), thj(nv), si2(*)
      integer k,l,m
      real*8 z,zdk
      z=0.d0
      m=1
      DO k=1,nv
         zdk=thi(k)-thj(k)
         if(k.gt.1) THEN
           DO l=1,k-1
              z=z+2.d0*(thi(l)-thj(l))*zdk*si2(m)
              m=m+1
           END DO
         ENDIF
         z=z+zdk*zdk*si2(m)
         m=m+1
      END DO
      KLdistsi=z
      RETURN
      END
