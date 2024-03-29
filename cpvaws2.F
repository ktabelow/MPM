C    differences to cvaws2.F
C
C   -     ndv is nv*(nv+1)/2
C   -     invcov replaces si2 and has dimension (ndv,nmask) instead of (nv,nv,nmask)
C              storage as  si_11, si_21, \dots, si_p1, si_22, \dots, si_2p, \dots si_pp
C   -     psize   patch extension  patch is a cube of (2*psize+1)^3
C
      subroutine cpvaws2(y,nv,nvd,n1,n2,n3,nmask,iii,jjj,hakt,lambda,
     1                theta,invcov,bi,bin,thnew,ncores,lwght,wght,
     2                swjy,psize)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      integer*4 nv,n1,n2,n3,ncores,nvd,nmask,psize
      logical aws
      integer*4 iii(nmask),jjj(*)
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,*),invcov(nvd,*),
     2  bin(*)
      integer*4 ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,
     2        dlw12,n12,k,thrednr
      double precision sij,swj,z1,z2,z3,wj,hakt2,hmax2,w1,w2,sijp
      integer*4 imask,jmask,ipmask,jpmask
      integer*4 ip1,ip2,ip3,ipind,jp1,jp2,jp3,jpind
      external lkern, KLdistsi
      double precision lkern, KLdistsi
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spmin=0.25d0
      spf=1.d0/(1.d0-spmin)
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
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,bin,nv,nvd,n1,n2,n3,hakt2,hmax2,theta,invcov,
C$OMP& ih3,lwght,wght,y,swjy,nmask,psize,iii,jjj)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,ip1,ip2,ip3,ipind,
C$OMP& jp1,jp2,jp3,jpind,sijp,l,m,imask,jmask,ipmask,jpmask)
C$OMP DO SCHEDULE(GUIDED)
      DO imask=1,nmask
!$         thrednr = omp_get_thread_num()+1
         iind = iii(imask)
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
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
                  if(jjj(jind).eq.0) CYCLE
                  jmask=jjj(jind)
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     sij=0.d0
                     DO ip1=i1-psize,i1+psize
                        if(ip1.le.0.or.ip1.gt.n1) CYCLE
                        jp1=ip1+jw1
                        DO ip2=i2-psize,i2+psize
                           if(ip2.le.0.or.ip2.gt.n2) CYCLE
                           jp2=ip2+jw2
                           DO ip3=i3-psize,i3+psize
                              if(sij.gt.1.d0) CYCLE
                              if(ip3.le.0.or.ip3.gt.n3) CYCLE
                              ipind=ip1+(ip2-1)*n1+(ip3-1)*n1*n2
                              if(jjj(ipind).eq.0) CYCLE
                              ipmask=jjj(ipind)
                              jp3=ip3+jw3
                              if(jp1.le.0.or.jp1.gt.n1) CYCLE
                              if(jp2.le.0.or.jp2.gt.n2) CYCLE
                              if(jp3.le.0.or.jp3.gt.n3) CYCLE
                              jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                              if(jjj(jpind).eq.0) CYCLE
                              jpmask=jjj(jpind)
C   need both ipind and jpind in mask
                              sijp=KLdistsi(theta(1,jpmask),
     1                                theta(1,ipmask),
     1                                invcov(1,ipmask),nv)
                              sij=max(sij,bi(ipmask)/lambda*sijp)
                           END DO
                        END DO
                     END DO
                     IF (sij.ge.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
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
      implicit none
      integer*4 kern
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
      integer*4 nv
      real*8 thi(nv), thj(nv), si2(*)
      integer*4 k,l,m
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
