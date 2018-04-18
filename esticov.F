      subroutine invcov3(theta,te,ind,sigma2i,nd,n,icov)
C   theta - array of estimates
C   te    - TE times
C   ind   - indicator for T1 (1) and PD (2) images
C   sigma2i - inverse of error variance
C   icov  - resulting inverse covariance matrix
      implicit none
      integer nd,n,ind(nd)
      double precision theta(3,n),sigma2i(n),icov(6,n),te(nd)
      integer i,j
      real*8 z11,z22,z13,z23,z33,aj,ajs,th1,th2,th1s,th2s,th3,tej,si
      do i=1,n
         z11=0.d0
         z22=0.d0
         z13=0.d0
         z23=0.d0
         z33=0.d0
         th3=theta(3,i)
         th1=theta(1,i)
         th2=theta(2,i)
         th1s=th1*th1
         th2s=th2*th2
         si=sigma2i(i)
         do j=1,nd
            tej=te(j)
            aj=exp(-th3*tej)
            ajs=aj*aj
            if(ind(j).eq.1) THEN
               z11=z11+ajs
               z13=z13+ajs*tej
               z33=z33+th1s*tej*tej*ajs
            END IF
            if(ind(j).eq.2) THEN
               z22=z22+ajs
               z23=z23+ajs*tej
               z33=z33+th2s*tej*tej*ajs
            END IF
         end do
         icov(1,i)=z11*si
         icov(2,i)=0.d0
         icov(3,i)=z22*si
         icov(4,i)=-th1*z13*si
         icov(5,i)=-th2*z23*si
         icov(6,i)=z33*si
      end do
      return
      end

      subroutine invcov4(theta,te,ind,sigma2i,nd,n,icov)
C   theta - array of estimates
C   te    - TE times
C   ind   - indicator for T1 (1) and PD (2) images
C   sigma2i - inverse of error variance
C   icov  - resulting inverse covariance matrix
      implicit none
      integer nd,n,ind(nd)
      double precision theta(4,n),sigma2i(n),icov(10,n),te(nd)
      integer i,j
      real*8 z11,z22,z33,z14,z24,z34,z44,aj,ajs,th1,th2,th3,
     1       th1s,th2s,th3s,th4,tej,si
      do i=1,n
         z11=0.d0
         z22=0.d0
         z33=0.d0
         z14=0.d0
         z24=0.d0
         z34=0.d0
         z44=0.d0
         th4=theta(4,i)
         th1=theta(1,i)
         th2=theta(2,i)
         th3=theta(3,i)
         th1s=th1*th1
         th2s=th2*th2
         th3s=th3*th3
         si=sigma2i(i)
         do j=1,nd
            tej=te(j)
            aj=exp(-th4*tej)
            ajs=aj*aj
            if(ind(j).eq.1) THEN
               z11=z11+ajs
               z14=z14+ajs*tej
               z44=z44+th1s*tej*tej*ajs
            END IF
            if(ind(j).eq.2) THEN
               z22=z22+ajs
               z24=z24+ajs*tej
               z44=z44+th2s*tej*tej*ajs
            END IF
            if(ind(j).eq.3) THEN
               z33=z33+ajs
               z34=z34+ajs*tej
               z44=z44+th3s*tej*tej*ajs
            END IF
         end do
         icov(1,i)=z11*si
         icov(2,i)=0.d0
         icov(3,i)=z22*si
         icov(4,i)=0.d0
         icov(5,i)=0.d0
         icov(6,i)=z33*si
         icov(7,i)=-th1*z14*si
         icov(8,i)=-th2*z24*si
         icov(9,i)=-th3*z34*si
         icov(10,i)=z44*si
      end do
      return
      end
