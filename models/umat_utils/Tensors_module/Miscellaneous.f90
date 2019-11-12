! jacobi_mekh       ADD DESCRIPTION!
! eigsrt            ADD DESCRIPTION!
! M_pow             Calculate C^m where C is 3x3 matrix
! eig_prob          ADD DESCRIPTION!
! us_solve          ADD DESCRIPTION!
! us_gausjor        ADD DESCRIPTION!

Module Miscellaneous_tensors
use Conversions
Implicit none
    CONTAINS
    
!==============================================================================  
recursive SUBROUTINE jacobi_mekh(a,n,np,d,v,nrot)
Implicit none

    INTEGER                     :: n,np,nrot,NMAX
    DOUBLE PRECISION            :: a(np,np),d(np),v(np,np)
    PARAMETER (NMAX=500)
    INTEGER i,ip,iq,j
    DOUBLE PRECISION            :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
    
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+&
      g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/dsqrt(1.d0+t**2.d0)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      write(*,*) 'too many iterations in jacobi_mekh'
      return
end subroutine
!  (C) Copr. 1986-92 Numerical Recipes Software #v([&.
!============================================================================== 

!============================================================================== 
recursive SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np
      DOUBLE PRECISION d(np),v(np,np)
      INTEGER i,j,k
      DOUBLE PRECISION p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
 END subroutine
!  (C) Copr. 1986-92 Numerical Recipes Software #v([&.
!============================================================================== 

!============================================================================== 
FUNCTION M_pow(a,m) RESULT (am)
! Old name: US_M_POW_GEN
Implicit none

    DOUBLE PRECISION, DIMENSION(:,:), INTENT (IN) ::a
    DOUBLE PRECISION, DIMENSION(size(a,1),size(a,2))::am
    DOUBLE PRECISION  ::av(9),amv(9),la(3),G(3,3),m1(9),m2(9),m3(9)
    DOUBLE PRECISION  :: m

    av(1)=a(1,1)
    av(2)=a(2,2)
    av(3)=a(3,3)
    av(4)=a(1,2)
    av(9)=a(3,2)
    av(5)=a(2,3)
    av(8)=a(2,1)
    av(7)=a(1,3)
    av(6)=a(3,1)

    CALL eig_prob(av,la,G,m1,m2,m3)

    amv=dabs(la(1))**m*m1+dabs(la(2))**m*m2+dabs(la(3))**m*m3
    
    am(1,1)=amv(1)
    am(2,2)=amv(2)
    am(3,3)=amv(3)
    am(1,2)=amv(4)
    am(2,1)=amv(8)
    am(2,3)=amv(5)
    am(3,2)=amv(9)
    am(1,3)=amv(7)
    am(3,1)=amv(6)

END FUNCTION 
!============================================================================== 

!============================================================================== 
recursive SUBROUTINE eig_prob(EPSE,lambdad,Gvec,M1V,M2V,M3V) 
! Purpose: Subroutine for solution of  eigenvalue problem, and 
!          some additional, where perturbations are used if neccessary.
!          This routine is a copy of EIGPROBNEW in the file USMATEIG.f
! 
! Written by:  Magnus Ekh  
! Modified by: Jim Brouzoulis
Implicit none
      
    DOUBLE PRECISION    :: EPSE,lambda,Gvec,M1V,M2V,M3V,I1,dfact
    DOUBLE PRECISION    :: EPSEM,slask,M1,M2,M3,lmax,TOL,DIST,SLASKM,lambdad,DELTA
    INTEGER             :: nrot,NM,NV,I,J
    PARAMETER (NM=3,NV=9)
    DIMENSION &
            EPSE(NV),lambda(NM),Gvec(NM,NM),M1V(NV),M2V(NV),M3V(NV),&
            dfact(NM),EPSEM(NM,NM),slask(NM),M1(NM,NM),M2(NM,NM),M3(NM,NM),&
            SLASKM(NM,NM),lambdad(NM),DELTA(NM,NM)
      
    EPSEM = V9_2_M(EPSE)

    CALL jacobi_mekh(EPSEM,NM,NM,lambda,Gvec,nrot)
    ! write(*,*) 'lambda ',lambda
    CALL EIGSRT(lambda,Gvec,NM,NM)
    lambdad=lambda
    ! write(*,*) 'Gvec1 ',Gvec(1:3,1)
    ! write(*,*) 'Gvec2 ',Gvec(1:3,2)
    ! write(*,*) 'Gvec3 ',Gvec(1:3,3)

    !C Find largest |lambda|
    !      lmax=-1.d9
    !      DO I=1,3
    !         IF (ABS(lambda(I)).GT.lmax) THEN
    !            lmax=ABS(lambda(I))
    !
    !         ENDIF
    !      ENDDO
    !
    !      lmax=max(lmax,1.d-12)
    !      TOL=1.d-7*lmax
    !      DIST=TOL/4.d0
    !
    !C Disturb if neccessary
    !      IF (ABS(lambda(1)-lambda(2))/lmax.LT.TOL) THEN
    !         lambdad(1)=lambda(1)*(1.d0+DIST)
    !         lambdad(2)=lambda(2)*(1.d0-DIST)
    !         lambdad(3)=lambda(3)/((1.d0+DIST)*(1.d0-DIST))
    !      ELSEIF (ABS(lambda(2)-lambda(3))/lmax.LT.TOL) THEN
    !         lambdad(2)=lambda(2)*(1.d0+DIST)
    !         lambdad(3)=lambda(3)*(1.d0-DIST)
    !         lambdad(1)=lambda(1)/((1.d0+DIST)*(1.d0-DIST))
    !      ELSE
    !         lambdad(1)=lambda(1)
    !         lambdad(2)=lambda(2)
    !         lambdad(3)=lambda(3)
    !      ENDIF

    !C If the eigenvalues are equal to zero
    !      IF (ABS(lambdad(1)).LT.1.d-12) THEN
    !         lambdad(1)=2.d-12
    !      ELSEIF (ABS(lambdad(2)).LT.1.d-12) THEN
    !         lambdad(2)=1.d-12
    !      ELSEIF (ABS(lambdad(3)).LT.1.d-12) THEN
    !         lambdad(3)=0.5d-12
    !      ENDIF
    !
    !      write(*,*) '???lambdad',lambdad(1:3)
    !C
    !
    !      I1=lambdad(1)+lambdad(2)+lambdad(3)
    !      detJ =lambdad(1)*lambdad(2)*lambdad(3)
    !
    !      DO I=1,NM
    !         dfact(I)=2.d0*lambdad(I)**2.d0-
    !     .    I1*lambdad(I)+detJ/lambdad(I)
    !      ENDDO
    !      write(*,*) '???dfact',dfact(1:3)
    !
    !      CALL USMMM(EPSEM,3,3,EPSEM,3,3,SLASKM)
    !      DO I=1,3
    !         DO J=1,3
    !
    !            DELTA(I,J)=0.d0
    !            IF (I.EQ.J) THEN
    !               DELTA(I,J)=1.d0
    !            ENDIF
    !         ENDDO
    !      ENDDO
    !C23456789012345678901234567890123456789012345678901234567890123456789012      
    DO I=1,3
        DO J=1,3
    !      M1(I,J)=1.d0/(dfact(1))*( 
    !     .   SLASKM(I,J)-(I1-lambdad(1))*EPSEM(I,J)
    !
    !     .        +detJ/(lambdad(1))*DELTA(I,J) )
           M1(I,J)=Gvec(I,1)*Gvec(J,1)

    !      M2(I,J)=1.d0/(dfact(2))*( 
    !     .  SLASKM(I,J)-(I1-lambdad(2))*EPSEM(I,J)
    !     .        +detJ/(lambdad(2))*DELTA(I,J) )
           M2(I,J)=Gvec(I,2)*Gvec(J,2)
    !      M3(I,J)=1.d0/(dfact(3))*( 
    !     .   SLASKM(I,J)-(I1-lambdad(3))*EPSEM(I,J)
    !     .        +detJ/(lambdad(3))*DELTA(I,J) )
           M3(I,J)=Gvec(I,3)*Gvec(J,3)

        ENDDO
    ENDDO

    M1V=M_2_V9(M1)
    M2V=M_2_V9(M2)
    M3V=M_2_V9(M3)

    RETURN 
END SUBROUTINE
!============================================================================== 

!============================================================================== 
recursive SUBROUTINE us_solve(LKONV,A,N,B,INVA,X)
      ! Matrix multiplication X=inv(A)*B
      ! Gauss-Jordan elimination with pivoting
      IMPLICIT NONE
      INTEGER N,LKONV
      DOUBLE PRECISION A(N,N),B(N),INVA(N,N),X(N)
      X=B;INVA=A;LKONV=1  
        CALL us_gausjor(LKONV,INVA,N,N,X,1,1)
      RETURN 

end subroutine us_solve
!============================================================================== 

!============================================================================== 
recursive SUBROUTINE us_gausjor(LKONV,A,N,NP,B,M,MP)
       !---------------------------------------------------
      ! SUBROUTINE FROM NUMERICAL RECIPES IN FORTRAN p.30
      ! Gauss-Jordan elimination with pivoting
      ! SOLVES A x= B
      ! INPUT: A(NP,NP), B(NP,N)
      ! OUTPUT: inverse matrix as A
      !         solution vector as B 
      !         LKOND=0 if singular matrix
      IMPLICIT LOGICAL (A-Z)
      INTEGER M,MP,N,NP,NMAX,LKONV
      DOUBLE PRECISION A(NP,NP),B(NP,MP)
      PARAMETER (NMAX=50)
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),IPIV(NMAX)
      DOUBLE PRECISION BIG,DUM,PIVINV      
      LKONV=1
      DO J=1,N
        IPIV(J)=0
      ENDDO
      DO I=1,N
        BIG=0.
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
              IF(IPIV(K).EQ.0)THEN
                IF (ABS(A(J,K)).GE.BIG) THEN
                  BIG=ABS(A(J,K))
                  IROW=J;ICOL=K
                ENDIF
              ELSEIF (IPIV(K).GT.1) THEN
                WRITE(*,*) '% singular matrix in us_gausjor'
                LKONV=0; GOTO 666
              ENDIF
            ENDDO
          ENDIF            
        ENDDO
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
          ENDDO
          DO L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
          ENDDO
        ENDIF           
        INDXR(I)=IROW
        INDXC(I)=ICOL          
        IF (A(ICOL,ICOL).EQ.0.d0) THEN          
          WRITE(*,*) '% singular matrix in us_gausjor'
          LKONV=0;  GOTO 666
        ENDIF           
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.d0         
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
        ENDDO
        DO L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
        ENDDO
        DO LL=1,N
          IF (LL.NE.ICOL) THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.d0
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            ENDDO
            DO L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
            ENDDO
          ENDIF
        ENDDO          
      ENDDO
      DO L=N,1,-1
        IF (INDXR(L).NE.INDXC(L)) THEN           
          DO K=1,N              
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
          ENDDO
        ENDIF
      ENDDO
 666  CONTINUE
      RETURN
end subroutine us_gausjor
!============================================================================== 
      
end module    