module gradientconjugue
  use matrices
  implicit none
  contains


        subroutine gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xn)
          implicit none
          double precision , dimension(:),allocatable ::b,xn, p, Ap, r,r1, produitmatriciel,nnz,colonnes,ic


          double precision :: alpha, rr, rr1, pAp,produitvectoriel
          integer :: j,i,k,statinfo,Np,me, NxNy,Nx,Ny,i1,iN,he
          allocate (r(NxNy))
          allocate (r1(NxNy))

          allocate (p(NxNy))
          allocate(xn(NxNy))

          xn=0
          r = b
          p = b
          call produitvec(Nx,Ny,r,r,produitvectoriel,statinfo,me,i1,iN)!calcul r_k'*r_k
          rr = produitvectoriel
          DO k = 1, NxNy
            call produitvec(Nx,Ny,r,r,produitvectoriel,statinfo,me,i1,iN)
            rr = produitvectoriel
            if (sqrt(rr)<=0.0001)then ! condition pour sortir de l'algortihme
              exit
            end if
            call produitmat(nnz,colonnes,ic,p,Nx,Ny,produitmatriciel,i1,iN,statinfo)
            Ap = produitmatriciel
            deallocate(produitmatriciel)


            call produitvec(Nx,Ny,p,Ap,produitvectoriel,statinfo,me,i1,iN)

            pAp=produitvectoriel
            alpha = rr / pAp! rr designe rk'*rk
            xn =xn + alpha * p

            r1 = r-alpha * Ap

            call produitvec(Nx,Ny,r1,r1,produitvectoriel,statinfo,me,i1,iN)
            rr1=produitvectoriel
            p=r1+p*(rr1/rr)
            r=r1
          end do



      end subroutine gc


    subroutine produitmat(nnz,colonnes,ic,xx,Nx,Ny,produitmatriciel,i1,iN,statinfo)
      implicit none
      integer::Nx,Ny,i,i1,iN,he,me,statinfo,Np,j1,jN,k1,k2,k3,j
      double precision,dimension(:,:),allocatable::A
      double precision,dimension(:),allocatable::xx,pm,produitmatriciel,nnz,colonnes,ic

      allocate(produitmatriciel(1:Nx*Ny))
      !produitmatriciel=0
      do i=1,Nx*Ny
        produitmatriciel(i)=0
        k1=ic(i)!indice du première éléments non nul de la ligne
        k2=ic(i+1)-1!indice du dernière éléement de la lgine non nul
        do j=k1,k2
          k3=(colonnes(j))!colonnes de la matrice qui correspond à l'indice dans x
          produitmatriciel(i)=produitmatriciel(i) +nnz(j)*xx(k3)!j correponds à l'indice de la colonnes considérée
        end do
      end do
    end subroutine produitmat

    subroutine produitvec(Nx,Ny,v,vv,produitvectoriel,statinfo,me,i1,iN)
      implicit none
      integer::Nx,Ny,	i1,iN,i,me,statinfo,Np
      double precision,dimension(nx*ny)::v,vv
      double precision:: pv,produitvectoriel

      produitvectoriel=0
      do i=1,Nx*Ny
        produitvectoriel= produitvectoriel + v(i)*vv(i)
      end do

    end subroutine produitvec




end module gradientconjugue
