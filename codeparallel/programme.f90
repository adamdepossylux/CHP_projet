			program programme
				use mpi
				use fonctions
				use matrices
				use gradientconjugue
				implicit none
				integer::D
				integer,parameter :: Nx=100,Ny=100,NxNy=(Nx+2)*(Ny+2)
				integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,Nt,n1,n2,n3
				double precision,parameter :: hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.,Ly=1.,dt=1.
				double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx
				double precision, dimension(:),allocatable :: solutionexacte,smf,solution,Uinitial
				integer, dimension( MPI_STATUS_SIZE) :: status
				double precision,dimension(:,:),allocatable::A,AA
				double precision:: produitscalaire,delta,temps_debut,temps_fin,temps_fin_max
				double precision,dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(100)::temps
				allocate(nnz((nx+2)*2+5*nx*ny+ny*2))
			  allocate(colonnes((nx+2)*2+5*nx*ny+ny*2))
!			  allocate(ic((nx+2)*2+nx*ny+ny*2))
				delta=1d0
				Nt=3

					call MPI_INIT(statinfo)

					temps_debut= MPI_WTIME()

					call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
					call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

					if(Np>=Nx+2) then
						if(Me==0) then
						print*,"Np = ",Np," est trop grand"
					end if

					else
					call charge(me,Np,(Nx+2)*(Ny+2),i1,iN)
					allocate(smf(i1:iN))
					allocate (b(i1:iN))
					allocate (solution(i1:iN))
					allocate (ic(1:iN-i1+1))


					do i=1,Nx+2
						x(i)=(i-1)*hx
					end do
					do i=1,Ny+2
						y(i)=(i-1)*hy
					end do
					do i=1,100
						temps(i)=(i-1)*delta
					end do

					call ajoutf(x,y,Nx,Ny,i1,iN,smf,temps(1),Ly,Lx)
					call sme(x,y,Nx,Ny,i1,iN,Uinitial)
					ic=0
					nnz=0
					colonnes=0
					!call laplaciencsr(hx,hy,dt,ic,colonnes,nnz,Nx+1,Ny+1,i1,iN)
					call matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)

!					n2=size(ic)
!					open(unit=2, file="ic.txt",form="formatted",access="sequential")
!					do i=1,n2
!						write(2,*)ic(i)
!					enddo
!					close(2)
					!deallocate(xn)
					b=Uinitial  + delta*smf
					b=Uinitial
					deallocate(smf)
					do i=1,Nt
						allocate(smf(i1:iN))
						allocate (xn(i1:iN))
						call gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xn,i1,iN)
						!call ajoutf(x,y,Nx,Ny,i1,iN,smf,temps(i),Ly,Lx)
						!b=xn+delta*smf
						b=xn
						deallocate(smf)
						solution=xn
						deallocate(xn)
					end do

!					do i=i1,iN
!						print*,solution(i)
!					enddo
					deallocate(ic)
					deallocate(solution)
					temps_fin= (MPI_WTIME()-temps_debut)

					call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
					if (me == 0) then
						print*,"Temps : ",temps_fin_max," secondes"

					end if
					end if
					call MPI_FINALIZE(statinfo)


end program programme
