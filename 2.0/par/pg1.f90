program programme
	use mpi
	use fonctions
	use matrices
	use gradientconjugue
	use parametres
	implicit none
	character*13 ::name
	integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,n2
	double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx,usol
	double precision, dimension(:),allocatable :: solutionexacte,smf,U,Uo,Uinitial
	integer, dimension( MPI_STATUS_SIZE) :: status
	double precision,dimension(:,:),allocatable::A,AA
	double precision:: produitscalaire,temps_debut,temps_fin,temps_fin_max
	double precision,dimension(Nx+2)::x
	double precision, dimension(Ny+2)::y
	double precision, dimension(100)::temps



	call MPI_INIT(statinfo)

	temps_debut= MPI_WTIME()

	call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
	call charge(me,Np,Nx*Ny,i1,iN)
	allocate (b(i1:iN))
	allocate (U(i1:iN))
	allocate (Uo(i1:iN))
	allocate(ssm(i1:iN))
	allocate(usol(i1:iN))
	! axe des abscisses
	do i=1,Nx+2
		x(i)=(i-1)*hx
	end do
	! axe des ordonnées
	do i=1,Ny+2
		y(i)=(i-1)*hy
	end do
	! échelle de temps
	do i=1,100
		temps(i)=(i-1)*dt
	end do
	b=0
	Uo=0
	!appel de la matrice du laplacien stockée en csr
  call laplaciencsr(hx,hy,dt,D,ic,colonnes,nnz,Nx,Ny,i1,iN)
	!appel du second membre pour chaque processus
			call sm2(x,y,Nx,Ny,ssm,dt,Uo,hx,hy,i1,iN)
			b=ssm
			do i=1,Nt
				!appel du gradient conjuguée
				call gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xn,i1,iN)
				U=xn
				call sm2(x,y,Nx,Ny,ssm,dt,U,hx,hy,i1,iN)
				b=ssm
				deallocate(xn)
			end do
			deallocate(ssm)
			call uex1m(x,y,Nx,Ny,i1,iN,usol)
			temps_fin= (MPI_WTIME()-temps_debut)
			!calcul de l'erreur
			print*,"erreur = ",abs(maxval(U-usol))
			call Rename(Me,name)
			!pour connaitre le temps de chaque processus et voir que le parallélisme distribue équitablement le travail à faire
			!print*,"Temps : ",me," ",temps_fin," secondes"
			open (unit =me, file=name)
			write(me,*) U
			call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
			if(me==0) then
				print*,"Temps : ",temps_fin_max," secondes"
			endif

		call MPI_FINALIZE(statinfo)

	contains

	subroutine Rename(Me,name)
		implicit none
		integer :: Me
		character*13 ::name
		character*3 :: tn
		integer :: i1,i2,i3
		i1 = Me/100
		i2 =( Me - 100*i1)/10
		i3 = Me - 100*i1 -10*i2
		tn = char(i1+48)//char(i2+48)//char(i3+48)
		name='sol'//tn//'.dat'
	end subroutine Rename


end program programme
