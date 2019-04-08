			program vr
				USE gradientconjugue
				use matrices
				implicit none
				integer,parameter :: Nx=1,Ny=1,NxNy=(Nx+2)*(Ny+2)
				integer::i,i1,iN,statinfo,Np,me
				double precision,parameter :: D=1.,hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.,Ly=1.,dt=1.
				
				real,dimension(:),allocatable::nnz,colonnes,ic,xn,ssm
				real,dimension(Nx+2)::x
				real, dimension(Ny+2)::y

			
				do i=1,Nx+2
					x(i)=(i-1)*hx
				end do
				do i=1,Ny+2
					y(i)=(i-1)*hy
				end do
				
				
				
				
				
				call matrice(Nx,Ny,hx,hy,nnz,colonnes,ic)
				call sm(x,y,Nx,Ny,i1,iN,ssm)
				
				print*,ssm
				call gc(nnz,colonnes,ic,ssm,NxNy,Nx,Ny,xn) 
				print*,xn
				print*,nnz
				
		
			
		
			end program
		




		