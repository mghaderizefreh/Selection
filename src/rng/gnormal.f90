subroutine gnormal(mean, cov, dim, N, output, seed)
  ! this subroutine produces normal distribution using boxMuller transformation


  ! Inputs:
  !     `dim`    : dimension of the distribution (integer)
  !     `mean`   : mean of the distribution (1D real array of shape `dim`)
  !     `cov`    : covariance matrix (2D real array of shape `dim` x `dim`)
  !     `N`   : number of samples (integer)
  !     `output` : (2D real array of shape `dim` x `N`) output
  !     `nseed`  : seed for random stream (integer)
  ! matrix `cov`. `N` is the length of the array and `output` is a real
  ! array of shape `dim` x `output`.
  !
  ! written by Masoud Ghaderi Zefreh
  ! first revision : 1 December 2020
  
  use constants

  implicit none
  integer, intent(in) :: dim, N
  real, dimension(1:dim), intent(in) :: mean
  real, dimension(1:dim, 1:dim), intent(in) :: cov
  real, dimension(1:N, 1:dim), intent(inout) :: output
  integer, intent(in), dimension(:), optional :: seed

  integer :: dim2, i, j
  real :: sd, val
  real, dimension(dim, dim) :: covcopy
  real, dimension(:,:), allocatable :: uniform

  covcopy = cov
  if (present(seed)) then
     call random_seed(put = seed)
  end if

  if (dim == 1) then
     sd = sqrt(cov(1,1))
     do i = 1, N
        call normdev(val)
        output(i,1) = mean(1) + sd * val
     end do
     return
  else

     dim2 = dim
     if (mod(dim, 2) .eq. 1) dim2 = dim2 + 1

     ! making uniform distribution
     allocate(uniform(N, dim2))
     call random_number(uniform)

     ! making independent pairs of normally distributed data
     do i = 1, dim, 2
        j = i + 1
        output(:,i) = sqrt(-2 * log(uniform(:,i))) * cos(2 * PI * uniform(:,j))
     end do
     do j = 2, dim, 2
        i = j - 1
        output(:,j) = sqrt(-2 * log(uniform(:,i))) * sin(2 * PI * uniform(:,j))
     end do

     if (dim .ne. dim2)then
        deallocate(uniform)
        allocate(uniform(N,dim))
     end if
     uniform = output

     ! cholesky factorisation
     call spotrf('L',dim,covcopy,dim,i)
     if (i.gt.0) then
        write(*,*) ' Cholesky factorisation failed!'
        write(*,*) ' Covariance matrix is not positive definite'
        write(*,*) ' Exiting...'
        stop
     elseif (i.lt.0) then
        write(*, *) "some other shit happened"
        stop
     end if
     do i = 1, dim - 1
        do j = i + 1, dim
           covcopy(i,j) = 0.0
        end do
     end do

     call sgemm('N','N',N,dim,dim,1.0,uniform,N,covcopy,dim,0.0,output,N)

     forall( i = 1:N)
        output(i,1:DIM) = mean(1:DIM) + output(i,1:DIM)
     end forall
  end if

end subroutine gnormal

subroutine normdev(val)
! making a random varaible form standard normal
! copy of ricardo's code (goto replaced by while)
! Routine 'gasdev' ,Numerical Recipes, page 203
  implicit none
  save
  real    ::val,v1,v2,fac, r = 2.
  real    ::gset =0.0
  integer ::iset =0

  if (iset .eq. 0) then
     do while(r.ge.1..or.r.eq.0.)
        call random_number(v1)
        call random_number(v2)
        v1 = 2. * v1 - 1.
        v2 = 2. * v2 - 1.
        r = v1 ** 2 + v2 ** 2
     end do
     fac = sqrt(-2. * log(r) / r)
     gset= v1 * fac
     val = v2 * fac
     iset= 1
  else
     val = gset
     iset = 0
  end if
end subroutine normdev

