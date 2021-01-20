subroutine inormal(mean, cov, dim, size, output, seed)
  ! this subroutine is a wrapper for vsrnggaussianmv which is the mkl
  ! function to produce an array of real numbers of given dimension which
  ! are drawn from the normal distribution with a given mean.
  ! Inputs:
  !     `dim`    : dimension of the distribution (integer)
  !     `mean`   : mean of the distribution (1D real array of shape `dim`)
  !     `cov`    : covariance matrix (2D real array of shape `dim` x `dim`)
  !     `size`   : number of samples (integer)
  !     `output` : (2D real array of shape `dim` x `size`) output
  !     `seed`   : seed for random stream (integer)
  ! matrix `cov`. `size` is the length of the array and `output` is a real
  ! array of shape `dim` x `output`.
  !
  ! written by Masoud Ghaderi Zefreh
  ! first revision : 1 November 2019

  use constants
  use mkl_vsl_type
  use mkl_vsl

  implicit none

  integer, intent(in) :: dim, size
  real, dimension(1:dim), intent(in) :: mean
  real, dimension(1:dim, 1:dim), intent(in) :: cov
  real, dimension(1:dim, 1:size), intent(out) :: output
  integer, intent(in), dimension(:), optional :: seed

  integer, dimension(:), allocatable :: seed2
  integer :: errcode, info, brng, method, mstore, seedi
  TYPE (VSL_STREAM_STATE) :: stream

  external :: spotrf

  if (present(seed)) then
     call random_seed(put = seed)
     seedi = seed(1)
  else
     call random_seed(size = info)
     allocate(seed2(info))
     call random_seed(get = seed2)
     seedi = seed2(1)
  end if

  brng = VSL_BRNG_MCG31
  method = VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
  mstore = VSL_MATRIX_STORAGE_FULL

  call spotrf('U',dim,cov,dim,info)
  if (info.lt.0) then
     write(STDERR, *) ' Cholesky factorisation failed!'
     write(STDERR, *) ' Covariance matrix is not positive definite'
     write(STDERR, *) ' Exiting...'
     stop 2
  end if

  errcode = vslnewstream(stream, brng, seedi)
  if (errcode.ne.0) then
     write(STDERR, *) ' Failed to create a new stream'
     write(STDERR, *) ' Exiting...'
     stop 2
  end if

  errcode = vsrnggaussianmv(method, stream, size, output, dim, mstore, mean, cov)
  if (errcode.ne.0) then
     write(STDERR, *), ' Failed to sample for uknown reasons'
     write(STDERR, *), ' Exiting...'
     stop 2
  end if

  errcode = vslDeleteStream(stream)
  if (errcode.ne.0) then
     write(STDERR, *) ' Failed to delete stream'
     write(STDERR, *) ' Exiting...'
     stop 2
  end if
end subroutine inormal
