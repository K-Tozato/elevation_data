module in_data
  double precision, allocatable :: xyz(:,:) 
  double precision :: xl(2), yl(2), dl, dx, dy
  integer :: n0, nx, ny, ndiv
  integer, allocatable :: ain(:,:), num(:,:)
end module

module out_data
  double precision, allocatable :: dat(:,:) 
  integer :: n1
  integer, allocatable :: aout(:,:)
end module


program maketopo
use in_data
use out_data
implicit none
integer :: i, j, k, l, n
double precision :: dum 

write(*,*) 'ndiv ?'
read(*,*) ndiv
 
call read_data

open(11,file='./input/area.txt')
read(11,*) xl(1), xl(2)
read(11,*) yl(1), yl(2)
read(11,*) dl, dum
close(11)

nx = nint(xl(2)-xl(1))/nint(dl)+1 ; ny = nint(yl(2)-yl(1))/nint(dl)+1
dx = xl(2) - xl(1) ; dy = yl(2) - yl(1)
write(*,'(a,f10.2,a)') 'dx = ', xl(2) - xl(1), '[m]'
write(*,'(a,f10.2,a)') 'dy = ', yl(2) - yl(1), '[m]'
write(*,*) 'nx & ny :', nx, ny

n1 = 0 
allocate(dat(nx*ny,3), ain(n0,ndiv*ndiv), aout(nx*ny, 2), num(ndiv, ndiv))
dat(:,:) = 0.d0
do j = 1, ny
  do i = 1, nx
    n1 = n1 + 1
    dat(n1,1) = xl(1) + dble(i-1)*dl
    dat(n1,2) = yl(1) + dble(j-1)*dl
    k = (i-1)/((nx-1)/ndiv) + 1
    l = (j-1)/((ny-1)/ndiv) + 1
    if(i==nx) k = (i-1)/((nx-1)/ndiv)
    if(j==ny) l = (j-1)/((ny-1)/ndiv)
    aout(n1,1) = k ; aout(n1,2) = l
  end do
end do
write(*,*) 'Number of out_data : ', n1

dum = 1.d0/dble(ndiv)
num(:,:) = 0
do n = 1,n0
  if(xyz(n,2)<=yl(1) + dum*dy) l = 1
  if(xyz(n,2)> yl(1) + (1.d0-dum)*dy) l = ndiv
  do j = 2,(ndiv-1)
     if(xyz(n,2)>yl(1)+dble(j-1)*dy*dum .and. xyz(n,2)<=yl(1)+dble(j)*dy*dum) then
       l = j
     end if
  end do
  if(xyz(n,1)<=xl(1) + dum*dx) k = 1
  if(xyz(n,1)> xl(1) + (1.d0-dum)*dx) k = ndiv
  do j = 2,(ndiv-1)
     if(xyz(n,1)>xl(1)+dble(j-1)*dx*dum .and. xyz(n,1)<=xl(1)+dble(j)*dx*dum) then
       k = j
     end if
  end do
  num(k,l) = num(k,l) + 1
  ain(num(k,l),(l-1)*ndiv+k) = n
end do
n = 0
do j = 1,ndiv
  do i = 1,ndiv
    n = n + num(i,j)
  end do
end do
write(*,*) n

do i = 1, n1
  call line_itp(i)
  if(mod(i,10000)==0) write(*,*) i
end do

open(12,file='./output/coordinate.txt',status='replace')
do i = 1,n1
  write(12,'(3f10.2)') dat(i,2), dat(i,1), dat(i,3)
end do
close(12)

end program

subroutine read_data
use in_data
implicit none
integer :: i, j, nn, nda
integer, allocatable :: n(:)
double precision :: dum
character(30) :: cdum
character(30), allocatable :: fname(:)

open(11,file='./input/fname_ndata.txt',status='old')
read(11,*) cdum, nda
allocate(n(nda), fname(nda))

n0 = 0
do i = 1,nda
  read(11,*) fname(i), n(i)
  n0 = n0 + n(i)
end do
close(11)

allocate(xyz(n0,3))
nn = 0

do j = 1, nda
  open(11,file=fname(j),status='old')
  do i = 1,n(j)
    nn = nn + 1
    read(11,*) dum, xyz(nn,1), xyz(nn,2), xyz(nn,3)
  end do
  write(*,*) j, fname(j), nn 
  close(11)
end do

end subroutine

subroutine line_itp(n)
use in_data 
use out_data
implicit none
integer,intent(in) :: n
integer :: i, j, k, l, ii, pp(4), ij(9,2)
double precision :: d0, dd(4), xx, yy, zz(2)

k=0
do j = -1,1
  do i = -1,1
    k = k + 1
    ij(k,1) = aout(n,1)+i
    ij(k,2) = aout(n,2)+j
  end do
end do

dd(:) = 1000.d0
pp(:) = 0
do k = 1,9
  i = ij(k,1) ; j = ij(k,2)
  if(mod(i,ndiv+1)==0 .or. mod(j,ndiv+1)==0) cycle
  
  do l = 1,num(i,j)
    ii = (j-1) * ndiv + i
    xx = xyz(ain(l,ii),1) - dat(n,1) ; yy = xyz(ain(l,ii),2) - dat(n,2)
    d0 = dsqrt(xx**2.d0 + yy**2.d0)
    !write(*,*) 'd0 = ', d0
    if(d0==0.d0) then
      dat(n,3) = xyz(ain(l,ii),3)
      exit
    elseif(d0>20.d0) then
      cycle
    else
      if(xx>=0.d0 .and. yy>=0.d0) then
        if(dd(1) > d0) dd(1) = d0 ; pp(1) = ain(l,ii)
      elseif(xx<=0.d0 .and. yy>=0.d0) then
        if(dd(2) > d0) dd(2) = d0 ; pp(2) = ain(l,ii)
      elseif(xx<=0.d0 .and. yy<=0.d0) then
        if(dd(3) > d0) dd(3) = d0 ; pp(3) = ain(l,ii)
      elseif(xx>=0.d0 .and. yy<=0.d0) then
        if(dd(4) > d0) dd(4) = d0 ; pp(4) = ain(l,ii)
      end if
    end if
  end do
 
  if(d0==0.d0) exit
end do

if(minval(pp)==0 .and. d0>0) then
  write(*,*) n
  write(*,*) 'Error !!', dat(n,1:2)
  write(*,*) xyz(pp(1),:)
  write(*,*) xyz(pp(2),:)
  write(*,*) xyz(pp(3),:)
  write(*,*) xyz(pp(4),:)
end if

if(d0 > 0.d0) then
  xx = (xyz(pp(1),2)-xyz(pp(2),2))/(xyz(pp(1),1)-xyz(pp(2),1))* &
       (dat(n,1) - xyz(pp(2),1)) + xyz(pp(2),2)
  yy = (xyz(pp(4),2)-xyz(pp(3),2))/(xyz(pp(4),1)-xyz(pp(3),1))* &
       (dat(n,1) - xyz(pp(3),1)) + xyz(pp(3),2)
  zz(1) = (xyz(pp(1),3)-xyz(pp(2),3))/(xyz(pp(1),1)-xyz(pp(2),1))* &
          (dat(n,1) - xyz(pp(2),1)) + xyz(pp(2),3)
  zz(2) = (xyz(pp(4),3)-xyz(pp(3),3))/(xyz(pp(4),1)-xyz(pp(3),1))* &
          (dat(n,1) - xyz(pp(3),1)) + xyz(pp(3),3)
  dat(n,3) = (zz(1)-zz(2))/(xx-yy)*(dat(n,2)-yy) + zz(2)
end if

!write(*,*) xyz(pp(1),:)
!write(*,*) xyz(pp(2),:)
!write(*,*) xyz(pp(3),:)
!write(*,*) xyz(pp(4),:)
!write(*,*) dat(n,:)

end subroutine


