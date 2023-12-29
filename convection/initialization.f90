module initialization
   use types_and_kinds
   use global
   implicit none
contains

   subroutine init_grid(grid, nx, ny, nz, nGhosts, init_value)
      real(rk), dimension(:, :, :), intent(out) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: init_value
      integer(ik) :: i, j, k
      ! to efficintly traverse the grid in the order of memory layout
      ! loop over 3d array with proper cache efficiency
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               grid(i, j, k) = init_value
            end do
         end do
      end do
   end subroutine init_grid

   ! similar to above subroutine, we now intialize the grid with a 3d gaussian centered in middle of grid
   real(rk) function gaussian3d(i, j, k, mu, sigma)
      real(rk), intent(in) :: mu, sigma
      integer(ik), intent(in) :: i, j, k
      real(rk) :: a = 1.0
      !gaussian3d = a*sin(  2*dacos(-1.d0)* ((i-3)/127.0) )
      gaussian3d = a*exp(-((i - mu)**2 + (j - mu)**2 + (k - mu)**2)/(2*sigma**2))
   end function gaussian3d

   ! use the above function to initialize the grid
   subroutine init_grid_gaussian(grid, nx, ny, nz, nGhosts, init_value_mu, init_value_sigma, magnitude, offset)
      real(rk), dimension(:, :, :), intent(out) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: init_value_mu, init_value_sigma, magnitude, offset
      integer(ik) :: i, j, k
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               grid(i, j, k) = magnitude*gaussian3d(i, j, k, init_value_mu, init_value_sigma) + offset
            end do
         end do
      end do
   end subroutine init_grid_gaussian
   elemental pure real(rk) function logic2dbl(a)
     logical, intent(in) :: a
     if (a) then
       logic2dbl = 1.d0
     else
       logic2dbl = 0.d0
     end if
   end function logic2dbl

   subroutine init_Kelvin_Helmholtz(rho,vx,vy,vz,p, nx, ny, nz, nGhosts,ds)
      real(rk), dimension(:, :, :), intent(inout) :: rho,vx,vy,vz,p
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) ::ds
      real(rk):: s,w,xs,zs
      integer(ik) :: i, k
      w = 0.1_rk
      s = 0.05 / sqrt(2.0_rk)
      do i = 1+nGhosts, nx-nGhosts
         do k = 1+nGhosts, nz-nGhosts
            xs=i*ds
            zs=k*ds
            rho(i,:, k) = 1.0 + logic2dbl(( abs(zs - ds*nz/2) < ds*nz/4 ))
            vx(i,:, k) = -0.5 + logic2dbl(abs(zs - ds*nz/2) < ds*nz/4)
            vz(i,:, k) = w * sin(4.0 * pi * xs) * (exp(- (zs - ds*nz/4)**2 / (2.0 * s**2)) + exp(- (zs - 3*ds*nz/4)**2 / (2.0 * s**2)))
            p(i, :,k) = 3
         end do
     end do
     vy=0
   end subroutine init_Kelvin_Helmholtz
   subroutine init_Rayleigh_Taylor(rho,vx,vy,vz,p, nx, ny, nz, nGhosts,ds)
      real(rk), dimension(:, :, :), intent(inout) :: rho,vx,vy,vz,p
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) ::ds
      real(rk):: s,w,xs,zs,v
      integer(ik) :: i, k
   end subroutine init_Rayleigh_Taylor
   subroutine init_Thermal_Rising_Bubble(rho,vx,vy,vz,p, nx, ny, nz, nGhosts,ds)
      real(rk), dimension(:, :, :), intent(inout) :: rho,vx,vy,vz,p
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) ::ds
      real(rk):: s,w,xs,zs,v
      integer(ik) :: i, k
   end subroutine init_Thermal_Rising_Bubble
   subroutine init_Equilibrium(rho,vx,vy,vz,p, nx, ny, nz, nGhosts,ds)
      real(rk), dimension(:, :, :), intent(inout) :: rho,vx,vy,vz,p
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) ::ds
      real(rk):: s,w,xs,zs,v,mag,rho_0,rho_1,t0,p0
      integer(ik) :: i,j,k
      if (abs(g)<1.0) then
         print*, "No gravity you stupid!!!"
      endif

      rho_1=1.5_rk
      rho_0=1.2_rk
      !rho=rho_
      p0=101000.0
      p=p0
      t0=20.0+273.15
      do k = 1, nz
         do j =1, ny
            do i = 1, nx 
               !How deep are we?
               zs=((nz-2.0*nGhosts)*ds+(0.0*ds)-(k-2.0)*ds)
               !zs=((k-2)*ds)
               !rho(i,j,k)=rho_0
               p(i,j,k)=p0+abs(g)*zs
               rho(i,j,k) = p(i,j,k)/(rs*t0) 
            end do
         end do
      end do

      print*, "Pressure Profile with Ghosts"
      p(:,:,2)=p(:,:,3)+abs(g)*ds
      p(:,:,1)=p(:,:,3)+2_rk*abs(g)*ds
      do k = 1, nz
               zs=((nz-2.0*nGhosts)*ds+(0.5*ds)-(k-2.0)*ds)
         print*,"Index = ",k,"Depth =", zs,"P=" ,p(nx/2,ny/2,k)
      enddo
   end subroutine init_Equilibrium
end module initialization

