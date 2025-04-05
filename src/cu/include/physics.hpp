/* File:   eulernv.cu
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 * */
#pragma once
#include "constants.h"
#include "include/matrix3d.hpp"
#include "kernels.h"
#include "spdlog/stopwatch.h"
#include <array>
#include <cmath>
#include <cuda_device_runtime_api.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace EULERCFD {
template <typename T>
void init_khi(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid) {

  EULERCFD::Grid<T,
                 GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                             EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                             EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                             EULERCFD::CONSTS::LZ},
                 BACKEND::HOST>
      hostgrid;
  hostgrid.p.fill(T(2.5));
  hostgrid.vy.fill(T(0));

  auto logic2dbl = [](auto condition) -> T { return condition ? T(1) : T(0); };

  T kappa = T(0.1);
  T s = T(0.05) / std::sqrt(T(2.0));

  // Iterate over the grid
  for (std::size_t i = EULERCFD::CONSTS::NGHOSTS;
       i < hostgrid.nx() - EULERCFD::CONSTS::NGHOSTS; ++i) {
    for (std::size_t j = EULERCFD::CONSTS::NGHOSTS;
         j < hostgrid.ny() - EULERCFD::CONSTS::NGHOSTS; ++j) {
      for (std::size_t k = EULERCFD::CONSTS::NGHOSTS;
           k < hostgrid.nz() - EULERCFD::CONSTS::NGHOSTS; ++k) {
        T xs = static_cast<T>(i) * hostgrid.dsx();
        T zs = static_cast<T>(k) * hostgrid.dsz();

        hostgrid.rho(i, j, k) =
            T(1) +
            logic2dbl((std::abs(zs - hostgrid.dsx() * hostgrid.nz() / T(2.0f)) <
                       hostgrid.dsx() * hostgrid.nz() / T(4.0)));

        hostgrid.vx(i, j, k) =
            T(-0.5f) +
            logic2dbl(std::abs(zs - hostgrid.dsx() * hostgrid.nz() / T(2)) <
                      hostgrid.dsx() * hostgrid.nz() / T(4));

        hostgrid.vz(i, j, k) =
            kappa * std::sin(T(4) * M_PI * xs) *
            (std::exp(
                 -std::pow(zs - hostgrid.dsx() * hostgrid.nz() / T(4), T(2)) /
                 (T(2) * std::pow(s, T(2)))) +
             std::exp(
                 -std::pow(zs - T(3) * hostgrid.dsx() * hostgrid.nz() / T(4),
                           T(2)) /
                 (T(2) * std::pow(s, T(2)))));
      }
    }
  }
  simgrid.rho = hostgrid.rho;
  simgrid.p = hostgrid.p;
  simgrid.vx = hostgrid.vx;
  simgrid.vy = hostgrid.vy;
  simgrid.vz = hostgrid.vz;
}

template <typename T>
void init_trb(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid) {

  EULERCFD::Grid<T,
                 GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                             EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                             EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                             EULERCFD::CONSTS::LZ},
                 BACKEND::HOST>
      hostgrid;

  std::size_t nx = hostgrid.nx();
  std::size_t ny = hostgrid.ny();
  std::size_t nz = hostgrid.nz();
  T ds = hostgrid.dsx(); // Assuming uniform spacing in all directions
  std::size_t ghosts = EULERCFD::CONSTS::NGHOSTS;

  hostgrid.vx.fill(T(0.0));
  hostgrid.vy.fill(T(0.0));
  hostgrid.vz.fill(T(0.0));

  const T p0 = T(101000.0); // Pressure at sea level
  const T t0 = T(293.0);    // Temperature in Kelvin
  constexpr T G = EULERCFD::CONSTS::G;
  constexpr T RS = EULERCFD::CONSTS::RS;

  // Build artificial equilibrium
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        T zs = static_cast<T>(nz * ds) - static_cast<T>(k) * ds;
        hostgrid.p(i, j, k) = p0 + std::abs(G) * zs;
        hostgrid.rho(i, j, k) = hostgrid.p(i, j, k) / (RS * t0);
      }
    }
  }

  // Add the bubble
  const T rad_outer = T(25.0);
  const T rad_inner = T(20.0);
  const T center_x = T(EULERCFD::CONSTS::NX) / 2.0;
  const T center_y = T(EULERCFD::CONSTS::NY) / 2.0;
  const T center_z = T(EULERCFD::CONSTS::NZ) / 2.0;

  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        T mag = std::sqrt(std::pow(static_cast<T>(i) - center_x, 2) +
                          std::pow(static_cast<T>(j) - center_y, 2) +
                          std::pow(static_cast<T>(k) - center_z, 2));

        if (mag < rad_outer) {
          hostgrid.rho(i, j, k) = T(1.1);
        }
        if (mag < rad_inner) {
          hostgrid.rho(i, j, k) = T(0.8);
        }
      }
    }
  }

  simgrid.rho = hostgrid.rho;
  simgrid.p = hostgrid.p;
  simgrid.vx = hostgrid.vx;
  simgrid.vy = hostgrid.vy;
  simgrid.vz = hostgrid.vz;
}

template <typename T>
void init_sod(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid) {

  EULERCFD::Grid<T,
                 GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                             EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                             EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                             EULERCFD::CONSTS::LZ},
                 BACKEND::HOST>
      hostgrid;

  std::size_t nx = hostgrid.nx();
  std::size_t ny = hostgrid.ny();
  std::size_t nz = hostgrid.nz();
  T ds = hostgrid.dsx();
  std::size_t ghosts = EULERCFD::CONSTS::NGHOSTS;

  hostgrid.vx.fill(T(0.0));
  hostgrid.vy.fill(T(0.0));
  hostgrid.vz.fill(T(0.0));
  assert(EULERCFD::CONSTS::GRAVITY == false && "Test does not support gravity");

  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        T xs = static_cast<T>(i * ds);
        if (xs < EULERCFD::CONSTS::LX / 2.0) {
          hostgrid.p(i, j, k) = 1.0;
          hostgrid.rho(i, j, k) = 1.0;
        } else {
          hostgrid.p(i, j, k) = 0.1;
          hostgrid.rho(i, j, k) = 0.125;
        }
      }
    }
  }

  simgrid.rho = hostgrid.rho;
  simgrid.p = hostgrid.p;
  simgrid.vx = hostgrid.vx;
  simgrid.vy = hostgrid.vy;
  simgrid.vz = hostgrid.vz;
}

template <typename T>
void init_greenhouse(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid) {

  EULERCFD::Grid<T,
                 GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                             EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                             EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                             EULERCFD::CONSTS::LZ},
                 BACKEND::HOST>
      hostgrid;

  std::size_t nx = hostgrid.nx();
  std::size_t ny = hostgrid.ny();
  std::size_t nz = hostgrid.nz();
  T ds = hostgrid.dsx(); // Assuming uniform spacing in all directions
  std::size_t ghosts = EULERCFD::CONSTS::NGHOSTS;

  hostgrid.vx.fill(T(0.0));
  hostgrid.vy.fill(T(0.0));
  hostgrid.vz.fill(T(0.0));

  const T p0 = T(101000.0); // Pressure at sea level
  const T t0 = T(293.0);    // Temperature in Kelvin
  constexpr T G = EULERCFD::CONSTS::G;
  constexpr T RS = EULERCFD::CONSTS::RS;

  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        T zs = static_cast<T>(nz * ds) - static_cast<T>(k) * ds;
        hostgrid.p(i, j, k) = p0 + std::abs(G) * zs;
        hostgrid.rho(i, j, k) = hostgrid.p(i, j, k) / (RS * t0);
      }
    }
  }

  simgrid.rho = hostgrid.rho;
  simgrid.p = hostgrid.p;
  simgrid.vx = hostgrid.vx;
  simgrid.vy = hostgrid.vy;
  simgrid.vz = hostgrid.vz;
}

template <typename T>
void init_square(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid) {

  EULERCFD::Grid<T,
                 GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                             EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                             EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                             EULERCFD::CONSTS::LZ},
                 BACKEND::HOST>
      hostgrid;

  constexpr std::size_t nx = EULERCFD::CONSTS::NX;
  constexpr std::size_t ny = EULERCFD::CONSTS::NY;
  constexpr std::size_t nz = EULERCFD::CONSTS::NZ;
  T ds = hostgrid.dsx(); // Assuming uniform spacing in all directions

  const float side = nx / 4;

  hostgrid.vx.fill(T(EULERCFD::CONSTS::NX));
  hostgrid.vz.fill(T(EULERCFD::CONSTS::NX));
  hostgrid.vy.fill(T(0.0));
  hostgrid.p.fill(T(101000.0));
  hostgrid.rho.fill(T(1.0));

  const T p0 = T(101000.0); // Pressure at sea level
  const T t0 = T(293.0);    // Temperature in Kelvin

  // Add the square
  const T center_x = T(EULERCFD::CONSTS::NX) / 2.0;
  const T center_y = T(EULERCFD::CONSTS::NY) / 2.0;
  const T center_z = T(EULERCFD::CONSTS::NZ) / 2.0;

  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        if (std::abs(center_x - i) < side && std::abs(center_x - j) < side &&
            std::abs(center_x - k) < side) {
          hostgrid.rho(i, j, k) *= T(2.0);
        }
      }
    }
  }

  simgrid.rho = hostgrid.rho;
  simgrid.p = hostgrid.p;
  simgrid.vx = hostgrid.vx;
  simgrid.vy = hostgrid.vy;
  simgrid.vz = hostgrid.vz;
}

template <typename T, GridInfo<T> Info, BACKEND Backend>
void build_sdf(Matrix3d<T, Info, Backend> &sdf_object) {
  static_assert(Backend == BACKEND::HOST &&
                "This can only be used on host grids tou idiot! So give me a "
                "copy of that and "
                "then copy back to device yourself!");

  constexpr T radious = EULERCFD::CONSTS::RADIOUS;
  constexpr T nx = EULERCFD::CONSTS::NX;
  constexpr T ny = EULERCFD::CONSTS::NY;
  constexpr T nz = EULERCFD::CONSTS::NZ;
  constexpr T cx = (nx / 2) * EULERCFD::CONSTS::DELTA;
  constexpr T cy = (ny / 2) * EULERCFD::CONSTS::DELTA;
  constexpr T cz = (nz / 2) * EULERCFD::CONSTS::DELTA;
  constexpr std::array<T, 3> c{cx, cy, cz};

  for (std::size_t i = 1; i < nx - 1; ++i) {
    for (std::size_t j = 1; j < ny - 1; ++j) {
      for (std::size_t k = 1; k < nz - 1; ++k) {

        const auto r = sim2real<T>(i, j, k);
        const auto x = r[0];
        const auto y = r[1];
        const auto z = r[2];
        const auto val = sdf({x, y, z}, {cx, cy, cz}, radious);
        if (val > 0) {
          sdf_object(i, j, k) = 100.0;
        }
        if (val < 0) {
          sdf_object(i, j, k) = -100.0;
        }
      }
    }
  }

  std::unordered_set<std::size_t> surface;
  for (std::size_t i = 1; i < nx - 1; ++i) {
    for (std::size_t j = 1; j < ny - 1; ++j) {
      for (std::size_t k = 1; k < nz - 1; ++k) {

        if (sdf_object(i,j,k)<=0.0 ){
          continue;
        }

        std::size_t ii=i;
        //X+
        while (ii < nx - 1) {
          ii++;
          if (sdf_object(ii, j, k) <0.0) {
            surface.insert(id_f(ii,j,k));
            break;
          }
        }
        
        ii=i;
        //X-
        while (ii >1) {
          ii--;
          if (sdf_object(ii, j, k) <0.0) {
            surface.insert(id_f(ii,j,k));
            break;
          }
        }
        std::size_t jj=j;
        //Y+
        while (jj < ny - 1) {
          jj++;
          if (sdf_object(i, jj, k) <0.0) {
            surface.insert(id_f(i,jj,k));
            break;
          }
        }
        
        jj=j;
        //Y-
        while (jj >1) {
          jj--;
          if (sdf_object(i, jj, k) <0.0) {
            surface.insert(id_f(i,jj,k));
            break;
          }
        }
        std::size_t kk = k;
        // Z+
        while (kk < nz - 1) {
          kk++;
          if (sdf_object(i, j, kk) <0.0) {
            surface.insert(id_f(i,j,kk));
            break;
          }
        }
        
        kk=k;
        //Z-
        while (kk >1) {
          kk--;
          if (sdf_object(i, j, kk) <0.0) {
            surface.insert(id_f(i,j,kk));
            break;
          }
        }
      }
    }
  }

  for (auto cand:surface){
    std::size_t i,j,k;
    id_f_t(cand,i,j,k);
    sdf_object(i,j,k)=0.0;
  }
  return;
}

template <typename T>
void init_tunnel(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid) {

  EULERCFD::Grid<T,
                 GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                             EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                             EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                             EULERCFD::CONSTS::LZ},
                 BACKEND::HOST>
      hostgrid;
  // Set on host
  hostgrid.rho.fill(T(1.0*1.225));
  hostgrid.p.fill(T(101325.0));
  hostgrid.vx.fill(T(0.0));
  hostgrid.vy.fill(T(0.0));
  hostgrid.vz.fill(T(0.0));
  hostgrid.sdf_object.fill(T(100));
  build_sdf(hostgrid.sdf_object);

  // Copy to device
  simgrid.rho = hostgrid.rho;
  simgrid.p = hostgrid.p;
  simgrid.vx = hostgrid.vx;
  simgrid.vy = hostgrid.vy;
  simgrid.vz = hostgrid.vz;
  simgrid.sdf_object = hostgrid.sdf_object;
}

template <typename T, std::size_t N>
void apply_boundaries(std::array<T *, N> src, T const *sdf_mask,
                      std::size_t len, std::array<std::size_t, 2> lp) {
  kernel_update_ghosts<T, 0, EULERCFD::CONSTS::bcs[0],
                       EULERCFD::CONSTS::NGHOSTS><<<lp[0], lp[1]>>>(src, len);
  kernel_update_ghosts<T, 1, EULERCFD::CONSTS::bcs[1],
                       EULERCFD::CONSTS::NGHOSTS><<<lp[0], lp[1]>>>(src, len);
  kernel_update_ghosts<T, 2, EULERCFD::CONSTS::bcs[2],
                       EULERCFD::CONSTS::NGHOSTS><<<lp[0], lp[1]>>>(src, len);
  kernel_update_ghosts<T, 3, EULERCFD::CONSTS::bcs[3],
                       EULERCFD::CONSTS::NGHOSTS><<<lp[0], lp[1]>>>(src, len);
  kernel_update_ghosts<T, 4, EULERCFD::CONSTS::bcs[4],
                       EULERCFD::CONSTS::NGHOSTS><<<lp[0], lp[1]>>>(src, len);
  kernel_update_ghosts<T, 5, EULERCFD::CONSTS::bcs[5],
                       EULERCFD::CONSTS::NGHOSTS><<<lp[0], lp[1]>>>(src, len);
  kernel_apply_sdf_object_bcs<T, 5, EULERCFD::CONSTS::bcs[5],
                              EULERCFD::CONSTS::NGHOSTS>
      <<<lp[0], lp[1]>>>(src, sdf_mask, len);
}

template <typename T, GridInfo<T> G>
std::array<T, 3> compute_step(
    EULERCFD::Grid<T,
                   GridInfo<T>{EULERCFD::CONSTS::NX, EULERCFD::CONSTS::NY,
                               EULERCFD::CONSTS::NZ, EULERCFD::CONSTS::NGHOSTS,
                               EULERCFD::CONSTS::LX, EULERCFD::CONSTS::LY,
                               EULERCFD::CONSTS::LZ},
                   BACKEND::DEVICE> &simgrid,
    T dt) {

  float max_tp = 0.0;
  float max_flops = 0.0;
  constexpr std::size_t KB = 1024ul;
  constexpr std::size_t TB = KB * KB * KB * KB;
  constexpr std::size_t primitives_size_bytes =
      (5ul) * sizeof(T) * EULERCFD::CONSTS::NX * EULERCFD::CONSTS::NY *
      EULERCFD::CONSTS::NZ;
  constexpr std::size_t gradient_size_bytes =
      (15ul) * sizeof(T) * EULERCFD::CONSTS::NX * EULERCFD::CONSTS::NY *
      EULERCFD::CONSTS::NZ;
  constexpr auto lp = launch_params(G.size());
  constexpr std::size_t workers = lp[0] * lp[1];
  std::array<cudaStream_t, 3> streams;
  cudaStreamCreate(&streams[0]);
  cudaStreamCreate(&streams[1]);
  cudaStreamCreate(&streams[2]);

  PROFILE_START("Compute Step");
  PROFILE_START("Calc Primitives");
  // Primitives
  spdlog::stopwatch sw0;
  kernel_calc_primitives<T, G.dv()><<<lp[0], lp[1]>>>(
      simgrid.get_conserved_pointers(), simgrid.get_primitive_pointers(),
      simgrid.sdf_object.data(), simgrid.size());
  cudaDeviceSynchronize();
  float tp = 2 * primitives_size_bytes / sw0.elapsed().count() / TB;
  float flops = 16 * workers / sw0.elapsed().count() / 1e9;
  max_tp = std::max(tp, max_tp);
  max_flops = std::max(flops, max_flops);
  spdlog::debug("KERNEL::calc_primitives [{0:d} x {1:d}] in {2:f} s ({3:f} "
                "TB/s | {4:f} GFLOPS) .",
                lp[0], lp[1], sw0, tp, flops);
  PROFILE_END();

  // Primitives BCs
  PROFILE_START("BCs Primitives");
  spdlog::stopwatch sw1;
  apply_boundaries<T>(simgrid.get_primitive_pointers(),
                      simgrid.sdf_object.data(), simgrid.size(), lp);
  // apply_boundaries<T>(simgrid.get_xtr, simgrid.size(), lp);
  spdlog::debug("KERNEL::apply_boundaries [{0:d} x {1:d}] in {2:f} s.", lp[0],
                lp[1], sw1);
  PROFILE_END();

  PROFILE_START("Calc timestep");
  // Calc and Reduce timestep
  spdlog::stopwatch sw2;
  dt = calc_timestep<T>(simgrid.get_primitive_pointers(),
                        simgrid.sdf_object.data(), simgrid.size(), lp);
  spdlog::debug("KERNEL::cal_timestep [{0:d} x {1:d}] in {2:f} s.", lp[0],
                lp[1], sw2);
  PROFILE_END();

  PROFILE_START("Calc Drag Coeff");
  // Calc and Redfuce drag coeff of immersed object
  spdlog::stopwatch sw21;
  T drag_coeff = calc_drag_coeff(simgrid.get_primitive_pointers(),
                       simgrid.sdf_object.data(), simgrid.size(), lp);
  spdlog::debug("KERNEL::cal_drag_coeff [{0:d} x {1:d}] in {2:f} s.", lp[0],
                lp[1], sw21);
  simgrid.scratch_space[0]=drag_coeff;
  PROFILE_END();
  
  PROFILE_START("Calc Gradients");
  // Gradients
  spdlog::stopwatch sw3;
  kernel_calc_gradients<T, G.dsx()><<<lp[0], lp[1]>>>(
      simgrid.get_primitive_pointers(), simgrid.get_gradients_pointers(),
      simgrid.sdf_object.data(), simgrid.size());
  cudaDeviceSynchronize();
  tp = (primitives_size_bytes + gradient_size_bytes) / sw3.elapsed().count() /
       TB;
  flops = 30 * workers / sw3.elapsed().count() / 1e9;
  max_tp = std::max(tp, max_tp);
  max_flops = std::max(flops, max_flops);
  spdlog::debug("KERNEL::cal_gradients [{0:d} x {1:d}] in {2:f} s ({3:f} TB/s "
                "| {4:f} GFLOPS) .",
                lp[0], lp[1], sw3, tp, flops);
  PROFILE_END();

  PROFILE_START("Extrapolate");
  // Extrapolation
  spdlog::stopwatch sw4;
  kernel_calc_xtr<T, G.dsx()><<<lp[0], lp[1]>>>(
      simgrid.get_primitive_pointers(), simgrid.get_primitive_xtr_pointers(),
      simgrid.get_gradients_pointers(), simgrid.sdf_object.data(),
      simgrid.size(), dt);
  cudaDeviceSynchronize();
  tp = (2 * primitives_size_bytes + gradient_size_bytes) /
       sw4.elapsed().count() / TB;
  flops = 61 * workers / sw4.elapsed().count() / 1e9;
  max_tp = std::max(tp, max_tp);
  max_flops = std::max(flops, max_flops);
  spdlog::debug("KERNEL::cal_xtr [{0:d} x {1:d}] in {2:f} s ({3:f} TB/s | "
                "{4:f} GFLOPS) .",
                lp[0], lp[1], sw4, tp, flops);
  PROFILE_END();

  PROFILE_START("BCs Primitives XTR");
  apply_boundaries<T>(simgrid.get_primitive_xtr_pointers(),
                      simgrid.sdf_object.data(), simgrid.size(), lp);
  PROFILE_END();

  PROFILE_START("Calc Fluxes");
  // Fluxes X
  spdlog::stopwatch sw5;
  kernel_calc_fluxes<T, G.dsx()><<<lp[0], lp[1], 0, streams[0]>>>(
      simgrid.fmass_x.data(), simgrid.fmomx_x.data(), simgrid.fmomy_x.data(),
      simgrid.fmomz_x.data(), simgrid.fe_x.data(), simgrid.drho_x.data(),
      simgrid.dvx_x.data(), simgrid.dvy_x.data(), simgrid.dvz_x.data(),
      simgrid.dp_x.data(), simgrid.rho_xtr.data(), simgrid.vx_xtr.data(),
      simgrid.vy_xtr.data(), simgrid.vz_xtr.data(), simgrid.p_xtr.data(),
      simgrid.sdf_object.data(), simgrid.size(),
      std::array<std::size_t, 3>{1, 0, 0}, 0);

  kernel_calc_fluxes<T, G.dsx()><<<lp[0], lp[1], 0, streams[1]>>>(
      simgrid.fmass_y.data(), simgrid.fmomy_y.data(), simgrid.fmomx_y.data(),
      simgrid.fmomz_y.data(), simgrid.fe_y.data(), simgrid.drho_y.data(),
      simgrid.dvy_y.data(), simgrid.dvx_y.data(), simgrid.dvz_y.data(),
      simgrid.dp_y.data(), simgrid.rho_xtr.data(), simgrid.vy_xtr.data(),
      simgrid.vx_xtr.data(), simgrid.vz_xtr.data(), simgrid.p_xtr.data(),
      simgrid.sdf_object.data(), simgrid.size(),
      std::array<std::size_t, 3>{0, 1, 0}, 1);

  kernel_calc_fluxes<T, G.dsx()><<<lp[0], lp[1], 0, streams[2]>>>(
      simgrid.fmass_z.data(), simgrid.fmomz_z.data(), simgrid.fmomx_z.data(),
      simgrid.fmomy_z.data(), simgrid.fe_z.data(), simgrid.drho_z.data(),
      simgrid.dvz_z.data(), simgrid.dvx_z.data(), simgrid.dvy_z.data(),
      simgrid.dp_z.data(), simgrid.rho_xtr.data(), simgrid.vz_xtr.data(),
      simgrid.vx_xtr.data(), simgrid.vy_xtr.data(), simgrid.p_xtr.data(),
      simgrid.sdf_object.data(), simgrid.size(),
      std::array<std::size_t, 3>{0, 0, 1}, 2);

  cudaStreamSynchronize(streams[0]);
  cudaStreamSynchronize(streams[1]);
  cudaStreamSynchronize(streams[2]);
  cudaStreamDestroy(streams[0]);
  cudaStreamDestroy(streams[1]);
  cudaStreamDestroy(streams[2]);
  tp = (3 * 3 * primitives_size_bytes) / sw5.elapsed().count() / TB;
  flops = 3 * 83 * workers / sw5.elapsed().count() / 1e9;
  max_tp = std::max(tp, max_tp);
  max_flops = std::max(flops, max_flops);
  spdlog::debug("KERNELS::calc_fluxes [{0:d} x {1:d}] in {2:f} s ({3:f} TB/s | "
                "{4:f} GFLOPS) .",
                lp[0], lp[1], sw5, tp, flops);
  PROFILE_END();

  PROFILE_START("Add Fluxes ");
  // Add Fluxes
  spdlog::stopwatch sw6;
  kernel_calc_addfluxes<T, G.dsx()><<<lp[0], lp[1]>>>(
      simgrid.get_fluxes_pointers(), simgrid.get_conserved_pointers(),
      simgrid.sdf_object.data(), simgrid.size(), dt);
  cudaDeviceSynchronize();
  tp = (4 * primitives_size_bytes) / sw6.elapsed().count() / TB;
  flops = 45 * workers / sw6.elapsed().count() / 1e9;
  max_tp = std::max(tp, max_tp);
  max_flops = std::max(flops, max_flops);
  spdlog::debug("KERNELS::add_fluxes [{0:d} x {1:d}] in {2:f} s ({3:f} TB/s | "
                "{4:f} GFLOPS) .",
                lp[0], lp[1], sw6, tp, flops);
  PROFILE_END();
  PROFILE_END();
  return {dt, max_tp, max_flops};
}

template <typename T, typename F, GridInfo<T> G>
void compute(Grid<T, G, BACKEND::DEVICE> &&simgrid, T total_time,
             std::size_t max_steps, T tout, F &&init_function) {

  auto ZeroPadNumber = [](int num) -> std::string {
    std::ostringstream ss;
    ss << std::setw(7) << std::setfill('0') << num;
    return ss.str();
  };

  constexpr float GB = 1024ul * 1024ul * 1024ul;
  constexpr std::size_t mem_needed_bytes =
      (45ul) * sizeof(T) * EULERCFD::CONSTS::NX * EULERCFD::CONSTS::NY *
      EULERCFD::CONSTS::NZ;
  spdlog::warn("MEM::Memory needed for this run = {0:f} GB.",
               mem_needed_bytes / GB);

  T time = 0;
  T wtime = 0;
  std::size_t wstep = 0;
  std::size_t tstep = 0;
  T dt = 0.5;
  init_function(simgrid);

  constexpr auto lp = launch_params(G.size());
  apply_boundaries<T>(simgrid.get_primitive_pointers(),
                      simgrid.sdf_object.data(), simgrid.size(), lp);
  spdlog::debug("Launching calc_conserved kernel [{0:d} x {1:d}] ", lp[0],
                lp[1]);
  kernel_calc_conserved<T, G.dv()><<<lp[0], lp[1]>>>(
      simgrid.get_primitive_pointers(), simgrid.get_conserved_pointers(),
      simgrid.sdf_object.data(), simgrid.size());
  cudaDeviceSynchronize();

  float max_tp = 0.0;
  float max_flops = 0.0;
  while (time < total_time && tstep < max_steps) {

    spdlog::stopwatch sw;
    auto retval = compute_step<T, G>(simgrid, dt);
    dt = retval[0];
    max_tp = retval[1];
    max_flops = retval[2];
    printf("Cd=%f\n",simgrid.scratch_space[0]);
    spdlog::info("Time,  tstep, dt = [{0:f},{1:d},{2:f}] in {3:f}seconds | "
                 "{4:f} TB/s | {5:f} GFLOPS]",
                 time, tstep, dt, sw, max_tp, max_flops);
    if (wtime > tout || tstep == 0) {
      spdlog::info("IO::Storing state file.");
      std::string fname = "state." + ZeroPadNumber(wstep) + ".h5";
      simgrid.store(fname.c_str());
      wstep++;
      wtime = 0.0;
    }
    tstep++;
    time += dt;
    wtime += dt;
  }
}
} // namespace EULERCFD
