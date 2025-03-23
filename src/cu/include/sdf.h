/* File: sdf.h
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
#ifndef dev_host
#define dev_host
#endif
#include <array>
namespace EULERCFD {

template <typename T>
constexpr dev_host T sdf_ellipsoid(T x, T y, T z, T cx, T cy, T cz, T rx, T ry,
                                   T rz, T radious) {

  T dx = x - cx;
  T dy = y - cy;
  T dz = z - cz;
  T k0 = std::sqrt((dx * dx) / (rx * rx) + (dy * dy) / (ry * ry) +
                   (dz * dz) / (rz * rz));
  T k1 = std::sqrt((dx * dx) / (rx * rx * rx * rx) +
                   (dy * dy) / (ry * ry * ry * ry) +
                   (dz * dz) / (rz * rz * rz * rz));
  return k0 * (k0 - T(1)) / k1;
}

template <typename T>
constexpr dev_host T sdf(T x, T y, T z, T cx, T cy, T cz, T radious) {
  return std::sqrt(std::pow((x - cx), T(2)) + std::pow((y - cy), T(2)) +
                   std::pow((z - cz), T(2))) -
         radious;
}

template <typename T>
constexpr dev_host T sdf(std::array<T, 3> p, std::array<T, 3> c, T radious) {
  return sdf(p[0], p[1], p[2], c[0], c[1], c[2], radious);
}
} // namespace EULERCFD
