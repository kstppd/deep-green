/* File:   neural_sdf.cu
* Authors: Kostis Papadakis and Adam Kit (2024)

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
USA.
* */
#include <iostream>
#include "genericTsPool.h"
#include "matrix.h"
#include "tinyAI.h"
#include "npy.hpp"
#include "spdlog/spdlog.h"

NumericMatrix::HostMatrix<float> parse_point_cloud(const char* filename) {

  npy::npy_data<float> data = npy::read_npy<float>(filename);
  const auto dims = data.shape;
  const std::size_t nrows = dims[0];
  const std::size_t ncols = dims[1];
  spdlog::info("Dims = {0:d},{1:d}", nrows, ncols);
  NumericMatrix::HostMatrix<float> x(nrows, ncols);
  for (std::size_t row = 0; row < x.nrows(); ++row) {
    for (std::size_t col = 0; col < x.ncols(); ++col) {
      auto val = data.data.at(row * x.ncols() + col);
      if (!std::isfinite(val)){
        // throw std::runtime_error("Not finite value detected");
        val=0.0;
      }
      x.set_value(row, col, val);
    }
  }
  return x;
}

void train_point_cloud(const NumericMatrix::HostMatrix<float>& data){
  
  NumericMatrix::HostMatrix<float> x(data.nrows(),3);
  NumericMatrix::HostMatrix<float> y(data.nrows(),1);
  // [:,0:3]-> spatial coords
  // [:,-1]-> spatial distance
  for (std::size_t row = 0; row < data.nrows(); ++row){
    x(row, 0)= data(row * data.ncols() + 0);
    x(row, 1)= data(row * data.ncols() + 1);
    x(row, 2)= data(row * data.ncols() + 2);
    y(row, 0)= data(row * data.ncols() + 3);
  }

  constexpr std::size_t bytes=4ull*1024ull*1024ull*1024ull;
  constexpr BACKEND HW = BACKEND::DEVICE;
  void*mem = nullptr;
  cudaMallocManaged(&mem,bytes);
  GENERIC_TS_POOL::MemPool p (mem,bytes);
  NumericMatrix::Matrix<float, HW> x_train(x.nrows(), x.ncols(), &p);
  NumericMatrix::Matrix<float, HW> y_train(y.nrows(), y.ncols(), &p);
  NumericMatrix::get_from_host(x_train, x);
  NumericMatrix::get_from_host(y_train, y);
   
  std::vector<int> arch{100,100,100,1};
  TINYAI::NeuralNetwork<float,HW,ACTIVATION::RELU,ACTIVATION::NONE,LOSSF::MSE> nn(arch, &p, x_train, y_train, 32);
  for (std::size_t i = 0; i < 10; i++) {
    auto error = nn.train(32, 1e-3);
    spdlog::info("-->Epoch [{0:d}] loss,patience=[{1:f}]", i, error);
  }

  //Build 3D mesh geometry and evaluate

}


int main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "ERROR: wrong usage!\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\t ./%s <sdf_file>\n", argv[0]);
  }

  const char *filename = argv[1];
  const NumericMatrix::HostMatrix<float> data = parse_point_cloud(filename);
  train_point_cloud(data);
  
  return 0;
}
