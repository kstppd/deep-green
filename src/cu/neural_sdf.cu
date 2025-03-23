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
#include "tinyAI.h"

void parse_point_cloud(){}
void train_point_cloud(){}
void export_sdf_mask(){}

int main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "ERROR: wrong usage!\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\t ./%s <sdf_file>\n", argv[0]);
  }

  const char *filename = argv[1];

  return 0;
}
