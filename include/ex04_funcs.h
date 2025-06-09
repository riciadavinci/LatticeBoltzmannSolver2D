#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iomanip>
#include "Globals.h"
#include "Grid.h"
#include "BoundingBox.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "FileReader.h"
#include <string>
#include <utility>

// class FileReader;

namespace lbm
{

void usage();

void save_to_vtk(const std::string& filename, const Grid<uint, 1>& Flags, const Grid<double, 1>& Density, const Grid<double, 2>& Velocity);

void collide_step(Grid<uint, 1>& Flags, Grid<double, 1>& Density, Grid<double, 2>& Velocity, Grid<double, 9>& Q_PDF, const Grid<int, 2>& Cq, const Grid<double, 1>& Wq, double viscosity);

void handle_boundary_conditions(Grid<uint, 1>& Flags, Grid<double, 1>& Density, Grid<double, 2>& Velocity, Grid<double, 9>& Q_PDF, const Grid<int, 2>& Cq, const Grid<double, 1>& Wq, BoundingBox box, const double& uin, const double& density_bc);

void streaming_step(Grid<uint, 1>& Flags, Grid<double, 9>& Q_PDF_src, Grid<double, 9>& Q_PDF_dst);

void calc_velocity_and_density(Grid<double, 1>& Density, Grid<double, 2>& Velocity, Grid<double, 9>& Q_PDF, const Grid<int, 2>& Cq);
//below function only for debug purposes
void initialize_bc_with_nonsense(Grid<double, 9>& Q_PDF);

void init_qpdf_with_wq(Grid<double, 9>& Q_PDF, const Grid<double, 1>& Wq);

void apply_bcs_to_flags(Grid<uint, 1>& Flags);

template<typename Type>
double get_distance(const Type x1, const Type y1, const Type x2, const Type y2);

void apply_circle_obstacle_to_flags(Grid<uint, 1>& Flags, const uint spherex, const uint spherey, const uint radius);

std::pair<BoundingBox, Grid<uint, 1>> get_bounding_box_and_flags_grid(FileReader fr);


}

#include "../src/ex04_funcs.cpp"
