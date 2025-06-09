#include "../include/ex04_funcs.h"

namespace lbm
{

void usage() {
    std::cout << "Incorrect usage of the program!\n\n"
              << "Correct usage: "
              << "    ./lbm <PARAMS>.dat \n\n"
              << "    => <PARAMS>.dat : parameters file \n";
}


void save_to_vtk(const std::string& filename, const Grid<uint, 1>& Flags, const Grid<double, 1>& Density, const Grid<double, 2>& Velocity) {
    auto ipfstream = std::ofstream(filename);
    // Writing file header data
    ipfstream << "# vtk DataFile Version 4.0\n";
    ipfstream << "SiWiRVisFile\n";
    ipfstream << "ASCII\n";
    ipfstream << "DATASET STRUCTURED_POINTS\n";
    ipfstream << "DIMENSIONS " << Flags.xsize() << " " << Flags.ysize() << " " << Flags.cellsize() << "\n";
    ipfstream << "ORIGIN 0 0 0\n";
    ipfstream << "SPACING 1 1 1\n";
    ipfstream << "POINT_DATA " << Flags.total_pts() << "\n\n";

    // Flags:
    ipfstream << "SCALARS flags unsigned_int 1\n";
    ipfstream << "LOOKUP_TABLE default\n";
    for(uint j = 0; j < Flags.ysize(); ++j) {
        for(uint i = 0; i < Flags.xsize(); ++i) {
            ipfstream << Flags(i, j, 0) << "\n";
        }
    }
    ipfstream << "\n";

    // Density:
    ipfstream << "SCALARS density double 1\n";
    ipfstream << "LOOKUP_TABLE default\n";
    for(uint j = 0; j < Density.ysize(); ++j) {
        for(uint i = 0; i < Density.xsize(); ++i) {
            ipfstream << Density(i, j, 0) << "\n";
        }
    }
    ipfstream << "\n";

    // Velocity:
    ipfstream << "VECTORS velocity double\n";
    for(uint j = 0; j < Velocity.ysize(); ++j) {
        for(uint i = 0; i < Velocity.xsize(); ++i) {
            ipfstream << std::setprecision(10)
                      << Velocity(i, j, 0) << " " 
                      << Velocity(i, j, 0) << " "
                      << 0 << "\n";
        }
    }
    ipfstream << "\n";
    ipfstream.close();
}


void collide_step(Grid<uint, 1>& Flags, Grid<double, 1>& Density, Grid<double, 2>& Velocity, Grid<double, 9>& Q_PDF, const Grid<int, 2>& Cq, const Grid<double, 1>& Wq, double viscosity) {
    const uint C_SIZE = 9;

    //find density and velocity of all cells
    uint xsize = Density.xsize()-1;
    uint ysize = Density.ysize()-1; //even "Flags" Grid can be used here, just to get the total number of cells
    
    double relax_time = (3.0 * viscosity) + 0.5;
    double inv_relax_time = 1.0/relax_time;    
    for(uint i = 1; i < xsize; ++i)
    {
        for(uint j = 1; j < ysize; ++j)
        {
            double cell_density = Density(i, j, 0);
            double cell_velx = Velocity(i, j, 0);
            double cell_vely = Velocity(i, j, 1);
            double vel_norm_sq =  (cell_velx * cell_velx) + 
                                  (cell_vely * cell_vely);

            for(uint q = 0; q < C_SIZE; ++q)
            {
                double termB = Cq[2*q]*cell_velx + Cq[2*q + 1]*cell_vely;
                double termC = termB * termB;
                double Q_PDF_Eq = Wq[q] * (cell_density + 
                                    (3.0 * termB) +
                                    (4.5 * termC) -
                                    (1.5 * vel_norm_sq));
                Q_PDF(i, j, q) = Q_PDF(i, j, q) - inv_relax_time*(Q_PDF(i, j, q) - Q_PDF_Eq);
            }
        }
    }
}



void handle_boundary_conditions(Grid<uint, 1>& Flags, Grid<double, 1>& Density, Grid<double, 2>& Velocity, Grid<double, 9>& Q_PDF, const Grid<int, 2>& Cq, const Grid<double, 1>& Wq, BoundingBox box, const double& uin, const double& density_bc) {
    // Iterating over all internal cells:
    const uint xsize = Q_PDF.xsize()-1;
    const uint ysize = Q_PDF.ysize()-1;

    const uint TOP_ROW_IDX = 0;
    const uint TOP_FLUID = TOP_ROW_IDX + 1;
    const uint BOTTOM_ROW_IDX = Q_PDF.ysize()-1;
    const uint BOTTOM_FLUID = BOTTOM_ROW_IDX - 1;

    const uint LEFT_COL_IDX = 0;
    const uint LEFT_FLUID = LEFT_COL_IDX + 1;
    const uint RIGHT_COL_IDX = Q_PDF.xsize()-1;
    const uint RIGHT_FLUID = RIGHT_COL_IDX - 1;
    
    // ###################################
    // Velocity Reflection for No slip cells:
    for(uint i = 1; i < xsize; ++i) {
        // Top-Boundary:
        Q_PDF(i-1, TOP_ROW_IDX, SE) = Q_PDF(i, TOP_FLUID, NW);
        Q_PDF(i, TOP_ROW_IDX, S) = Q_PDF(i, TOP_FLUID, N);
        Q_PDF(i+1, TOP_ROW_IDX, SW) = Q_PDF(i, TOP_FLUID, NE);

        // Bottom-Boundary:
        Q_PDF(i-1, BOTTOM_ROW_IDX, NE) = Q_PDF(i, BOTTOM_FLUID, SW);
        Q_PDF(i, BOTTOM_ROW_IDX, N) = Q_PDF(i, BOTTOM_FLUID, S);
        Q_PDF(i+1, BOTTOM_ROW_IDX, NW) = Q_PDF(i, BOTTOM_FLUID, SE);
    }
    // ###################################


    // ###################################
    // Velocity reflection for velocity and density type cells:
    for(uint j = 1; j < ysize; ++j) {
        // Left-Boundary Handling - Velocity Boundary:
        // Update Q_PDF values on the LEFT_FLUID Cells & 
        //push the updates values in the LEFT_FLUID Cells into the velocity boundary pseudo layer  
        Q_PDF(LEFT_COL_IDX, j-1, SE) = Q_PDF(LEFT_FLUID, j, NW); 
        Q_PDF(LEFT_COL_IDX, j, E) = Q_PDF(LEFT_FLUID, j, W);
        Q_PDF(LEFT_COL_IDX, j+1, NE) = Q_PDF(LEFT_FLUID, j, SW); 

        Q_PDF(LEFT_COL_IDX, j-1, SE) -= (6.0 * Wq[NW] * ((-1.0 * uin)));
        Q_PDF(LEFT_COL_IDX, j, E) -= (6.0 * Wq[W] * ((-1.0 * uin)));
        Q_PDF(LEFT_COL_IDX, j+1, NE) -= (6.0 * Wq[SW] * ((-1.0 * uin)));

        // Right-Boundary Handling - Density Boundary:
        // Update Q_PDF values on the RIGHT_FLUID Cells
        double VelX = Velocity(RIGHT_FLUID, j, 0);
        double VelY = Velocity(RIGHT_FLUID, j, 1);
        double VelX2 = (VelX * VelX);
        double VelY2 = (VelY * VelY);
        double termA = 1.5 * (VelX2 + VelY2);
        Q_PDF(RIGHT_COL_IDX, j-1, SW) = Q_PDF(RIGHT_FLUID, j, NE);
        Q_PDF(RIGHT_COL_IDX, j, W) = Q_PDF(RIGHT_FLUID, j, E);
        Q_PDF(RIGHT_COL_IDX, j+1, NW) = Q_PDF(RIGHT_FLUID, j, SE);

        Q_PDF(RIGHT_COL_IDX, j-1, SW) = (-1.0* Q_PDF(RIGHT_COL_IDX, j-1, SW)) + (2*Wq[NE]*(density_bc + 
                                                                      (4.5*std::pow((Cq(2,0,0)*VelX + Cq(2,0,1)*VelY), 2)) - 
                                                                      termA));
        Q_PDF(RIGHT_COL_IDX, j, W) = (-1.0*Q_PDF(RIGHT_COL_IDX, j, W)) + (2*Wq[E]*(density_bc + 
                                                                      (4.5*std::pow((Cq(2,1,0)*VelX + Cq(2,1,1)*VelY), 2)) - 
                                                                      termA));
        Q_PDF(RIGHT_COL_IDX, j+1, NW) = (-1.0*Q_PDF(RIGHT_COL_IDX, j+1, NW)) + (2*Wq[SE]*(density_bc + 
                                                                      (4.5*std::pow((Cq(2,2,0)*VelX + Cq(2,2,1)*VelY), 2)) - 
                                                                      termA));
    }
    // ###################################

    // ###################################
    // Handling for obstacle boundary (only circle for now):
    const uint LEFT_BOUNDARY = box.left();
    const uint RIGHT_BOUNDARY = box.right();
    const uint TOP_BOUNDARY = box.top();
    const uint BOTTOM_BOUNDARY = box.bottom();

    for(uint i = LEFT_BOUNDARY; i <= RIGHT_BOUNDARY; ++i) {
        for(uint j = TOP_BOUNDARY; j <= BOTTOM_BOUNDARY; ++j) {
            if(Flags(i, j, 0) == NO_SLIP) {
                // Top-left:
                if(Flags(i-1, j-1, 0) == FLUID_CELL) {
                    Q_PDF(i, j, NW) = Q_PDF(i-1, j-1, SE);
                }

                // Top:
                if(Flags(i, j-1, 0) == FLUID_CELL) {
                    Q_PDF(i, j, N) = Q_PDF(i-1, j, S);
                }

                // Top-right:
                if(Flags(i+1, j-1, 0) == FLUID_CELL) {
                    Q_PDF(i, j, NE) = Q_PDF(i+1, j-1, SW);
                }

                // Left:
                if(Flags(i-1, j, 0) == FLUID_CELL) {
                    Q_PDF(i, j, W) = Q_PDF(i-1, j, E);
                }

                // Right:
                if(Flags(i+1, j, 0) == FLUID_CELL) {
                    Q_PDF(i, j, E) = Q_PDF(i+1, j, W);
                }

                // Bottom-left:
                if(Flags(i-1, j+1, 0) == FLUID_CELL) {
                    Q_PDF(i, j, SW) = Q_PDF(i-1, j+1, NE);
                }

                // Bottom:
                if(Flags(i, j+1, 0) == FLUID_CELL) {
                    Q_PDF(i, j, S) = Q_PDF(i, j+1, N);
                }

                // Bottom-right:
                if(Flags(i+1, j+1, 0) == FLUID_CELL) {
                    Q_PDF(i, j, SE) = Q_PDF(i+1, j+1, NW);
                }
            }
        }
    }
    // ###################################
}


void streaming_step(Grid<uint, 1>& Flags, Grid<double, 9>& Q_PDF_src, Grid<double, 9>& Q_PDF_dst) {
    const uint X_BEGIN = 1;
    const uint X_END = Q_PDF_src.xsize()-1;
    const uint Y_BEGIN = 1;
    const uint Y_END = Q_PDF_src.ysize()-1;

    for(uint i = X_BEGIN; i < X_END; ++i) {
        for(uint j = Y_BEGIN; j < Y_END; ++j) {
            if(Flags(i, j, 0) == FLUID_CELL) {
                // Top row:
                Q_PDF_dst(i, j, NW) = Q_PDF_src(i+1, j+1, NW);
                Q_PDF_dst(i, j, N) = Q_PDF_src(i, j+1, N);
                Q_PDF_dst(i, j, NE) = Q_PDF_src(i-1, j+1, NE);

                // Middle row:
                Q_PDF_dst(i, j, W) = Q_PDF_src(i+1, j, W);
                Q_PDF_dst(i, j, C) = Q_PDF_src(i, j, C);
                Q_PDF_dst(i, j, E) = Q_PDF_src(i-1, j, E);

                // Bottom row:
                Q_PDF_dst(i, j, SW) = Q_PDF_src(i+1, j-1, SW);
                Q_PDF_dst(i, j, S) = Q_PDF_src(i, j-1, S);
                Q_PDF_dst(i, j, SE) = Q_PDF_src(i-1, j-1, SE);
            }
        }
    }

}

void calc_velocity_and_density(Grid<double, 1>& Density, Grid<double, 2>& Velocity, Grid<double, 9>& Q_PDF, const Grid<int, 2>& Cq)
{
    const uint C_SIZE = 9;
    uint xsize = Density.xsize()-1;
    uint ysize = Density.ysize()-1;
    
    for(uint i = 1; i < xsize; ++i)
    {
        for(uint j = 1; j < ysize; ++j)
        {
            //sum up all QPDF values inside the cell
            double d_sum = 0.0;
            double vx_sum = 0.0;
            double vy_sum = 0.0;
            for(uint q = 0; q < C_SIZE; ++q)
            {
                d_sum += Q_PDF(i, j, q);
                vx_sum += (Q_PDF(i, j, q) * Cq[2*q]);
                vy_sum += (Q_PDF(i, j, q) * Cq[2*q + 1]);
            }
            if(d_sum == 0) {
                std::cout<<"density becoming zero!!\n"; 
                std::cout<<"i : "<<i<<"\n";
                std::cout<<"j : "<<j<<"\n";
                return;}
            Density(i, j, 0) = d_sum;
            Velocity(i, j, 0) = vx_sum;
            Velocity(i, j, 1) = vy_sum;
        }
    }
}


//below function only for debug purpose
void initialize_bc_with_nonsense(Grid<double, 9>& Q_PDF)
{
    srand(time(NULL));

    //std::cout <<  rand()%100 << "\n";
    const uint TOP_ROW_IDX = 0;
    const uint BOTTOM_ROW_IDX = Q_PDF.ysize()-1;

    const uint LEFT_COL_IDX = 0;
    const uint RIGHT_COL_IDX = Q_PDF.xsize()-1;
    
    //top and bottom no-slip bc
    for(uint i = 0; i <= RIGHT_COL_IDX; ++i)
    {
        Q_PDF(i, TOP_ROW_IDX, SW) = 10.0;
        Q_PDF(i, TOP_ROW_IDX, S) = 20.0;
        Q_PDF(i, TOP_ROW_IDX, SE) = 30.0;

        Q_PDF(i, BOTTOM_ROW_IDX, NW) = 40.0;
        Q_PDF(i, BOTTOM_ROW_IDX, N) = 50.0;
        Q_PDF(i, BOTTOM_ROW_IDX, NE) = 60.0;
    }
    
    //left and right bc
    for(uint j = 1; j < BOTTOM_ROW_IDX; ++j)
    {
        Q_PDF(LEFT_COL_IDX, j, NE) = 1.0;
        Q_PDF(LEFT_COL_IDX, j, E) = 1.0;
        Q_PDF(LEFT_COL_IDX, j, SE) = 1.0;

        Q_PDF(RIGHT_COL_IDX, j, NW) = 1.0;
        Q_PDF(RIGHT_COL_IDX, j, W) = 1.0;
        Q_PDF(RIGHT_COL_IDX, j, SW) = 1.0;
    }
}


void init_qpdf_with_wq(Grid<double, 9>& Q_PDF, const Grid<double, 1>& Wq) {
    const uint xsize = Q_PDF.xsize();
    const uint ysize = Q_PDF.ysize();
    const uint cellsize = Q_PDF.cellsize();

    for (uint i = 0; i < xsize; ++i) {
        for (uint j = 0; j < ysize; ++j) {
            for(uint q = 0; q < cellsize; ++q)
            {
                Q_PDF(i, j, q) = Wq[q];
            }
        }
    }
}


void apply_bcs_to_flags(Grid<uint, 1>& Flags) {
    // Apply top and bottom no-slip bcs/cell types:
    const uint NUM_COLS = Flags.xsize();
    const uint FIRST_COL = 0;
    const uint LAST_COL = NUM_COLS - 1;

    const uint NUM_ROWS = Flags.ysize();
    const uint FIRST_ROW = 0;
    const uint SECOND_ROW = FIRST_ROW + 1;
    const uint LAST_ROW = NUM_ROWS - 1;
    
    for(uint i = FIRST_COL; i < NUM_COLS; ++i) {
        Flags(i, FIRST_ROW, 0) = static_cast<uint>(NO_SLIP);
        Flags(i, LAST_ROW, 0) = static_cast<uint>(NO_SLIP);
    }

    // Apply left and right velocity and densiry bcs/cell types:
    for (uint j = SECOND_ROW; j < LAST_ROW; ++j) {
        Flags(FIRST_COL, j, 0) = static_cast<uint>(VEL_BCS);
        Flags(LAST_COL, j, 0) = static_cast<uint>(DENS_BCS);
    }
}

template<typename Type>
double get_distance(const Type x1, const Type y1, const Type x2, const Type y2) {
    return static_cast<double>(std::pow(std::abs(x1-x2), 2)) +  \
           static_cast<double>(std::pow(std::abs(y1-y2), 2));
}

void apply_circle_obstacle_to_flags(Grid<uint, 1>& Flags, const uint spherex, const uint spherey, const uint radius) {
    // -1 coz we are iterating over internal fluid cells
    const uint LEFT_BOUND = spherex - radius;
    const uint RIGHT_BOUND = spherex + radius;
    const uint UPPER_BOUND = spherey - radius;
    const uint LOWER_BOUND = spherey + radius;

    const double r2 = static_cast<double>(radius * radius);

    for (uint i = LEFT_BOUND; i <= RIGHT_BOUND; ++i) {
        for (uint j = UPPER_BOUND; j <= LOWER_BOUND; ++j) {
            if (get_distance(static_cast<int>(i), static_cast<int>(j), static_cast<int>(spherex), static_cast<int>(spherey)) <= r2) {
                Flags(i, j, 0) = NO_SLIP;
            }
        }
    }
}

std::pair<BoundingBox, Grid<uint, 1>> get_bounding_box_and_flags_grid(FileReader fr) {
    auto geometry = fr.getParam<std::string>("geometry");
    if (geometry != "") {
        // Arbitrary geometry!
        const std::string ip_geom_file = "files/" + geometry;        
        auto Flags = Grid<uint, 1>::from_pgm_file(ip_geom_file);
        auto box = BoundingBox::from_flags_grid(Flags);
        return std::make_pair(box, Flags);
    }
    // Circle obstacle:
    auto sizex = fr.getParam<uint>("sizex")+2;
    auto sizey = fr.getParam<uint>("sizey")+2;
    auto spherex = fr.getParam<uint>("spherex");
    auto spherey = fr.getParam<uint>("spherey");
    auto radius = fr.getParam<uint>("diameter")/2;
    auto box = BoundingBox::from_circle(radius, spherex, spherey);
    auto Flags = Grid<uint, 1>(sizex, sizey);
    apply_circle_obstacle_to_flags(Flags, spherex, spherey, radius);
    return std::make_pair(box, Flags);
}

}
