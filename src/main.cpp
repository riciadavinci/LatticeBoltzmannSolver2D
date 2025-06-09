#include <string>
#include <iostream>
#include <iomanip>
#include <memory>


#include "../include/Globals.h"
#include "../include/Grid.h"
#include "../include/ex04_funcs.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        lbm::usage();
        exit(-1);
    }
    const auto filename = std::string(argv[1]);

    std::cout << "Parameter file: " << filename << "\n";

    // Reading file params
    FileReader fr(filename);
    std::string outfile = fr.getParam<std::string>("vtk_file");
    uint timesteps = fr.getParam<uint>("timesteps");
    uint vtk_step = fr.getParam<uint>("vtk_step");
    double uin = fr.getParam<double>("uin");
    uint Re = fr.getParam<uint>("Re");

    auto box_and_flags_grid = lbm::get_bounding_box_and_flags_grid(fr);
    auto box = box_and_flags_grid.first;
    auto Flags = box_and_flags_grid.second;
    auto sizex = Flags.xsize();
    auto sizey = Flags.ysize();

    double viscosity = (uin * static_cast<double>(sizey)) / static_cast<double>(Re);
    double density_bc = 1.0;

    lbm::apply_bcs_to_flags(Flags);
    
    // Constant matrices:
    const auto Wq = lbm::Grid<double, 1>(3, 3, {
        1.0/36.0, 1.0/9.0, 1.0/36.0,
        1.0/9.0, 4.0/9.0, 1.0/9.0,
        1.0/36.0, 1.0/9.0, 1.0/36.0
    });

    // PDF Directions:
    const auto Cq = lbm::Grid<int, 2>(3, 3, {
        -1, 1, 0, 1, 1, 1,
        -1, 0, 0, 0, 1, 0,
        -1, -1, 0, -1, 1, -1 
    });

    // Init for PDF Q:
    auto Q_PDF_src = lbm::Grid<double, 9>(sizex, sizey);
    lbm::init_qpdf_with_wq(Q_PDF_src, Wq);

    auto Q_PDF_dst = lbm::Grid<double, 9>(sizex, sizey);
    lbm::init_qpdf_with_wq(Q_PDF_dst, Wq);


    auto Density = lbm::Grid<double, 1>(sizex, sizey);
    // Init density to 1 for incompressible flow:
    Density.fill(1.0);

    auto Velocity = lbm::Grid<double, 2>(sizex, sizey);
    
    // Collide Step:
    for (size_t tstep = 0; tstep <= timesteps; ++tstep) {
        lbm::collide_step(Flags, Density, Velocity, Q_PDF_src, Cq, Wq, viscosity);

        //The below function is only for debugging purposes
        //lbm::initialize_bc_with_nonsense(Q_PDF_src);
        
        lbm::handle_boundary_conditions(Flags, Density, Velocity, Q_PDF_src, Cq, Wq, box, uin, density_bc);
        lbm::streaming_step(Flags, Q_PDF_src, Q_PDF_dst);
        lbm::calc_velocity_and_density(Density, Velocity, Q_PDF_dst, Cq);

        lbm::swap(Q_PDF_src, Q_PDF_dst);
        if (tstep%vtk_step == 0) {
            const std::string out_filename = outfile + std::to_string(tstep) + ".vtk";
            lbm::save_to_vtk(out_filename, Flags, Density, Velocity);
            std::cout << "Saved VTK file for timestep: " << tstep << "\n";
        }
    }
    std::cout << "Done!\n";
    return 0;
}
