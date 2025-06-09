#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <fstream>
#include <string>
#include "Globals.h"

namespace lbm
{

template<typename Type, uint CellSize>
class Grid
{
public:
    Grid() = default;
    Grid(uint xsize, uint ysize);
    Grid(uint xsize, uint ysize, std::initializer_list<Type> l);
    
    inline Type operator()(uint x, uint y, uint pos) const;
    inline Type& operator()(uint x, uint y, uint pos);

    inline uint xsize() const;
    inline uint ysize() const;
    inline uint cellsize() const;

    void swap(Grid& other);

    inline uint total_pts() const;
    void fill(Type initVal);

    Type operator[](size_t idx) const;
    Type& operator[](size_t idx);

    Type get_max_value() const;
    Type get_min_value() const;

    static Grid<Type, CellSize> from_pgm_file(const std::string& filename);

private:
    uint xsize_;
    uint ysize_;
    std::vector<Type> data_;
};

template<typename Type, uint CellSize>
void swap(Grid<Type, CellSize>& grid1, Grid<Type, CellSize>& grid2);

}

#include "../src/Grid.cpp"
