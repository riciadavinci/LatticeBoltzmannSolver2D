#include "Grid.h"

namespace lbm
{

template<typename Type, uint CellSize>
Grid<Type, CellSize>::Grid(uint xsize, uint ysize) 
    : xsize_ {xsize},
      ysize_ {ysize},
      data_ {std::vector<Type>(xsize * ysize * CellSize)}
{}


template<typename Type, uint CellSize>
Grid<Type, CellSize>::Grid(uint xsize, uint ysize, std::initializer_list<Type> l) 
    : xsize_ {xsize},
      ysize_ {ysize},
      data_ {l}
{}


template<typename Type, uint CellSize>
inline Type Grid<Type, CellSize>::operator()(uint x, uint y, uint pos) const {
    assert((x < xsize_) && (y < ysize_) && (pos < CellSize));
    return data_[pos + x*CellSize + y*CellSize*xsize_];
}


template<typename Type, uint CellSize>
inline Type& Grid<Type, CellSize>::operator()(uint x, uint y, uint pos) {
    assert((x < xsize_) && (y < ysize_) && (pos < CellSize));
    return data_[pos + x*CellSize + y*CellSize*xsize_];
}


template<typename Type, uint CellSize>
inline uint Grid<Type, CellSize>::xsize() const {
    return xsize_;
}


template<typename Type, uint CellSize>
inline uint Grid<Type, CellSize>::ysize() const {
    return ysize_;
}


template<typename Type, uint CellSize>
inline uint Grid<Type, CellSize>::cellsize() const {
    return CellSize;
}



template<typename Type, uint CellSize>
void Grid<Type, CellSize>::swap(Grid& other) {
    if(&other != this) {
        std::swap(xsize_, other.xsize_);
        std::swap(ysize_, other.ysize_);
        std::swap(data_, other.data_);
    }
}

template<typename Type, uint CellSize>
inline uint Grid<Type, CellSize>::total_pts() const {
    return xsize_ * ysize_ * CellSize;
}


template<typename Type, uint CellSize>
void Grid<Type, CellSize>::fill(Type initVal) {
    std::fill(data_.begin(), data_.end(), initVal);
}


template<typename Type, uint CellSize>
Type Grid<Type, CellSize>::operator[](size_t idx) const {
    assert(idx < this->total_pts());
    return data_[idx];
}

template<typename Type, uint CellSize>
Type& Grid<Type, CellSize>::operator[](size_t idx) {
    assert(idx < this->total_pts());
    return data_[idx];
}


template<typename Type, uint CellSize>
void swap(Grid<Type, CellSize>& grid1, Grid<Type, CellSize>& grid2) {
    grid1.swap(grid2);
}


template<typename Type, uint CellSize>
Type Grid<Type, CellSize>::get_max_value() const {
    auto max_elem = std::max_element(data_.begin(), data_.end());
    return (*max_elem);
}

template<typename Type, uint CellSize>
Type Grid<Type, CellSize>::get_min_value() const {
    auto min_elem = std::min_element(data_.begin(), data_.end());
    return (*min_elem);
}

template<typename Type, uint CellSize>
Grid<Type, CellSize> Grid<Type, CellSize>::from_pgm_file(const std::string& filename) {
    auto ipfstream = std::ifstream(filename);
    std::string line;
    // Skipping header lines
    const uint header_lines_to_skip = 2;
    for (uint i = 0; i < header_lines_to_skip; ++i) { std::getline(ipfstream, line); }
    // Reading the dimensions!
    std::getline(ipfstream, line);
    uint xsize, ysize;
    auto ss = std::stringstream(line);
    ss >> xsize >> ysize;
    auto grid = Grid<Type, 1>(xsize+2, ysize+2);
    // auto grid = Grid<Type, 1>(xsize, ysize);
    uint value;

    grid.fill(static_cast<uint>(FLUID_CELL));

    for (uint j = grid.ysize()-2; j > 0; --j) {
        for (uint i = 1; i < grid.xsize()-1; ++i) {
            std::getline(ipfstream, line);
            value = std::stoul(line);
            if (value == 0) {
                grid(i, j, 0) = static_cast<uint>(NO_SLIP);
            }
        }
    }
    ipfstream.close();
    return grid;
}


}


