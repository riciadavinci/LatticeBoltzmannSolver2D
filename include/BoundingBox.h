#pragma once

#include <fstream>
#include "Globals.h"
#include "Grid.h"

namespace lbm
{

class BoundingBox
{
public:

    uint left() const;
    uint right() const;
    uint top() const;
    uint bottom() const;
    static BoundingBox from_circle(const uint radius, const uint spherex, const uint spherey);
    static BoundingBox from_flags_grid(const Grid<uint, 1>& grid);

private:
    const uint LEFT_;
    const uint RIGHT_;
    const uint TOP_;
    const uint BOTTOM_;
    BoundingBox(const uint left_bound, const uint right_bound, const uint top_bound, const uint bottom_bound);
};

}

