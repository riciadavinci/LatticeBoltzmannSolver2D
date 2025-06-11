#include "BoundingBox.h"

namespace lbm
{

BoundingBox::BoundingBox(const uint left_bound, const uint right_bound, const uint top_bound, const uint bottom_bound)
    : LEFT_ {left_bound},
      RIGHT_ {right_bound},
      TOP_ {top_bound},
      BOTTOM_ {bottom_bound}
{}

uint BoundingBox::left() const {return LEFT_;}

uint BoundingBox::right() const {return RIGHT_;}

uint BoundingBox::top() const {return TOP_;}

uint BoundingBox::bottom() const {return BOTTOM_;}

BoundingBox BoundingBox::from_circle(uint radius, uint spherex, uint spherey) {
    return BoundingBox(spherex - radius, spherex + radius, spherey - radius, spherey + radius);
}

BoundingBox BoundingBox::from_flags_grid(const Grid<uint, 1>& grid) {
    // TODO: Implement the function
    uint upper_bound = 0, lower_bound = 1;
    uint left_bound = 0, right_bound = 1;
    const uint xbegin = 1;
    const uint ybegin = 1;
    const uint xsize = grid.xsize() - 1;
    const uint ysize = grid.ysize() - 1;
    for (uint j = ybegin; j < ysize; ++j) {
        for (uint i = xbegin; i < xsize; ++i) {
            if (grid(i, j, 0) == NO_SLIP) {
                if (upper_bound == 0) {upper_bound = j;}
                if ((left_bound == 0) || (i < left_bound)) {left_bound = i;}
                if (j > lower_bound) {lower_bound = j;}
                if (i > right_bound) {right_bound = i;}
            }
        }
    }
    return BoundingBox(left_bound, right_bound, upper_bound, lower_bound);
}

}

