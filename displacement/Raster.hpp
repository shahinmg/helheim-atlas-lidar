/*****************************************************************************
 *   Copyright (c) 2024, Hobu, Inc. (info@hobu.co)                           *
 *                                                                           *
 *   All rights reserved.                                                    *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/

#pragma once

#include <cmath>
#include <limits>
#include <vector>

namespace AtlasProcessor
{

// Splatted surface raster: z is the weighted mean value per cell, w the
// total splat weight — a measure of how much real data backs the cell
// (0 = empty). Correlation code weights cells by w rather than treating
// interpolated halo cells as equal to data-bearing ones.
struct Raster
{
    Raster() : width(0), height(0) {}
    Raster(int width, int height)
        : width(width), height(height),
          z(width * height, 0.f), w(width * height, 0.f)
    {}

    float& at(int x, int y)
        { return z[y * width + x]; }
    float at(int x, int y) const
        { return z[y * width + x]; }
    float& weight(int x, int y)
        { return w[y * width + x]; }
    float weight(int x, int y) const
    {
        if (x < 0 || x >= width || y < 0 || y >= height)
            return 0.f;
        return w[y * width + x];
    }

    int width;
    int height;
    std::vector<float> z;
    std::vector<float> w;
};

} // namespace AtlasProcessor
