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

struct Raster
{
    Raster() : width(0), height(0) {}
    Raster(int width, int height)
        : width(width), height(height),
          z(width * height, std::numeric_limits<float>::quiet_NaN())
    {}

    float& at(int x, int y)
        { return z[y * width + x]; }
    float at(int x, int y) const
        { return z[y * width + x]; }
    bool valid(int x, int y) const
    {
        return x >= 0 && x < width && y >= 0 && y < height
            && !std::isnan(z[y * width + x]);
    }

    int width;
    int height;
    std::vector<float> z;
};

} // namespace AtlasProcessor
