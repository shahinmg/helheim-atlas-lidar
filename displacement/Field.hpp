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

#include "Types.hpp"

namespace AtlasProcessor
{

class Field
{
public:
    Field(size_t width, size_t height) : m_width(width), m_height(height),
        m_x(width * height, -9999.f), m_y(width * height, -9999.f), m_z(width * height, -9999.f),
        m_medianX(width * height, -9999.f), m_medianY(width * height, -9999.f),
        m_medianZ(width * height, -9999.f),
        m_madX(width * height, -9999.f), m_madY(width * height, -9999.f),
        m_madZ(width * height, -9999.f),
        m_before(width * height, -9999.f), m_after(width * height, -9999.f),
        m_valid(width * height),
        m_matchCount(width * height, -9999.f), m_rmsResidual(width * height, -9999.f),
        m_rejected(width * height, -9999.f),
        m_maxLen2(std::numeric_limits<double>::lowest())
    {}

    const float *xdata() const
    { return m_x.data(); }

    const float *ydata() const
    { return m_y.data(); }

    const float *zdata() const
    { return m_z.data(); }

    const float *medianXdata() const
    { return m_medianX.data(); }

    const float *medianYdata() const
    { return m_medianY.data(); }

    const float *medianZdata() const
    { return m_medianZ.data(); }

    const float *madXdata() const
    { return m_madX.data(); }

    const float *madYdata() const
    { return m_madY.data(); }

    const float *madZdata() const
    { return m_madZ.data(); }

    const float *bdata() const
    { return m_before.data(); }

    const float *adata() const
    { return m_after.data(); }

    const float *matchCountData() const
    { return m_matchCount.data(); }

    const float *rmsResidualData() const
    { return m_rmsResidual.data(); }

    const float *rejectedData() const
    { return m_rejected.data(); }

    size_t width() const
    { return m_width; }

    size_t height() const
    { return m_height; }

    double maxLen2() const
    { return m_maxLen2; }

    void setOffset(Coord c, Point displacement)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        double len2 = displacement.x * displacement.x + displacement.y * displacement.y;
        m_maxLen2 = std::max(len2, m_maxLen2);
        m_x[idx] = displacement.x;
        m_y[idx] = displacement.y;
        m_z[idx] = displacement.z;
        m_valid[idx] = true;
    }

    void setMatchQuality(Coord c, int count, float rms)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        m_matchCount[idx] = static_cast<float>(count);
        m_rmsResidual[idx] = rms;
    }

    void setRejected(Coord c, int count)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        m_rejected[idx] = static_cast<float>(count);
    }

    void setMedianOffset(Coord c, Point displacement)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        m_medianX[idx] = displacement.x;
        m_medianY[idx] = displacement.y;
        m_medianZ[idx] = displacement.z;
    }

    void setMadOffset(Coord c, Point mad)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        m_madX[idx] = mad.x;
        m_madY[idx] = mad.y;
        m_madZ[idx] = mad.z;
    }

    void setBeforeCount(Coord c, float count)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        m_before[idx] = count;
    }

    void setAfterCount(Coord c, float count)
    {
        int idx = pos(c);
        if (idx < 0)
            return;

        m_after[idx] = count;
    }

    Point offset(Coord c)
    {
        int idx = pos(c);
        if (idx < 0)
            return Point();
        return Point(m_x[idx], m_y[idx], m_z[idx]);
    }

    bool valid(Coord c)
    { return pos(c) >= 0; }

    // Returns {valid, estimated_offset, rms_spread_of_neighbors}.
    // Spread is the weighted RMS deviation of neighbor offsets from the mean;
    // used to set an adaptive match threshold in matchShapes.
    // With neighborsOnly, the cell's own value is ignored and the estimate
    // always comes from the neighbor average (used by the second grid pass).
    std::tuple<bool, Point, double> initialOffset(Coord c, bool neighborsOnly = false)
    {
        const double Sqrt2Recip = 0.70710678118;
        double x = 0;
        double y = 0;
        double sum = 0;

        auto addit = [&x, &y, &sum, this](Coord c, double rel)
        {
            int idx = pos(c);
            if (idx >= 0 && m_valid[idx])
            {
                x += (m_x[idx] * rel);
                y += (m_y[idx] * rel);
                sum += rel;
            }
        };

        int idx = pos(c);
        if (idx >= 0)
        {
            if (m_valid[idx] && !neighborsOnly)
                return { true, { m_x[idx], m_y[idx] }, 0.0 };

            // We weight each full neighbor equally and each corner neighbor
            // by 1/sqrt(2)
            addit(Coord(c.first + 1, c.second), 1);
            addit(Coord(c.first - 1, c.second), 1);
            addit(Coord(c.first, c.second + 1), 1);
            addit(Coord(c.first, c.second - 1), 1);
            addit(Coord(c.first + 1, c.second + 1), Sqrt2Recip);
            addit(Coord(c.first - 1, c.second + 1), Sqrt2Recip);
            addit(Coord(c.first + 1, c.second - 1), Sqrt2Recip);
            addit(Coord(c.first - 1, c.second - 1), Sqrt2Recip);
            if (sum)
            {
                x /= sum;
                y /= sum;

                // Weighted RMS deviation of neighbors from the computed mean.
                double variance = 0;
                auto addVariance = [&](Coord nc, double rel)
                {
                    int nidx = pos(nc);
                    if (nidx >= 0 && m_valid[nidx])
                    {
                        double dx = m_x[nidx] - x;
                        double dy = m_y[nidx] - y;
                        variance += rel * (dx * dx + dy * dy);
                    }
                };
                addVariance(Coord(c.first + 1, c.second), 1);
                addVariance(Coord(c.first - 1, c.second), 1);
                addVariance(Coord(c.first, c.second + 1), 1);
                addVariance(Coord(c.first, c.second - 1), 1);
                addVariance(Coord(c.first + 1, c.second + 1), Sqrt2Recip);
                addVariance(Coord(c.first - 1, c.second + 1), Sqrt2Recip);
                addVariance(Coord(c.first + 1, c.second - 1), Sqrt2Recip);
                addVariance(Coord(c.first - 1, c.second - 1), Sqrt2Recip);
                double spread = std::sqrt(variance / sum);

                return { true, { x, y }, spread };
            }
        }
        return { false, { 0, 0 }, 3.0 };
    }

private:
    size_t m_width;
    size_t m_height;

    // Easy output to GDAL
    std::vector<float> m_x;
    std::vector<float> m_y;
    std::vector<float> m_z;
    std::vector<float> m_medianX;
    std::vector<float> m_medianY;
    std::vector<float> m_medianZ;
    std::vector<float> m_madX;
    std::vector<float> m_madY;
    std::vector<float> m_madZ;
    std::vector<float> m_before;
    std::vector<float> m_after;
    std::vector<float> m_matchCount;
    std::vector<float> m_rmsResidual;
    std::vector<float> m_rejected;
    std::vector<bool> m_valid;
    double m_maxLen2;

    int pos(Coord c)
    {
        int x = c.first;
        int y = c.second;

        if (x < 0 || x >= m_width || y < 0 || y >= m_height)
            return -1;
        return (y * m_width) + x;
    }
};

}
