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

#include <pdal/PipelineManager.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include "Field.hpp"
#include "Grid.hpp"
#include "Raster.hpp"

namespace AtlasProcessor
{

// Least-squares plane z ~= mz + bx*(x - mx) + by*(y - my) fitted to a tile.
struct PlaneFit
{
    double mx, my, mz;
    double bx, by;

    double residual(double x, double y, double z) const
    { return z - (mz + bx * (x - mx) + by * (y - my)); }
};

struct ShapeRecord
{
    std::vector<GridIndex> indices;
    Point origin;
    bool isBefore;
    Coord tile;
    size_t id;
    bool matched;
    size_t matchId;  // SIZE_MAX if unmatched
};

// Deferred NCC work captured during the (serial, seeding-dependent) shape-match
// spiral so the pure-compute refinement can run in parallel afterward. Holds
// the flier-removed views (PDAL views are cheap PointId index vectors over a
// shared table, so capturing them does not copy point data), the fitted planes,
// the core-tile bounds, and the converged shape-match offset.
struct NccWork
{
    Coord coord;
    pdal::PointViewPtr before;
    pdal::PointViewPtr after;
    PlaneFit bPlane;
    PlaneFit aPlane;
    pdal::BOX2D box;
    Point offset;
};

class Atlas
{
public:
    Atlas();
    void run(const pdal::StringList& s);

private:
    Histogram histogram(pdal::PointViewPtr v, pdal::Dimension::Id dim,
        pdal::PointViewPtr debugView, const std::string& label);
    bool removeFliers(pdal::PointViewPtr& v, const Histogram& hist);
    PlaneFit fitPlane(pdal::PointViewPtr v);
    pdal::PointViewPtr surfaceSlice(pdal::PointViewPtr v, const PlaneFit& plane);
    std::vector<pdal::PointId> rankByResidual(pdal::PointViewPtr v,
        const PlaneFit& plane);
    pdal::PointViewPtr sliceFromRank(pdal::PointViewPtr v,
        const std::vector<pdal::PointId>& rank, double rankLo, double rankHi);
    Raster buildRaster(pdal::PointViewPtr v, Point origin, int width,
        int height, const PlaneFit& plane, int kernel);
    std::tuple<Point, double, bool> nccOffset(const Raster& br,
        const Raster& ar, int searchRadius);
    GridPtr buildGrid(pdal::PointViewPtr v, Point origin);
    void sortShapes(GridPtr& g);
    void dumpShapes(GridPtr& g);
    void dumpSurrounding();
    Coord splitterCoord(const Coord& c) const;
    std::vector<ShapePair> matchShapes(GridPtr& bg, GridPtr& ag, double threshold,
        const PlaneFit& bPlane, const PlaneFit& aPlane);
    void collectPairDisplacements(GridPtr& bg, GridPtr& ag,
        const std::vector<ShapePair>& shapes,
        std::vector<double>& pairX, std::vector<double>& pairY,
        std::vector<double>& pairZ, std::vector<Point>& centers,
        size_t dedupCount);
    std::tuple<Point, Point, double, Point> aggregateOffset(
        const std::vector<double>& pairX, const std::vector<double>& pairY,
        const std::vector<double>& pairZ);
    void addArgs();
    void load();
    void parse(const pdal::StringList& s);
    bool process(Coord c, Point& offset, double spread);
    void processGrid();
    void nccPass();
    void throwError(const std::string& s);
    void writeSimple();
    void writeTiff(const std::string& filename);
    void writeSvg(const std::string& filename, const pdal::BOX2D& extent);
    void read(const std::string& filename);
    bool shouldDump() const;
    void recordShapes(GridPtr& g, bool isBefore, Coord tile, Point origin,
        const std::vector<ShapePair>& matches);
    void writeShapeGeoJSON(const std::string& filename);

    pdal::ProgramArgs m_args;
    std::string m_beforeFilename;
    std::string m_afterFilename;

    // Base shift
    Point m_shift;

    // Splitter Origin
    Point m_origin;

    double m_xShift;
    double m_yShift;
    std::unique_ptr<Grid> m_grid;
    Field m_field;
    Coord m_coord;
    Coord m_dumpij;
    Point m_dumpxy;
    double m_dumpfrac;
    int m_minShape;
    int m_passes;
    int m_slices;
    double m_gridLen;
    bool m_ncc;
    int m_nccRadius;
    double m_nccLen;
    double m_zGate;
    bool m_zGatePlane;
    double m_matchMin;
    double m_matchMax;

    pdal::PipelineManager m_beforeMgr;
    pdal::PipelineManager m_afterMgr;
    const double m_len = 100.0;
    const double m_overlap = 20.0;
    std::vector<ShapeRecord> m_shapeRecords;
    std::vector<NccWork> m_nccWork;
    std::string m_tiffDir;
    std::string m_geojsonDir;
    std::string m_svgDir;
};

} // namespace
