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

namespace AtlasProcessor
{

// Per-cell displacement statistics computed from matched shape pairs that
// survive outlier rejection.
struct OffsetStats
{
    Point mean;
    Point median;
    Point mad;
    double rmsResidual = 0;
    size_t used = 0;      // pairs kept after outlier rejection
    size_t rejected = 0;  // pairs dropped by outlier rejection
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

class Atlas
{
public:
    Atlas();
    void run(const pdal::StringList& s);

private:
    Histogram histogram(pdal::PointViewPtr v, pdal::Dimension::Id dim,
        pdal::PointViewPtr debugView, const std::string& label);
    bool removeFliers(pdal::PointViewPtr& v, const Histogram& hist);
    pdal::PointViewPtr surfaceSlice(pdal::PointViewPtr v);
    GridPtr buildGrid(pdal::PointViewPtr v, Point origin);
    void sortShapes(GridPtr& g);
    void dumpShapes(GridPtr& g);
    void dumpSurrounding();
    Coord splitterCoord(const Coord& c) const;
    std::vector<ShapePair> matchShapes(GridPtr& bg, GridPtr& ag, double threshold);
    OffsetStats calculateOffset(GridPtr& bg, GridPtr& ag,
        const std::vector<ShapePair>& shapes);
    void addArgs();
    void load();
    void parse(const pdal::StringList& s);
    bool process(Coord c, Point& offset, double spread);
    void processGrid();
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
    int m_minMatch;
    int m_passes;
    double m_gridLen;

    pdal::PipelineManager m_beforeMgr;
    pdal::PipelineManager m_afterMgr;
    const double m_len = 100.0;
    const double m_overlap = 20.0;
    std::vector<ShapeRecord> m_shapeRecords;
    std::string m_tiffDir;
    std::string m_geojsonDir;
    std::string m_svgDir;
};

} // namespace
