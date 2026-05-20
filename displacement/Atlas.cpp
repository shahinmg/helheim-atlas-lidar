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

#include "Atlas.hpp"
#include "Draw.hpp"

#include <algorithm>
#include <filesystem>
#include <math.h>
#include <numeric>

#include <pdal/private/gdal/GDALUtils.hpp>
#include <pdal/private/gdal/Raster.hpp>
#include <pdal/private/MathUtils.hpp>
#include <pdal/filters/SplitterFilter.hpp>
#include <pdal/util/FileUtils.hpp>
#include <pdal/io/GDALReader.hpp>


namespace pdal
{
namespace Utils
{

template<>
std::string toString(const AtlasProcessor::Coord& c)
{
    return std::to_string(c.first) + ", " + std::to_string(c.second);
}

template<>
StatusWithReason fromString(const std::string& from, AtlasProcessor::Coord& c)
{
    std::string comma;
    auto pos = from.find(',');
    if (pos == std::string::npos || pos + 1 == from.size())
        return false;
    c.first = std::stoi(from.substr(0, pos));
    c.second = std::stoi(from.substr(pos + 1));
    return true;
}

}
}


namespace AtlasProcessor
{

Coord COORD;

void fatal(const std::string& err)
{
    std::cerr << "atlas: " << err << "\n";
    exit(-1);
}

//Origin is set to lower left corner.
Atlas::Atlas() : m_origin({525500, 7356000}), m_field(XCellCount, YCellCount)
{}

void Atlas::throwError(const std::string& s)
{
    throw std::runtime_error(s);
}

void Atlas::addArgs()
{
    m_args.add("before", "Filename of scene at time 't'",
        m_beforeFilename).setPositional();
    m_args.add("after", "Filename of scene at time 't + n'",
        m_afterFilename).setPositional();
    m_args.add("xshift", "X distance shift", m_shift.x).setPositional();
    m_args.add("yshift", "Y distance shift", m_shift.y).setPositional();
    m_args.add("dumpxy", "XY pos in UTM meters to dump", m_dumpxy, Point(-1000, -1000));
    m_args.add("dumpij", "IJ pos in grid index to dump", m_dumpij, Coord(-1000, -1000));
    m_args.add("dumpfrac", "Top Fraction to be used for dump", m_dumpfrac, .1);
    m_args.add("tiff", "Output directory for GeoTIFF", m_tiffDir, std::string());
    m_args.add("geojson", "Output directory for shapes GeoJSON (omit to skip)",
        m_geojsonDir, std::string());
    m_args.add("svg", "Output directory for vector SVG (omit to skip)",
        m_svgDir, std::string());
}

void Atlas::parse(const pdal::StringList& slist)
{
    try
    {
        m_args.parse(slist);
    }
    catch (const pdal::arg_error& err)
    {
        fatal(err.what());
    }
    if (m_dumpxy != Point(-1000, -1000) && m_dumpij != Coord(-1000, -1000))
        fatal("Can't specifify both 'dumpxy' and 'dumpij'.");
    if (m_dumpfrac < 0 || m_dumpfrac > 1)
        fatal("'dumpfrac' must be a value in the range (0,1].");
    m_dumpfrac = 1 - m_dumpfrac;
}

void Atlas::run(const pdal::StringList& s)
{
    addArgs();
    parse(s);
    try
    {
        load();
        processGrid();
        namespace fs = std::filesystem;
        std::string stem = pdal::FileUtils::stem(m_beforeFilename) + "_" +
            pdal::FileUtils::stem(m_afterFilename);

        if (!m_geojsonDir.empty())
        {
            fs::path geojsonPath = fs::path(m_geojsonDir) / (stem + ".geojson");
            writeShapeGeoJSON(geojsonPath.string());
        }

        if (!m_svgDir.empty())
        {
            pdal::SplitterFilter *splitter =
                dynamic_cast<pdal::SplitterFilter *>(m_beforeMgr.getStage());
            fs::path svgPath = fs::path(m_svgDir) / (stem + ".svg");
            writeSvg(svgPath.string(), splitter->extent());
        }

        fs::path tiffPath = m_tiffDir.empty()
            ? fs::path(stem + ".tif")
            : fs::path(m_tiffDir) / (stem + ".tif");
        writeTiff(tiffPath.string());
//        read(filename);
    }
    catch (const pdal::pdal_error& err)
    {
        fatal(err.what());
    }
}


// Splitter uses Y-up coordinates.
Coord Atlas::splitterCoord(const Coord& c) const
{
    return {c.first, YCellCount - c.second - 1};
}


bool Atlas::shouldDump() const
{
    return m_dumpij != Coord(-1000, -1000);
}


void Atlas::load()
{
    using namespace pdal;

    if (m_dumpxy != Point(-1000, -1000))
    {
        double xlen = XCellCount * m_len;
        double ylen = YCellCount * m_len;
        m_dumpij.first = (m_dumpxy.x - m_origin.x) / m_len;
        m_dumpij.second = (m_dumpxy.y - (m_origin.y - ylen)) / m_len;
    }

    Options splitterOpts;
    splitterOpts.add("length", m_len);
    splitterOpts.add("origin_x", m_origin.x);
    splitterOpts.add("origin_y", m_origin.y);
    splitterOpts.add("overlap", m_overlap);

    StageCreationOptions bOps { m_beforeFilename };
    Stage& beforeReader = m_beforeMgr.makeReader(bOps);
    Stage& beforeSplitter = m_beforeMgr.makeFilter("filters.splitter",
        beforeReader, splitterOpts);
    SplitterFilter& s = dynamic_cast<SplitterFilter&>(beforeSplitter);
    m_beforeMgr.execute(ExecMode::Standard);

    Options afterSplitterOpts;
    afterSplitterOpts.add("length", m_len);
    afterSplitterOpts.add("origin_x", m_origin.x + m_shift.x);
    afterSplitterOpts.add("origin_y", m_origin.y + m_shift.y);
    afterSplitterOpts.add("overlap", m_overlap);

    StageCreationOptions aOps { m_afterFilename };
    Stage& afterReader = m_afterMgr.makeReader(aOps);
    Stage& afterSplitter = m_afterMgr.makeFilter("filters.splitter",
        afterReader, afterSplitterOpts);
    m_afterMgr.execute(ExecMode::Standard);
}


void Atlas::processGrid()
{
    Coord start {XCellCount / 2, YCellCount / 2};

    // Spiral must reach the farthest cell from start; process() filters
    // out-of-bounds coords via m_field.valid().
    int maxRadius = std::max({start.first, start.second,
        XCellCount - start.first - 1, YCellCount - start.second - 1});

    for (int i = 0; i <= maxRadius; ++i)
    {
        int y = -i;
        int x = -i;
        while (true)
        {
            Coord c{start.first + x, start.second + y};
            m_coord = c;

            // If this cell has an initial offset, use it, otherwise use
            // the base one.
            auto [valid, pos, spread] = m_field.initialOffset(c);
            if (!valid)
            {
                pos = m_shift;
                spread = 3.0;
            }

            // The offset is updated if processing works.
            if (process(c, pos, spread))
                m_field.setOffset(c, pos);

            // This starts in the lower left corner and goes around in
            // a "circle".  Once we get back to the beginning, we break.
            if (y == -i && x != i)
                x++;
            else if (x == i && y != i)
                y++;
            else if (y == i && x != -i)
                x--;
            else if (x == -i && y != -i)
                y--;
            if (x == -i && y == -i)
                break;
        }
    }

    if (shouldDump())
        dumpSurrounding();
}

void Atlas::dumpSurrounding()
{
    int repI = m_dumpij.first + 100;
    int repJ = m_dumpij.second + 50;

    std::cout << "\t" << (m_dumpij.first - 1) << "\t" << (m_dumpij.first) << "\t" <<
        (m_dumpij.first + 1) << "\n";
    for (int j = -1; j <= 1; ++j)
    {
        int jj = m_dumpij.second + j;

        std::cout << jj << "\t";
        for (int i = -1; i <= 1; ++i)
        {
            int ii = m_dumpij.first + i;
            Point p = m_field.offset(Coord(ii, jj));
            double disp = (std::sqrt(p.x * p.x + p.y * p.y));
            std::cout << disp << "\t";
        }
        std::cout << "\n";
    }
}

bool Atlas::process(Coord coord, Point& offset, double spread)
{
    using namespace pdal;

    if (! m_field.valid(coord))
        return false;

    SplitterFilter *splitter = dynamic_cast<SplitterFilter *>(
        m_beforeMgr.getStage());

    PointViewPtr v = splitter->view(splitterCoord(coord));
    if (!v)
        return false;
    m_field.setBeforeCount(coord, v->size());

    pdal::PointViewPtr debugView;
    if (coord == m_dumpij)
        debugView = v->makeNew();
    Histogram hist = histogram(v, Dimension::Id::Z, debugView, "Before");
    for (int pass = 0; pass < 3; ++pass)
        if (!removeFliers(v, hist))
            break;
        else
            hist = histogram(v, Dimension::Id::Z, debugView, "Before");

    PointViewPtr slice = hist[hist.size() - 1];
    if (slice->empty())
        return false;

    // box is the nominal core tile bounds (without overlap).
    // Shift the grid origin back by m_overlap so buffer-zone points bin into
    // negative grid indices and flood-fill can cross the tile boundary.
    BOX2D box = splitter->bounds(splitterCoord(coord));
    Point bgOrigin { box.minx - m_overlap, box.miny - m_overlap };
    if (coord == m_dumpij)
        slice = debugView;
    GridPtr bg = buildGrid(slice, bgOrigin);
    bg->findShapes(2);

    // Remove shapes whose center is outside the core tile. These are pure
    // buffer-zone blobs that belong to an adjacent cell's core area.
    const double gridLen = 2.0;
    auto inCoreBox = [&](const Shape& s, const Point& origin, const BOX2D& core) {
        double utmX = origin.x + (s.exactCenter().x + 0.5) * gridLen;
        double utmY = origin.y + (s.exactCenter().y + 0.5) * gridLen;
        return utmX >= core.minx && utmX < core.maxx &&
               utmY >= core.miny && utmY < core.maxy;
    };

    auto& bgShapes = bg->shapes();
    bgShapes.erase(
        std::remove_if(bgShapes.begin(), bgShapes.end(),
            [&](const Shape& s){ return !inCoreBox(s, bgOrigin, box); }),
        bgShapes.end());

    splitter = dynamic_cast<SplitterFilter *>(m_afterMgr.getStage());
    v = splitter->view(splitterCoord(coord));
    if (!v)
        return false;
    m_field.setAfterCount(coord, v->size());

    if (coord == m_dumpij)
        debugView = v->makeNew();
    hist = histogram(v, Dimension::Id::Z, debugView, "After");
    for (int pass = 0; pass < 3; ++pass)
        if (!removeFliers(v, hist))
            break;
        else
            hist = histogram(v, Dimension::Id::Z, debugView, "After");
    slice = hist[hist.size() - 1];
    if (slice->empty())
        return false;

    // The splitter for the "after" tiles is offset by m_shift.
    // If we want the grid to line up with the tile, the offset should be
    // m_shift, but we use an offset here that tries to account for what
    // we've learned since we started processing.
    //
    // Note here that "box" is the origin of the BEFORE points.
    BOX2D afterCore { box.minx + offset.x, box.miny + offset.y,
                      box.maxx + offset.x, box.maxy + offset.y };
    Point agOrigin { box.minx + offset.x - m_overlap, box.miny + offset.y - m_overlap };
    if (coord == m_dumpij)
        slice = debugView;
    GridPtr ag = buildGrid(slice, agOrigin);
    ag->findShapes(2);

    auto& agShapes = ag->shapes();
    agShapes.erase(
        std::remove_if(agShapes.begin(), agShapes.end(),
            [&](const Shape& s){ return !inCoreBox(s, agOrigin, afterCore); }),
        agShapes.end());



    sortShapes(bg);
    sortShapes(ag);

    std::vector<ShapePair> shapes = matchShapes(bg, ag, spread);
    recordShapes(bg, true,  coord, bgOrigin, shapes);
    recordShapes(ag, false, coord, agOrigin, shapes);
    if (shapes.empty())
        return false;
    auto [mean, median, rmsResidual] = calculateOffset(bg, ag, shapes);
    offset = mean;
    m_field.setMedianOffset(coord, median);
    m_field.setMatchQuality(coord, shapes.size(), rmsResidual);



    if (m_dumpij == coord)
    {
        /**
        dumpShapes(bg);
        dumpShapes(ag);
        **/

        // Draw with point counts.
        bg->draw(0, 0, 49, 49, false);
        ag->draw(0, 0, 49, 49, false);

        // Draw with shape numbers.
        bg->draw(0, 0, 49, 49, true);
        ag->draw(0, 0, 49, 49, true);
    }

    return true;
}


void Atlas::sortShapes(GridPtr& g)
{
    std::vector<Shape>& shapes = g->shapes();
    sort(shapes.begin(), shapes.end(),
        [](const Shape& s1, const Shape& s2) { return s1.size() > s2.size(); });
}


void Atlas::dumpShapes(GridPtr& g)
{
    std::vector<Shape>& shapes = g->shapes();
    for (const Shape& s : shapes)
        std::cout << "Size = " << s.size() << "/" <<
            s.center().first << "/" << s.center().second << "!\n";
    std::cout << "\n";
}


std::vector<ShapePair> Atlas::matchShapes(GridPtr& bg, GridPtr& ag, double spread)
{
    std::vector<ShapePair> matches;
    double threshold = std::clamp(spread, 1.0, 3.0);

    // Build a list of shape pointers
    std::list<Shape *> asp;
    for (Shape& s : ag->shapes())
        asp.push_back(&s);

    int matchCount = 0;
    for (const Shape& s : bg->shapes())
    {
        std::list<Shape *>::iterator exactMatch;
        double minExactDist = (std::numeric_limits<double>::max)();
        for (auto it = asp.begin(); it != asp.end(); ++it)
        {
            Shape *ts = *it;

            double sizeRatio = (double)s.size() / ts->size();
            if (sizeRatio > 3.0 || sizeRatio < (1.0 / 3.0))
                continue;

            Point bp = s.exactCenter();
            Point ap = ts->exactCenter();
            double exactDist = std::sqrt(std::pow(bp.x - ap.x, 2) + std::pow(bp.y - ap.y, 2));
            if (exactDist < minExactDist)
            {
                minExactDist = exactDist;
                exactMatch = it;
            }
        }
        if (minExactDist > threshold)
            continue;

        if (m_dumpij == m_coord)
            std::cout << "Match: " << s.id() << " = " << (*exactMatch)->id() << "\n";
        ShapePair sp {&s, *exactMatch};
        matches.push_back(sp);
        asp.erase(exactMatch);
        matchCount++;
        if (asp.empty())
            break;
    }
/**
    std::cerr << "Matched " << matchCount << " of " << bg->shapes().size() <<
        "!\n\n";
**/
    return matches;
}


std::tuple<Point, Point, double> Atlas::calculateOffset(GridPtr& bg, GridPtr& ag,
    const std::vector<ShapePair>& shapes)
{
    std::vector<double> pairX, pairY, pairZ;

    for (const ShapePair& sp : shapes)
    {
        Point bCenter, bHigh;
        Point aCenter, aHigh;

        pdal::BOX2D bExtent = bg->location(sp.first, bCenter, bHigh);
        pdal::BOX2D aExtent = ag->location(sp.second, aCenter, aHigh);

        pairX.push_back(((aCenter.x - bCenter.x) +
                         (aExtent.minx - bExtent.minx) +
                         (aExtent.maxx - bExtent.maxx)) / 3.0);
        pairY.push_back(((aCenter.y - bCenter.y) +
                         (aExtent.miny - bExtent.miny) +
                         (aExtent.maxy - bExtent.maxy)) / 3.0);
        pairZ.push_back(aCenter.z - bCenter.z);
    }

    auto calcMedian = [](std::vector<double> v) -> double {
        std::sort(v.begin(), v.end());
        size_t n = v.size();
        return (n % 2) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
    };

    double sumX = std::accumulate(pairX.begin(), pairX.end(), 0.0);
    double sumY = std::accumulate(pairY.begin(), pairY.end(), 0.0);
    double sumZ = std::accumulate(pairZ.begin(), pairZ.end(), 0.0);
    size_t n = shapes.size();

    Point mean { sumX / n, sumY / n, sumZ / n };
    Point median { calcMedian(pairX), calcMedian(pairY), calcMedian(pairZ) };

    double rmsResidual = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double dx = pairX[i] - mean.x;
        double dy = pairY[i] - mean.y;
        rmsResidual += dx * dx + dy * dy;
    }
    rmsResidual = std::sqrt(rmsResidual / n);

    return { mean, median, rmsResidual };
}


Histogram Atlas::histogram(pdal::PointViewPtr v, pdal::Dimension::Id dim,
    pdal::PointViewPtr debugView, const std::string& label)
{
    using namespace pdal;

    double mx = std::numeric_limits<double>::lowest();
    double mn = (std::numeric_limits<double>::max)();
    for (size_t id = 0; id < v->size(); ++id)
    {
        double val = v->getFieldAs<double>(dim, id);
        mx = (std::max)(mx, val);
        mn = (std::min)(mn, val);
    }
    double rangeInc = (mx - mn) / 10;

    Histogram splits;
    for (size_t idx = 0; idx < splits.size(); ++idx)
        splits[idx] = v->makeNew();

    for (size_t id = 0; id < v->size(); ++id)
    {
        double val = v->getFieldAs<double>(dim, id);
        size_t idx =  (val - mn) / rangeInc;
        if (idx > 9)
            idx = 9;
        splits[idx]->appendPoint(*v, id);

        if (debugView && (val - mn) / (mx - mn) > m_dumpfrac)
            debugView->appendPoint(*v, id);
    }

    // Dump histogram.
    if (shouldDump() && debugView)
    {
        size_t total = 0;
        std::cout << label << " - ";
        for (size_t i = 0; i < splits.size(); i++)
        {
            total += splits[i]->size();
            std::cout << splits[i]->size() << "    ";
        }
        std::cout << " - Total = " << total << "\n";
    }
    return splits;
}


bool Atlas::removeFliers(pdal::PointViewPtr& v, const Histogram& hist)
{
    double ev = 0;
    for (size_t i = 1; i <= hist.size(); ++i)
        ev += (hist[i - 1]->size() * i) / (double)v->size();
    if (ev < 2)
    {
        v = v->makeNew();
        for (size_t i = 0; i < hist.size() && hist[i]->size(); ++i)
            v->append(*hist[i]);
        return true;
    }
    return false;
}


GridPtr Atlas::buildGrid(pdal::PointViewPtr v, Point origin)
{
    using namespace pdal;

    auto cmp = [](const PointRef& p1, const PointRef& p2)
    {
        return p1.compare(Dimension::Id::Z, p2);
    };

    std::sort(v->begin(), v->end(), cmp);

    return GridPtr(new Grid(v, 2, origin));
}


void Atlas::writeSimple()
{
    char sep;
    const float *xp = m_field.xdata();
    const float *yp = m_field.ydata();

    for (size_t ypos = 0; ypos < m_field.height(); ++ypos)
    {
        for (size_t xpos = 0; xpos < m_field.width(); ++xpos)
        {
            int x = std::round(*xp * 10);
            int y = std::round(*yp * 10);
            if (y < 0)
            {
                sep = '-';
                y = -y;
            }
            else
                sep = '+';
            std::cout << std::setfill('0') << std::setw(2) << x << sep <<
                std::setfill('0') << std::setw(2) << y << " ";
            xp++; yp++;
        }
        std::cout << "\n";
    }
}


void Atlas::writeSvg(const std::string& filename, const pdal::BOX2D& extent)
{
    size_t xsize = extent.maxx - extent.minx + 1;
    size_t ysize = extent.maxy - extent.miny + 1;

    Draw d(2000, 1000, xsize, ysize, std::sqrt(m_field.maxLen2()));

    const float *xp = m_field.xdata();
    const float *yp = m_field.ydata();

    for (size_t ypos = 0; ypos < m_field.height(); ++ypos)
        for (size_t xpos = 0; xpos < m_field.width(); ++xpos)
            d.doVector({xpos, ypos}, {*xp++, *yp++});
}


void Atlas::writeTiff(const std::string& filename)
{
    using namespace pdal;

    size_t xsize = m_field.width();
    size_t ysize = m_field.height();

    std::array<double, 6> pixelToPos;
    pixelToPos[0] = m_origin.x;
    pixelToPos[1] = m_len;
    pixelToPos[2] = 0;
    pixelToPos[3] = m_origin.y + (m_len * ysize);
    pixelToPos[4] = 0;
    pixelToPos[5] = -m_len;

    gdal::registerDrivers();
    gdal::Raster raster(filename, "GTiff", "EPSG:32624", pixelToPos);
    gdal::GDALError err = raster.open(xsize, ysize, 10,
        Dimension::Type::Float, -9999, pdal::StringList());

    if (err != gdal::GDALError::None)
        throwError(raster.errorMsg());
    {
        raster.writeBand(m_field.xdata(),            -9999.0, 1,  "X");
        raster.writeBand(m_field.ydata(),            -9999.0, 2,  "Y");
        raster.writeBand(m_field.zdata(),            -9999.0, 3,  "Z");
        raster.writeBand(m_field.bdata(),            -9999.0, 4,  "BEFORE");
        raster.writeBand(m_field.adata(),            -9999.0, 5,  "AFTER");
        raster.writeBand(m_field.medianXdata(),      -9999.0, 6,  "MEDIAN_X");
        raster.writeBand(m_field.medianYdata(),      -9999.0, 7,  "MEDIAN_Y");
        raster.writeBand(m_field.medianZdata(),      -9999.0, 8,  "MEDIAN_Z");
        raster.writeBand(m_field.matchCountData(),   -9999.0, 9,  "MATCH_COUNT");
        raster.writeBand(m_field.rmsResidualData(),  -9999.0, 10, "RMS_RESIDUAL");
    }
}


void Atlas::read(const std::string& filename)
{
    using namespace pdal;

    GDALReader gr;

    Options opts;
    opts.add("filename", filename);
    opts.add("header", "VX,VY");
    gr.setOptions(opts);

    PointTable t;
    gr.prepare(t);
    PointViewSet s = gr.execute(t);
    PointViewPtr v = *s.begin();
    PointLayoutPtr l = t.layout();
    Dimension::Id Vx = l->findDim("VX");
    Dimension::Id Vy = l->findDim("VY");

    for (size_t i = 0; i < v->size(); ++i)
    {
        double x = v->getFieldAs<double>(Dimension::Id::X, i);
        double y = v->getFieldAs<double>(Dimension::Id::Y, i);
        double vx = v->getFieldAs<double>(Vx, i);
        double vy = v->getFieldAs<double>(Vy, i);

        std::cout << "x/y/vx/vy/total " <<
            x << "/" << y << "/" << vx << "/" << vy << "/" << std::sqrt(vx * vx + vy * vy) << "!\n";
    }
}

void Atlas::recordShapes(GridPtr& g, bool isBefore, Coord tile, Point origin,
    const std::vector<ShapePair>& matches)
{
    for (const Shape& s : g->shapes())
    {
        bool matched = false;
        size_t matchId = std::numeric_limits<size_t>::max();
        for (const ShapePair& sp : matches)
        {
            const Shape* check   = isBefore ? sp.first  : sp.second;
            const Shape* partner = isBefore ? sp.second : sp.first;
            if (check == &s)
            {
                matched = true;
                matchId = partner->id();
                break;
            }
        }
        m_shapeRecords.push_back({s.indices(), origin, isBefore, tile, s.id(), matched, matchId});
    }
}


void Atlas::writeShapeGeoJSON(const std::string& filename)
{
    std::ofstream out(filename);
    out << std::fixed << std::setprecision(2);
    out << "{\"type\":\"FeatureCollection\","
           "\"crs\":{\"type\":\"name\",\"properties\":{\"name\":\"EPSG:32624\"}},"
           "\"features\":[\n";

    bool firstFeature = true;
    for (const ShapeRecord& rec : m_shapeRecords)
    {
        if (!firstFeature) out << ",\n";
        firstFeature = false;

        out << "{\"type\":\"Feature\",\"geometry\":"
               "{\"type\":\"MultiPolygon\",\"coordinates\":[";
        bool firstCell = true;
        for (const GridIndex& gi : rec.indices)
        {
            if (!firstCell) out << ",";
            firstCell = false;
            double minx = rec.origin.x + gi.x() * 2.0;
            double miny = rec.origin.y + gi.y() * 2.0;
            double maxx = minx + 2.0;
            double maxy = miny + 2.0;
            out << "[["
                << "[" << minx << "," << miny << "],"
                << "[" << maxx << "," << miny << "],"
                << "[" << maxx << "," << maxy << "],"
                << "[" << minx << "," << maxy << "],"
                << "[" << minx << "," << miny << "]"
                << "]]";
        }
        // Flip tile_j from internal Y-down to Y-up so it aligns with the
        // map orientation in QGIS (large tile_j = north).
        int tileJ = YCellCount - rec.tile.second - 1;
        out << "]},\"properties\":{"
            << "\"scan\":\"" << (rec.isBefore ? "before" : "after") << "\","
            << "\"tile_i\":" << rec.tile.first << ","
            << "\"tile_j\":" << tileJ << ","
            << "\"shape_id\":" << rec.id << ","
            << "\"size\":" << rec.indices.size() << ","
            << "\"matched\":" << (rec.matched ? "true" : "false") << ",";
        if (rec.matched)
            out << "\"match_id\":" << rec.matchId;
        else
            out << "\"match_id\":null";
        out << "}}";
    }
    out << "\n]}\n";
}

} // namespace
