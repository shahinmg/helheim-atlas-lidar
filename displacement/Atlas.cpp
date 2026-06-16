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

#include <gdal_priv.h>

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
    m_args.add("topfrac", "Top fraction of points by height above the tile's fitted "
        "plane used for surface slice (default 0.5)", m_dumpfrac, 0.5);
    m_args.add("minshape", "Minimum shape size in grid cells; smaller shapes are dropped "
        "before matching (default 1, disables filter)", m_minShape, 1);
    m_args.add("passes", "Number of grid passes; passes after the first reprocess "
        "every cell seeded from its converged neighborhood (default 1)", m_passes, 1);
    m_args.add("slices", "Number of contiguous height bands the top 'topfrac' is "
        "split into for shape detection; shapes from all bands are matched and "
        "their displacements pooled (default 1 = single top slice)", m_slices, 1);
    m_args.add("gridlen", "Grid cell size in meters for flood-fill shape detection "
        "(default 1.0)", m_gridLen, 1.0);
    m_args.add("ncc", "Refine each cell's shape-match displacement by normalized "
        "cross-correlation of detrended mean-height rasters, with subpixel peak fit. "
        "Adds bands NCC_X, NCC_Y, NCC_PEAK.", m_ncc);
    m_args.add("nccradius", "NCC search radius in raster cells around the shape-match "
        "displacement (default 5)", m_nccRadius, 5);
    m_args.add("ncclen", "NCC raster cell size in meters, independent of 'gridlen' "
        "(default 2.0)", m_nccLen, 2.0);
    m_args.add("zgate", "Max mean-height difference in meters for a shape match; "
        "rejects elevation-inconsistent pairs so the distance gate can be loosened "
        "without admitting false matches. 0 disables (default 0)", m_zGate, 0.0);
    m_args.add("zgate-plane", "Gate 'zgate' on mean height above the tile's fitted "
        "plane instead of absolute Z, removing the bulk slope term so one zgate "
        "value works across baselines (requires zgate > 0)", m_zGatePlane);
    m_args.add("matchmin", "Lower bound in meters for the centroid-distance match "
        "threshold, applied at well-converged interior cells; lower to tighten "
        "matching for short baselines (default 2.5)", m_matchMin, 2.5);
    m_args.add("matchmax", "Upper bound in meters for the centroid-distance match "
        "threshold; raise to admit more matches (default 4.0)", m_matchMax, 4.0);
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
    if (m_dumpfrac <= 0 || m_dumpfrac > 1)
        fatal("'topfrac' must be a value in the range (0,1].");
    if (m_minShape < 1)
        fatal("'minshape' must be >= 1.");
    if (m_passes < 1)
        fatal("'passes' must be >= 1.");
    if (m_slices < 1)
        fatal("'slices' must be >= 1.");
    if (m_gridLen <= 0)
        fatal("'gridlen' must be > 0.");
    if (m_zGate < 0)
        fatal("'zgate' must be >= 0.");
    if (m_zGatePlane && m_zGate <= 0)
        fatal("'zgate-plane' requires 'zgate' > 0.");
    if (m_matchMin <= 0)
        fatal("'matchmin' must be > 0.");
    if (m_matchMax < m_matchMin)
        fatal("'matchmax' must be >= 'matchmin'.");
    if (m_nccRadius < 1)
        fatal("'nccradius' must be >= 1.");
    if (m_nccLen <= 0)
        fatal("'ncclen' must be > 0.");
    if (m_nccRadius * m_nccLen > m_overlap)
        fatal("'nccradius' x 'ncclen' must not exceed the tile overlap (" +
            std::to_string((int)m_overlap) + "m).");
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
        auto stripCopc = [](std::string s) -> std::string {
            auto pos = s.rfind(".copc");
            if (pos != std::string::npos)
                s.erase(pos, 5);
            return s;
        };
        std::string stem = stripCopc(pdal::FileUtils::stem(m_beforeFilename)) + "_" +
            stripCopc(pdal::FileUtils::stem(m_afterFilename));

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
    splitterOpts.add("buffer", m_overlap);

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
    afterSplitterOpts.add("buffer", m_overlap);

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

    // The first pass solves the field spiraling out from the center, seeding
    // each cell from already-solved neighbors. Optional extra passes
    // (--passes) reprocess every cell seeded from its full converged
    // neighborhood, which removes the one-sided frontier bias of the spiral
    // order and gives failed cells another chance — at the cost of a full
    // grid traversal each.
    for (int pass = 0; pass < m_passes; ++pass)
    {
        // Shape records from earlier passes are superseded.
        if (pass > 0)
            m_shapeRecords.clear();

        for (int i = 0; i <= maxRadius; ++i)
        {
            int y = -i;
            int x = -i;
            while (true)
            {
                Coord c{start.first + x, start.second + y};
                m_coord = c;

                // If this cell has an initial offset, use it, otherwise use
                // the base one. On later passes, seed from the neighbor
                // average even when the cell itself was solved.
                auto [valid, pos, spread] = m_field.initialOffset(c, pass > 0);
                if (!valid && pass > 0)
                {
                    // No valid neighbors; fall back to the cell's own value.
                    auto own = m_field.initialOffset(c, false);
                    valid = std::get<0>(own);
                    pos = std::get<1>(own);
                    spread = std::get<2>(own);
                }
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

    PlaneFit bPlane = fitPlane(v);
    PointViewPtr beforeView = v;

    // box is the nominal core tile bounds (without overlap), from the BEFORE
    // splitter (still 'splitter' here). Shift the grid origin back by m_overlap
    // so buffer-zone points bin into negative grid indices and flood-fill can
    // cross the tile boundary.
    BOX2D box = splitter->bounds(splitterCoord(coord));
    Point bgOrigin { box.minx - m_overlap, box.miny - m_overlap };

    // Drop shapes whose center is outside the core tile (pure buffer-zone blobs
    // belonging to an adjacent cell) or smaller than the minimum size.
    auto inCoreBox = [&](const Shape& s, const Point& origin, const BOX2D& core) {
        double utmX = origin.x + (s.exactCenter().x + 0.5) * m_gridLen;
        double utmY = origin.y + (s.exactCenter().y + 0.5) * m_gridLen;
        return utmX >= core.minx && utmX < core.maxx &&
               utmY >= core.miny && utmY < core.maxy;
    };
    auto pruneShapes = [&](GridPtr& g, const Point& origin, const BOX2D& core) {
        auto& sh = g->shapes();
        sh.erase(std::remove_if(sh.begin(), sh.end(),
            [&](const Shape& s){
                return !inCoreBox(s, origin, core) || s.size() < (size_t)m_minShape;
            }), sh.end());
    };

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
    PlaneFit aPlane = fitPlane(v);
    PointViewPtr afterView = v;

    // The "after" tile grid is placed at the BEFORE origin shifted by the
    // current displacement estimate ("box" is the origin of the BEFORE points).
    BOX2D afterCore { box.minx + offset.x, box.miny + offset.y,
                      box.maxx + offset.x, box.maxy + offset.y };
    Point agOrigin { box.minx + offset.x - m_overlap, box.miny + offset.y - m_overlap };

    // Sort each cloud by height above its plane ONCE; bands are sliced from the
    // rank below rather than re-sorting per band.
    std::vector<PointId> bRank = rankByResidual(beforeView, bPlane);
    std::vector<PointId> aRank = rankByResidual(afterView, aPlane);

    // Detect and match shapes over m_slices contiguous residual bands of the top
    // 'topfrac', pooling the per-pair displacements. With m_slices == 1 the band
    // is [1-topfrac, 1] -- the original top slice -- so the result is identical
    // to the single-slice path. Splitting the top fraction into height bands lets
    // features that merge into one blob at a single threshold be detected
    // separately (more candidate shapes -> more matches). Bands run top-first and
    // a pair is dropped if its before-centroid coincides with one already pooled
    // from a higher band, so a feature spanning bands is counted once -- keeping
    // MATCH_COUNT a count of distinct matches for the downstream sigma=RMS/sqrt(N).
    std::vector<double> pairX, pairY, pairZ;
    std::vector<Point> centers;   // before-centroid (UTM) of each pooled pair
    for (int k = m_slices - 1; k >= 0; --k)
    {
        double lo = (1.0 - m_dumpfrac) + (double)k * (m_dumpfrac / m_slices);
        // Pin the top band's upper edge to 1.0 so floating-point error in
        // k*topfrac/N can't drop the single highest-residual point.
        double hi = (k == m_slices - 1) ? 1.0
                  : (1.0 - m_dumpfrac) + (double)(k + 1) * (m_dumpfrac / m_slices);

        PointViewPtr bSlice = sliceFromRank(beforeView, bRank, lo, hi);
        PointViewPtr aSlice = sliceFromRank(afterView, aRank, lo, hi);
        if (bSlice->empty() || aSlice->empty())
            continue;

        GridPtr bg = buildGrid(bSlice, bgOrigin);
        bg->findShapes(2);
        pruneShapes(bg, bgOrigin, box);

        GridPtr ag = buildGrid(aSlice, agOrigin);
        ag->findShapes(2);
        pruneShapes(ag, agOrigin, afterCore);

        sortShapes(bg);
        sortShapes(ag);

        if (m_dumpij == coord && k == 0)
        {
            bg->draw(0, 0, 49, 49, false);
            ag->draw(0, 0, 49, 49, false);
            bg->draw(0, 0, 49, 49, true);
            ag->draw(0, 0, 49, 49, true);
        }

        std::vector<ShapePair> shapes = matchShapes(bg, ag, spread, bPlane, aPlane);
        recordShapes(bg, true,  coord, bgOrigin, shapes);
        recordShapes(ag, false, coord, agOrigin, shapes);
        size_t priorBands = centers.size();
        collectPairDisplacements(bg, ag, shapes, pairX, pairY, pairZ,
            centers, priorBands);
    }

    if (pairX.empty())
        return false;

    auto [mean, median, rmsResidual, mad] = aggregateOffset(pairX, pairY, pairZ);
    offset = mean;
    m_field.setMedianOffset(coord, median);
    m_field.setMadOffset(coord, mad);
    m_field.setMatchQuality(coord, pairX.size(), rmsResidual);

    if (m_ncc)
    {
        // Rasterize the full flier-removed clouds (not the top slice — NCC
        // wants the whole surface, the slice exists only for blob detection)
        // over the core tile plus the overlap buffer. The after raster
        // covers the same area shifted by the converged shape-match offset,
        // with a search margin; cells outside the after view's coverage stay
        // NaN and drop out of the correlation. NCC then measures the
        // residual displacement the shape match missed.
        int rw = (int)std::ceil((m_len + 2 * m_overlap) / m_nccLen);
        int rh = rw;
        int R = m_nccRadius;
        double margin = R * m_nccLen;

        // Splat kernel adapts to the sparser of the two tiles: dense tiles
        // get sharp single-cell binning (best subpixel precision), sparse
        // tiles get progressively wider Gaussian splats so the before/after
        // rasters still overlap. Both rasters use the same kernel to keep
        // the smoothing symmetric.
        double density = (double)std::min(beforeView->size(), v->size()) /
            ((double)rw * rh);
        int kernel = density >= 0.15 ? 0 : (density >= 0.05 ? 1 : 2);

        Raster br = buildRaster(beforeView,
            {box.minx - m_overlap, box.miny - m_overlap}, rw, rh, bPlane,
            kernel);
        Raster ar = buildRaster(v,
            {box.minx - m_overlap + offset.x - margin,
             box.miny - m_overlap + offset.y - margin},
            rw + 2 * R, rh + 2 * R, aPlane, kernel);

        auto [residual, peak, onEdge] = nccOffset(br, ar, R);
        if (peak >= -1.0)
            m_field.setNccPeak(coord, (float)peak);

        // Accept the refinement only with a confident, interior peak.
        // The refined total goes to the NCC bands; bands 1-13 and the
        // neighbor-seeding field keep the pure shape-match result so the
        // two estimates stay directly comparable per cell.
        const double minPeak = 0.5;
        if (peak >= minPeak && !onEdge)
            m_field.setNccOffset(coord,
                { offset.x + residual.x, offset.y + residual.y });
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


std::vector<ShapePair> Atlas::matchShapes(GridPtr& bg, GridPtr& ag, double spread,
    const PlaneFit& bPlane, const PlaneFit& aPlane)
{
    std::vector<ShapePair> matches;
    double threshold = std::clamp(spread, m_matchMin, m_matchMax);

    // Mean height per shape — a descriptor that lets us admit matches at larger
    // horizontal distance without false pairings: a genuine pair has nearly
    // equal mean height, while a spurious match to a horizontally-distant shape
    // on sloped ground does not. With --zgate-plane the height is taken above
    // the tile's fitted plane instead of absolute Z, which cancels the bulk
    // flow-down-slope term (= displacement x tan(slope)) so one zgate works
    // across baselines. The plane is affine, so the mean height-above-plane of
    // a shape is just plane.residual at its centroid (no per-point pass).
    // Computed once per shape (reusing Grid::location) so the O(n^2) mutual-NN
    // search below stays cheap. Only populated when the gate is on.
    std::unordered_map<const Shape *, double> meanZ;
    auto loadHeights = [&](GridPtr& g, const PlaneFit& plane) {
        Point c, hi;
        for (Shape& s : g->shapes())
        {
            g->location(&s, c, hi);
            meanZ[&s] = m_zGatePlane ? plane.residual(c.x, c.y, c.z) : c.z;
        }
    };
    if (m_zGate > 0)
    {
        loadHeights(bg, bPlane);
        loadHeights(ag, aPlane);
    }

    // Closest candidate passing the 3:1 size-ratio gate and, when enabled, the
    // mean-height gate. Returns the candidate and its horizontal distance.
    auto nearest = [&](const Shape& s, std::vector<Shape>& candidates)
        -> std::pair<const Shape *, double>
    {
        const Shape *best = nullptr;
        double minDist = (std::numeric_limits<double>::max)();
        for (const Shape& ts : candidates)
        {
            double sizeRatio = (double)s.size() / ts.size();
            if (sizeRatio > 3.0 || sizeRatio < (1.0 / 3.0))
                continue;
            if (m_zGate > 0 && std::abs(meanZ[&s] - meanZ[&ts]) > m_zGate)
                continue;

            Point sp = s.exactCenter();
            Point tp = ts.exactCenter();
            double dist = std::sqrt(std::pow(sp.x - tp.x, 2) +
                                    std::pow(sp.y - tp.y, 2));
            if (dist < minDist)
            {
                minDist = dist;
                best = &ts;
            }
        }
        return { best, minDist };
    };

    // Mutual nearest neighbor: accept a pair only if each shape is the other's
    // closest eligible candidate. One-to-one by construction and independent of
    // visit order, so a large shape can't claim an after-shape that better
    // matches a different before-shape (the old greedy matcher's failure mode).
    for (const Shape& s : bg->shapes())
    {
        auto [after, dist] = nearest(s, ag->shapes());
        if (!after || dist > threshold)
            continue;
        if (nearest(*after, bg->shapes()).first != &s)
            continue;

        if (m_dumpij == m_coord)
            std::cout << "Match: " << s.id() << " = " << after->id() << "\n";
        matches.push_back({ &s, after });
    }
    return matches;
}


// Append the per-pair (x, y, z) displacement of each matched shape pair to the
// running vectors. Split out from the aggregation so multiple slice bands can
// pool their pairs before the mean/median/MAD/RMS are computed once.
void Atlas::collectPairDisplacements(GridPtr& bg, GridPtr& ag,
    const std::vector<ShapePair>& shapes,
    std::vector<double>& pairX, std::vector<double>& pairY,
    std::vector<double>& pairZ, std::vector<Point>& centers,
    size_t dedupCount)
{
    for (const ShapePair& sp : shapes)
    {
        Point bCenter, bHigh;
        Point aCenter, aHigh;

        pdal::BOX2D bExtent = bg->location(sp.first, bCenter, bHigh);
        pdal::BOX2D aExtent = ag->location(sp.second, aCenter, aHigh);

        // Drop a pair whose before-centroid coincides (within one grid cell)
        // with a pair already pooled from an EARLIER slice band -- the same
        // feature re-detected at a different height. Only the first dedupCount
        // centers (prior bands) are checked, so pairs within the current band
        // are never removed and --slices 1 (dedupCount 0) is unchanged.
        bool dup = false;
        for (size_t j = 0; j < dedupCount && j < centers.size(); ++j)
            if (std::hypot(bCenter.x - centers[j].x,
                           bCenter.y - centers[j].y) < m_gridLen)
            {
                dup = true;
                break;
            }
        if (dup)
            continue;

        pairX.push_back(((aCenter.x - bCenter.x) +
                         (aExtent.minx - bExtent.minx) +
                         (aExtent.maxx - bExtent.maxx)) / 3.0);
        pairY.push_back(((aCenter.y - bCenter.y) +
                         (aExtent.miny - bExtent.miny) +
                         (aExtent.maxy - bExtent.maxy)) / 3.0);
        pairZ.push_back(aCenter.z - bCenter.z);
        centers.push_back(bCenter);
    }
}


// Mean / median / RMS-residual / MAD over the pooled per-pair displacements.
std::tuple<Point, Point, double, Point> Atlas::aggregateOffset(
    const std::vector<double>& pairX, const std::vector<double>& pairY,
    const std::vector<double>& pairZ)
{
    auto calcMedian = [](std::vector<double> v) -> double {
        std::sort(v.begin(), v.end());
        size_t n = v.size();
        return (n % 2) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
    };

    auto calcMAD = [&calcMedian](const std::vector<double>& v, double med) -> double {
        std::vector<double> absdev(v.size());
        for (size_t i = 0; i < v.size(); ++i)
            absdev[i] = std::abs(v[i] - med);
        return calcMedian(absdev);
    };

    double sumX = std::accumulate(pairX.begin(), pairX.end(), 0.0);
    double sumY = std::accumulate(pairY.begin(), pairY.end(), 0.0);
    double sumZ = std::accumulate(pairZ.begin(), pairZ.end(), 0.0);
    size_t n = pairX.size();

    Point mean { sumX / n, sumY / n, sumZ / n };
    Point median { calcMedian(pairX), calcMedian(pairY), calcMedian(pairZ) };
    Point mad { calcMAD(pairX, median.x), calcMAD(pairY, median.y), calcMAD(pairZ, median.z) };

    double rmsResidual = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double dx = pairX[i] - mean.x;
        double dy = pairY[i] - mean.y;
        rmsResidual += dx * dx + dy * dy;
    }
    rmsResidual = std::sqrt(rmsResidual / n);

    return { mean, median, rmsResidual, mad };
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

        if (debugView && (val - mn) / (mx - mn) > (1.0 - m_dumpfrac))
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


PlaneFit Atlas::fitPlane(pdal::PointViewPtr v)
{
    using namespace pdal;

    PlaneFit p { 0, 0, 0, 0, 0 };
    size_t n = v->size();
    if (n == 0)
        return p;

    for (PointId i = 0; i < n; ++i)
    {
        p.mx += v->getFieldAs<double>(Dimension::Id::X, i);
        p.my += v->getFieldAs<double>(Dimension::Id::Y, i);
        p.mz += v->getFieldAs<double>(Dimension::Id::Z, i);
    }
    p.mx /= n;
    p.my /= n;
    p.mz /= n;

    double sxx = 0, sxy = 0, syy = 0, sxz = 0, syz = 0;
    for (PointId i = 0; i < n; ++i)
    {
        double dx = v->getFieldAs<double>(Dimension::Id::X, i) - p.mx;
        double dy = v->getFieldAs<double>(Dimension::Id::Y, i) - p.my;
        double dz = v->getFieldAs<double>(Dimension::Id::Z, i) - p.mz;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
        sxz += dx * dz;
        syz += dy * dz;
    }

    // Plane slopes from the centered normal equations. Degenerate point
    // geometry (e.g. a single scan line) gets a flat plane, which reduces
    // to raw-Z ranking.
    double det = sxx * syy - sxy * sxy;
    if (det > 1e-10 * sxx * syy)
    {
        p.bx = (sxz * syy - syz * sxy) / det;
        p.by = (syz * sxx - sxz * sxy) / det;
    }
    return p;
}


// Return the top 'topfrac' of points ranked by height above the tile's
// least-squares plane. Ranking on plane residuals rather than raw Z keeps
// the slice from collapsing to the uphill portion of sloped tiles, so the
// same locally-high features (serac tops, crevasse ridges) are selected
// across the whole tile in both epochs. The points themselves are not
// modified; only slice membership depends on the plane.
pdal::PointViewPtr Atlas::surfaceSlice(pdal::PointViewPtr v, const PlaneFit& plane)
{
    // The top 'topfrac' of points = the residual-rank band [1 - topfrac, 1].
    return sliceFromRank(v, rankByResidual(v, plane), 1.0 - m_dumpfrac, 1.0);
}


// PointIds of v sorted ascending by height above the plane (0 = lowest above-
// plane point). Computing the rank once and slicing bands from it keeps
// multi-slice O(n log n) per cloud instead of re-sorting for every band.
std::vector<pdal::PointId> Atlas::rankByResidual(pdal::PointViewPtr v,
    const PlaneFit& plane)
{
    using namespace pdal;

    std::vector<std::pair<double, PointId>> residuals;
    residuals.reserve(v->size());
    for (PointId i = 0; i < v->size(); ++i)
    {
        double x = v->getFieldAs<double>(Dimension::Id::X, i);
        double y = v->getFieldAs<double>(Dimension::Id::Y, i);
        double z = v->getFieldAs<double>(Dimension::Id::Z, i);
        residuals.push_back({plane.residual(x, y, z), i});
    }
    std::sort(residuals.begin(), residuals.end());

    std::vector<PointId> rank;
    rank.reserve(residuals.size());
    for (const auto& pr : residuals)
        rank.push_back(pr.second);
    return rank;
}


// Points whose residual rank falls in [rankLo, rankHi) (0 = lowest above-plane,
// 1 = highest), from a precomputed rank. surfaceSlice is the band [1-topfrac, 1];
// multi-slice (--slices N) splits that into N contiguous bands so locally-high
// features that merge at one threshold are detected separately. Points are not
// modified; only band membership uses the plane.
pdal::PointViewPtr Atlas::sliceFromRank(pdal::PointViewPtr v,
    const std::vector<pdal::PointId>& rank, double rankLo, double rankHi)
{
    pdal::PointViewPtr slice = v->makeNew();
    size_t n = rank.size();
    if (n == 0)
        return slice;

    size_t loIdx = (size_t)(n * rankLo);
    size_t hiIdx = (size_t)(n * rankHi);
    if (hiIdx > n)
        hiIdx = n;
    for (size_t i = loIdx; i < hiIdx; ++i)
        slice->appendPoint(*v, rank[i]);
    return slice;
}


// Splat raster of height-above-plane at ncclen resolution. 'kernel' is the
// splat radius in cells: 0 bins each point into its containing cell (sharp,
// for dense tiles); 1 or 2 spreads each point over a Gaussian neighborhood
// (sigma = kernel/2) so sparse tiles produce a continuous weighted surface
// instead of isolated single-point cells whose joint before/after occupancy
// is near zero. Detrending by the tile's fitted plane removes slope, so the
// correlation peak reflects feature structure rather than surface tilt.
Raster Atlas::buildRaster(pdal::PointViewPtr v, Point origin, int width,
    int height, const PlaneFit& plane, int kernel)
{
    using namespace pdal;

    int K = kernel;
    double invTwoSigma2 = K ? 1.0 / (2 * (0.5 * K) * (0.5 * K)) : 0;

    Raster r(width, height);
    for (PointId i = 0; i < v->size(); ++i)
    {
        double x = v->getFieldAs<double>(Dimension::Id::X, i);
        double y = v->getFieldAs<double>(Dimension::Id::Y, i);
        double z = v->getFieldAs<double>(Dimension::Id::Z, i);

        // Fractional position in cell units.
        double fx = (x - origin.x) / m_nccLen;
        double fy = (y - origin.y) / m_nccLen;
        int cx = (int)std::floor(fx);
        int cy = (int)std::floor(fy);
        if (cx < -K || cx >= width + K || cy < -K || cy >= height + K)
            continue;

        float res = (float)plane.residual(x, y, z);
        for (int yi = cy - K; yi <= cy + K; ++yi)
        {
            if (yi < 0 || yi >= height)
                continue;
            for (int xi = cx - K; xi <= cx + K; ++xi)
            {
                if (xi < 0 || xi >= width)
                    continue;
                double dx = xi + 0.5 - fx;
                double dy = yi + 0.5 - fy;
                float wgt = K ?
                    (float)std::exp(-(dx * dx + dy * dy) * invTwoSigma2) :
                    1.f;
                r.at(xi, yi) += wgt * res;
                r.weight(xi, yi) += wgt;
            }
        }
    }
    for (int i = 0; i < width * height; ++i)
        if (r.w[i] > 0)
            r.z[i] /= r.w[i];
    return r;
}


// Weighted ZNCC search over integer cell shifts with subpixel peak fit.
// Cells are weighted by the product of the two splat weights, so halo-only
// cells contribute in proportion to the data behind them. 'ar' must extend
// 'searchRadius' cells beyond 'br' on every side so each candidate shift
// has full support. Returns {residual displacement in meters, peak
// correlation, peak-on-search-edge}; a peak of -2 means no shift had
// enough effective support to score.
std::tuple<Point, double, bool> Atlas::nccOffset(const Raster& br,
    const Raster& ar, int searchRadius)
{
    // Effective sample size (sum(w)^2 / sum(w^2)) required before a
    // correlation score is meaningful.
    const double minEffective = 50;

    int R = searchRadius;
    int diameter = 2 * R + 1;
    std::vector<double> surface(diameter * diameter, -2.0);
    auto nccAt = [&](int dx, int dy) -> double& {
        return surface[(dy + R) * diameter + (dx + R)];
    };

    double bestNcc = -2.0;
    int bestDx = 0, bestDy = 0;

    for (int dy = -R; dy <= R; ++dy)
        for (int dx = -R; dx <= R; ++dx)
        {
            double sw = 0, sw2 = 0, sumB = 0, sumA = 0;
            for (int y = 0; y < br.height; ++y)
                for (int x = 0; x < br.width; ++x)
                {
                    double w = br.weight(x, y) *
                        ar.weight(x + dx + R, y + dy + R);
                    if (w <= 0)
                        continue;
                    sw += w;
                    sw2 += w * w;
                    sumB += w * br.at(x, y);
                    sumA += w * ar.at(x + dx + R, y + dy + R);
                }
            if (sw <= 0 || (sw * sw) / sw2 < minEffective)
                continue;

            double meanB = sumB / sw;
            double meanA = sumA / sw;
            double num = 0, varB = 0, varA = 0;
            for (int y = 0; y < br.height; ++y)
                for (int x = 0; x < br.width; ++x)
                {
                    double w = br.weight(x, y) *
                        ar.weight(x + dx + R, y + dy + R);
                    if (w <= 0)
                        continue;
                    double db = br.at(x, y) - meanB;
                    double da = ar.at(x + dx + R, y + dy + R) - meanA;
                    num += w * db * da;
                    varB += w * db * db;
                    varA += w * da * da;
                }
            double denom = std::sqrt(varB * varA);
            if (denom < 1e-10)
                continue;

            double ncc = num / denom;
            nccAt(dx, dy) = ncc;
            if (ncc > bestNcc)
            {
                bestNcc = ncc;
                bestDx = dx;
                bestDy = dy;
            }
        }

    if (bestNcc < -1.0)
        return { Point(0, 0), -2.0, true };

    bool onEdge = std::abs(bestDx) == R || std::abs(bestDy) == R;

    // Subpixel refinement, one axis at a time. Gaussian (log-parabola) fit:
    // correlation peaks are near-Gaussian, and fitting the parabola to the
    // log of the samples removes most of the pixel-locking bias a plain
    // parabola shows when the peak is narrower than the cell spacing.
    // Non-positive samples fall back to the plain parabolic fit. The second
    // difference must be positive (a genuine maximum). The vertex of the
    // parabola through (-1,lo), (0,c), (1,hi) is at +0.5*(hi-lo)/(2c-lo-hi).
    auto subpixel = [](double lo, double c, double hi) -> double
    {
        if (lo > 0 && c > 0 && hi > 0)
        {
            double ll = std::log(lo), lc = std::log(c), lh = std::log(hi);
            double d = 2 * lc - ll - lh;
            if (d > 1e-10)
                return 0.5 * (lh - ll) / d;
        }
        double d = 2 * c - lo - hi;
        if (d > 1e-10)
            return 0.5 * (hi - lo) / d;
        return 0.0;
    };

    double subDx = bestDx;
    double subDy = bestDy;
    if (!onEdge)
    {
        double l = nccAt(bestDx - 1, bestDy);
        double r = nccAt(bestDx + 1, bestDy);
        if (l > -1.5 && r > -1.5)
            subDx = bestDx + subpixel(l, bestNcc, r);
        double dn = nccAt(bestDx, bestDy - 1);
        double up = nccAt(bestDx, bestDy + 1);
        if (dn > -1.5 && up > -1.5)
            subDy = bestDy + subpixel(dn, bestNcc, up);
    }

    return { Point(subDx * m_nccLen, subDy * m_nccLen), bestNcc, onEdge };
}


GridPtr Atlas::buildGrid(pdal::PointViewPtr v, Point origin)
{
    using namespace pdal;

    auto cmp = [](const PointRef& p1, const PointRef& p2)
    {
        return p1.compare(Dimension::Id::Z, p2);
    };

    std::sort(v->begin(), v->end(), cmp);

    return GridPtr(new Grid(v, m_gridLen, origin));
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
        {
            float x = *xp++;
            float y = *yp++;
            if (x != -9999.f)
                d.doVector({xpos, ypos}, {x, y});
        }
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

    pdal::gdal::registerDrivers();

    // Write to in-memory GTiff, then copy to COG so the output is cloud-optimized.
    std::string memPath = "/vsimem/atlas_temp.tif";
    {
        pdal::gdal::Raster raster(memPath, "GTiff", "EPSG:32624", pixelToPos);
        pdal::gdal::GDALError err = raster.open(xsize, ysize, m_ncc ? 16 : 13,
            Dimension::Type::Float, -9999, pdal::StringList());

        if (err != pdal::gdal::GDALError::None)
            throwError(raster.errorMsg());

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
        raster.writeBand(m_field.madXdata(),         -9999.0, 11, "MAD_X");
        raster.writeBand(m_field.madYdata(),         -9999.0, 12, "MAD_Y");
        raster.writeBand(m_field.madZdata(),         -9999.0, 13, "MAD_Z");
        if (m_ncc)
        {
            raster.writeBand(m_field.nccXdata(),     -9999.0, 14, "NCC_X");
            raster.writeBand(m_field.nccYdata(),     -9999.0, 15, "NCC_Y");
            raster.writeBand(m_field.nccPeakData(),  -9999.0, 16, "NCC_PEAK");
        }
    }

    GDALDataset* srcDs = static_cast<GDALDataset*>(GDALOpen(memPath.c_str(), GA_ReadOnly));
    GDALDriver* cogDriver = GetGDALDriverManager()->GetDriverByName("COG");
    const char* cogOptions[] = { "COMPRESS=DEFLATE", nullptr };
    GDALDataset* cogDs = cogDriver->CreateCopy(
        filename.c_str(), srcDs, FALSE, const_cast<char**>(cogOptions), nullptr, nullptr);
    GDALClose(cogDs);
    GDALClose(srcDs);
    VSIUnlink(memPath.c_str());
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
            double minx = rec.origin.x + gi.x() * m_gridLen;
            double miny = rec.origin.y + gi.y() * m_gridLen;
            double maxx = minx + m_gridLen;
            double maxy = miny + m_gridLen;
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
