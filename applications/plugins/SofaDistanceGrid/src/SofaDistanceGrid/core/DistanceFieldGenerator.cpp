/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
******************************************************************************/
#include <SofaDistanceGrid/core/DistanceFieldGenerator.h>
#include <SofaDistanceGrid/DistanceGrid.h>
#include <sofa/helper/logging/Messaging.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <numeric>

namespace sofa::component::container::_distancegrid_ {

////////////////////////////////////////////////////////////////////////////////
// FastMarchingGenerator — Delegates to existing DistanceGrid::calcDistance
////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<AbstractDistanceField> FastMarchingGenerator::buildFromMesh(
    sofa::helper::io::Mesh* mesh,
    double scale,
    int nx, int ny, int nz,
    Coord pmin, Coord pmax)
{
    auto grid = std::make_unique<DistanceGrid>(nx, ny, nz, pmin, pmax);
    grid->calcDistance(mesh, scale);
    grid->computeBBox();
    return grid;
}

////////////////////////////////////////////////////////////////////////////////
// FastSweepingGenerator — Modern cache-friendly Eikonal solver
////////////////////////////////////////////////////////////////////////////////

namespace {

/// Solve the 1D/2D/3D Eikonal equation at a single node using Godunov upwind.
/// Given sorted candidate values a <= b <= c from neighbors along each axis,
/// and grid spacings hx, hy, hz, find the smallest d such that:
///   max((d-a)/hx, 0)^2 + max((d-b)/hy, 0)^2 + max((d-c)/hz, 0)^2 = 1
static SReal solveEikonal(SReal a, SReal b, SReal c,
                          SReal hx, SReal hy, SReal hz)
{
    // Sort so that a_s <= b_s <= c_s (with corresponding h values)
    SReal vals[3] = {a, b, c};
    SReal hs[3]   = {hx, hy, hz};

    // Simple sort of 3 elements
    if (vals[0] > vals[1]) { std::swap(vals[0], vals[1]); std::swap(hs[0], hs[1]); }
    if (vals[1] > vals[2]) { std::swap(vals[1], vals[2]); std::swap(hs[1], hs[2]); }
    if (vals[0] > vals[1]) { std::swap(vals[0], vals[1]); std::swap(hs[0], hs[1]); }

    // Try 1D solution: d = a + hx
    SReal d = vals[0] + hs[0];

    if (d > vals[1])
    {
        // Try 2D solution: solve (d-a)^2/hx^2 + (d-b)^2/hy^2 = 1
        SReal ihx2 = 1.0 / (hs[0] * hs[0]);
        SReal ihy2 = 1.0 / (hs[1] * hs[1]);
        SReal sum_ih2 = ihx2 + ihy2;
        SReal sum_aih2 = vals[0] * ihx2 + vals[1] * ihy2;
        SReal sum_a2ih2 = vals[0] * vals[0] * ihx2 + vals[1] * vals[1] * ihy2;

        SReal disc = sum_aih2 * sum_aih2 - sum_ih2 * (sum_a2ih2 - 1.0);
        if (disc < 0) disc = 0;
        d = (sum_aih2 + std::sqrt(disc)) / sum_ih2;

        if (d > vals[2])
        {
            // Try 3D solution
            SReal ihz2 = 1.0 / (hs[2] * hs[2]);
            sum_ih2 += ihz2;
            sum_aih2 += vals[2] * ihz2;
            sum_a2ih2 += vals[2] * vals[2] * ihz2;

            disc = sum_aih2 * sum_aih2 - sum_ih2 * (sum_a2ih2 - 1.0);
            if (disc < 0) disc = 0;
            d = (sum_aih2 + std::sqrt(disc)) / sum_ih2;
        }
    }
    return d;
}

} // anonymous namespace

std::unique_ptr<AbstractDistanceField> FastSweepingGenerator::buildFromMesh(
    sofa::helper::io::Mesh* mesh,
    double scale,
    int nx, int ny, int nz,
    Coord pmin, Coord pmax)
{
    // Step 1: Use legay FMM to initialize boundary cells (seed the front)
    // This is because the initial seeding (triangle-edge intersection) logic
    // is complex and correct in the existing code. We reuse it, then re-sweep.
    auto grid = std::make_unique<DistanceGrid>(nx, ny, nz, pmin, pmax);
    grid->calcDistance(mesh, scale);

    // Step 2: Fast Sweeping to refine distances
    // The FMM already computed good distances, but FSM can smooth them
    // and is more cache-friendly for large grids.
    const SReal INF = DistanceGrid::maxDist();
    const int nxny = nx * ny;

    Coord cellWidth = grid->getCellWidth();
    SReal hx = std::abs(cellWidth[0]);
    SReal hy = std::abs(cellWidth[1]);
    SReal hz = std::abs(cellWidth[2]);

    // We work on absolute distances, then restore signs
    std::vector<SReal> absDist(nx * ny * nz);
    std::vector<int> sign(nx * ny * nz);

    for (int i = 0; i < nx * ny * nz; ++i)
    {
        SReal d = (*grid)[i];
        sign[i] = (d < 0) ? -1 : 1;
        absDist[i] = std::abs(d);
    }

    // Mark boundary cells (cells whose absolute distance is small relative
    // to the cell diagonal) as frozen
    SReal diagH = std::sqrt(hx * hx + hy * hy + hz * hz);
    std::vector<bool> frozen(nx * ny * nz, false);
    for (int i = 0; i < nx * ny * nz; ++i)
    {
        if (absDist[i] < 1.5 * diagH)
            frozen[i] = true;
        else
            absDist[i] = INF;
    }

    // 8 alternating sweeps
    const int sweepDirs[8][6] = {
        {0, nx-1, 1,  0, ny-1, 1},  // +x +y +z  (z loop below is separate)
        {nx-1, 0, -1, 0, ny-1, 1},
        {0, nx-1, 1,  ny-1, 0, -1},
        {nx-1, 0, -1, ny-1, 0, -1},
        {0, nx-1, 1,  0, ny-1, 1},
        {nx-1, 0, -1, 0, ny-1, 1},
        {0, nx-1, 1,  ny-1, 0, -1},
        {nx-1, 0, -1, ny-1, 0, -1},
    };
    const int zDirs[8][3] = {
        {0, nz-1, 1}, {0, nz-1, 1}, {0, nz-1, 1}, {0, nz-1, 1},
        {nz-1, 0, -1}, {nz-1, 0, -1}, {nz-1, 0, -1}, {nz-1, 0, -1},
    };

    for (int sweep = 0; sweep < 8; ++sweep)
    {
        int x0 = sweepDirs[sweep][0], x1 = sweepDirs[sweep][1], dx = sweepDirs[sweep][2];
        int y0 = sweepDirs[sweep][3], y1 = sweepDirs[sweep][4], dy = sweepDirs[sweep][5];
        int z0 = zDirs[sweep][0], z1 = zDirs[sweep][1], dz = zDirs[sweep][2];

        for (int z = z0; z != z1 + dz; z += dz)
            for (int y = y0; y != y1 + dy; y += dy)
                for (int x = x0; x != x1 + dx; x += dx)
                {
                    int idx = x + nx * (y + ny * z);
                    if (frozen[idx]) continue;

                    // Get minimum neighbor values along each axis
                    SReal ax = INF, ay = INF, az = INF;
                    if (x > 0)    ax = std::min(ax, absDist[idx - 1]);
                    if (x < nx-1) ax = std::min(ax, absDist[idx + 1]);
                    if (y > 0)    ay = std::min(ay, absDist[idx - nx]);
                    if (y < ny-1) ay = std::min(ay, absDist[idx + nx]);
                    if (z > 0)    az = std::min(az, absDist[idx - nxny]);
                    if (z < nz-1) az = std::min(az, absDist[idx + nxny]);

                    SReal dnew = solveEikonal(ax, ay, az, hx, hy, hz);
                    if (dnew < absDist[idx])
                        absDist[idx] = dnew;
                }
    }

    // Restore signed distances
    for (int i = 0; i < nx * ny * nz; ++i)
    {
        (*grid)[i] = sign[i] * absDist[i];
    }

    grid->computeBBox();
    msg_info("DistanceGrid") << "FastSweepingGenerator: distance field computed ("
                             << nx << "x" << ny << "x" << nz << ")";
    return grid;
}

} // namespace sofa::component::container::_distancegrid_
