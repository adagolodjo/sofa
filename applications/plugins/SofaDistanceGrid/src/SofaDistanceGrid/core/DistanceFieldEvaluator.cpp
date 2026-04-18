/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
******************************************************************************/
#include <SofaDistanceGrid/core/DistanceFieldEvaluator.h>
#include <SofaDistanceGrid/DistanceGrid.h>
#include <cmath>

namespace sofa::component::container::_distancegrid_ {

////////////////////////////////////////////////////////////////////////////////
// TrilinearEvaluator — Direct delegation to DistanceGrid::interp / grad
////////////////////////////////////////////////////////////////////////////////

SReal TrilinearEvaluator::evaluateDistance(const AbstractDistanceField& field, const Coord& p) const
{
    // Dynamic cast to the concrete type we know how to evaluate
    const auto* grid = dynamic_cast<const DistanceGrid*>(&field);
    if (!grid) return 0;
    return grid->eval(p);
}

Coord TrilinearEvaluator::evaluateGradient(const AbstractDistanceField& field, const Coord& p) const
{
    const auto* grid = dynamic_cast<const DistanceGrid*>(&field);
    if (!grid) return Coord();
    return grid->grad(p);
}

////////////////////////////////////////////////////////////////////////////////
// TricubicEvaluator — Catmull-Rom spline interpolation (C1 continuous)
////////////////////////////////////////////////////////////////////////////////

namespace {

/// Catmull-Rom basis functions and their derivatives.
/// For parameter t in [0,1]:
///   h0(t) = -0.5t^3 + t^2 - 0.5t        (weight for P_{i-1})
///   h1(t) =  1.5t^3 - 2.5t^2 + 1         (weight for P_i)
///   h2(t) = -1.5t^3 + 2t^2 + 0.5t        (weight for P_{i+1})
///   h3(t) =  0.5t^3 - 0.5t^2             (weight for P_{i+2})
struct CatmullRom
{
    SReal w[4];   // basis weights
    SReal dw[4];  // derivative weights

    void compute(SReal t)
    {
        SReal t2 = t * t;
        SReal t3 = t2 * t;

        w[0] = -0.5*t3 +     t2 - 0.5*t;
        w[1] =  1.5*t3 - 2.5*t2 + 1.0;
        w[2] = -1.5*t3 + 2.0*t2 + 0.5*t;
        w[3] =  0.5*t3 - 0.5*t2;

        dw[0] = -1.5*t2 + 2.0*t - 0.5;
        dw[1] =  4.5*t2 - 5.0*t;
        dw[2] = -4.5*t2 + 4.0*t + 0.5;
        dw[3] =  1.5*t2 - 1.0*t;
    }
};

/// Safely get distance value with clamped indices
inline SReal getDistClamped(const DistanceGrid* grid, int x, int y, int z)
{
    int nx = grid->getNx();
    int ny = grid->getNy();
    int nz = grid->getNz();
    x = std::max(0, std::min(x, nx - 1));
    y = std::max(0, std::min(y, ny - 1));
    z = std::max(0, std::min(z, nz - 1));
    return (*grid)[grid->index(x, y, z)];
}

} // anonymous namespace

SReal TricubicEvaluator::evaluateDistance(const AbstractDistanceField& field, const Coord& p) const
{
    const auto* grid = dynamic_cast<const DistanceGrid*>(&field);
    if (!grid) return 0;
    if (!grid->inGrid(p)) return grid->eval(p); // fallback for out-of-grid

    Coord cellW = grid->getCellWidth();
    Coord pmin = grid->getPMin();
    int nx = grid->getNx(), ny = grid->getNy();

    // Continuous grid coordinates
    SReal gx = (p[0] - pmin[0]) / cellW[0];
    SReal gy = (p[1] - pmin[1]) / cellW[1];
    SReal gz = (p[2] - pmin[2]) / cellW[2];

    int ix = (int)std::floor(gx);
    int iy = (int)std::floor(gy);
    int iz = (int)std::floor(gz);

    SReal tx = gx - ix;
    SReal ty = gy - iy;
    SReal tz = gz - iz;

    CatmullRom crx, cry, crz;
    crx.compute(tx);
    cry.compute(ty);
    crz.compute(tz);

    SReal result = 0;
    for (int dz = -1; dz <= 2; ++dz)
    {
        SReal valYZ = 0;
        for (int dy = -1; dy <= 2; ++dy)
        {
            SReal valX = 0;
            for (int dx = -1; dx <= 2; ++dx)
            {
                valX += crx.w[dx + 1] * getDistClamped(grid, ix + dx, iy + dy, iz + dz);
            }
            valYZ += cry.w[dy + 1] * valX;
        }
        result += crz.w[dz + 1] * valYZ;
    }
    return result;
}

Coord TricubicEvaluator::evaluateGradient(const AbstractDistanceField& field, const Coord& p) const
{
    const auto* grid = dynamic_cast<const DistanceGrid*>(&field);
    if (!grid) return Coord();
    if (!grid->inGrid(p)) return grid->grad(p); // fallback for out-of-grid

    Coord cellW = grid->getCellWidth();
    Coord invCellW(1.0 / cellW[0], 1.0 / cellW[1], 1.0 / cellW[2]);
    Coord pmin = grid->getPMin();

    SReal gx = (p[0] - pmin[0]) / cellW[0];
    SReal gy = (p[1] - pmin[1]) / cellW[1];
    SReal gz = (p[2] - pmin[2]) / cellW[2];

    int ix = (int)std::floor(gx);
    int iy = (int)std::floor(gy);
    int iz = (int)std::floor(gz);

    SReal tx = gx - ix;
    SReal ty = gy - iy;
    SReal tz = gz - iz;

    CatmullRom crx, cry, crz;
    crx.compute(tx);
    cry.compute(ty);
    crz.compute(tz);

    // Compute gradient by differentiating the Catmull-Rom in each axis
    Coord result(0, 0, 0);

    for (int dz = -1; dz <= 2; ++dz)
    {
        for (int dy = -1; dy <= 2; ++dy)
        {
            for (int dx = -1; dx <= 2; ++dx)
            {
                SReal val = getDistClamped(grid, ix + dx, iy + dy, iz + dz);

                // d/dx: differentiate x basis, keep y and z normal
                result[0] += crx.dw[dx + 1] * cry.w[dy + 1] * crz.w[dz + 1] * val;
                // d/dy: differentiate y basis
                result[1] += crx.w[dx + 1] * cry.dw[dy + 1] * crz.w[dz + 1] * val;
                // d/dz: differentiate z basis
                result[2] += crx.w[dx + 1] * cry.w[dy + 1] * crz.dw[dz + 1] * val;
            }
        }
    }

    // The derivatives are in grid-space, convert to world-space
    result[0] *= invCellW[0];
    result[1] *= invCellW[1];
    result[2] *= invCellW[2];

    return result;
}

} // namespace sofa::component::container::_distancegrid_
