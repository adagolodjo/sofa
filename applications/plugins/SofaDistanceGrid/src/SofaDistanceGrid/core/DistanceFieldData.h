#pragma once

#include <SofaDistanceGrid/config.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa::component::container::_distancegrid_ {

using sofa::type::Vec3;
typedef Vec3 Coord;

/**
 * @brief Abstract base interface for any Distance Field.
 *
 * This interface decouples the evaluation of a signed distance field
 * from its concrete storage (dense grid, sparse grid, neural SDF, etc.).
 * All existing code that depends on DistanceGrid can progressively
 * migrate to this interface.
 */
class SOFA_SOFADISTANCEGRID_API AbstractDistanceField
{
public:
    typedef type::Vec3 Coord;

    virtual ~AbstractDistanceField() = default;

    /// Evaluate the signed distance at point x (with clamping for out-of-grid points)
    virtual SReal eval(const Coord& x) const = 0;

    /// Compute the gradient of the distance field at point p
    virtual Coord grad(const Coord& p) const = 0;

    /// Check if point p is inside the grid domain
    virtual bool inGrid(const Coord& p) const = 0;

    /// Get the bounding box of the object (tighter than the grid bounds)
    virtual const Coord& getBBMin() const = 0;
    virtual const Coord& getBBMax() const = 0;

    /// Get the grid domain bounds
    virtual const Coord& getPMin() const = 0;
    virtual const Coord& getPMax() const = 0;
};

} // namespace sofa::component::container::_distancegrid_
