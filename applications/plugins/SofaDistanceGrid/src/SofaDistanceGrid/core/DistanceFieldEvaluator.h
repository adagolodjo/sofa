#pragma once

#include <SofaDistanceGrid/core/DistanceFieldData.h>

namespace sofa::component::container::_distancegrid_ {

/**
 * @brief Strategy Interface for SDF Evaluation (Interpolation)
 */
class SOFA_SOFADISTANCEGRID_API AbstractDistanceEvaluator 
{
public:
    using Coord = AbstractDistanceField::Coord;
    using SReal = Coord::value_type;

    virtual ~AbstractDistanceEvaluator() = default;

    virtual SReal evaluateDistance(const AbstractDistanceField& field, const Coord& p) const = 0;
    virtual Coord evaluateGradient(const AbstractDistanceField& field, const Coord& p) const = 0;
};

/**
 * @brief Standard Trilinear Interpolation (C0 continuity)
 */
class SOFA_SOFADISTANCEGRID_API TrilinearEvaluator : public AbstractDistanceEvaluator
{
public:
    SReal evaluateDistance(const AbstractDistanceField& field, const Coord& p) const override;
    Coord evaluateGradient(const AbstractDistanceField& field, const Coord& p) const override;
};

/**
 * @brief Advanced Tricubic Interpolation (C1 continuity)
 */
class SOFA_SOFADISTANCEGRID_API TricubicEvaluator : public AbstractDistanceEvaluator
{
public:
    SReal evaluateDistance(const AbstractDistanceField& field, const Coord& p) const override;
    Coord evaluateGradient(const AbstractDistanceField& field, const Coord& p) const override;
};

} // namespace sofa::component::container::_distancegrid_
