#pragma once

#include <SofaDistanceGrid/config.h>
#include <SofaDistanceGrid/DistanceGrid.h>
#include <SofaDistanceGrid/core/DistanceFieldEvaluator.h>
#include <SofaDistanceGrid/core/DistanceFieldGenerator.h>

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/core/objectmodel/Data.h>

namespace sofa::component::container::_distancegrid_ {

/**
 * @brief SOFA Engine that generates and exposes a Distance Grid.
 *
 * This component replaces the hidden loadShared() cache mechanism
 * with a proper SOFA scene graph element. Other components can link
 * to this engine's output data via standard SOFA data bindings:
 *
 *   <DistanceGridEngine name="mySDF" filename="mesh.obj" />
 *   <DistanceGridForceField distanceGrid="@mySDF" />
 */
class SOFA_SOFADISTANCEGRID_API DistanceGridEngine : public core::DataEngine
{
public:
    SOFA_CLASS(DistanceGridEngine, core::DataEngine);

    DistanceGridEngine();
    ~DistanceGridEngine() override = default;

    void init() override;
    void doUpdate() override;

    /// Access the computed distance grid
    DistanceGrid* getGrid() { return m_grid.get(); }
    const DistanceGrid* getGrid() const { return m_grid.get(); }

    /// Access the evaluator
    const AbstractDistanceEvaluator* getEvaluator() const { return m_evaluator.get(); }

    // ---- Data fields ----
    sofa::core::objectmodel::DataFileName d_filename;
    Data<double> d_scale;
    Data<double> d_sampling;
    Data<int> d_nx;
    Data<int> d_ny;
    Data<int> d_nz;
    Data<type::Vec3> d_pmin;
    Data<type::Vec3> d_pmax;
    Data<std::string> d_generatorType; ///< "fmm" or "fsm"
    Data<std::string> d_evaluatorType; ///< "trilinear" or "tricubic"

private:
    std::shared_ptr<DistanceGrid> m_grid;
    std::unique_ptr<AbstractDistanceEvaluator> m_evaluator;
};

} // namespace sofa::component::container::_distancegrid_

namespace sofa::component::container {
    using _distancegrid_::DistanceGridEngine;
} // namespace sofa::component::container
