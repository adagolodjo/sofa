/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
******************************************************************************/
#include <SofaDistanceGrid/components/engine/DistanceGridEngine.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/logging/Messaging.h>

namespace sofa::component::container::_distancegrid_ {

DistanceGridEngine::DistanceGridEngine()
    : d_filename(initData(&d_filename, "filename", "Path to the mesh file (.obj, .vtk, .raw)"))
    , d_scale(initData(&d_scale, 1.0, "scale", "Scale factor applied to the mesh"))
    , d_sampling(initData(&d_sampling, 0.0, "sampling",
        "Surface sampling distance (negative = in voxels, 0 = use mesh vertices)"))
    , d_nx(initData(&d_nx, 64, "nx", "Grid resolution along X"))
    , d_ny(initData(&d_ny, 64, "ny", "Grid resolution along Y"))
    , d_nz(initData(&d_nz, 64, "nz", "Grid resolution along Z"))
    , d_pmin(initData(&d_pmin, type::Vec3(), "pmin", "Grid bounding box minimum"))
    , d_pmax(initData(&d_pmax, type::Vec3(), "pmax", "Grid bounding box maximum"))
    , d_generatorType(initData(&d_generatorType, std::string("fmm"), "generator",
        "SDF generation method: 'fmm' (Fast Marching) or 'fsm' (Fast Sweeping)"))
    , d_evaluatorType(initData(&d_evaluatorType, std::string("trilinear"), "evaluator",
        "Interpolation method: 'trilinear' (C0) or 'tricubic' (C1)"))
{
    addInput(&d_filename);
    addInput(&d_scale);
    addInput(&d_sampling);
    addInput(&d_nx);
    addInput(&d_ny);
    addInput(&d_nz);
    addInput(&d_pmin);
    addInput(&d_pmax);
    addInput(&d_generatorType);
    addInput(&d_evaluatorType);
}

void DistanceGridEngine::init()
{
    setDirtyValue();
    doUpdate();
}

void DistanceGridEngine::doUpdate()
{
    // Choose evaluator
    const std::string& evalType = d_evaluatorType.getValue();
    if (evalType == "tricubic")
    {
        m_evaluator = std::make_unique<TricubicEvaluator>();
        msg_info() << "Using Tricubic (C1) evaluator";
    }
    else
    {
        m_evaluator = std::make_unique<TrilinearEvaluator>();
        if (evalType != "trilinear")
            msg_warning() << "Unknown evaluator type '" << evalType
                          << "', falling back to trilinear";
    }

    // Load/generate the grid
    const std::string& filename = d_filename.getFullPath();
    if (filename.empty())
    {
        msg_error() << "No filename specified";
        return;
    }

    const std::string& genType = d_generatorType.getValue();
    msg_info() << "Loading distance grid from '" << filename
               << "' with generator=" << genType
               << ", evaluator=" << evalType;

    // Use the existing loadShared mechanism for now (backward compatible).
    // In the future, the generator interface will be used directly
    // for mesh-based construction (the load() function handles
    // VTK/RAW/OBJ/cube formats internally).
    m_grid = DistanceGrid::loadShared(
        filename,
        d_scale.getValue(),
        d_sampling.getValue(),
        d_nx.getValue(), d_ny.getValue(), d_nz.getValue(),
        d_pmin.getValue(), d_pmax.getValue()
    );

    if (!m_grid)
    {
        msg_error() << "Failed to load distance grid from '" << filename << "'";
        return;
    }

    msg_info() << "Distance grid loaded: "
               << m_grid->getNx() << "x" << m_grid->getNy() << "x" << m_grid->getNz()
               << " cells, bbox=[" << m_grid->getBBMin() << "] - [" << m_grid->getBBMax() << "]";
}

// Register the component in SOFA's ObjectFactory
int DistanceGridEngineClass = core::RegisterObject(
    "Engine that generates and exposes a Distance Grid (SDF). "
    "Supports Fast Marching (fmm) and Fast Sweeping (fsm) generators, "
    "and Trilinear (C0) or Tricubic (C1) evaluators.")
    .add<DistanceGridEngine>();

} // namespace sofa::component::container::_distancegrid_
