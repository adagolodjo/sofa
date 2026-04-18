#pragma once

#include <SofaDistanceGrid/core/DistanceFieldData.h>
#include <sofa/helper/io/Mesh.h>
#include <memory>
#include <string>

namespace sofa::component::container::_distancegrid_ {

/**
 * @brief Strategy Interface for SDF Generation / Initialization
 */
class SOFA_SOFADISTANCEGRID_API AbstractDistanceGenerator 
{
public:
    using Coord = AbstractDistanceField::Coord;

    virtual ~AbstractDistanceGenerator() = default;

    /**
     * @brief Build a SDF from an explicit surface mesh.
     */
    virtual std::unique_ptr<AbstractDistanceField> buildFromMesh(
        sofa::helper::io::Mesh* mesh, 
        double scale, 
        int nx, int ny, int nz, 
        Coord pmin, Coord pmax) = 0;
};

/**
 * @brief Legacy Fast Marching Method (FMM) Generator.
 */
class SOFA_SOFADISTANCEGRID_API FastMarchingGenerator : public AbstractDistanceGenerator
{
public:
    std::unique_ptr<AbstractDistanceField> buildFromMesh(
        sofa::helper::io::Mesh* mesh, 
        double scale, 
        int nx, int ny, int nz, 
        Coord pmin, Coord pmax) override;
};

/**
 * @brief Modern Fast Sweeping Method (FSM) Generator - Cache friendly.
 */
class SOFA_SOFADISTANCEGRID_API FastSweepingGenerator : public AbstractDistanceGenerator
{
public:
    std::unique_ptr<AbstractDistanceField> buildFromMesh(
        sofa::helper::io::Mesh* mesh, 
        double scale, 
        int nx, int ny, int nz, 
        Coord pmin, Coord pmax) override;
};

} // namespace sofa::component::container::_distancegrid_
