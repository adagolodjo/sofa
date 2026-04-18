/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_INL
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_INL

#include <sofa/core/visual/VisualParams.h>
#include "DistanceGridForceField.h"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <cassert>
#include <iostream>

#include <sofa/core/behavior/MultiMatrixAccessor.h>


namespace sofa
{

namespace component
{

namespace forcefield
{


template<class DataTypes>
void DistanceGridForceField<DataTypes>::init()
{
    Inherit::init();


    if (fileDistanceGrid.getValue().empty())
    {
        if (grid == nullptr)
            msg_error() << "DistanceGridForceField requires an input filename." ;
        /// the grid has already been set
        return;
    }
    msg_info() << " creating "<<nx.getValue()<<"x"<<ny.getValue()<<"x"<<nz.getValue()<< msgendl
               << " DistanceGrid from file '"<< fileDistanceGrid.getValue() <<"'.";

    msg_info_when(scale.getValue()!=1.0) << " scale="<<scale.getValue();

    msg_info_when(box.getValue()[0][0]<box.getValue()[1][0])
            <<" bbox=<"<<box.getValue()[0]<<">-<"<<box.getValue()[1]<<">";

    grid = DistanceGrid::loadShared(fileDistanceGrid.getFullPath(), scale.getValue(), 0.0,
                                    nx.getValue(),ny.getValue(),nz.getValue(),
                                    box.getValue()[0],box.getValue()[1]);

    if (grid == nullptr)
    {
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        msg_error() << "Failed to initialize: Invalid distance grid";
        return;
    }

    if (this->stiffnessArea.getValue() != 0 && this->mstate)
    {
        core::topology::BaseMeshTopology* topology = this->getContext()->getMeshTopology();
        if (topology && topology->getNbTriangles() > 0)
        {
            const core::topology::BaseMeshTopology::SeqTriangles& triangles = topology->getTriangles();
            Real sumArea = 0;
            Real sumSArea = 0;
            const VecCoord& p1 = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
            pOnBorder.resize(p1.size(), false);
            for (unsigned int ti = 0; ti < triangles.size(); ++ti)
            {
                const auto& t = triangles[ti];
                Coord B = p1[t[1]]-p1[t[0]];
                Coord C = p1[t[2]]-p1[t[0]];
                Coord tN = cross(B, C);
                Real area = tN.norm()/2;
                Coord sN = grid->grad((p1[t[0]]+p1[t[1]]+p1[t[2]])*(1.0/3.0));
                sN.normalize();
                sumSArea += (sN*tN)*0.5f;
                sumArea += area;
                pOnBorder[t[0]] = true;
                pOnBorder[t[1]] = true;
                pOnBorder[t[2]] = true;
            }
            msg_info() << "Surface area : " << sumArea << " ( mean " << sumArea / triangles.size() << " per triangle )" ;
            flipNormals = (sumSArea < 0);
        }
        else
        {
            msg_error() << "No triangles found in topology";
        }
    }

    if (this->stiffnessVolume.getValue() != 0 && this->mstate)
    {
        core::topology::BaseMeshTopology* topology = this->mstate->getContext()->getMeshTopology();
        if (topology && topology->getNbTetrahedra() > 0)
        {
            const core::topology::BaseMeshTopology::SeqTetrahedra& tetrahedra = topology->getTetrahedra();
            Real sumVolume = 0;
            const VecCoord& p1 = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
            for (unsigned int ti = 0; ti < tetrahedra.size(); ++ti)
            {
                const auto & t = tetrahedra[ti];
                Coord A = p1[t[1]]-p1[t[0]];
                Coord B = p1[t[2]]-p1[t[0]];
                Coord C = p1[t[3]]-p1[t[0]];
                Real volume = (A*cross(B, C))/6.0f;
                sumVolume += volume;
            }
            msg_info() << "Volume : " << sumVolume << " ( mean " << sumVolume / tetrahedra.size() << " per tetra )" ;
        }
        else
        {
            msg_error() << "No tetrahedra found in topology";
        }
    }

    sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void DistanceGridForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV )
{
    VecDeriv& f1 = *(dataF.beginEdit());
    const VecCoord& p1=dataX.getValue();
    const VecDeriv& v1=dataV.getValue();


    if (!grid) return;
    //this->dfdd.resize(p1.size());
    f1.resize(p1.size());

    sofa::type::vector<Contact>& contacts = *this->contacts.beginEdit();
    contacts.clear();

    unsigned int ibegin = 0;
    unsigned int iend = p1.size();

    if (localRange.getValue()[0] >= 0)
        ibegin = localRange.getValue()[0];

    if (localRange.getValue()[1] >= 0 && (unsigned int)localRange.getValue()[1]+1 < iend)
        iend = localRange.getValue()[1]+1;

    const Real stiffIn = stiffnessIn.getValue();
    const Real stiffOut = stiffnessOut.getValue();
    const Real damp = damping.getValue();
    const Real maxdist = maxDist.getValue();
    unsigned int nbIn = 0;
    for (unsigned int i=ibegin; i<iend; i++)
    {
        if (i < pOnBorder.size() && !pOnBorder[i]) continue;
        Real d = grid->teval(p1[i]);
        if(d > 0)
        {
            if (d >= maxdist || stiffOut == 0) continue;
            Deriv grad = grid->tgrad(p1[i]);
            grad.normalize();
            Real forceIntensity = -stiffOut * d;
            Real dampingIntensity = forceIntensity * damp;
            Deriv force = grad * forceIntensity - v1[i]*dampingIntensity;
            f1[i]+=force;
            Contact c;
            c.index = i;
            c.normal = grad;
            c.fact = forceIntensity;
            contacts.push_back(c);
        }
        else if (d < 0)
        {
            if (-d >= maxdist || stiffIn == 0) continue;
            Deriv grad = grid->tgrad(p1[i]);
            grad.normalize();
            Real forceIntensity = -stiffIn * d;
            Real dampingIntensity = forceIntensity * damp;
            Deriv force = grad * forceIntensity - v1[i]*dampingIntensity;
            f1[i]+=force;
            Contact c;
            c.index = i;
            c.normal = grad;
            c.fact = forceIntensity;
            contacts.push_back(c);
            nbIn++;
        }
    }

    dmsg_info() << " number of points " << nbIn ;

    this->contacts.endEdit();

    sofa::type::vector<TContact>& tcontacts = *this->tcontacts.beginEdit();
    tcontacts.clear();

    const Real stiffA = stiffnessArea.getValue();
    const Real minA = minArea.getValue();
    if (stiffA != 0)
    {
        core::topology::BaseMeshTopology* topology = this->getContext()->getMeshTopology();
        if (topology && topology->getNbTriangles() > 0)
        {
            const core::topology::BaseMeshTopology::SeqTriangles& triangles = topology->getTriangles();
            for (unsigned int ti = 0; ti < triangles.size(); ++ti)
            {
                const auto& t = triangles[ti];
                Coord B = p1[t[1]]-p1[t[0]];
                Coord C = p1[t[2]]-p1[t[0]];
                Coord tN = cross(B, C);
                Coord tcenter = (p1[t[0]]+p1[t[1]]+p1[t[2]])*(1.0/3.0);
                Coord sN = grid->tgrad(tcenter);
                if (flipNormals) sN = -sN;
                sN.normalize();
                Real area = (tN * sN) * 0.5f;

                if (area < minA)
                {

                    // lets consider A = (0,0,0)
                    // area = 0.5*(sNx*tNx + sNy * tNy + sNz * tNz)
                    //      = 0.5*(sNx*(By*Cz-Bz*Cy) + sNy * (Bz*Cx-Bx*Cz) + sNz * (Bx*Cy-By*Cx))

                    // d(area) / dBx = 0.5*(sNz*Cy-sNy*Cz)
                    // d(area) / dBy = 0.5*(sNx*Cz-sNz*Cx)
                    // d(area) / dBz = 0.5*(sNy*Cx-sNx*Cy)
                    // d(area) / dB = 0.5*cross(C,sN)
                    // d(area) / dC = 0.5*cross(sN,B)

                    Real forceIntensity = (stiffA * (minA-area));
                    Coord fB = cross(C,sN)*forceIntensity; f1[t[1]] += fB;
                    Coord fC = cross(sN,B)*forceIntensity; f1[t[2]] += fC;
                    Coord fA = -(fB+fC);                   f1[t[0]] += fA;

                    TContact c;
                    c.index = t.array();
                    c.fact = minA-area;
                    c.normal = sN;
                    c.B = B;
                    c.C = C;
                    tcontacts.push_back(c);
                }
            }
        }
    }
    this->tcontacts.endEdit();

    sofa::type::vector<VContact>& vcontacts = *this->vcontacts.beginEdit();
    vcontacts.clear();

    const Real stiffV = stiffnessVolume.getValue();
    const Real minV = minVolume.getValue();
    if (stiffV != 0)
    {
        core::topology::BaseMeshTopology* topology = this->mstate->getContext()->getMeshTopology();
        if (topology && topology->getNbTetrahedra() > 0)
        {
            const core::topology::BaseMeshTopology::SeqTetrahedra& tetrahedra = topology->getTetrahedra();
            const Real v1_6 = (Real)(1.0/6.0);
            for (unsigned int ti = 0; ti < tetrahedra.size(); ++ti)
            {
                const auto& t = tetrahedra[ti];
                Coord A = p1[t[1]]-p1[t[0]];
                Coord B = p1[t[2]]-p1[t[0]];
                Coord C = p1[t[3]]-p1[t[0]];
                Real volume = (A*cross(B, C))*v1_6;

                if (volume < minV)
                {
                    // vol = 1/6*(A(BxC))
                    //      = 1/6*(Ax*(By*Cz-Bz*Cy) + Ay * (Bz*Cx-Bx*Cz) + Az * (Bx*Cy-By*Cx))

                    // d(vol) / dBx = 1/6*(Az*Cy-Ay*Cz)
                    // d(vol) / dBy = 1/6*(Ax*Cz-Az*Cx)
                    // d(vol) / dBz = 1/6*(Ay*Cx-Ax*Cy)
                    // d(vol) / dA = 1/6*cross(B,C)
                    // d(vol) / dB = 1/6*cross(C,A)
                    // d(vol) / dC = 1/6*cross(A,B)

                    Real forceIntensity = v1_6*(stiffV * (minV-volume));
                    Coord fA = cross(B,C)*forceIntensity; f1[t[1]] += fA;
                    Coord fB = cross(C,A)*forceIntensity; f1[t[2]] += fB;
                    Coord fC = cross(A,B)*forceIntensity; f1[t[3]] += fC;
                    Coord f0 = -(fA+fB+fC);               f1[t[0]] += f0;

                    VContact c;
                    c.index = t.array();
                    c.fact = minV-volume;
                    c.A = A;
                    c.B = B;
                    c.C = C;
                    vcontacts.push_back(c);
                }
            }
        }
    }
    this->vcontacts.endEdit();

    dataF.endEdit();

}

template<class DataTypes>
void DistanceGridForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df1      = *(datadF.beginEdit());
    const VecCoord& dx1=   datadX.getValue()  ;
    Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());

    if (!grid)
        return;

    const sofa::type::vector<Contact>& contacts = this->contacts.getValue();
    const sofa::type::vector<TContact>& tcontacts = this->tcontacts.getValue();
    const sofa::type::vector<VContact>& vcontacts = this->vcontacts.getValue();

    if (contacts.empty() && tcontacts.empty() && vcontacts.empty())
        return;

    df1.resize(dx1.size());
    const Real fact = (Real)(kFactor);

    for (unsigned int i=0; i<contacts.size(); i++)
    {
        const Contact& c = (this->contacts.getValue())[i];
        Coord du = dx1[c.index];

        Real dd = du * c.normal;
        Deriv dforce = c.normal * (dd * c.fact * fact);
        df1[c.index] += dforce;
    }

    const Real factA = (Real)( -this->stiffnessArea.getValue()* kFactor );
    for (unsigned int i=0; i<tcontacts.size(); i++)
    {
        const TContact& c = (this->tcontacts.getValue())[i];
        const type::fixed_array<unsigned int,3>& t = c.index;
        Coord dB = dx1[t[1]]-dx1[t[0]];
        Coord dC = dx1[t[2]]-dx1[t[0]];

        // d(area) = dot(0.5*cross(C,sN), dB) + dot(0.5*cross(sN,B), dC)
        Real darea = 0.5f*(cross(c.C,c.normal)*dB + cross(c.normal,c.B)*dC);

        // fB = (C x sN) * (stiffA * (minA-area));
        // dfB = (C x sN) * (-stiffA * d(area)) + (dC x sN) * (stiffA * (minA-area));
        Coord dfB = cross(c.C, c.normal) * (factA * darea) - cross (dC, c.normal) * (factA*c.fact);
        Coord dfC = cross(c.normal, c.B) * (factA * darea) - cross (c.normal, dB) * (factA*c.fact);
        Coord dfA = -(dfB+dfC);
        df1[t[0]] += dfA;
        df1[t[1]] += dfB;
        df1[t[2]] += dfC;
    }

    const Real factV = (Real)( - this->stiffnessVolume.getValue()*(1.0/6.0)* kFactor );
    for (unsigned int i=0; i<vcontacts.size(); i++)
    {
        const Real v1_6 = (Real)(1.0/6.0);
        const VContact& c = (this->vcontacts.getValue())[i];
        const type::fixed_array<unsigned int,4>& t = c.index;
        Coord dA = dx1[t[1]]-dx1[t[0]];
        Coord dB = dx1[t[2]]-dx1[t[0]];
        Coord dC = dx1[t[3]]-dx1[t[0]];
        // d(vol) = 1/6*(dot(cross(B,C), dA) + dot(cross(C,A), dB) + dot(cross(A,B), dC))
        Real dvolume = v1_6*(dA*cross(c.B,c.C) + dB*cross(c.C,c.A) + dC*cross(c.A,c.B));

        // fA = (1/6)*(B x C)*(stiffV * (minV-volume))
        // dfA = (stiffV*1/6) * ((dB x C + B x dC) * (minV-volume) - (B x C) * dvol)
        // dfB = (stiffV*1/6) * ((dC x A + C x dA) * (minV-volume) - (C x A) * dvol)
        // dfC = (stiffV*1/6) * ((dA x B + A x dB) * (minV-volume) - (A x B) * dvol)
        Coord dfA = cross(c.B, c.C) * (factV * dvolume) - (cross(dB, c.C) + cross(c.B, dC)) * (factV*c.fact);
        Coord dfB = cross(c.C, c.A) * (factV * dvolume) - (cross(dC, c.A) + cross(c.C, dA)) * (factV*c.fact);
        Coord dfC = cross(c.A, c.B) * (factV * dvolume) - (cross(dA, c.B) + cross(c.A, dB)) * (factV*c.fact);
        Coord df0 = -(dfA+dfB+dfC);
        df1[t[1]] += dfA;
        df1[t[2]] += dfB;
        df1[t[3]] += dfC;
        df1[t[0]] += df0;
    }

    datadF.endEdit();
}

template<class DataTypes>
void DistanceGridForceField<DataTypes>::addKToMatrix(const sofa::core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());
    unsigned int &offset = r.offset;
    sofa::linearalgebra::BaseMatrix* mat = r.matrix;

    if (!r) return;
    if (!grid) return;

    // Helper: add a 3x3 outer-product block  a ⊗ b * factor  at matrix position (pi, pj)
    auto addOuterProduct = [&](int pi, int pj, const Coord& a, const Coord& b, Real factor)
    {
        for (sofa::Size l = 0; l < Deriv::total_size; ++l)
            for (sofa::Size c = 0; c < Deriv::total_size; ++c)
                mat->add(offset + pi * Deriv::total_size + l,
                         offset + pj * Deriv::total_size + c,
                         a[l] * b[c] * factor);
    };

    // Helper: add a 3x3 skew-symmetric block [n]× * factor at matrix position (pi, pj)
    // [n]× · v = cross(n, v), so the matrix is:
    //   [  0   -n2   n1 ]
    //   [  n2   0   -n0 ]
    //   [ -n1   n0   0  ]
    auto addSkewBlock = [&](int pi, int pj, const Coord& n, Real factor)
    {
        mat->add(offset + pi*3 + 0, offset + pj*3 + 1, -n[2] * factor);
        mat->add(offset + pi*3 + 0, offset + pj*3 + 2,  n[1] * factor);
        mat->add(offset + pi*3 + 1, offset + pj*3 + 0,  n[2] * factor);
        mat->add(offset + pi*3 + 1, offset + pj*3 + 2, -n[0] * factor);
        mat->add(offset + pi*3 + 2, offset + pj*3 + 0, -n[1] * factor);
        mat->add(offset + pi*3 + 2, offset + pj*3 + 1,  n[0] * factor);
    };

    // ===== 1) Point contacts (existing code) =====
    const sofa::type::vector<Contact>& contacts = this->contacts.getValue();
    for (unsigned int i = 0; i < contacts.size(); i++)
    {
        const Contact& c = contacts[i];
        const int p = c.index;
        const Deriv& normal = c.normal;
        addOuterProduct(p, p, normal, normal, c.fact * -kFactor);
    }

    // ===== 2) Triangle contacts (area penalty) =====
    // From addDForce:
    //   a = cross(C, N),  b = cross(N, B)
    //   darea = 0.5*(a·dB + b·dC)    where dB=dx1-dx0, dC=dx2-dx0
    //   dfB = a*(factA*darea) + [N]×*dC*(factA*fact)
    //   dfC = b*(factA*darea) - [N]×*dB*(factA*fact)   (cross(N,dB) = -[N]×·dB... checking sign)
    //   dfA = -(dfB+dfC)
    // The derivative blocks K[ti,tj] are assembled from these relationships.
    const sofa::type::vector<TContact>& tcontacts = this->tcontacts.getValue();
    const Real factA = -this->stiffnessArea.getValue() * kFactor;
    for (unsigned int i = 0; i < tcontacts.size(); i++)
    {
        const TContact& tc = tcontacts[i];
        const int t0 = tc.index[0], t1 = tc.index[1], t2 = tc.index[2];
        const Coord& N = tc.normal;
        const Coord& B = tc.B;
        const Coord& C = tc.C;
        const Real cFact = tc.fact;

        Coord a = cross(C, N); // d(area)/d(B) direction * 2
        Coord b = cross(N, B); // d(area)/d(C) direction * 2

        // dfB w.r.t. dB: factA * 0.5 * a ⊗ a     (from darea term)
        // dfB w.r.t. dC: factA * 0.5 * a ⊗ b + factA*cFact * [N]×  (from darea + cross term)
        // Since dB = dx1-dx0 and dC = dx2-dx0, chain rule gives blocks for t0,t1,t2

        // K[t1,t1] += factA * 0.5 * a ⊗ a
        addOuterProduct(t1, t1, a, a, factA * 0.5);

        // K[t1,t0] -= factA * 0.5 * a ⊗ a
        addOuterProduct(t1, t0, a, a, -factA * 0.5);

        // K[t1,t2] += factA * 0.5 * a ⊗ b + factA*cFact * [N]×
        addOuterProduct(t1, t2, a, b, factA * 0.5);
        addSkewBlock(t1, t2, N, factA * cFact);

        // K[t1,t0] -= factA * 0.5 * a ⊗ b + factA*cFact * [N]× (from -dC w.r.t. dx0)
        addOuterProduct(t1, t0, a, b, -factA * 0.5);
        addSkewBlock(t1, t0, N, -factA * cFact);

        // dfC w.r.t. dC: factA * 0.5 * b ⊗ b
        addOuterProduct(t2, t2, b, b, factA * 0.5);
        addOuterProduct(t2, t0, b, b, -factA * 0.5);

        // dfC w.r.t. dB: factA * 0.5 * b ⊗ a - factA*cFact * [N]×
        addOuterProduct(t2, t1, b, a, factA * 0.5);
        addSkewBlock(t2, t1, N, -factA * cFact);

        addOuterProduct(t2, t0, b, a, -factA * 0.5);
        addSkewBlock(t2, t0, N, factA * cFact);

        // dfA = -(dfB+dfC): assemble by adding negated rows to t0
        // K[t0,tj] = -K[t1,tj] - K[t2,tj]  for each tj
        // This is equivalent to: for each block we added for t1 and t2,
        // add the negative to t0.
        // Rather than re-derive, we use the constraint that sum of rows = 0
        // by going through t0 explicitly:
        Coord apb = a + b; // combined derivative factor
        addOuterProduct(t0, t1, apb, a, -factA * 0.5);
        addOuterProduct(t0, t2, apb, b, -factA * 0.5);
        addOuterProduct(t0, t0, apb, a, factA * 0.5);
        addOuterProduct(t0, t0, apb, b, factA * 0.5);
    }

    // ===== 3) Tetrahedron contacts (volume penalty) =====
    // From addDForce:
    //   ga = cross(B,C), gb = cross(C,A), gc = cross(A,B)
    //   dvol = (1/6)*(ga·dA + gb·dB + gc·dC)
    //   dfA = factV * (ga * dvol - (cross(dB,C)+cross(B,dC)) * cFact)
    //   dfB = factV * (gb * dvol - (cross(dC,A)+cross(C,dA)) * cFact)
    //   dfC = factV * (gc * dvol - (cross(dA,B)+cross(A,dB)) * cFact)
    //   df0 = -(dfA+dfB+dfC)
    const sofa::type::vector<VContact>& vcontacts = this->vcontacts.getValue();
    const Real factV = -this->stiffnessVolume.getValue() * (1.0 / 6.0) * kFactor;
    for (unsigned int i = 0; i < vcontacts.size(); i++)
    {
        const VContact& vc = vcontacts[i];
        const int v0 = vc.index[0], v1 = vc.index[1], v2 = vc.index[2], v3 = vc.index[3];
        const Coord& A = vc.A;
        const Coord& B = vc.B;
        const Coord& C = vc.C;
        const Real cFact = vc.fact;

        Coord ga = cross(B, C);
        Coord gb = cross(C, A);
        Coord gc = cross(A, B);

        // Volume-derivative outer products: K_vol[vi,vj] = factV * (1/6) * g_i ⊗ g_j
        // where dA=dx1-dx0, dB=dx2-dx0, dC=dx3-dx0
        const Real v16 = (Real)(1.0 / 6.0);

        // dfA (at v1) w.r.t. dA (from v1): factV * v16 * ga ⊗ ga
        addOuterProduct(v1, v1, ga, ga, factV * v16);
        addOuterProduct(v1, v0, ga, ga, -factV * v16);
        addOuterProduct(v1, v2, ga, gb, factV * v16);
        addOuterProduct(v1, v0, ga, gb, -factV * v16);
        addOuterProduct(v1, v3, ga, gc, factV * v16);
        addOuterProduct(v1, v0, ga, gc, -factV * v16);

        // Cross-product terms for dfA: -factV*cFact*(cross(dB,C)+cross(B,dC))
        // cross(dB,C) = -[C]× dB   and  cross(B,dC) = [B]× dC (skew matrices)
        addSkewBlock(v1, v2, C, -factV * cFact);    // -[C]× · (dx2-dx0)
        addSkewBlock(v1, v0, C, factV * cFact);
        addSkewBlock(v1, v3, B, factV * cFact);     // [B]× · (dx3-dx0)
        addSkewBlock(v1, v0, B, -factV * cFact);

        // dfB (at v2)
        addOuterProduct(v2, v1, gb, ga, factV * v16);
        addOuterProduct(v2, v0, gb, ga, -factV * v16);
        addOuterProduct(v2, v2, gb, gb, factV * v16);
        addOuterProduct(v2, v0, gb, gb, -factV * v16);
        addOuterProduct(v2, v3, gb, gc, factV * v16);
        addOuterProduct(v2, v0, gb, gc, -factV * v16);

        addSkewBlock(v2, v3, A, -factV * cFact);
        addSkewBlock(v2, v0, A, factV * cFact);
        addSkewBlock(v2, v1, C, factV * cFact);
        addSkewBlock(v2, v0, C, -factV * cFact);

        // dfC (at v3)
        addOuterProduct(v3, v1, gc, ga, factV * v16);
        addOuterProduct(v3, v0, gc, ga, -factV * v16);
        addOuterProduct(v3, v2, gc, gb, factV * v16);
        addOuterProduct(v3, v0, gc, gb, -factV * v16);
        addOuterProduct(v3, v3, gc, gc, factV * v16);
        addOuterProduct(v3, v0, gc, gc, -factV * v16);

        addSkewBlock(v3, v1, B, -factV * cFact);
        addSkewBlock(v3, v0, B, factV * cFact);
        addSkewBlock(v3, v2, A, factV * cFact);
        addSkewBlock(v3, v0, A, -factV * cFact);

        // df0 (at v0) = -(dfA+dfB+dfC)
        // Row v0 gets the negative sum of rows v1, v2, v3
        Coord gsum = ga + gb + gc;
        addOuterProduct(v0, v1, gsum, ga, -factV * v16);
        addOuterProduct(v0, v2, gsum, gb, -factV * v16);
        addOuterProduct(v0, v3, gsum, gc, -factV * v16);
        addOuterProduct(v0, v0, gsum, ga, factV * v16);
        addOuterProduct(v0, v0, gsum, gb, factV * v16);
        addOuterProduct(v0, v0, gsum, gc, factV * v16);
    }
}

template<class DataTypes>
void DistanceGridForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    drawDistanceGrid(vparams);
}

template<class DataTypes>
void DistanceGridForceField<DataTypes>::drawDistanceGrid(const core::visual::VisualParams* vparams,float size)
{
    if (!grid) return;
    if (size == 0.0f) size = (float)drawSize.getValue();

    const VecCoord& p1 = this->mstate->read(core::vec_id::read_access::position)->getValue();

    std::vector< type::Vec3 > pointsLineIn;
    std::vector< type::Vec3 > pointsLineOut;
    // lines for points penetrating the distancegrid

    unsigned int ibegin = 0;
    unsigned int iend = p1.size();

    if (localRange.getValue()[0] >= 0)
        ibegin = localRange.getValue()[0];

    if (localRange.getValue()[1] >= 0 && (unsigned int)localRange.getValue()[1]+1 < iend)
        iend = localRange.getValue()[1]+1;

    const Real stiffIn = stiffnessIn.getValue();
    const Real stiffOut = stiffnessOut.getValue();
    const Real maxdist = maxDist.getValue();

    type::Vec3 point1,point2;
    for (unsigned int i=ibegin; i<iend; i++)
    {
        if (i < pOnBorder.size() && !pOnBorder[i]) continue;
        Real d = grid->teval(p1[i]);
        if (d > 0)
        {
            if (d >= maxdist || stiffOut == 0) continue;
        }
        else if (d < 0)
        {
            if (-d >= maxdist || stiffIn == 0) continue;
        }
        else continue;
        Coord p2 = p1[i];
        Deriv normal = grid->tgrad(p1[i]); //normal.normalize();
        p2 += normal*(-d);
        point1 = DataTypes::getCPos(p1[i]);
        point2 = DataTypes::getCPos(p2);
        if (d > 0)
        {
            pointsLineOut.push_back(point1);
            pointsLineOut.push_back(point2);
        }
        else //if (d < 0)
        {
            pointsLineIn.push_back(point1);
            pointsLineIn.push_back(point2);
        }
    }
    if (!pointsLineIn.empty())
        vparams->drawTool()->drawLines(pointsLineIn, 1, sofa::type::RGBAColor(1.,0.,0.,1.));
    if (!pointsLineOut.empty())
        vparams->drawTool()->drawLines(pointsLineOut, 1, sofa::type::RGBAColor(1.,0.,1.,1.));

    const sofa::type::vector<TContact>& tcontacts = this->tcontacts.getValue();
    if (!tcontacts.empty())
    {
        std::vector< type::Vec3 > pointsTri;
        for (unsigned int i=0; i<tcontacts.size(); i++)
        {
            const TContact& c = (this->tcontacts.getValue())[i];
            type::Vec3 p;
            for (int j=0; j<3; ++j)
            {
                p = DataTypes::getCPos(p1[c.index[j]]);
                pointsTri.push_back(p);
            }
        }
        vparams->drawTool()->drawTriangles(pointsTri, sofa::type::RGBAColor{ 1.0f,0.2f,0.2f,0.5f });
    }
    const sofa::type::vector<VContact>& vcontacts = this->vcontacts.getValue();
    if (!vcontacts.empty())
    {
        std::vector< type::Vec3 > pointsTet;
        for (unsigned int i=0; i<vcontacts.size(); i++)
        {
            const VContact& c = (this->vcontacts.getValue())[i];
            const type::fixed_array<unsigned int,4>& t = c.index;
            type::Vec3 p[4];
            Coord pc = (p1[t[0]]+p1[t[1]]+p1[t[2]]+p1[t[3]])*0.25f;
            for (int j=0; j<4; ++j)
            {
                Coord pj = p1[c.index[j]];
                pj += (pc-pj)*0.2f;
                p[j] = DataTypes::getCPos(pj);
            }
            pointsTet.push_back(p[0]);
            pointsTet.push_back(p[1]);
            pointsTet.push_back(p[2]);
            pointsTet.push_back(p[0]);
            pointsTet.push_back(p[2]);
            pointsTet.push_back(p[3]);
            pointsTet.push_back(p[0]);
            pointsTet.push_back(p[3]);
            pointsTet.push_back(p[1]);
            pointsTet.push_back(p[1]);
            pointsTet.push_back(p[3]);
            pointsTet.push_back(p[2]);
        }
        vparams->drawTool()->drawTriangles(pointsTet, sofa::type::RGBAColor{ 0.8f,0.8f,0.0f,0.25f });
    }

    if (drawPoints.getValue())
    {
        std::vector< type::Vec3 > distancePointsIn;
        std::vector< type::Vec3 > distancePointsOut;

        for (int i=0; i < grid->getNx(); i++)
            for (int j=0; j < grid->getNy(); j++)
                for (int k=0; k < grid->getNz(); k++)
                {
                    Coord cellCoord = grid->coord(i,j,k);
                    if (grid->teval(cellCoord) < 0.0)
                        distancePointsIn.push_back(cellCoord);
                    else
                        distancePointsOut.push_back(cellCoord);
                }

        if (distancePointsIn.size())
            vparams->drawTool()->drawPoints(distancePointsIn, (float)drawSize.getValue(), sofa::type::RGBAColor{ 0.8f,0.2f,0.2f,1.0f });
        if (distancePointsOut.size())
            vparams->drawTool()->drawPoints(distancePointsOut, (float)drawSize.getValue() * 1.2f, sofa::type::RGBAColor{ 0.2f,0.8f,0.2f,1.0f });
    }

}
} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_INTERACTIONFORCEFIELD_DISTANCEGRIDFORCEFIELD_INL
