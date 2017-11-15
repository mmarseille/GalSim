/* -*- c++ -*-
 * Copyright (c) 2012-2016 by the GalSim developers team on GitHub
 * https://github.com/GalSim-developers
 *
 * This file is part of GalSim: The modular galaxy image simulation toolkit.
 * https://github.com/GalSim-developers/GalSim
 *
 * GalSim is free software: redistribution and use in source and binary forms,
 * with or without modification, are permitted provided that the following
 * conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions, and the disclaimer given in the accompanying LICENSE
 *    file.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions, and the disclaimer given in the documentation
 *    and/or other materials provided with the distribution.
 */

#ifndef GalSim_SBSecondKickImpl_H
#define GalSim_SBSecondKickImpl_H

#include "SBProfileImpl.h"
#include "SBSecondKick.h"
#include "LRUCache.h"
#include "OneDimensionalDeviate.h"
#include "Table.h"

namespace galsim {

    //
    //
    //
    //SecondKickInfo
    //
    //
    //

    class SecondKickInfo
    {
    public:
        SecondKickInfo(double lam, double r0, double L0, double kcrit, bool doDelta,
                       const GSParamsPtr& gsparams);

        ~SecondKickInfo() {}

        double stepK() const { return _stepk; }
        double maxK() const { return _maxk; }
        double getDeltaAmplitude() const { return _deltaAmplitude; }
        double getHalfLightRadius() const { return _hlr; }

        double kValue(double) const;
        double xValue(double) const;
        double structureFunction(double) const;
        double structureFunctionDirect(double) const;
        double phasePower(double) const;

    private:
        SecondKickInfo(const SecondKickInfo& rhs); ///<Hide the copy constructor
        void operator=(const SecondKickInfo& rhs); ///<Hide the assignment operator

        double vKStructureFunction(double) const;
        double complementaryStructureFunction(double) const;
        double kValueNoTrunc(double) const;
        void computeDeltaAmplitude();

        double _lam; // Wavelength in meters
        double _r0; // Fried parameter in meters
        double _r0m53; // r0^(-5./3)
        double _L0; // Outer scale in meters
        double _2piL02; // (2pi/L0)^(2)
        double _r0L0m53; // (r0/L0)^(-5/3)
        double _kcrit;
        double _stepk;
        double _maxk;
        double _deltaAmplitude;
        bool _doDelta;
        double _hlr; // half-light-radius

        // Magic constants that we can compute once and store.
        const static double magic1;
        const static double magic2;
        const static double magic3;
        const static double magic4;
        const static double magic5;

        const GSParamsPtr _gsparams;
    };

    class SBSecondKick::SBSecondKickImpl : public SBProfileImpl
    {
    public:
        SBSecondKickImpl(double lam, double r0, double L0, double kcrit,
                         double flux, double scale, bool doDelta, const GSParamsPtr& gsparams);
        ~SBSecondKickImpl() {}

        bool isAxisymmetric() const { return true; }
        bool hasHardEdges() const { return false; }
        bool isAnalyticX() const { return false; }
        bool isAnalyticK() const { return true; }

        double maxK() const;
        double stepK() const;
        double getDeltaAmplitude() const;
        double getHalfLightRadius() const;

        Position<double> centroid() const { return Position<double>(0., 0.); }

        double getFlux() const { return _flux; }
        double getLam() const { return _lam; }
        double getR0() const { return _r0; }
        double getL0() const { return _L0; }
        double getKCrit() const { return _kcrit; }
        double getScale() const { return _scale; }
        bool getDoDelta() const {return _doDelta; }
        double maxSB() const { return 1.0; }  // no idea how right/wrong this is.

        /**
         * @brief SBSecondKick photon-shooting is done numerically with `OneDimensionalDeviate`
         * class.
         *
         * @param[in] N Total number of photons to produce.
         * @param[in] ud UniformDeviate that will be used to draw photons from distribution.
         * @returns PhotonArray containing all the photons' info.
         */
        boost::shared_ptr<PhotonArray> shoot(int N, UniformDeviate ud) const
        { throw SBError("SBSecondKick::shoot() is not implemented"); }

        double xValue(double) const;
        double xValue(const Position<double>& p) const;
        double kValue(double) const;
        std::complex<double> kValue(const Position<double>& p) const;

        double phasePower(double kappa) const;
        double structureFunction(double rho) const;
        double structureFunctionDirect(double rho) const;

        std::string serialize() const;

    private:
        double _lam;
        double _r0;
        double _L0;
        double _kcrit;
        double _flux;
        double _scale;
        bool _doDelta;

        boost::shared_ptr<SecondKickInfo> _info;

        // Copy constructor and op= are undefined.
        SBSecondKickImpl(const SBSecondKickImpl& rhs);
        void operator=(const SBSecondKickImpl& rhs);

        static LRUCache<boost::tuple<double,double,double,double,bool,GSParamsPtr>,SecondKickInfo>
            cache;
    };
}

#endif
