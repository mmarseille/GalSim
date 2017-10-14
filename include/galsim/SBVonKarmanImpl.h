/* -*- c++ -*-
 * Copyright (c) 2012-2017 by the GalSim developers team on GitHub
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

#ifndef GalSim_SBVonKarmanImpl_H
#define GalSim_SBVonKarmanImpl_H

#include "SBProfileImpl.h"
#include "SBVonKarman.h"
#include "LRUCache.h"
#include "OneDimensionalDeviate.h"
#include "Table.h"

namespace galsim {

    //
    //
    //
    //VonKarmanInfo
    //
    //
    //

    // I'm not totally sure if it's worth having this class.  It *may* be that we get some
    // reuse out of it as we change wavelengths, but keep L0, gsparams, and kcrit the
    // same.  But I'm not sure if we'd want to keep kcrit the same across wavelengths, or
    // if it should scale like r0^-1, i.e., lambda^(-6/5).  Keeping it here for now
    // though.
    class VonKarmanInfo
    {
    public:
        VonKarmanInfo(double L0, const GSParamsPtr& gsparams, double kcrit=0.0);

        ~VonKarmanInfo() {}

        double structureFunction(double rho) const;
        double structureFunctionTab(double rho) const;
        void buildSFTab();
        double maxRho(double r0m53) const;
        // boost::shared_ptr<PhotonArray> shoot(int N, UniformDeviate ud) const;

    private:
        VonKarmanInfo(const VonKarmanInfo& rhs); ///<Hide the copy constructor
        void operator=(const VonKarmanInfo& rhs); ///<Hide the assignment operator

        const double _L0invsq;
        const GSParamsPtr _gsparams;
        const double _kcrit; ///<Lower turbulence mode cutoff

        TableDD _sfTab;
    };

    //
    //
    //
    //SBVonKarmanImpl
    //
    //
    //

    class SBVonKarman::SBVonKarmanImpl : public SBProfileImpl
    {
    public:
        SBVonKarmanImpl(double lam, double r0, double L0, double kcrit, double flux, double scale,
                        const GSParamsPtr& gsparams);
        ~SBVonKarmanImpl() {}

        bool isAxisymmetric() const { return true; }
        bool hasHardEdges() const { return false; }
        bool isAnalyticX() const { return false; }
        bool isAnalyticK() const { return true; }

        double maxK() const;
        double stepK() const;

        Position<double> centroid() const { return Position<double>(0., 0.); }

        double getFlux() const { return _flux; }
        double getLam() const { return _lam; }
        double getR0() const { return _r0; }
        double getL0() const { return _L0; }
        double getKCrit() const { return _kcrit; }
        double getScale() const { return _scale; }
        // double maxSB();// const { return _xnorm * _info->xValue(0.); }
        double maxSB() const { return 1.0; }  // no idea how right/wrong this is.

        /**
         * @brief SBVonKarman photon-shooting is done numerically with `OneDimensionalDeviate`
         * class.
         *
         * @param[in] N Total number of photons to produce.
         * @param[in] ud UniformDeviate that will be used to draw photons from distribution.
         * @returns PhotonArray containing all the photons' info.
         */
        boost::shared_ptr<PhotonArray> shoot(int N, UniformDeviate ud) const;

        double xValue(const Position<double>& p) const;
        double xValue(double r) const;
        std::complex<double> kValue(const Position<double>& p) const;
        std::complex<double> kValue(double k) const;

        double structureFunction(double rho) const;

        std::string serialize() const;

    private:

        double _lam;
        double _r0;
        double _r0m53;
        double _L0;
        double _kcrit;
        double _flux;
        double _scale;
        double _kappa_min;

        boost::shared_ptr<VonKarmanInfo> _info; ///< Points to info structure for this
                                                ///  beta_min = L0*kcrit

        // Copy constructor and op= are undefined.
        SBVonKarmanImpl(const SBVonKarmanImpl& rhs);
        void operator=(const SBVonKarmanImpl& rhs);

        mutable boost::shared_ptr<FluxDensity> _radial;
        mutable boost::shared_ptr<OneDimensionalDeviate> _sampler;
        // mutable Table<double,double> _structure_fn;
        // mutable Table<double,double> _PSF;

        static LRUCache<boost::tuple<double,GSParamsPtr,double>,VonKarmanInfo> cache;
    };
}

#endif
