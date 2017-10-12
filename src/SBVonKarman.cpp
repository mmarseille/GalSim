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

// #define DEBUGLOGGING

#include "galsim/IgnoreWarnings.h"

#define BOOST_NO_CXX11_SMART_PTR
#include <boost/math/special_functions/bessel.hpp>

#include "SBVonKarman.h"
#include "SBVonKarmanImpl.h"
#include "fmath/fmath.hpp"

namespace galsim {

    inline double fast_pow(double x, double y)
    { return fmath::expd(y * std::log(x)); }

    //
    //
    //
    //SBVonKarman
    //
    //
    //

    SBVonKarman::SBVonKarman(double lam, double r0, double L0, double kcrit, double flux,
                             double scale, const GSParamsPtr& gsparams) :
        SBProfile(new SBVonKarmanImpl(lam, r0, L0, kcrit, flux, scale, gsparams)) {}

    SBVonKarman::SBVonKarman(const SBVonKarman &rhs) : SBProfile(rhs) {}

    SBVonKarman::~SBVonKarman() {}

    double SBVonKarman::getLam() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getLam();
    }

    double SBVonKarman::getR0() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getR0();
    }

    double SBVonKarman::getL0() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getL0();
    }

    double SBVonKarman::getKCrit() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getKCrit();
    }

    double SBVonKarman::getScale() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getScale();
    }

    double SBVonKarman::structureFunction(double rho) const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).structureFunction(rho);
    }

    //
    //
    //
    //VonKarmanInfo
    //
    //
    //

    VonKarmanInfo::VonKarmanInfo(double L0, const GSParamsPtr& gsparams, double kcrit) :
        _L0invsq(1./L0/L0), _gsparams(gsparams), _kcrit(kcrit), _sfTab(TableDD::spline)
    {
        buildSFTab();
    }

    class VonKarmanInfoIntegrand : public std::unary_function<double,double>
    {
    public:
        VonKarmanInfoIntegrand(double L0invsq, double rho) :
            _L0invsq(L0invsq), _rho(rho) {}
        double operator()(double kappa) const
        { return fast_pow(kappa*kappa+_L0invsq, -11./6)*(1.-j0(_rho*kappa))*kappa; }
    private:
        double _rho;
        double _L0invsq;
    };

    double VonKarmanInfo::structureFunction(double rho) const {
        VonKarmanInfoIntegrand I(_L0invsq, rho);
        const double magic = 0.033/0.423*8*M_PI*M_PI;
        return magic*integ::int1d(I, _kcrit, integ::MOCK_INF,
                                  _gsparams->integration_relerr,
                                  1e-4*_gsparams->integration_abserr);
    }

    void VonKarmanInfo::buildSFTab() {
        // 1 mm corresponds to a realspace folding scale of 30 arcsec at a wavelength of
        // 300nm.  I don't expect we'd ever need anything more demanding than that.  Plus,
        // the turbulence inner scale starts to matter around here, and we're not modeling
        // that.  So a mm seems sufficient.
        double log_rho_min = log(1e-3);
        // LSST is probably the largest aperture in play for GalSim, so using 10 m should
        // be sufficient.
        double log_rho_max = log(10.0);
        double dlogrho = _gsparams->table_spacing * sqrt(sqrt(_gsparams->kvalue_accuracy / 10.0));

        for (double logrho = log_rho_min; logrho < log_rho_max; logrho += dlogrho)
            _sfTab.addEntry(logrho, structureFunction(exp(logrho)));
    }

    double VonKarmanInfo::structureFunctionTab(double rho) const {
        return _sfTab(log(rho));
    }

    double VonKarmanInfo::maxRho(double r0m53) const {
        // param[in] r0m53 = r0^(-5/3)
        //
        // Want to find kmax s.t. MTF(k>kmax) < maxk_threshold
        // MTF = exp(-0.5*r0^(-5/3)*structureFunction(rho))
        // so really want to find rhomax s.t. rho > rhomax implies
        // structureFunction(rho) > -2*r0^(5/3)*log(maxk_threshold)
        double thresh = -2.0/r0m53*log(_gsparams->maxk_threshold);
        double rhomax = _sfTab.argMax();
        const std::vector<double> args = _sfTab.getArgs();
        const std::vector<double> vals = _sfTab.getVals();
        typedef std::vector<double>::const_iterator CIter;
        CIter a, v;
        for (a=args.begin(), v=vals.begin();
             (a+1) != args.end();
             a++, v++)
        {
            if (*v < thresh)
                rhomax = exp(*(a+1));
        }
        return rhomax;
    }

    LRUCache<boost::tuple<double,GSParamsPtr,double>,VonKarmanInfo>
        SBVonKarman::SBVonKarmanImpl::cache(sbp::max_vonKarman_cache);

    //
    //
    //
    //SBVonKarmanImpl
    //
    //
    //

    SBVonKarman::SBVonKarmanImpl::SBVonKarmanImpl(double lam, double r0, double L0, double kcrit,
                                                  double flux, double scale,
                                                  const GSParamsPtr& gsparams) :
        SBProfileImpl(gsparams),
        _lam(lam),
        _r0(r0),
        _r0m53(fast_pow(r0, -5./3)),
        _L0(L0),
        _kcrit(kcrit),
        _flux(flux),
        _scale(scale),
        _info(cache.get(boost::make_tuple(L0, this->gsparams.duplicate(), kcrit)))
    {
        dbg<<"SBVonKarmanImpl constructor: gsparams = "<<gsparams.get()<<std::endl;
        dbg<<"this->gsparams = "<<this->gsparams.get()<<std::endl;
    }

    double SBVonKarman::SBVonKarmanImpl::maxK() const
    { return _info->maxRho(_r0m53)/_lam*(2*M_PI)*_scale; }

    double SBVonKarman::SBVonKarmanImpl::stepK() const
    { return 0.1; }

    std::string SBVonKarman::SBVonKarmanImpl::serialize() const
    {
        std::ostringstream oss(" ");
        oss.precision(std::numeric_limits<double>::digits10 + 4);
        oss << "galsim._galsim.SBVonKarman("
            <<getLam()<<", "
            <<getR0()<<", "
            <<getL0()<<", "
            <<getKCrit()<<", "
            <<getFlux()<<", galsim.GSParams("<<*gsparams<<"))";
        return oss.str();
    }

    class VonKarmanRadialFunction : public FluxDensity
    {
    public:
        VonKarmanRadialFunction() {}
        double operator()(double r) const { return 0.; }
    private:
    };

    boost::shared_ptr<PhotonArray>
    SBVonKarman::SBVonKarmanImpl::shoot(int N, UniformDeviate ud) const
    {
        // Placeholder.  Nothing here...
        dbg<<"VonKarman::shoot: N = "<<N<<std::endl;
        dbg<<"Target flux = 1.0"<<std::endl;

        if (!_sampler) {
            std::vector<double> range(2, 0.);
            range[1] = 1.0; //
            _radial.reset(new VonKarmanRadialFunction());
            _sampler.reset(new OneDimensionalDeviate(*_radial, range, true, gsparams));
        }

        assert(_sampler.get());
        boost::shared_ptr<PhotonArray> result = _sampler->shoot(N,ud);
        dbg<<"VonKarman Realized flux = "<<result->getTotalFlux()<<std::endl;
        return result;
    }

    double SBVonKarman::SBVonKarmanImpl::structureFunction(double rho) const
    {
        dbg<<"rho = "<<rho<<'\n';
        return _r0m53*_info->structureFunction(rho);
    }

    std::complex<double> SBVonKarman::SBVonKarmanImpl::kValue(const Position<double>& p) const
    {
        double k = sqrt(p.x*p.x+p.y*p.y)/_scale;
        dbg<<"k = "<<k<<'\n';
        return fmath::expd(-0.5*_r0m53*_info->structureFunctionTab(k*_lam/(2.*M_PI)));
    }
}
