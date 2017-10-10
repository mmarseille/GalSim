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

    //SBVonKarman

    SBVonKarman::SBVonKarman(double lam, double r0, double L0, double kcrit, double flux,
                             double maxk, const GSParamsPtr& gsparams) :
        SBProfile(new SBVonKarmanImpl(lam, r0, L0, kcrit, flux, maxk, gsparams)) {}

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

    double SBVonKarman::structureFunction(double rho) const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).structureFunction(rho);
    }

    //VonKarmanInfo

    VonKarmanInfo::VonKarmanInfo(const GSParamsPtr& gsparams, const double beta_min) :
        _beta_min(beta_min), _gsparams(gsparams), _sfTab(TableDD::spline)
    {
        
    }

    class VonKarmanInfoIntegrand : public std::unary_function<double,double>
    {
    public:
        VonKarmanInfoIntegrand(double nu) : _nu(nu) {}
        double operator()(double beta) const
        { return fast_pow((1.0+beta*beta), -11.0/6.0)*(1.0-j0(_nu*beta))*beta; }
    private:
        double _nu;
    };

    double VonKarmanInfo::structureFunction(double nu) const {
        VonKarmanInfoIntegrand I(nu);
        return integ::int1d(I, 0., integ::MOCK_INF,
                            _gsparams->integration_relerr,
                            _gsparams->integration_abserr);
    }

    LRUCache<boost::tuple<GSParamsPtr,double>,VonKarmanInfo>
        SBVonKarman::SBVonKarmanImpl::cache(sbp::max_vonKarman_cache);

    //SBVonKarmanImpl

    SBVonKarman::SBVonKarmanImpl::SBVonKarmanImpl(double lam, double r0, double L0, double kcrit,
                                                  double flux, double maxk,
                                                  const GSParamsPtr& gsparams) :
        SBProfileImpl(gsparams),
        _lam(lam),
        _r0(r0),
        _L0(L0),
        _kcrit(kcrit),
        _flux(flux),
        _maxk(maxk),
        _info(cache.get(boost::make_tuple(this->gsparams.duplicate(), L0*kcrit)))
    {
        dbg<<"SBVonKarmanImpl constructor: gsparams = "<<gsparams.get()<<std::endl;
        dbg<<"this->gsparams = "<<this->gsparams.get()<<std::endl;
    }

    double SBVonKarman::SBVonKarmanImpl::maxK() const
    { return 25.0; }

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
            <<getFlux()<<", "
            <<maxK()<<", galsim.GSParams("<<*gsparams<<"))";
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
        return 0.033/0.423*8*M_PI*M_PI*_info->structureFunction(rho/_L0)*fast_pow(_r0/_L0, -5.0/3.0);
    }

    std::complex<double> SBVonKarman::SBVonKarmanImpl::kValue(const Position<double>& p) const
    {
        double k = sqrt(p.x*p.x+p.y*p.y);
        return _flux*fmath::expd(-0.5*structureFunction(k*_lam/(2.*M_PI)));
    }
}
