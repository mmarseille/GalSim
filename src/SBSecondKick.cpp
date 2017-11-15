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

// #define DEBUGLOGGING

#include "galsim/IgnoreWarnings.h"

#define BOOST_NO_CXX11_SMART_PTR
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include "SBSecondKick.h"
#include "SBSecondKickImpl.h"
#include "fmath/fmath.hpp"


namespace galsim {

    const double ARCSEC2RAD = 180.*60*60/M_PI;  // ~206265
    const double MOCK_INF = 1.e300;

    inline double fast_pow(double x, double y)
    { return fmath::expd(y * std::log(x)); }

    //
    //
    //
    //SBSecondKick
    //
    //
    //

    SBSecondKick::SBSecondKick(double lam, double r0, double L0, double kcrit, double flux,
                               double scale, bool doDelta, const GSParamsPtr& gsparams) :
        SBProfile(new SBSecondKickImpl(lam, r0, L0, kcrit, flux, scale, doDelta, gsparams)) {}

    SBSecondKick::SBSecondKick(const SBSecondKick &rhs) : SBProfile(rhs) {}

    SBSecondKick::~SBSecondKick() {}

    double SBSecondKick::getLam() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getLam();
    }

    double SBSecondKick::getR0() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getR0();
    }

    double SBSecondKick::getL0() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getL0();
    }

    double SBSecondKick::getKCrit() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getKCrit();
    }

    double SBSecondKick::getScale() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getScale();
    }

    bool SBSecondKick::getDoDelta() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getDoDelta();
    }

    double SBSecondKick::getDeltaAmplitude() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getDeltaAmplitude();
    }

    double SBSecondKick::getHalfLightRadius() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getHalfLightRadius();
    }

    double SBSecondKick::phasePower(double kappa) const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).phasePower(kappa);
    }

    double SBSecondKick::structureFunction(double rho) const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).structureFunction(rho);
    }

    double SBSecondKick::structureFunctionDirect(double rho) const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).structureFunctionDirect(rho);
    }

    //
    //
    //
    //SecondKickInfo
    //
    //
    //

    const double SecondKickInfo::magic1 = 2*boost::math::tgamma(11./6)/(pow(2, 5./6)*pow(M_PI, 8./3))
                                          * pow(24/5.*boost::math::tgamma(6./5), 5./6);
    const double SecondKickInfo::magic2 = boost::math::tgamma(5./6)/pow(2., 1./6);
    const double SecondKickInfo::magic3 = SecondKickInfo::magic1*boost::math::tgamma(-5./6)/pow(2., 11./6);
    const double SecondKickInfo::magic4 = boost::math::tgamma(11./6)*boost::math::tgamma(5./6)
                                          / pow(M_PI,8./3)
                                          * pow(24./5*boost::math::tgamma(6./5),5./6);
    const double SecondKickInfo::magic5 = boost::math::tgamma(11./6)*boost::math::tgamma(11./6)
                                          * pow(2., 8./3)
                                          * pow(24./5*boost::math::tgamma(6./5),5./6);

    SecondKickInfo::SecondKickInfo(double lam, double r0, double L0, double kcrit, bool doDelta,
                                   const GSParamsPtr& gsparams) :
        _lam(lam), _r0(r0), _r0m53(fast_pow(r0, -5./3)), _L0(L0), _2piL02(4*M_PI*M_PI/L0/L0),
        _r0L0m53(fast_pow(r0/L0, -5./3)), _kcrit(kcrit), _doDelta(doDelta), _gsparams(gsparams)
    {
        computeDeltaAmplitude();
    }

    double SecondKickInfo::phasePower(double kappa) const {
        return magic5*_r0m53*fast_pow(kappa*kappa+_2piL02, -11./6);
    }

    double SecondKickInfo::vKStructureFunction(double rho) const {
    // rho in meters
        double rhoL0 = rho/_L0;
        if (rhoL0 < 1e-6) {
            return -magic3*fast_pow(2*M_PI*rho/_r0, 5./3);
        } else {
            double x = 2*M_PI*rhoL0;
            return magic1*_r0L0m53*(magic2-fast_pow(x, 5./6)*boost::math::cyl_bessel_k(5./6, x));
        }
    }

    class SFIntegrand : public std::unary_function<double,double>
    {
    public:
        SFIntegrand(double rho, const SecondKickInfo& ski) : _rho(rho), _ski(ski) {}
        double operator()(double kappa) const {
            return _ski.phasePower(kappa)*(1-j0(kappa*_rho))*kappa;
        }
    private:
        const double _rho;
        const SecondKickInfo& _ski;
    };

    double SecondKickInfo::structureFunctionDirect(double rho) const {
        SFIntegrand I(rho, *this);
        return integ::int1d(I, 0.0, integ::MOCK_INF,
                            _gsparams->integration_relerr, _gsparams->integration_abserr)/M_PI;
    }

    double SecondKickInfo::complementaryStructureFunction(double rho) const {
        SFIntegrand I(rho, *this);
        return integ::int1d(I, 0.0, _kcrit,
                            _gsparams->integration_relerr, _gsparams->integration_abserr)/M_PI;
    }

    double SecondKickInfo::structureFunction(double rho) const {
        return vKStructureFunction(rho) - complementaryStructureFunction(rho);
    }

    class DeltaFunctionAmplitudeIntegrand : public std::unary_function<double,double>
    {
    public:
        DeltaFunctionAmplitudeIntegrand(const SecondKickInfo& ski) : _ski(ski) {}
        double operator()(double kappa) const { return _ski.phasePower(kappa)*kappa; }
    private:
        const SecondKickInfo& _ski;
    };

    void SecondKickInfo::computeDeltaAmplitude() {
        DeltaFunctionAmplitudeIntegrand I(*this);
        double integral = integ::int1d(
            I, 0.0, _kcrit, _gsparams->integration_relerr, _gsparams->integration_abserr)/M_PI;
        _deltaAmplitude = exp(-0.5*(magic4*_r0L0m53 - integral));
    }

    LRUCache<boost::tuple<double,double,double,double,bool,GSParamsPtr>,SecondKickInfo>
        SBSecondKick::SBSecondKickImpl::cache(sbp::max_secondKick_cache);

    //
    //
    //
    //SBSecondKickImpl
    //
    //
    //

    SBSecondKick::SBSecondKickImpl::SBSecondKickImpl(double lam, double r0, double L0, double kcrit,
                                                     double flux, double scale, bool doDelta,
                                                     const GSParamsPtr& gsparams) :
        SBProfileImpl(gsparams),
        _lam(lam),
        _r0(r0),
        _L0(L0),
        _kcrit(kcrit),
        _flux(flux),
        _scale(scale),
        _doDelta(doDelta),
        _info(cache.get(boost::make_tuple(1e-9*lam, r0, L0, kcrit, doDelta,
                                          this->gsparams.duplicate())))
    { }

    double SBSecondKick::SBSecondKickImpl::maxK() const
    { return _info->maxK()*_scale; }

    double SBSecondKick::SBSecondKickImpl::stepK() const
    { return _info->stepK()*_scale; }

    double SBSecondKick::SBSecondKickImpl::getDeltaAmplitude() const
    { return _info->getDeltaAmplitude()*_flux; }

    double SBSecondKick::SBSecondKickImpl::getHalfLightRadius() const
    { return _info->getHalfLightRadius()*_scale; }

    std::string SBSecondKick::SBSecondKickImpl::serialize() const
    {
        std::ostringstream oss (" ");
        oss.precision(std::numeric_limits<double>::digits10 + 4);
        oss << "galsim._galsim.SBSecondKick("
            <<getLam()<<", "
            <<getR0()<<", "
            <<getL0()<<", "
            <<getKCrit()<<", "
            <<getFlux()<<", "
            <<getScale()<<", "
            <<getDoDelta()<<", "
            <<"galsim.GSParams("<<*gsparams<<"))";
        return oss.str();
    }

    double SBSecondKick::SBSecondKickImpl::phasePower(double kappa) const
    { return _info->phasePower(kappa); }

    double SBSecondKick::SBSecondKickImpl::structureFunction(double rho) const
    { return _info->structureFunction(rho); }

    double SBSecondKick::SBSecondKickImpl::structureFunctionDirect(double rho) const
    { return _info->structureFunctionDirect(rho); }
}
