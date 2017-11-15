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
#include "Solve.h"


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

    class SKIkValueResid {
    public:
        SKIkValueResid(const SecondKickInfo& ski, double mkt) : _ski(ski), _mkt(mkt) {}
        double operator()(double k) const {
            double val = _ski.kValue(k)-_mkt;
            xdbg<<"resid(k="<<k<<")="<<val<<'\n';
            return val;
         }
    private:
        const double _mkt;
        const SecondKickInfo& _ski;
    };

    SecondKickInfo::SecondKickInfo(double lam, double r0, double L0, double kcrit, bool doDelta,
                                   const GSParamsPtr& gsparams) :
        _lam(lam), _r0(r0), _r0m53(fast_pow(r0, -5./3)), _L0(L0), _2piL02(4*M_PI*M_PI/L0/L0),
        _r0L0m53(fast_pow(r0/L0, -5./3)), _kcrit(kcrit), _doDelta(doDelta), _gsparams(gsparams)
    {
        computeDeltaAmplitude();
        // determine maxK
        // want abs(kValue(k>maxK))/kValue(0.0) <= _gsparams->maxk_threshold;
        // note that kValue(0.0) = 1.
        // Harder this time compared to VonKarman, since kValue may oscillate back above 0.0...
        double mkt = _gsparams->maxk_threshold;
        if (_doDelta) {
            if (mkt < _deltaAmplitude) {
                // If the delta function amplitude is too large, then no matter how far out in k we
                // go, kValue never drops below that amplitude.
                // _maxk = std::numeric_limits<double>::infinity();
                _maxk = MOCK_INF;
            } else {
                mkt = mkt*(1-_deltaAmplitude)+_deltaAmplitude;
            }
        }
        // if (_maxk != std::numeric_limits<double>::infinity()) {
        if (_maxk != MOCK_INF) {
            SKIkValueResid skikvr(*this, mkt);
            Solve<SKIkValueResid> solver(skikvr, 0.1, 1);
            solver.bracket();
            solver.setMethod(Brent);
            _maxk = solver.root();
        }
        dbg<<"_maxk = "<<_maxk<<" arcsec^-1\n";
        dbg<<"SB(maxk) = "<<kValue(_maxk)<<'\n';
        dbg<<"_deltaAmplitude = "<<_deltaAmplitude<<'\n';
        _stepk = 0.1;
    }

    double SecondKickInfo::kValueNoTrunc(double k) const {
    // k in inverse arcsec
        return fmath::expd(-0.5*structureFunction(_lam*k*ARCSEC2RAD/(2*M_PI)));
    }

    double SecondKickInfo::kValue(double k) const {
    // k in inverse arcsec
    // We're subtracting the asymptotic kValue limit here so that kValue->0 as k->inf.
    // This means we should also rescale by (1-_deltaAmplitude) though, so we still retain
    // kValue(0)=1.
        double val = (kValueNoTrunc(k) - _deltaAmplitude)/(1-_deltaAmplitude);
        if (std::abs(val) < std::numeric_limits<double>::epsilon())
            return 0.0;
        return val;
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

    class SKIXIntegrand : public std::unary_function<double,double>
    {
    public:
        SKIXIntegrand(double r, const SecondKickInfo& ski) : _r(r), _ski(ski) {}
        double operator()(double k) const { return _ski.kValue(k)*j0(k*_r)*k; }
    private:
        const double _r;  //arcsec
        const SecondKickInfo& _ski;
    };

    double SecondKickInfo::xValue(double r) const {
    // r in arcsec
        SKIXIntegrand I(r, *this);
        integ::IntRegion<double> reg(0, integ::MOCK_INF);
        return integ::int1d(I, reg,
                            _gsparams->integration_relerr,
                            _gsparams->integration_abserr)/(2.*M_PI);
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

    double SBSecondKick::SBSecondKickImpl::kValue(double k) const
    // this kValue assumes k is in inverse arcsec
    {
        return _info->kValue(k)*_flux;
    }

    std::complex<double> SBSecondKick::SBSecondKickImpl::kValue(const Position<double>& p) const
    // k in units of _scale.
    {
        return kValue(sqrt(p.x*p.x+p.y*p.y)/_scale);
    }

    double SBSecondKick::SBSecondKickImpl::xValue(double r) const
    {
        return _info->xValue(r)*_flux;
    }

    double SBSecondKick::SBSecondKickImpl::xValue(const Position<double>& p) const
    {
        return xValue(sqrt(p.x*p.x+p.y*p.y)*_scale);
    }

}
