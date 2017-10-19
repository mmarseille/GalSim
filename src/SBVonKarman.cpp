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
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include "SBVonKarman.h"
#include "SBVonKarmanImpl.h"
#include "fmath/fmath.hpp"
#include "Solve.h"
#include "bessel/Roots.h"

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

    SBVonKarman::SBVonKarman(double lam, double r0, double L0, double flux,
                             double scale, const GSParamsPtr& gsparams) :
        SBProfile(new SBVonKarmanImpl(lam, r0, L0, flux, scale, gsparams)) {}

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

    double SBVonKarman::getScale() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getScale();
    }

    double SBVonKarman::getDeltaAmplitude() const
    {
        assert(dynamic_cast<const SBVonKarmanImpl*>(_pimpl.get()));
        return static_cast<const SBVonKarmanImpl&>(*_pimpl).getDeltaAmplitude();
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

    class VKIkValueResid {
    public:
        VKIkValueResid(const VonKarmanInfo& vki, double mkt) : _vki(vki), _mkt(mkt) {}
        double operator()(double k) const {
            double val = _vki.kValue(k)-_mkt;
            xdbg<<"resid(k="<<k<<")="<<val<<'\n';
            return val;
         }
    private:
        double _mkt;
        const VonKarmanInfo& _vki;
    };

    double VonKarmanInfo::magic1 = 2*boost::math::tgamma(11./6)/(pow(2, 5./6)*pow(M_PI, 8./3))
                                    *pow(24/5.*boost::math::tgamma(6./5), 5./6);
    double VonKarmanInfo::magic2 = boost::math::tgamma(5./6)/pow(2., 1./6);
    double VonKarmanInfo::magic3 = VonKarmanInfo::magic1*boost::math::tgamma(-5./6)/pow(2., 11./6);
    double VonKarmanInfo::magic4 = boost::math::tgamma(11./6)*boost::math::tgamma(5./6)
                                   / pow(M_PI,8./3)
                                   * pow(24./5*boost::math::tgamma(6./5),5./6);

    VonKarmanInfo::VonKarmanInfo(double lam, double r0, double L0, const GSParamsPtr& gsparams) :
        _lam(lam), _r0(r0), _L0(L0), _r0L0m53(pow(r0/L0, -5./3)), _gsparams(gsparams),
        _deltaAmplitude(exp(-0.5*magic4*_r0L0m53)),
        _radial(TableDD::spline)
    {
        dbg << "magic1 = " << magic1 << '\n';
        dbg << "magic2 = " << magic2 << '\n';
        dbg << "magic3 = " << magic3 << '\n';
        dbg << "magic4 = " << magic4 << '\n';
        dbg << "r0 = " << r0 << '\n';
        dbg << "L0 = " << L0 << '\n';
        dbg << "r0L0m53 = " << _r0L0m53 << '\n';
        dbg << "deltaAmplitude = " << exp(-0.5*magic4*_r0L0m53) << '\n';
        dbg << "deltaAmplitudeInit = " << _deltaAmplitude << '\n';

        // determine maxK
        // want kValue(maxK)/kValue(0.0) = _gsparams->maxk_threshold;
        // note that kValue(0.0) = 1.
        VKIkValueResid vkikvr(*this, _gsparams->maxk_threshold);
        Solve<VKIkValueResid> solver(vkikvr, 1e5, 1e6);
        solver.bracket();
        solver.setMethod(Brent);
        _maxk = solver.root();
        dbg<<"_maxk = "<<_maxk<<" rad^-1 = "<<_maxk*M_PI/180/60/60<<" arcsec^-1\n";

        // build the radial function, and along the way, set _stepk
        _buildRadialFunc();
    }

    double VonKarmanInfo::structureFunction(double rho) const {
        double rhoL0 = rho/_L0;
        if (rhoL0 < 1e-6) {
            return -magic3*fast_pow(2*M_PI*rho/_r0, 5./3);
        } else {
            double x = 2*M_PI*rhoL0;
            return magic1*_r0L0m53*(magic2-fast_pow(x, 5./6)*boost::math::cyl_bessel_k(5./6, x));
        }
    }

    double VonKarmanInfo::kValue(double k) const {
        return fmath::expd(-0.5*structureFunction(_lam*k/(2*M_PI)));
    }

    class VKIXIntegrand : public std::unary_function<double,double>
    {
    public:
        VKIXIntegrand(double r, const VonKarmanInfo& vki) : _r(r), _vki(vki) {}
        double operator()(double k) const { return _vki.kValue(k)*j0(k*_r)*k; }
    private:
        double _r;
        const VonKarmanInfo& _vki;
    };

    double VonKarmanInfo::xValue(double r) const {
        VKIXIntegrand I(r, *this);
        integ::IntRegion<double> reg(0, integ::MOCK_INF);
        return integ::int1d(I, reg,
                            _gsparams->integration_relerr,
                            _gsparams->integration_abserr)/(2.*M_PI);
    }

    void VonKarmanInfo::_buildRadialFunc() {
        set_verbose(2);
        double r = 0.0;
        double val = xValue(0.0);
        _radial.addEntry(r, val);
        xdbg<<"f(0) = "<<val<<" rad^-2 = "<<val*M_PI*M_PI/180/180/60/60/60/60<<" arcsec^-2\n";

        double dr = _gsparams->table_spacing * sqrt(sqrt(_gsparams->xvalue_accuracy / 10.));
        dr *= M_PI/180/60/60;  // arcsec -> rad
        xdbg<<"dr = "<<dr<<" rad = "<<dr*180*60*60/M_PI<<" arcsec\n";

        double sum = 0.0;
        double thresh0 = 0.5 / (2.*M_PI*dr);
        double thresh1 = (1.-_gsparams->folding_threshold) / (2.*M_PI*dr);
        double thresh2 = (1.-_gsparams->folding_threshold/2) / (2.*M_PI*dr);
        double R = 1e10, hlr = 1e10;
        for(r = dr;
            (r < _gsparams->stepk_minimum_hlr*hlr) || (r < R) || (sum < thresh2);
            r += dr)
        {
            val = xValue(r);
            xdbg<<"f("<<r<<") = "<<val<<'\n';
            _radial.addEntry(r, val);
            sum += r*val;
            xdbg<<"sum = "<<sum<<"  thresh0 = "<<thresh0<<"  thresh1 = "<<thresh1<<'\n';
            xdbg<<"sum*2*pi*dr = "<<sum*2.*M_PI*dr<<'\n';
            if (hlr == 1e10 && sum > thresh0) hlr = r;
            if (R == 1e10 && sum > thresh1) R = r;
        }
        dbg<<"Done loop to build radial function.\n";
        dbg<<"R = "<<R<<" rad = "<<R*180*60*60/M_PI<<" arcsec\n";
        dbg<<"hlr = "<<hlr<<" rad = "<<hlr*180*60*60/M_PI<<" arcsec\n";
        R = std::max(R, _gsparams->stepk_minimum_hlr*hlr);
        _stepk = M_PI / R;
        dbg<<"stepk = "<<_stepk<<" rad^-1 = "<<_stepk*M_PI/180/60/60<<" arcsec^-1\n";
        dbg<<"sum*2*pi*dr = "<<sum*2.*M_PI*dr<<"   (should ~= 0.999)\n";

        std::vector<double> range(2, 0.);
        range[1] = _radial.argMax();
        _sampler.reset(new OneDimensionalDeviate(_radial, range, true, _gsparams));
    }

    boost::shared_ptr<PhotonArray> VonKarmanInfo::shoot(int N, UniformDeviate ud) const
    {
        dbg<<"VonKarmanInfo shoot: N = "<<N<<std::endl;
        dbg<<"Target flux = 1.0\n";
        assert(_sampler.get());
        boost::shared_ptr<PhotonArray> result = _sampler->shoot(N,ud);
        dbg<<"VonKarmanInfo Realized flux = "<<result->getTotalFlux()<<std::endl;
        return result;
    }

    LRUCache<boost::tuple<double,double,double,GSParamsPtr>,VonKarmanInfo>
        SBVonKarman::SBVonKarmanImpl::cache(sbp::max_vonKarman_cache);

    //
    //
    //
    //SBVonKarmanImpl
    //
    //
    //

    SBVonKarman::SBVonKarmanImpl::SBVonKarmanImpl(double lam, double r0, double L0, double flux,
                                                  double scale, const GSParamsPtr& gsparams) :
        SBProfileImpl(gsparams),
        _lam(lam),
        _r0(r0),
        _L0(L0),
        _flux(flux),
        _scale(scale),
        _info(cache.get(boost::make_tuple(lam, r0, L0, this->gsparams.duplicate())))
    {
        dbg<<"SBVonKarmanImpl constructor: gsparams = "<<gsparams.get()<<std::endl;
        dbg<<"this->gsparams = "<<this->gsparams.get()<<std::endl;
    }

    double SBVonKarman::SBVonKarmanImpl::maxK() const
    { return _info->maxK()*_scale; }

    double SBVonKarman::SBVonKarmanImpl::stepK() const
    { return _info->stepK()*_scale; }

    double SBVonKarman::SBVonKarmanImpl::getDeltaAmplitude() const
    { return _info->getDeltaAmplitude()*_flux; }

    std::string SBVonKarman::SBVonKarmanImpl::serialize() const
    {
        std::ostringstream oss(" ");
        oss.precision(std::numeric_limits<double>::digits10 + 4);
        oss << "galsim._galsim.SBVonKarman("
            <<getLam()<<", "
            <<getR0()<<", "
            <<getL0()<<", "
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

    double SBVonKarman::SBVonKarmanImpl::structureFunction(double rho) const
    {
        dbg<<"rho = "<<rho<<'\n';
        return _info->structureFunction(rho);
    }

    double SBVonKarman::SBVonKarmanImpl::kValue(double k) const
    // this kValue assumes k is in inverse radians
    {
        return _info->kValue(k);
    }

    std::complex<double> SBVonKarman::SBVonKarmanImpl::kValue(const Position<double>& p) const
    // k in units of _scale.
    {
        return kValue(sqrt(p.x*p.x+p.y*p.y)/_scale);
    }

    class VKXIntegrand : public std::unary_function<double,double>
    {
    public:
        VKXIntegrand(double r, const SBVonKarman::SBVonKarmanImpl& sbvki) :
            _r(r), _sbvki(sbvki)
        {}

        double operator()(double k) const { return std::real(_sbvki.kValue(k))*k*j0(k*_r); }
    private:
        double _r;
        const SBVonKarman::SBVonKarmanImpl& _sbvki;
    };

    double SBVonKarman::SBVonKarmanImpl::xValue(double r) const {
        VKXIntegrand I(r, *this);
        return integ::int1d(I, 0.0, integ::MOCK_INF,
                            gsparams->integration_relerr, gsparams->integration_abserr);
    }

    double SBVonKarman::SBVonKarmanImpl::xValue(const Position<double>& p) const
    // r in units of _scale
    {
        return xValue(sqrt(p.x*p.x+p.y*p.y)*_scale);
    }

    boost::shared_ptr<PhotonArray> SBVonKarman::SBVonKarmanImpl::shoot(
        int N, UniformDeviate ud) const
    {
        dbg<<"VonKarman shoot: N = "<<N<<std::endl;
        dbg<<"Target flux = "<<getFlux()<<std::endl;
        // Get photons from the VonKarmanInfo structure, rescale flux and size for this instance
        boost::shared_ptr<PhotonArray> result = _info->shoot(N,ud);
        result->scaleFlux(_flux);
        result->scaleXY(_scale);
        dbg<<"VonKarman Realized flux = "<<result->getTotalFlux()<<std::endl;
        return result;
    }

}
