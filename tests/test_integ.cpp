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

#include "galsim/integ/Int.h"

#define BOOST_TEST_DYN_LINK

#include "galsim/IgnoreWarnings.h"

#define BOOST_NO_CXX11_SMART_PTR
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>

const double test_sigma = 7.;           // test value of Gaussian sigma for integral tests
const double test_rel_err = 1.e-7;      // the relative accuracy at which to test
const double test_abs_err = 1.e-13;     // the absolute accuracy at which to test
const double test_mock_inf = 2.e10;     // number large enough to get interpreted as infinity by
                                        // integration routines
// Note: all "true" answers in the below tests are found using Wolfram Alpha.

// A simple Gaussian that works as a functional object
class Gauss : public std::unary_function<double,double>
{
public :
    Gauss(double sig) : _sig(sig) {}

    double operator()(double x) const
    { return exp(-0.5*pow(x/_sig,2)); }

private :
    double _sig;
};

// A simple power law
class Power : public std::unary_function<double,double>
{
public :
    Power(double expon) : _expon(expon) {}

    double operator()(double x) const
    { return pow(x,_expon); }

private :
    double _expon;
};

// A straight function, rather than a functional class:
double osc_func(double x)
{ return sin(pow(x,2)) * exp(-std::abs(x)); }

// A simple function:
// f(x,y) = x*(3*x+y) + y
double twod_func(double x, double y)
{ return x * (3.*x + y) + y; }

BOOST_AUTO_TEST_SUITE(integ_tests);

BOOST_AUTO_TEST_CASE( TestGaussian )
{
    Gauss gauss(test_sigma);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(gauss, -1., 1., test_rel_err, test_abs_err),
        1.99321805307377285009,
        // Note: it seems like this should be test_rel_err, but boost has this last
        // parameter as a _percentage_ error.  (!?!)  So need to multiply by 100.
        100 * test_rel_err); 

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(gauss, 0., 20., test_rel_err, test_abs_err),
        8.73569586966967345835,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(gauss, -50., -40., test_rel_err, test_abs_err),
        9.66426031085587421984e-8,
        // Difference is allowed to be test_abs_err, which corresponds to a relative error of
        // test_abs_err / 9.66e-8
        100 * test_abs_err / 9.66e-8);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(gauss, 0., test_mock_inf, test_rel_err, test_abs_err),
        8.77319896120850210849,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(gauss, -test_mock_inf, 5.4, test_rel_err, test_abs_err),
        13.68221660030048620971,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(gauss, -test_mock_inf, test_mock_inf, test_rel_err, test_abs_err),
        17.54639792241700421699,
        100 * test_rel_err);
}

BOOST_AUTO_TEST_CASE( TestOscillatory )
{
    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(std::ptr_fun(osc_func), -1., 1., test_rel_err, test_abs_err),
        0.30182513444548879567,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(std::ptr_fun(osc_func), 0., 20., test_rel_err, test_abs_err),
        0.27051358019041255485,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(std::ptr_fun(osc_func), -15., -14., test_rel_err, test_abs_err),
        7.81648378350593176887e-9,
        100 * test_abs_err / 7.82e-8);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(std::ptr_fun(osc_func), 0., test_mock_inf, test_rel_err,
                             test_abs_err),
        0.27051358016221414426,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(std::ptr_fun(osc_func), -test_mock_inf, 5.4, test_rel_err,
                             test_abs_err),
        0.5413229824941895221,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(std::ptr_fun(osc_func), -test_mock_inf, test_mock_inf, test_rel_err,
                             test_abs_err),
        0.54102716032442828852,
        100 * test_rel_err);
}

BOOST_AUTO_TEST_CASE( TestPole )
{
    Power powm05(-0.5);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(powm05, 0., 1., test_rel_err, test_abs_err),
        2.,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(powm05, 0., 300., test_rel_err, test_abs_err),
        34.64101615137754587055,
        100 * test_rel_err);

    Power powm2(-2.);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(powm2, 1., 2., test_rel_err, test_abs_err),
        0.5,
        100 * test_rel_err);

    BOOST_CHECK_CLOSE(
        galsim::integ::int1d(powm2, 1., test_mock_inf, test_rel_err, test_abs_err),
        1.,
        100 * test_rel_err);

#if 0
    // This works if everything uses the same compiler.  But boost testing may have
    // been installed with a different compiler, in which case it is unable to catch
    // exceptions thrown from GalSim.  So we skip this test normally.
    // (Developers working on the integrator should reenable this during development.)
    BOOST_CHECK_THROW(
        galsim::integ::int1d(powm2, 0., 1., test_rel_err, test_abs_err),
        galsim::integ::IntFailure);
#endif
}

BOOST_AUTO_TEST_CASE( Test2d )
{
    BOOST_CHECK_CLOSE(
        galsim::integ::int2d(std::ptr_fun(twod_func),0.,1.,0.,1., test_rel_err, test_abs_err),
        1.75,
        100 * test_rel_err);

    // This next one is an irrelevant test, but it suppresses a warning with gcc 4.2.1 that 
    // for some reason doesn't get suppressed with the pragma at the top of this file.
    BOOST_CHECK_SMALL(0., 0.1);
}

BOOST_AUTO_TEST_SUITE_END();

