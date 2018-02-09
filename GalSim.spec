Name:		GalSim
Version:	1.5.1
Release:	1%{?dist}
Summary:	GalSim Euclid specific RPM deployment package with Python 3 full support
Group:		Development/Libraries
License:	GPLv2
Source0:	%{name}-%{version}.tar.gz

BuildRequires:	scons
BuildRequires:	TMV-devel >= 0.72
BuildREquires:	fftw-devel
BuildRequires:	gcc-c++
BuildRequires:	boost
Requires:	boost-python
Requires:	python-astropy
Requires:	numpy

Patch1:  GalSim-1.5.1-1-euclid.patch

# Some macros disabled
%define debug_package %{nil}
%define __arch_install_post %{nil}
%define __os_install_post %{nil}


%description
GalSim is open-source software for simulating images of astronomical objects (stars, galaxies) in a variety of ways. The bulk of the calculations are carried out in C++, and the user interface is in python. In addition, the code can operate directly on "config" files, for those users who prefer not to work in python. The impetus for the software package was a weak lensing community data challenge, called GREAT3:

http://great3challenge.info/

However, the code has numerous additional capabilities beyond those needed for the challenge, and has been useful for a number of projects that needed to simulate high-fidelity galaxy images with accurate sizes and shears.

%prep
%setup -q
# Applying patch
cd ..
%patch1 -p0
cd %{name}-%{version}
# Must be done for clean GalSim compilation
scons -c
rm -rf .scons*

%build
# With python 3 support from cvmfs
scons PREFIX=%{_prefix} PYTHON=%{_prefix}/bin/python3 PYPREFIX=%{_prefix}/lib/python3.6/site-packages/

%install
scons install PREFIX=%{buildroot}%{_prefix}
# prefix not used for python script install : moving them manually to BUILDROOT+PREFIX path
mkdir -p %{buildroot}%{_prefix}/lib/python3.6/site-packages
cp -r %{_prefix}/lib/python3.6/site-packages/galsim %{buildroot}%{_prefix}/lib/python3.6/site-packages/.
# Necessary for check build without BUILDROOT in binaries
sed -i "s|%{buildroot}||g" %{buildroot}%{_prefix}/lib/python3.6/site-packages/galsim/meta_data.py

%check
# scons tests PREFIX=%{_prefix}

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

# BUILDROOT is automatically prepended to files path
%files
%{_prefix}/bin/galsim
%{_prefix}/bin/galsim_download_cosmos
%{_bindir}/galsim_yaml
%{_prefix}/bin/galsim_json
%{_prefix}/lib/python3.6/site-packages/galsim
%{_prefix}/lib/libgalsim.so*
%{_prefix}/share/galsim

%clean
rm -rf $RPM_BUILD_ROOT

%changelog
* Fri Feb 9 2018 Matthieu Marseille <matthieu.marseille@thales-services.fr> 1.5.1-1
- Initial release
