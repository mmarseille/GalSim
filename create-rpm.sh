#!/bin/bash

appliance=GalSim
version=1.5.1
patch_version=1

# Patch name
appliance_name=${appliance}-${version}

# Making the DIFF for the patch
echo "Creating diff ..."
diff -ru --exclude '.*' ${appliance_name} ${appliance_name}-patch > $HOME/rpmbuild/SOURCES/${appliance_name}-${patch_version}-euclid.patch

# Creating source tarball file
echo "Creating tarball files ..."
source_tarball=${appliance_name}.tar.gz
cd ${appliance_name}
rm -rf .git
cd ..
tar czf ${source_tarball} ${appliance_name}

echo "Moving tarball to RPM SOURCE folder ..."
cp ${source_tarball} ~/rpmbuild/SOURCES/.

# Going for rpmbuild
echo "Building SRPM file"
rpmbuild -bs --define '_prefix %{getenv:EUCLID_CUSTOM_PREFIX}' GalSim.spec

echo "Building RPM from specfile"
rpmbuild --define '_prefix %{getenv:EUCLID_CUSTOM_PREFIX}' --rebuild ~/rpmbuild/SRPMS/GalSim-1.5.1-1.el7.centos.src.rpm
