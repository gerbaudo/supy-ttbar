
# GCCTAG=x86_64-slc5-gcc43-opt
# TAG=i686-slc5-gcc43-opt
# BASEDIR=/software/atlas/ifae/prod/releases
# REL=rel_16-4.9
# 
# source ${BASEDIR}/${REL}/AtlasSite/gcc-links/${GCCTAG}/setup.sh
# source ${BASEDIR}/${REL}/LCGCMT/LCGCMT_59a/InstallArea/${TAG}/bin/thisroot.sh
# 
# export PATH=${BASEDIR}/${REL}/AtlasSite/gcc-links/${GCCTAG}/bin:${PATH}
# export LD_LIBRARY_PATH=${BASEDIR}/${REL}/AtlasSite/gcc-links/${GCCTAG}/lib64:${LD_LIBRARY_PATH}
# 
# export PATH=${BASEDIR}/${REL}/sw/lcg/external/Python/2.6.5/${TAG}/bin:${PATH}
# export LD_LIBRARY_PATH=${BASEDIR}/${REL}/sw/lcg/external/Python/2.6.5/${TAG}/lib:${LD_LIBRARY_PATH}

# The instructions above for some reason do not work.
# Use the ones below in the meantime
. /software/at3/etc/profile.d/grid-env.sh
asetup 16.7.0