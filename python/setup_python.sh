#!/bin/sh

NTHREADS=20

setup_python() {
    URL=$1
    base=`basename $URL`
    pkg=`echo $base | sed 's/\.tgz$//'`
    PYTHON=$2
    env_nm=$3
    if [ ! -d $pkg ] ; then
        wget -O $base $URL
        tar xvfz $base
        rm -f $base
    fi

    if [ ! -f $pkg/install-root/bin/$PYTHON ] ; then
        cd $pkg
        if [ ! -f Makefile ] ; then
            ./configure --prefix=`pwd`/install-root
        fi
        if [ ! -f install-root/bin/$PYTHON ] ; then
            make -j $NTHREADS
            make install
        fi
        cd ..
    fi
    
    VE=$pkg/install-root/bin/pyvenv
    if [ ! -f $env_nm ] ; then
        $VE $env_nm
    fi
    
    source $env_nm/bin/activate
    
    if ! pip list | grep numpy ; then
        git clone git@bitbucket.org:pypy/numpy.git
        cd numpy && ../$pkg/install-root/bin/$PYTHON setup.py install
        rm -rf numpy
    fi
    
    pip list
    deactivate
}

setup_python https://www.python.org/ftp/python/3.4.3/Python-3.4.3.tgz python3.4 python3_qsim_env

