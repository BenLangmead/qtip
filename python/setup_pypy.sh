#!/bin/sh

setup_pypy() {
    URL=$1
    base=`basename $URL`
    pkg=`echo $URL | sed 's/\.tar\.bz2$//'`
    env_nm=$2
    if [ ! -d $pkg ] ; then
        wget -O $base $URL
        tar xvfj $base
        rm -rf $base
    fi
    
    VE=$pkg/bin/virtualenv-pypy
    if [ ! -f $env_nm ] ; then
        $VE $env_nm
    fi
    
    source $env_nm/bin/activate
    
    if ! pip list | grep numpy ; then
        git clone git@bitbucket.org:pypy/numpy.git
        cd numpy && pypy setup.py install
    fi
    
    pip list
    deactivate
}

setup_pypy https://bitbucket.org/squeaky/portable-pypy/downloads/pypy-2.6-linux_x86_64-portable.tar.bz2 pypy2_qsim_env
setup_pypy https://bitbucket.org/squeaky/portable-pypy/downloads/pypy3-2.4-linux_x86_64-portable.tar.bz2 pypy3_qsim_env
