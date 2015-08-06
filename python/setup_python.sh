#!/bin/sh

NTHREADS=20

setup_python_3() {
    URL=$1
    base=`basename $URL`
    pkg=`echo $base | sed 's/\.tgz$//'`
    PYTHON=$2
    env_nm=$3

    if [ ! -d $pkg ] ; then
        echo
        echo "=== Getting package ==="
        echo
        wget -O $base $URL
        tar xvfz $base
        rm -f $base
    fi

    if [ ! -f $pkg/install-root/bin/$PYTHON ] ; then
        cd $pkg
        if [ ! -f Makefile ] ; then
            echo
            echo "=== Configuring Python build ==="
            echo
            ./configure --prefix=`pwd`/install-root
        fi
        if [ ! -f install-root/bin/$PYTHON ] ; then
            echo
            echo "=== Making Python ==="
            echo
            make -j $NTHREADS

            echo
            echo "=== Installing Python ==="
            echo
            make install
        fi
        cd ..
    fi
    
    VE=$pkg/install-root/bin/pyvenv
    if [ ! -d $env_nm ] ; then
        echo
        echo "=== Making virtual environment ==="
        echo
        $VE $env_nm
    fi
    
    source $env_nm/bin/activate
    
    for pk in numpy scipy scikit-learn ; do
        if ! pip list | grep $pk ; then
            echo
            echo "=== Installing $pk with pip ==="
            echo
            pip --trusted-host pypi.python.org install $pk
        fi
    done

    pip list
    deactivate
}

setup_python_2() {
    URL=$1
    base=`basename $URL`
    pkg=`echo $base | sed 's/\.tgz$//'`
    PYTHON=$2
    env_nm=$3
    if [ ! -d $pkg ] ; then
        echo
        echo "=== Getting package ==="
        echo
        wget -O $base $URL
        tar xvfz $base
        rm -f $base
    fi

    if [ ! -f $pkg/install-root/bin/$PYTHON ] ; then
        cd $pkg
        if [ ! -f Makefile ] ; then
            echo
            echo "=== Configuring Python build ==="
            echo
            ./configure --prefix=`pwd`/install-root
        fi
        if [ ! -f install-root/bin/$PYTHON ] ; then
            echo
            echo "=== Making Python ==="
            echo
            make -j $NTHREADS

            echo
            echo "=== Installing Python ==="
            echo
            make install
        fi
        cd ..
    fi
    
    # install pip
    if [ ! -f $pkg/install-root/bin/pip ] ; then
        echo
        echo "=== Installing pip ==="
        echo
        wget --no-check https://bootstrap.pypa.io/get-pip.py
        $pkg/install-root/bin/$PYTHON get-pip.py
        rm -f get-pip.py
    fi

    # install virtualenv
    if [ ! -f $pkg/install-root/bin/virtualenv ] ; then
        echo
        echo "=== Installing virtualenv ==="
        echo
        $pkg/install-root/bin/pip install virtualenv
    fi

    VE=$pkg/install-root/bin/virtualenv
    if [ ! -d $env_nm ] ; then
        echo
        echo "=== Making virtual environment ==="
        echo
        $VE $env_nm
    fi
    
    source $env_nm/bin/activate
    
    for pk in numpy scipy scikit-learn ; do
        if ! pip list | grep $pk ; then
            echo
            echo "=== Installing $pk with pip ==="
            echo
            pip --trusted-host pypi.python.org install $pk
        fi
    done
    
    pip list
    deactivate
}

setup_python_3 https://www.python.org/ftp/python/3.4.3/Python-3.4.3.tgz python3.4 python3_qsim_env
setup_python_2 https://www.python.org/ftp/python/2.7.10/Python-2.7.10.tgz python2.7 python2_qsim_env
