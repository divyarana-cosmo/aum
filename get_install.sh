#!/bin/zsh
export LDFLAGS="$(gsl-config --libs)"
export CPPFLAGS="$(gsl-config --cflags)"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}$(gsl-config --prefix)/lib"
pip install --prefix=/Users/divyarana/github/aum/install .
pip install --prefix=/Users/divyarana/github/aum/install .


#export LDFLAGS="-L`gsl-config --prefix`/lib"
#export CPPFLAGS="-I`gsl-config --prefix`/include"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:`gsl-config --prefix`/lib"
