[![Build Status](https://travis-ci.org/tpaviot/smesh.png?branch=master)](https://travis-ci.org/tpaviot/smesh)
[![Build status](https://ci.appveyor.com/api/projects/status/github/tpaviot/smesh?branch=master&svg=true)](https://ci.appveyor.com/project/tpaviot/smesh)

Description
-----------

A complete OpenCascade based MESH framework. Note this is not the original SALOME SMESH project but an effort to create a standalone mesh framework based on the existing one from SALOME project. This is a fork of the salomesmesh project available at:
http://sourceforge.net/projects/salomesmesh/ (Original project by Fotis Soutis). This fork is intended for a use in the pythonocc project.

We use the following online resources:
  * Sources
       https://github.com/tpaviot/smesh
  * Bug tracker
       https://github.com/tpaviot/smesh/issues
  * Mailing list
       http://groups.google.com/group/smesh-dev/about
  * Travic-CI
       https://travis-ci.org/tpaviot/smesh
  * Appveyor
       https://ci.appveyor.com/project/tpaviot/smesh

Just ask @tpaviot (tpaviot@gmail.com) for a request regarding write access to the repository.

How to create a local copy of the repository?
---------------------------------------------

    git clone git://github.com/tpaviot/smesh.git


How to stay up to date with latest developements?
-------------------------------------------------

    cd smesh
    git pull

Install with conda
------------------

    conda install -c pythonocc -c dlr-sc -c oce smesh=6.7.6

Build - Install
---------------

For both OSX, Linux and Windows, the instructions are the same.

Requirements
------------
  * a c++ and a fortran compiler 

  * cmake 2.8 or higher

  * oce 0.18.x

Build
-----

    $ mkdir cmake-build
    $ cd cmake-build
    $ cmake ..
    $ make
    $ make install

Build with conda
----------------

    $ conda build ci/conda
