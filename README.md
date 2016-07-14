
[![License](https://binstar.org/jf/smesh/badges/license.svg)](https://github.com/tpaviot/smesh/blob/master/LICENCE.lgpl.txt)
[![Build Status](https://travis-ci.org/tpaviot/smesh.png?branch=master)](https://travis-ci.org/tpaviot/smesh)
[![Conda installer](https://binstar.org/jf/smesh/badges/installer/conda.svg)](https://binstar.org/jf/smesh/)

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

Just ask @tpaviot (tpaviot@gmail.com) for a request regarding write access to the repository.

How to create a local copy of the repository?
---------------------------------------------

    git clone git://github.com/tpaviot/smesh.git


How to stay up to date with latest developements?
-------------------------------------------------

    cd smesh
    git pull

Build - Install
---------------

For both OSX, Linux and Windows (please use Mingw on Windows), the instructions are the same.

Requirements
------------
  * a c++ and a fortran compiler 

  * cmake 2.8 or higher

  * oce 0.17.1 or oce-0.17.2

Build
-----

    $ mkdir cmake-build
    $ cd cmake-build
    $ make
    $ make install
