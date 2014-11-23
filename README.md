
[![Build Status](https://travis-ci.org/tpaviot/smesh.png?branch=master)](https://travis-ci.org/tpaviot/smesh)

Description
-----------

A complete OpenCascade based MESH framework.Note this is not the original SALOME SMESH project but an effort to create a standalone mesh framework based on the existing one from SALOME project.

This is a fork of the salomesmesh project available at:
http://sourceforge.net/projects/salomesmesh/

Original project by Fotis Soutis.

This fork is intended for a use in the pythonocc project.

Build - Install
---------------

For both OSX, Linux and Windows (please use Mingw on Windows), the instructions are the same.

Requirements
------------
  * a c++ and a fortran compiler 

  * cmake 2.6 or higher

  * oce 0.16

Build
-----

  $ mkdir cmake-build

  $ cd cmake-build

  $ make

  $ make install
