#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Process previously analysed dependency database. For fun and profit!

from __future__ import print_function;

import os.path

###############################################################################
# Process dependency database.
#
class Processor():
    ###########################################################################
    # Constructor.
    #
    # Arguments:
    #   logger    - A Logger object to write output to.
    #   database  - The database containing dependency information.
    #   objectdir - The directory which holds .o files.
    #   moduledir - The directory which holds .mod files.
    #
    def __init__( self, logger, database, objectDirectory, moduleDirectory):
        self._logger          = logger
        self._database        = database
        self._objectDirectory = objectDirectory
        self._moduleDirectory = moduleDirectory

    ###########################################################################
    # Generate a Makefile snippet containing dependency rules for all objects.
    #
    # Arguments:
    #   output - File object into which to write the snippet.
    #
    def dependencyRules( self, output ):
        print( '# Object dependencies', file=output )

        for dependor, dependees in self._database.getAllFileDependencies():
            objectFile = os.path.join( self._objectDirectory, \
                                       self._replaceExtension( dependor, 'o' ) )
            print('#', file=output)

            moduleFiles = [os.path.join( self._moduleDirectory, \
                                         self._replaceExtension( dependee, \
                                                                 'mod' ) ) \
                           for dependee in dependees]
            print( '{} {}: {} {}'.format( objectFile, \
                                          self._replaceExtension( objectFile, \
                                                                  'mod' ), \
                                          dependor, \
                                          ' '.join( moduleFiles ) ), \
                   file=output )

    ###########################################################################
    # Generate a Makefile snippet containing a macro listing all objects
    # required to build each program.
    #
    # Arguments:
    #   output - File object into which to write the snippet.
    #
    def programObjects( self, output ):
        print('# Program object list', file=output)
        print('#', file=output)

        for (program, sources) in self._database.getProgramSourceList():
            objectPaths = [self._replaceExtension( os.path.join( self._objectDirectory, \
                                                                 source), \
                                                   'o' ) \
                           for source in sources]
            print( '{}_OBJS = {}'.format( program.upper(), \
                                          ' '.join( objectPaths ) ), \
                   file=output )

        programObject = []
        for name, filename in self._database.getProgramSources().iteritems():
            objectFilename = os.path.join( self._objectDirectory, \
                                           self._replaceExtension( filename, \
                                                                   'o' ) )
            programObject.append( objectFilename )

        print( 'PROG_OBJS = {}'.format( ' '.join( programObject ) ), \
               file=output )

    ###########################################################################
    # Recursively descend the dependency tree building up a list of all the
    # source files the specified source file depends on.
    #
    # Arguments:
    #   filename    - The file to find dependencies for.
    #   sourceFiles - A list to receive the files depended upon.
    #
    def _descendDependencyTree(self, filename, sourceFiles):
        sourceFiles.append(filename)
        if filename in self._filenameToDependencies:
            for module in self._filenameToDependencies[filename]:
                if self._moduleToFilename[module] not in sourceFiles:
                    self._descendDependencyTree(self._moduleToFilename[module], \
                                               sourceFiles)

    ###########################################################################
    # Replace a filename's extension with another.
    #
    # Arguments:
    #   filename  - The filename to replace the extension of.
    #   extension - The new extension. Sans '.'.
    #
    def _replaceExtension( self, filename, extension ):
        (base, rubbish) = os.path.splitext( filename )
        return '{}.{}'.format( base, extension)
