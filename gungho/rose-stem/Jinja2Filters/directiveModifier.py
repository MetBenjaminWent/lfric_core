#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from __future__ import print_function

import math
import re

    

def directiveModifier( directive, cores, walltime ):

    def choose_replacement(name, arguments):
        if name == 'nodes':
            nodeSize = int( arguments[0] )
            middlebit = str(int(math.ceil( float(cores) / float(nodeSize) )))
        elif name == 'cores':
            middlebit = str(cores)
        elif name == 'time_hhmmss':
            middlebit = str(walltime)
        else:
            raise Exception( 'Unrecognised function name "{}"'.format( name ) )
        return middlebit

    if type(directive) == type('string'):
        directive=[directive]

        
    pattern = re.compile( '<<(\S+?)\((\S*?)\)>>' )

    plist=[]
    for directive_i in directive:
        processed = directive_i
        for match in pattern.finditer( directive_i ):
            name = match.group(1)
            arguments = match.group(2).split( ',' )
            match_string= '<<%s(%s)>>' % (match.group(1), match.group(2))
            processed = processed.replace(match_string, choose_replacement(name, arguments))
        plist.append(processed)

    return plist

if __name__=='__main__':
    TARGET_RUN_DIRECTIVES = '-l=select=<<nodes(36)>>,walltime=<<time_hhmmss()>>'
    print( directiveModifier(TARGET_RUN_DIRECTIVES, 122, '00:05:00'))
