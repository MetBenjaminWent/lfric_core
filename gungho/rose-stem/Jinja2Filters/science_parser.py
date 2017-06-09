#!/usr/bin/env python
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a parser for the science mission command line
'''
def science_parser(s):
    """
    Split `s` by commas where those in parentheses, i.e. (),[],{}, are ignored.
    """

    # Parse the string tracking whether the current character is within
    # parentheses.
    balance = 0
    parts = []
    part = ''

    def c_comp(s):
        ' return complementary string to s '
        r=None
        if s=='(':r=')'
        if s=='{':r='}'
        if s=='[':r=']'
        return r

    for c in s:
        part += c
        if balance == 0: bc=None
        if c in ['(','{','['] and bc in [None,c]:
            balance += 1
            bc=c
        elif c == c_comp(bc):
            balance -= 1
        elif c == ',' and balance == 0:
            parts.append(part[:-1].strip())
            part = ''

    # Capture last part
    if len(part):
        parts.append(part.strip())

    return parts
