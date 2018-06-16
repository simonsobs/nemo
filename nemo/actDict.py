# -*- coding: iso-8859-1 -*-
import os
import string

def ask_for( key ):
    s = input( "ACTDict: enter value for '%s': " % key )
    try:
        val = eval(s)
    except NameError:
        # allow people to enter unquoted strings
        val = s
    return val

class ACTDict( dict ):

    def __getitem__( self, key ):
        if key not in self:
            print("ACTDict: parameter '%s' not found" % key)
            val = ask_for( key )
            print("ACTDict: setting '%s' = %s" % (key,repr(val)))
            dict.__setitem__( self, key, val )
        return dict.__getitem__( self, key )

    def read_from_file( self, filename ):
        f = open( filename )
        old = ''
        for line in f:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            s = line.split('#')
            line = s[0]
            s = line.split('\\')
            if len(s) > 1:
                old = string.join([old, s[0]])
                continue
            else:
                line = string.join([old, s[0]])
                old = ''
            s = line.split('=')
            if len(s) != 2:
                print("Error parsing line:")
                print(line)
                continue
            try:
                key = s[0].strip()
                val = eval(s[1].strip()) # XXX:make safer
            except:
                raise Exception("can't parse line: %s" % (line))
            self[key] = val
        f.close()

    def write_to_file( self, filename, mode = 'w' ):
        f = open( filename, mode )
        keys = list(self.keys())
        keys.sort()
        for key in keys:
            f.write( "%s = %s\n" % (key,repr(self[key])) )
        f.close()

    def cmp( self, otherDict ):
        
        diff = []
        ks = list(self.keys())
        for k in ks:
            try:
                if otherDict[k] == self.params[k]:
                    continue
                diff += [k]
                break
            except KeyError:
                diff += [k]
        return otherDict
