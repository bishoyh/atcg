#!/usr/bin/env python
# Copyright (c) 2007 Shuzhao Li.
#
# Licensed under the MIT License.
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
Compute DNA/DNA and RNA/DNA binding energy.
by Nearest Neighbor Model (see SantaLucia 1998).

Three sets of parameters are included from:

DNA/DNA parameters, SantaLucia Jr. 1998. PNAS 95:1460.
RNA/RNA parameters, Wu et al., 2002. Eur. J. Biochem. 269:2821.
Sugimoto et al. 1995. Biochemistry 34:11211.


####
Table Wu et al., 2002. Eur. J. Biochem. 269:2821.
seq delta_G@37 (kcal/mol)
rAA/dTT    -1    -0.4    
rAC/dTG    -2.1    -1.6
rAG/dTC    -1.8    -1.4
rAU/dTA    -0.9    -1    
rCA/dGT    -0.9    -1    
rCC/dGG    -2.1    -1.5
rCG/dGC    -1.7    -1.2
rCU/dGA    -0.9    -0.9
rGA/dCT    -1.3    -1.4
rGC/dCG    -2.7    -2.4
rGG/dCC    -2.9    -2.2
rGU/dCA    -1.1    -1.5
rUA/dAT    -0.6    -0.3
rUC/dAG    -1.5    -0.8
rUG/dAC    -1.6    -1
rUU/dAA    -0.2    -0.2
Initiation    3.1    1
"""

#
# rd95, rd02 RNA/DNA NN stacking energies;
# Wu et al., 2002
#
rd95={'TT':-1.0,
      'TG':-2.1,
      'TC':-1.8,
      'TA':-0.9,
      'GT':-0.9,
      'GG':-2.1,
      'GC':-1.7,
      'GA':-0.9,
      'CT':-1.3,
      'CG':-2.7,
      'CC':-2.9,
      'CA':-1.1,
      'AT':-0.6,
      'AG':-1.5,
      'AC':-1.6,
      'AA':-0.2,
      'init':3.1
    }

rd02={'TT':-0.4,
      'TG':-1.6,
      'TC':-1.4,
      'TA':-1.0,
      'GT':-1.0,
      'GG':-1.5,
      'GC':-1.2,
      'GA':-0.9,
      'CT':-1.4,
      'CG':-2.4,
      'CC':-2.2,
      'CA':-1.5,
      'AT':-0.3,
      'AG':-0.8,
      'AC':-1.0,
      'AA':-0.2,
      'init':1.0
    }


# 
# DNA/DNA parameters, SantaLucia Jr. 1998.
# init = 1, not distinguishing 0.98/1.03
#
dd98={'TT':-1.00,
      'TG':-1.45,
      'TC':-1.30,
      'TA':-0.58,
      'GT':-1.44,
      'GG':-1.84,
      'GC':-2.24,
      'GA':-1.30,
      'CT':-1.28,
      'CG':-2.17,
      'CC':-1.84,
      'CA':-1.45,
      'AT':-0.88,
      'AG':-1.28,
      'AC':-1.44,
      'AA':-1.00,
      'init':1.0
    }



def seqcheck(s):
    q = 1
    for ch in s:
        if ch not in ['A', 'T', 'C', 'G']:
            q -= 1
    return q

def seq2g(s, rd):
    """
    take a 5'-3' DNA probe sequence, return delta_G_37 (kcal/mol)
    """
    if seqcheck(s) == 1:
        g = rd['init']
        for i in range(len(s)-1):
            g += rd[s[-1-i]+s[-2-i]]
        return g
    else:
        return 'Illegal BASE!'

