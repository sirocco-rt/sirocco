#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Routine to estimate the size of the CommBuffer


Command line usage (if any):

    usage: GetCommBufferSize.py filename

    where filename is the name of a text file that contains only the  MPI_PACK commands one wished to size 

Description:  

Primary routines:

    doit

Notes:
                                       
History:

250927 ksl Coding begun

'''


# # Get the Combuffer Size for Broadcasting a set of data
# 
# This is a routine to extract the CommBuffer size from a subroutine intended to share variables with into other threads using MPI routines.  The routine is specifically designed to handle size the comm_buffer needed for plasmamain.
# 
# Before using this routine one needs to create a text file that contains the MPI_PAC

# In[21]:



from astropy.table import Table,join
import numpy as np


import re
from pathlib import Path

_control_keywords = {'if', 'for', 'while', 'switch', 'else', 'do', 'return',
                     'sizeof', 'case', 'default', 'catch', 'goto'}

def remove_comments_preserve_strings(code: str) -> str:
    """Remove // and /* */ comments but preserve string/char literals."""
    out = []
    i = 0
    n = len(code)
    while i < n:
        c = code[i]
        if c == '"' or c == "'":
            quote = c
            out.append(c)
            i += 1
            while i < n:
                out.append(code[i])
                if code[i] == '\\':
                    i += 1
                    if i < n:
                        out.append(code[i])
                elif code[i] == quote:
                    i += 1
                    break
                i += 1
            continue
        if c == '/' and i + 1 < n and code[i+1] == '/':
            i += 2
            while i < n and code[i] != '\n':
                i += 1
            if i < n and code[i] == '\n':
                out.append('\n')
                i += 1
            continue
        if c == '/' and i + 1 < n and code[i+1] == '*':
            i += 2
            while i + 1 < n and not (code[i] == '*' and code[i+1] == '/'):
                if code[i] == '\n':
                    out.append('\n')
                i += 1
            if i + 1 < n:
                i += 2
            continue
        out.append(c)
        i += 1
    return ''.join(out)

def build_paren_brace_maps(code: str):
    """Return maps for parentheses and braces and brace depth for each opening brace."""
    paren_open_to_close = {}
    paren_close_to_open = {}
    brace_open_to_close = {}
    brace_close_to_open = {}
    brace_depth = {}

    paren_stack = []
    brace_stack = []

    i = 0
    n = len(code)
    while i < n:
        c = code[i]
        if c == '"' or c == "'":
            quote = c
            i += 1
            while i < n:
                if code[i] == '\\':
                    i += 2
                    continue
                if code[i] == quote:
                    i += 1
                    break
                i += 1
            continue
        if c == '(':
            paren_stack.append(i)
        elif c == ')':
            if paren_stack:
                o = paren_stack.pop()
                paren_open_to_close[o] = i
                paren_close_to_open[i] = o
        elif c == '{':
            depth = len(brace_stack)
            brace_stack.append(i)
            brace_depth[i] = depth
        elif c == '}':
            if brace_stack:
                o = brace_stack.pop()
                brace_open_to_close[o] = i
                brace_close_to_open[i] = o
        i += 1

    return paren_open_to_close, paren_close_to_open, brace_open_to_close, brace_close_to_open, brace_depth

def extract_name_before_paren(code: str, open_paren_pos: int):
    """Find identifier (function name) immediately before '('; return None if none."""
    j = open_paren_pos - 1
    while j >= 0 and code[j].isspace():
        j -= 1
    end = j
    while j >= 0 and (code[j].isalnum() or code[j] == '_'):
        j -= 1
    start = j + 1
    name = code[start:end+1].strip()
    if not name:
        return None
    if name in _control_keywords:
        return None
    if name.startswith('__attribute') or name == 'attribute' or name == '__declspec':
        return None
    return name

def separate(filename: str):
    """
    Return (functions, names) for filename.
     - functions: list of function source strings (comments removed)
     - names: list of function names (parallel to functions)
    """
    txt = Path(filename).read_text(encoding='utf-8', errors='ignore')
    code = remove_comments_preserve_strings(txt)
    po2c, pc2o, bo2c, bc2o, brace_depth = build_paren_brace_maps(code)

    funcs = []
    names = []

    for b_open in sorted(bo2c.keys()):
        if brace_depth.get(b_open, 0) != 0:
            continue
        b_close = bo2c[b_open]

        # candidates: open-parens whose close is before this brace
        candidates = [(o, po2c[o]) for o in po2c.keys() if po2c[o] < b_open]
        candidates.sort(key=lambda x: x[0], reverse=True)  # rightmost first

        found = False
        for open_paren_pos, close_paren_pos in candidates:
            name = extract_name_before_paren(code, open_paren_pos)
            if not name:
                continue
            if close_paren_pos < b_open:
                # capture from start-of-signature-line to matching brace
                sig_start = code.rfind('\n', 0, open_paren_pos)
                if sig_start == -1:
                    sig_start = 0
                func_text = code[sig_start:b_close+1].strip()
                funcs.append(func_text)
                names.append(name)
                found = True
                break

        if not found:
            # fallback: capture block and try a regex
            sig_start = code.rfind('\n', 0, b_open)
            if sig_start == -1:
                sig_start = 0
            func_text = code[sig_start:b_close+1].strip()
            funcs.append(func_text)
            m = re.search(r'([A-Za-z_]\w*)\s*\([^)]{0,200}\)\s*$', code[:b_open+1], flags=re.DOTALL)
            names.append(m.group(1) if m else None)

    return funcs, names



def parse(lines):
# Read the file containg the Pack commands
    xvar=[]
    xtype=[]
    xsize=[]
    for one in lines:
        if one.count('MPI_Pack'):
            word=one.split(',')
            # print(word)
            xxvar=word[0].split('>')[-1]
            xxsize=word[1]
            xxtype=word[2]
            # print(xxvar,xxsize,xxtype)
            xvar.append(xxvar.strip())
            xtype.append(xxtype.strip())
            xsize.append(xxsize.strip())
    ctab=Table([xvar,xtype,xsize],names=['Variable','Type','Dimension'])
    return ctab



def xmake_string(xsum):
    name,count=np.unique(xsum['Dimension'],return_counts=True)
    xstring=''
    i=0
    while i<len(name):
        # print('%20s %s' % (name[i],count[i]))
        if name[i]=='1':
            ntot=count[i]
        else:
            xstring+=' + %d * (%s)' % (count[i],name[i])
        i+=1
    xxstring='%d % s' % (ntot,xstring)
    return xxstring



    
warning='''
*******************WARNING********************

These are the nominal answers that were for commitn 864b30fdc4ad02a843240991cf04396d913a455

integers

(1 + 20 + nphot_total + nions + NXBANDS + 2 * N_PHOT_PROC),

doubles
 
(73 + 11 * nions + nlte_levels + 2 * nphot_total + n_inner_tot +11 * NXBANDS + NBINS_IN_CELL_SPEC + 6 * NFLUX_ANGLES +N_DMO_DT_DIRECTIONS + 12 * NFORCE_DIRECTIONS));


For this the e current routine produced (on 250927):

Space for integers:

 24  + 1 * (NXBANDS) + 2 * (N_PHOT_PROC) + 1 * (nions) + 1 * (nphot_total)

Space for doubles:

 73  + 1 * (NBINS_IN_CELL_SPEC) + 6 * (NFLUX_ANGLES) + 12 * (NFORCE_DIRECTIONS) + 11 * (NXBANDS) + 1 * (N_DMO_DT_DIRECTIONS) + 1 * (n_inner_tot) + 11 * (nions) + 1 * (nlte_levels) + 2 * (nphot_total)

 These are identical except for the number of integers, where this program procuces a slightly larger integer space, due to a few more single integers., 

 One should be very carful about using the routine blindly. and it depends in detail on the exact way lines for MPI_PAC are written. 

Proably the best apprach is to use it in situtations where one is adding or subtracting variables from plasmamain, which one wishes to compare. In that case, comparing
the results from two different versions of Sirocco, one which one knows is correct, and one which might not be correct is best.

For Plasmamain one has to be carefule to note that the size given here is for only element of Plasmamain


 Note also this routine is currently only set up to work with  Plasmamain

*******************

'''


def doit(filename='plasma.c'):
    print(warning)
    funcs,names=separate(filename)
    i=0
    while i<len(names):
        print('\n\nStarting %s in %s' % (names[i],filename))
        lines=funcs[i].split('\n')
        ctab=parse(lines)
        print('For %s found %s MPI_Pack lines' % (names[i],len(ctab)))
        if len(ctab)>0:
              

            ctab.write('xall_%s.txt' % names[i], format='ascii.fixed_width_two_line',overwrite=True)
            xint=ctab[ctab['Type']=='MPI_INT']
            xint.write('xint_%s.txt' % names[i], format='ascii.fixed_width_two_line',overwrite=True)
            xdoub=ctab[ctab['Type']=='MPI_DOUBLE']
            xdoub.write('xdoub_%s.txt' % names[i], format='ascii.fixed_width_two_line',overwrite=True)
            int_string= xmake_string(xint)
            print('Space for integers:\n', int_string)
            doub_string=xmake_string(xdoub)
            print('Space for doubles:\n', doub_string)
        else:
            print('There is nothing to do')
        print('All done for %s' % names[i])
        i+=1
    




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print (__doc__)
