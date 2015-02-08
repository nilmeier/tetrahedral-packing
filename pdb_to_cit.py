# this will convert a pdb format to the citrine format as given here:
# http://www.citrine.io/tech-challenge
#
import sys
import re

pdbname=sys.argv[1]
outfile=re.sub('.pdb','.cit',pdbname) 

# basis1x basis1y basis1z
# basis2x basis2y basis2z
# basis3x basis3y basis3z
# tet1vert1x tet1vert1y tet1vert1z tet1vert2x tet1vert2y tet1vert2z ... tet1vert4z
# tet2vert1x ...

fin=open(pdbname,'r')
fout=open(outfile,'w')
outline=''
for inline in fin:
    if 'CRYST1' in inline:
        bx = inline.split()[1]; by=inline.split()[2]; bz=inline.split()[3]
        #bx = '%9.3f'%(float(bx)); #first column not padded 
        by='%9.3f'%(float(by)); 
        bz='%9.3f'%(float(bz)); 
        z   = '%3.3f'%(0.0) #not padded
        zp  = '%9.3f'%(0.0)  #padded
        fout.write(bx+ zp + zp +'\n')
        fout.write(z + by + zp +'\n')
        fout.write(z + zp + bz +'\n')
 
    if 'HETATM' in inline:
        if 'H1' in inline:
            outline+=inline.split()[5]+' '+inline.split()[6]+' '+inline.split()[7]+' '
        if '4H1' in inline:
            outline+='\n'
            fout.write(outline)
            outline=''
fout.close()
        

#JN DEBUG LINE ================================================
#import rlcompleter; import readline; readline.parse_and_bind("tab: complete")
#print "jndb\n"; import code; code.interact(local=locals())
#JN DEBUG LINE =================================================

