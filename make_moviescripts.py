import sys
import re

print "jn don't forget to add options"
print "wanna print volume too?"
dirname = sys.argv[1]
clustnum =  sys.argv[2] 
trialnum =  sys.argv[3]
stepnum  =  sys.argv[4]

class ChimeraScript:
    def __init__(self, name, lines=[]):
        self.name=name
        self.lines=lines
    def loadfile(self):
        fin=open('metascripts/'+self.name+'.meta.cmd','r')
        self.lines=fin.read()
        fin.close()
    def loadfile_byline(self):
        fin=open('metascripts/'+self.name+'.meta.cmd','r')
        for inline in  fin:
            self.lines.append( inline )
        print self.lines 
        fin.close()
    def writefile(self):
        fout=open(self.name+'.cmd','w')
        fout.write(self.lines)
        fout.close()

#===

def parse_runmovie():
    header=[]; s1_sn=[]; stepping=[]; footer=[]
    header_flag=False; s1_sn_flag=False; stepping_flag=False; footer_flag=False
    for line in runmovie.lines:
        print "==>"+ line
        if 'start initialization' in line:
            header_flag = True
        if header_flag:
            header.append(line)
        if 'end initialization' in line:
            header_flag = False
        if 'start showing steps' in line:
            s1_sn_flag = True
        if s1_sn_flag:
           s1_sn.append(line)
        if 'end showing steps' in line:
            s1_sn_flag = False

        if 'start stepping' in line:
            stepping_flag=True  
        if stepping_flag:
            stepping.append(line)
        if 'end stepping' in line:
            stepping_flag=False  
            
        if 'start footer' in line:
            footer_flag=True  
        if footer_flag:
            footer.append(line)
        if 'end footer' in line:
            footer_flag=True  
    return header,s1_sn,stepping,footer 

moviemain       = ChimeraScript('moviemain')
assignpositions = ChimeraScript('assignpositions')
assigntriangles = ChimeraScript('assigntriangles')
runmovie        = ChimeraScript('runmovie')

# moviemain processing ========
moviemain.loadfile()
# processing .meta.cmd file explicitly here
moviemain.lines=re.sub('CLUSTNUM', clustnum,moviemain.lines)
moviemain.lines=re.sub('STEPNUM' , stepnum,moviemain.lines)
moviemain.lines=re.sub('DIRNAME' , dirname,moviemain.lines)
moviemain.lines=re.sub('TRIALNUM', trialnum,moviemain.lines)

moviemain.writefile()

# assigntriangles processing ==
assigntriangles.loadfile()

trisnippet=assigntriangles.lines
assigntriangles.lines=re.sub('Sp50','51',trisnippet)
assigntriangles.lines=re.sub('STEPNUM','1',trisnippet)
assigntriangles.lines=''
for i in range(0 ,int(stepnum)+1 ):
    I=str(i+1)  #starting at model #50
    Ip50=str(i+1+50)
    tmp=re.sub('Sp50',Ip50,trisnippet)

#    assigntriangles.lines+=re.sub('Sp50',Ip50,trisnippet)
    assigntriangles.lines+=re.sub('STEPNUM',I,tmp)
assigntriangles.writefile()

# assignpositions processing ==
assignpositions.loadfile()
possnippet=assignpositions.lines
assignpositions.lines=re.sub('STEPNUM','1',possnippet)
for i in range(1 ,int(stepnum)+1 ):
    I=str(i+1)
    assignpositions.lines+=re.sub('STEPNUM',I,possnippet)
assignpositions.writefile()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# runmovie processing 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#    wait times
time1=str(40); time2=str(10); time3=str(30); time4=str(10)
n_range=5;

runmovie.loadfile_byline()
header,s1_sn,stepping,footer=parse_runmovie() #pulls out headers, etc.
runmovie.lines=header[:]
runmovie.lines.append('#\n')
        
for i in range(1,int(stepnum)+1,n_range):
    for snippetline in s1_sn:
        if 'S1 to SN' in snippetline: 
            outline = re.sub('S1 to SN',str(i)+' to '+str(i+n_range-1) ,snippetline)
            runmovie.lines.append(outline)

        if 'S1-SN' in snippetline:
            for k in range(n_range):
                outline = re.sub('S1-SN',str(i+k),snippetline)
                runmovie.lines.append(outline)
    
    runmovie.lines.append('\n')

    for k in range(n_range):
        for snippetline in stepping:
            if 'S1 to SN' in snippetline: 
                outline = re.sub('S1 to SN',str(i+k)+' to '+str(i+n_range-1) ,snippetline)
                runmovie.lines.append(outline)
            else:
                outline = re.sub('S_P0',str(i+k),snippetline)
                outline = re.sub('S_P1',str(i+k+1), outline)
                outline = re.sub('TIME_1',time1, outline)
                outline = re.sub('TIME_2',time2, outline)
                outline = re.sub('TIME_3',time3, outline)
                outline = re.sub('TIME_4',time4, outline)
                runmovie.lines.append(outline)

        runmovie.lines.append('\n')
    runmovie.lines.append('\n')


runmovie.lines.append('\n')
runmovie.lines.extend(footer) # adding footer 

tmp=runmovie.lines[:]; runmovie.lines=''
for outline in tmp:
    runmovie.lines+=outline
runmovie.lines=re.sub('DIRNAME',dirname,runmovie.lines)
filename='clust-'+str(clustnum)+'-trial'+str(trialnum)+'-n-'+str(stepnum)
runmovie.lines=re.sub('FILENAME',filename,runmovie.lines)
    
runmovie.writefile()

# JN DEBUG===============================
#from rlcompleter import readline; readline.parse_and_bind("tab:complete")
#import code; code.interact(local=locals())
# JN DEBUG===============================
