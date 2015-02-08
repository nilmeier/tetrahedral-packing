"""Cluster Module:  Base classes for cluster of tetrahedrons"""
import numpy as np
import ray_triangle
import polyhedra_volume 

db=False
edge  = 2.0
root2 = 2**(0.5)
root3 = 3**(0.5)
atom_labels = [' C1', '1H1','2H1','3H1','4H1']
pdbstring="HETATM %4d %s  METH %4d     %7.3f %7.3f %7.3f  1.00  0.00" 
#  to be used with <F3>

# following some conventions listed in 
# http://en.wikipedia.org/wiki/Tetrahedron

# tetrahedral geometry globals 
#  interior angle (H1-C-H2):
theta = np.arctan(root2) #not used
#  complementary angle
phi = np.pi-theta        #not used

tetrahedron_volume = edge**3/6./root2
# distances along tetrahedral axis (h-c)
base_to_h3_distance = edge*root2/root3
base_to_c1_distance = edge/root3*(1./( np.sin(phi)**2 ) - 1.0 )**(0.5)  # actually c1-c1 distance
dsq_centroid_min = edge*4.0/3.0  # tetrahedra never overlap at this distance
eps = 0.2
dsq_centroid_tooclose = (base_to_c1_distance-eps)**2 # tetrahedra always overlap at this distance

#c1-h3 distance
bond_length = ( base_to_c1_distance/2.**2 + edge**2/3.0 )**0.5



class Cluster(object):
    def __init__(self,size):
        """Initializes Cluster"""
        self.size=size #number of tetrahedra in cluster
        self.openfaces=[]      
        ## List of labels. They can also strings if needed.
        self.tetrahedra_labels = [1] #labels added only after tetrahedron is appended
        self.clustersize=len(self.tetrahedra_labels)
        ## Array of states for each tetrahedron.
        ## 1st element in tuple is tetrahedraon index, (i_tet) 
#       ## and second is face index (starts at 1)
        self.openfaces.extend([(0,1),(0,2),(0,3),(0,4)])
        self.X  = np.zeros((5*self.size,3))  # coordinates at origin
 
#       # building a cubic lattice that will hold specified number of tetrahedra at 
#       # 85% packing fraction
        self.cube_edge=(self.size*tetrahedron_volume/0.85)**(1.0/3.0)
        self.cube_edge=(self.size*tetrahedron_volume/0.65)**(1.0/3.0)
        
        # initializing rectangular coordinates as a box
        self.box=[0.0 for i in range(4)];  #box dimensions array  
        self.box[0]=self.cube_edge;  self.box[1] =self.cube_edge; self.box[2]=self.cube_edge 
        self.box[3]=self.box[0]*self.box[1]*self.box[2]
        #building cubic symmetries (27 reflections needed!)     
        self.DXB = []; self.printlabel=[]
        for i in range(27): self.DXB.append( np.zeros((1,3)) )
        count=0; sft=[0,1.0,-1.0]; sftstr=['0','p','m']
        self.n_tet=0 
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.DXB[count][0,:]=\
                      [ sft[i]*self.box[0],sft[j]*self.box[1],sft[k]*self.box[2] ]
                    self.printlabel.append( 'x' + sftstr[i] + 'y' + sftstr[j] + 'z'+sftstr[k] )
                    count+=1
  
        self.place_first_tetrahedron()
        self.filename=''
    
    def place_first_tetrahedron(self):
        ## placing coordinates at origin
        self.X  = np.zeros((5*self.size,3))  # coordinates at origin
        #assigning fixed coordinates to first state   
        self.X[0,:]=[    0,    0,          0 ]
        self.X[1,:]=[  0.5,    0, -0.5/root2 ]
        self.X[2,:]=[ -0.5,    0, -0.5/root2 ]
        self.X[3,:]=[  0.0,  0.5,  0.5/root2 ]
        self.X[4,:]=[  0.0, -0.5,  0.5/root2 ]
        self.X=np.multiply(self.X,edge) #scale to edge value  
        ##global face_normal array that is updated frequently
        self.face_normal=np.zeros((1,3))  
        self.face_origin=np.zeros((1,3))
        A=self.random_rotation_matrix()
        self.X[range(5),:]=np.dot( A, self.X[range(5),:].T ).T
        # centering for xtal symmetry
        self.X[range(5),:]+=self.cube_edge/2.0 
        
    def random_rotation_matrix(self):
        th_x=np.random.random()*2.0*np.pi
        th_y=np.random.random()*2.0*np.pi
        th_z=np.random.random()*2.0*np.pi
        cx=np.cos(th_x); sx=np.sin(th_x)
        cy=np.cos(th_y); sy=np.sin(th_y)
        cz=np.cos(th_z); sz=np.sin(th_z)

        A = np.zeros((3,3))
        # setting A to Rx
        A[0,0] =1.0
        A[1,1] = cx;  A[1,2] = -sx
        A[2,1] = sx;  A[2,2] = cx
        # defining Ry
        R = np.zeros((3,3))
        R[0,0] = cy;  R[0,2] = sy
        R[1,1] = 1.0
        R[2,0] = -sy; R[2,2]=cy 
        A = np.dot(A,R)
        # defining Rz
        R = np.zeros((3,3))
        R[0,0] = cz;  R[0,1] = -sz
        R[1,0] = sz;  R[1,1] = cz
        R[2,2] = 1.0 
        A = np.dot(A,R)
        
        return(A)

    def print_xtal_copies(self,filename):
        for xflag in range(1,27):  #excludes 0th model 
            self.print_pdb(filename,xflag)
    
    def print_pdb(self,filename,xflag):
        self.filename=filename
        if xflag==0:
            f = open(self.filename+'.pdb','w')
            dX=self.DXB[0]
        else:
            f = open(self.filename+"."+self.printlabel[xflag]+'.pdb','w')

        dX=self.DXB[xflag]

        xtal_line='CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 2 3\n'\
           %(self.box[0], self.box[1],self.box[2],90.,90.,90.)
        f.write(xtal_line)   
        for i_tet in range( len(self.tetrahedra_labels) ):
            tetrahedron_index=i_tet+1 #tetrahedron count starts at 1
        #    print "printing tet #" + str(i_tet) + " to "+ self.filename
            X_index_start=5*i_tet
            for j in range(5):
                j_X=j+X_index_start  # locating j_X in np coordinate array
                atom_index=j_X+1     # indexing atoms starting with 1
                atom_label=atom_labels[j]  #atom labels from list
                x=self.X[j_X,0]; y=self.X[j_X,1]; z=self.X[j_X,2] #xyz coordinates
                x+=dX[0,0]; y+=dX[0,1]; z+=dX[0,2] 
                outline = pdbstring%( atom_index, atom_label, tetrahedron_index,x,y,z )
                #prin/ct outline
                f.write(outline+'\n')
        
        f.write('END')
        f.close()

    def print_packing_fraction(self,filename):
        
        self.filename=filename
        f = open(self.filename+'.vol.txt','w')
        for i_tet in range(len(self.tetrahedra_labels)):
            tetrahedron_index=i_tet+1 #tetrahedron count starts at 1
        #    print "printing tet #" + str(i_tet) + " to "+ self.filename
            X_index=range(5*tetrahedron_index)
            X_vol=self.X[X_index,:]
            hull_volume=polyhedra_volume.convex_hull_volume_bis(X_vol)
            packing_fraction=tetrahedron_index*tetrahedron_volume/hull_volume
            #cubic_packing_fraction=tetrahedron_index*tetrahedron_volume/self.cube_edge**3.0
            prism_packing_fraction=tetrahedron_index*tetrahedron_volume/self.box[3]
            f.write( "%5d  %10.3f %10.3f\n" %(tetrahedron_index, packing_fraction,prism_packing_fraction) )
        
        f.close()

    def check_box_for_overlap(self):
         # checks to make sure that: 
         # 1)centroids are in cell
         # 2) vertices don't overlap adjacent cell
 
         for i_tet in range(self.n_tet):
            # checking centroids manually:
            for j_xyz in range(3):
                dbe=self.X[(i_tet)*5,j_xyz]  - self.box[j_xyz]
                #print "%9.3f %9.3f"%(db,self.box[j_xyz])
                if ( dbe > 0.0 ): 
                    self.box[j_xyz]+=dbe
                    print "dbe:"+str(dbe)  + "box"+str(self.box[j_xyz])
                dbe=self.X[(i_tet)*5,j_xyz] 
                if ( dbe < 0  ): 
                    self.box[j_xyz]-=dbe; 
                    print "b:"+str(self.X[5*i_tet,:])
                    print "dbe:"+str(dbe)
                    self.X[range(5*self.n_tet),j_xyz]-=dbe                  
                    print "a:"+str(self.X[5*i_tet,:])
#
            self.box[3]=self.box[0]*self.box[1]*self.box[2]
            self.update_symmetry_vectors()
            
    def gen_bigbox(self):  #creates a prism according to present cluster shape:
        # this routine places boxes at centroid boundaries and 
        # adjusts as needed to accomodate vertex overlaps
        maxlist=[0 for i in range(3)]; minlist=[0 for i in range(3)]
        #centroids=range(0,5*(self.n_tet),5)
        all_points=range(5*self.n_tet)
        for i_tet in range(self.n_tet):
            for j_xyz in range(3):
                minlist[j_xyz]=self.X[all_points,j_xyz].min()
                maxlist[j_xyz]=self.X[all_points,j_xyz].max()
                self.box[j_xyz]=maxlist[j_xyz]-minlist[j_xyz]       
                self.box[j_xyz]+=2.0*edge  #give room for 2 tetrahedra in all directions
         
        self.box[3]=self.box[0]*self.box[1]*self.box[2]
        self.X[all_points,:] -= minlist 
        self.update_symmetry_vectors()   

    def scale_box_by_overlap(self):
        # this routine places boxes at centroid boundaries and 
        # adjusts as needed to accomodate vertex overlaps
        maxlist=[0 for i in range(3)]; minlist=[0 for i in range(3)]
        #centroids=range(0,5*(self.n_tet),5)
        all_points=range(5*self.n_tet)
        for i_tet in range(self.n_tet):
            for j_xyz in range(3):
                minlist[j_xyz]=self.X[all_points,j_xyz].min()
                maxlist[j_xyz]=self.X[all_points,j_xyz].max()
                self.box[j_xyz]=maxlist[j_xyz]-minlist[j_xyz]       
        
        self.box[3]=self.box[0]*self.box[1]*self.box[2]
        self.X[all_points,:] -= minlist 
        self.update_symmetry_vectors()             

        print "box: %9.3f %9.3f %9.3f %9.3f"%(self.box[0],self.box[1],self.box[2],self.box[3] )
        # starting from box that is tightly packed around centroids...
        # gradually increase box size to no more than 2x the bond length while screening for 
        # overlaps
        nsteps=50 #nominally hard coded nsteps
        #boxrange=np.array(range(0,nsteps+1))*(2.*bond_length)/nsteps 
        dbe=2.*bond_length/nsteps 


         #'shrinkwrap' loop: this allows us to optimize each side individually
        for j_xyz in range(3):
            doublebreak=False #have to jump out of two loops when an overlap is found
            for k in range(1,nsteps+1):
                if not(doublebreak):
                    # update box before looping through cluster
                    self.box[j_xyz]-=dbe
                    self.X[all_points,j_xyz]-=dbe/2.0
                    self.update_symmetry_vectors()


                    for i_tet in range(self.n_tet):
                        if db: print str(i_tet)+")=======dim "+str(j_xyz)
                        if not(doublebreak):
                            overlap = self.multiple_screens_new_tetrahedron(i_tet-1,self_exclude=True) 
                            if db: print "overlap "+str(overlap)
                            if overlap: 
                                if db: print "breaking and restoring"
                                self.box[j_xyz]+=dbe
                                self.X[all_points,j_xyz]+=dbe/2.0
                                self.update_symmetry_vectors()
                                overlap = self.multiple_screens_new_tetrahedron(i_tet-1,self_exclude=True) 
                                if db: print "overlap after restore?"+str(overlap)
                                doublebreak=True
                                print "ending compression of dim# "+str(j_xyz)+ " at k= "+str(k)+" tetra# "+str(i_tet)
                                #Xprint=self.X[(5*i_tet),:]
                                #print "  X restore: %9.3f %9.3f %9.3f"%(Xprint[0],Xprint[1],Xprint[2] )
                                #print "DXB restore: %9.3f %9.3f %9.3f"%(self.DXB[-1][0,0],self.DXB[-1][0,1],self.DXB[-1][0,2] )
                                #print "box restore: %9.3f %9.3f %9.3f %9.3f"%(self.box[0],self.box[1],self.box[2],self.box[3] )
                                if overlap:  
                                    print "restored coordinates stil overlap?" 
                                    self.box[j_xyz]+=dbe
                                    self.X[all_points,j_xyz]+=dbe/2.0
                                    self.update_symmetry_vectors()
                                    overlap = self.multiple_screens_new_tetrahedron(i_tet-1,self_exclude=True) 
                                    print "second restore:  overlap?"+ str(overlap)

   
    # repeated function in rigid_body!
    def update_symmetry_vectors(self):
        count=0; sft=[0,1.0,-1.0] 
        self.box[3]=1.0
        for i in range(3):
            self.box[3]*=self.box[i]
            for j in range(3):
                for k in range(3):
                    self.DXB[count][0,:]=\
                      [ sft[i]*self.box[0],sft[j]*self.box[1],sft[k]*self.box[2] ]
                    count+=1
   
    def set_all_openfaces(self):
       self.openfaces=[] #sets every face to open...should get filtered rapidly
       for i_tet in range(self.n_tet-1):
           #print str(i_tet) +" of " + str(self.n_tet)   
           self.openfaces.extend( [(i_tet,1),(i_tet,2),(i_tet,3)] )


    def multiple_screens_new_tetrahedron(self,i_tet,self_exclude):
        # in the future, a full sweep should not be required, but 
        # this will work for now.  Only nearest ne
        # also, a linked list would be faster here
        # finally, a check against openfaces would help
        closelist=[]
        for k_cell in range(27): closelist.append([])
        overlap=False

        # checking if new tetrahedron centroid is within cubic lattice
        if self_exclude:
            if ( self.X[(i_tet+1)*5,0]  > self.box[0] ) or \
               ( self.X[(i_tet+1)*5,1]  > self.box[1] ) or \
               ( self.X[(i_tet+1)*5,2]  > self.box[2] ) or \
               ( self.X[(i_tet+1)*5,0]  < 0.0  ) or \
               ( self.X[(i_tet+1)*5,1]  < 0.0  ) or \
               ( self.X[(i_tet+1)*5,2]  < 0.0  ):
                   print "centroid-box overlap" 
                   Xprint=self.X[5*(i_tet+1),:]
                   print "  X inline: %9.3f %9.3f %9.3f"%(Xprint[0],Xprint[1],Xprint[2] )
                   return True
        else:
            # used to add a slight offset here...
            if ( self.X[(i_tet+1)*5,0]  > self.box[0] ) or \
               ( self.X[(i_tet+1)*5,1]  > self.box[1] ) or \
               ( self.X[(i_tet+1)*5,2]  > self.box[2] ) or \
               ( self.X[(i_tet+1)*5,0]  < 0.0  ) or \
               ( self.X[(i_tet+1)*5,1]  < 0.0  ) or \
               ( self.X[(i_tet+1)*5,2]  < 0.0  ):
                   return True 
        
        # for crystal packing, we exclude the 'self' search when screening for 
        # overlap
        if self_exclude: cells_to_search=range(1,27)
        else: cells_to_search=range(27)

        # checking if centroids are too close by fast sphere overlap checking

        for j_check in range(i_tet+1):
            for k_cell in cells_to_search:  #looping through crystal copies 
                X_check=self.X[j_check*5,:]+self.DXB[k_cell]
                dsq_check=np.sum( ( self.X[(i_tet+1)*5,:]-X_check )**2 )

                if dsq_check<dsq_centroid_tooclose: #early screen 
                    if (self_exclude): print "centroid too close" 
                    return True
      
                if dsq_check<dsq_centroid_min:
                #if True:
                    #print "closelist"
                    closelist[k_cell].append( j_check )
            
                    #also checking if vertices get too close to centroid:
                    for l_vertex in range(1,5):
                        L=j_check*5+l_vertex
                        X_vertex_check=self.X[L,:]+self.DXB[k_cell]
                        dsq_vertex_check=np.sum((self.X[(i_tet+1)*5,:]-X_vertex_check)**2)
                        if dsq_vertex_check<(0.9*bond_length)**2:
                            if (self_exclude): 
                                print "vertex in centroid sphere: "+ str(dsq_check)
                            return True
                
        # this could be done in situ at some point
        # X_new contains the protruding (4th) vertex 
        X_new=self.X[range((i_tet+1)*5,(i_tet+1)*5+5),:]
        if not(self_exclude):  #standard branch
            for check in self.openfaces:
                for k_cell in cells_to_search:
                    if check[0] in closelist[k_cell]:
                        # face indexing in ray_triangle.py: 0:3
                        # face indexing in cluster.py: 1:4
                        # vertex 4 is the protruding vertex, and is 
                        # the only one that needs to be checked.
                        # passed in as an array [3] 

                        # X_ref contains the face to be checked
                        # doing a fast distance screen:
                        X_ref=self.X[range(check[0]*5,check[0]*5+5),:]
                        X_ref+=self.DXB[k_cell]

                        if db: print "checking tet: "+ str(check[0]) + " face: " + str(check[1]) + " "+str(overlap)
                        face_ref=check[1]-1  
                        if check[0] in closelist[k_cell]:
                            overlap=ray_triangle.check_overlap( X_new,3, X_ref, [face_ref] ) 
                            if overlap: return True


           
        if self_exclude:  # cell optimization branch
            for k_cell in cells_to_search:
                for cc in closelist[k_cell]:   #running through list of known close centroids (cc)
                    openfaces_cc =[(cc,1),(cc,2),(cc,3),(cc,4)] #making all faces of cc searchable 
                    for check in openfaces_cc:
                        X_ref=self.X[range(check[0]*5,check[0]*5+5),:]
                        X_ref+=self.DXB[k_cell]
                        face_ref=check[1]-1  
                        for vertex in range(3):  #checking all vertices for xtal packing
                                overlap=ray_triangle.check_overlap( X_new,vertex, X_ref, [face_ref] ) 
                                if overlap: 
                                    print "failed ray trace"
                                    return True    

              
        return overlap                                     

    def place_next_tetrahedron(self,i_tet):
        """adds a tetrahedron to the ith tetrahedron at specified face_index"""

        success=False
        if len(self.openfaces)==0:
            #"no more faces!" 
            return False 

        while success==False:
            # a dense notation for picking the face_index randomly from the openfaces list:
            (j_dock,face_index) = self.openfaces[np.random.randint(0,len(self.openfaces)) ] 
            # j_dock is the location of the site for growing the (i_tet+1)th tetrahedron

            if db: print str(i_tet)+") trying" + str((j_dock,face_index))   
            #do the easy part..coordinates H1,H2,H3, get same vertices as face_index                    
            face_coordinates_list=range(1,1+4)
            face_coordinates_list.remove(face_index)  
            self.X[(i_tet+1)*5+1,:] = self.X[ j_dock*5 + face_coordinates_list[0] , : ]
            self.X[(i_tet+1)*5+2,:] = self.X[ j_dock*5 + face_coordinates_list[1] , : ]
            self.X[(i_tet+1)*5+3,:] = self.X[ j_dock*5 + face_coordinates_list[2] , : ]
            #jndb:  ^ there is redundancy here...see get_face_normal()

            #getting center and outward normal of triangluar plane
            self.get_face_normal(j_dock,face_index)
            # C1 and H3 are on the same axis as the norma:e classl to the triangular plane
            # adding C1 coordinate
            norm=np.linalg.norm(self.face_normal )

            self.X[(i_tet+1)*5+0,:] = self.face_origin + self.face_normal/norm*base_to_c1_distance

            # addding H4 coordinate with C1 as origin   
            self.X[(i_tet+1)*5+4,:] = self.X[(i_tet+1)*5+0,:]+ self.face_normal   

            # screening for next round of addition------------------
            overlap=self.multiple_screens_new_tetrahedron(i_tet,self_exclude=False)
            if overlap:
                self.openfaces.remove((j_dock,face_index))

                if db: 
                    print "overlap!  now what?  rinse and repeat" 
                    print "removing both faces"
                    print self.openfaces
            else:
                self.tetrahedra_labels.append(self.tetrahedra_labels[-1]+1)
                # remove  (i_tet, face_index) tuple from openfaces array
                self.openfaces.remove((j_dock,face_index))
                # add c/reate a new openfaces array for i_tet+1.  The first 3 faces are opposite 
                # h1,h2,h3.  The face opposite h4 is always excluded, because it is constructed on the 
                # interface to the i_tet'th shape
                self.openfaces.extend( [(i_tet+1,1),(i_tet+1,2),(i_tet+1,3)] )
                self.n_tet=self.tetrahedra_labels[-1] #internal tet count
                success=True

            if len(self.openfaces)==0:
                #print "no more faces!"
                self.X[range( (i_tet+1)*5+0,(i_tet+1)*5+5 ),: ]=0.0 
                return False 
        return success        

    def get_face_normal(self,j_tet,face_index):
        """ returns a face center and outward normal of tetrahedron face """
        # face_index goes from 1 to 4...this matches coordinates indexing (C1 is at 0 index)
        
        # computing vectors from the first vertex in the face to the remaining 
        # two vertices...the sum of these vectors points to the center of the triangular face       
        x0=self.X[ j_tet*5+face_index, : ] # distal hydrogen
        x1=self.X[ j_tet*5,: ]             # central carbon
        
        self.face_normal = x1-x0 
        self.face_origin = x1 + self.face_normal \
                           /base_to_h3_distance*base_to_c1_distance
        
        self.face_origin = x1  
