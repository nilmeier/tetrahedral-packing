import cluster 
import polyhedra_volume
import copy
from optparse import OptionParser
import os
import rigid_body

db=False
parser = OptionParser()

parser.add_option("-t","--n_tet", 
               dest="n_t",default=10,
               help="number of trials per output cluster")
parser.add_option("-c","--n_clusters", 
               dest="n_c",default=7,
               help="number of output clusters in set")
parser.add_option("-b","--n_buffer", 
               dest="n_b",default=6,
               help="number of trial (buffer) clusters")
parser.add_option("-l","--l_trial", 
               dest="l_t",default=5,
               help="number of tetrahedrons per trial")
parser.add_option("-s","--random_seed",default=1,
               help="random seed")
parser.add_option("-o","--n_opt_steps",
                  dest="n_o",default=100,
                  help="# optimization steps")
(options,args) = parser.parse_args()


def main():
    a=[];b=[]; volume=[];fraction=[]
    n_t=int(options.n_t)
    n_c=int(options.n_c)
    n_b=int(options.n_b)  
    l_t=int(options.l_t)
    n_o=int(options.n_o)
    extra_steps=0
    cluster.np.random.seed(options.random_seed)
    output_directory = "out_t-"  + str(n_t) +\
                       "_c-" + str(n_c) + \
                       "_b-" +str(n_b) + \
                       "_l-" + str(l_t) + \
                       "_s-" + str(options.random_seed)+\
                       "_o-"  + str(options.n_o)
    print "printing to " + output_directory
    print "total # of tetrahedra per output file: " + str(n_t*l_t)
    if os.path.exists(output_directory):
        print "overwriting data in " + output_directory
    else:
        os.mkdir(output_directory)
            

# may need to think about memory mgmt in next revision
    for i_clust in range(n_c):
        a.append( cluster.Cluster(n_t*l_t+1))
    for i_buff in range(n_c*n_b):
        b.append(cluster.Cluster(n_t*l_t+1))        
        volume.append(0)
        fraction.append(0)

    for j_trial in range(n_t):  # number of trial runs of length=l_t
        # need to add some extra bookkeepping for when 
        # the cell fills up...this will apply to the entire buffer

        fraction=[0.0 for i in range(len(fraction)) ]
        volume=[0.0 for i in range(len(volume)) ]
        buffercount=[]
        for i_clust in range(n_c): # running through each cluster
            # buffer contains n_b copies of a[i_clust]A
            buffercount.append(0) #every cluster has at least one buffer
            if len(a[i_clust].openfaces)==0:  #no faces to populate...skip bufferlist altogether
                #print "cell is full...no buffers generated!"
                buffercount[i_clust]=1
                for k_buffer in range(n_b):  
                    l_buffer_full = k_buffer+i_clust*n_b  # each cluster gets n_b buffers
                    if k_buffer==0:
                        b[l_buffer_full] = copy.deepcopy( a[i_clust] )
                    fraction[l_buffer_full]=0.0
            else: # create buffers if source cluster has openfaces
                success=True
                for k_buffer in range(n_b):  # looping through temporary (buffer) sets
                    if success:
                        buffercount[i_clust]+=1  # keeping track of buffers
                        l_buffer_full = k_buffer+i_clust*n_b  # each cluster gets n_b buffers
                        # making a copy of i_clust to build off of
                        b[l_buffer_full] = copy.deepcopy( a[i_clust] )

                        # building out l_t tetrahedrons
                        if db: 
                             print "working on cluster" + str(i_clust)
                             print "working on buffer"+str(l_buffer_full)

                        for m_trial in range(l_t):
                            # each completed trial contains l_t additional tetrahedra    
                            n_tet = m_trial + j_trial*l_t 
                            
                            if db: print "i_clust: "+str(i_clust)+ "j_trial: "+str(j_trial)+"l_buff: "+ str(l_buffer_full)+" n_tet = "+ str(n_tet)+ " m= "+str(m_trial)
                            search_flag=True
                            while search_flag:
                              success=b[l_buffer_full].place_next_tetrahedron(n_tet)
                              if success or len(b[l_buffer_full].openfaces)==0: search_flag=False


        # populating arrays with successful trials:
        for i_clust in range(n_c): # running through each cluster
                # print "f:" + str(fraction)
                for k_buffer in range(buffercount[i_clust]): #looping through successful buffers  
                    l_buffer_full = k_buffer+i_clust*n_b  # each cluster gets n_b buffersi

                    if j_trial==n_t-1:  #sorting on prism packing for last step
                        volume[ l_buffer_full ]=b[i_clust].box[3]
                    else:               # all other steps sort on hull volume
                        volume[ l_buffer_full ]=polyhedra_volume.convex_hull_volume_bis(b[l_buffer_full].X)
                    
                    fraction[ l_buffer_full ]=cluster.tetrahedron_volume*(b[l_buffer_full].n_tet+1)/volume[l_buffer_full] 
                    
                      
                if db: 
                    print "vol [ " + str(l_buffer_full) + " ]= " + str( volume[l_buffer_full] )
                    print "fraction = " + str(cluster.tetrahedron_volume*(n_tet+1)/volume[l_buffer_full])

        volume_indices=[ i[0] for i in sorted(enumerate(fraction), key=lambda x:x[1], reverse=True) ]
        out_string = output_directory + ', step-'+ str(n_tet+1)  
        print "working on "+ out_string 

        for i_clust in range(n_c):
            # taking top n_c candidates (lowest volume first) 
            if db: print str(i_clust)+" )taking "+str(volume_indices[i_clust])+"th buffer ...f= %9.3f"%(fraction[volume_indices[i_clust]]) 
            # only update cluster if fraction is unique (artifacts have 0.0 fraction)
            if fraction[volume_indices[i_clust]]!=0:
                a[i_clust]=copy.deepcopy(b[volume_indices[i_clust]])                

            #printing intermediates
            filename = \
            output_directory+'/clust-'+str(i_clust)+'-trial-'+str( j_trial )+"-n-"+str(a[i_clust].n_tet)   

            # optimizing and printing for comparison:
            if (j_trial%5==0 and j_trial>0) or j_trial==n_t-1 :
                filename = output_directory+'/clust-'+str(i_clust)+'-trial-'+str( j_trial )+"-n-"+str(a[i_clust].n_tet)   
                print "\nprinting "+filename + " \n   with " + str(a[i_clust].n_tet) + " tetrahedra"
                a[i_clust].print_pdb(filename,xflag=0)
                a[i_clust].print_xtal_copies(filename) 
                a[i_clust].print_packing_fraction(filename)
                #optimizing and printing
                extra_steps  = j_trial/10 * n_o/2 
                extra_steps += 0 #100
                active_coords = range( 5*a[i_clust].n_tet) # for nonuniform cluster lengths
                Xin=copy.deepcopy(a[i_clust].X[active_coords,:])
                boxin=copy.deepcopy(a[i_clust].box)
                print "optimizing "+ str(j_trial) +"th step for "+str(n_o+extra_steps)+" steps" 
                box_opt=False #optimize box edges flag
                Xout,boxout=rigid_body.optimize_packing(Xin,boxin, n_o,'condense_cluster')
                a[i_clust].box=copy.deepcopy(boxout)
                a[i_clust].X[active_coords,:] = copy.deepcopy(Xout)

                if j_trial==n_t-1:  #final optimization branch
                    box_opt=True
                    fin_factor=3    #running for extra steps in final optimization
                    Xout,boxout=rigid_body.optimize_packing(Xin,boxin, fin_factor*n_o,'condense_cluster')
                    a[i_clust].box=copy.deepcopy(boxout)
                    a[i_clust].X[active_coords,:] = copy.deepcopy(Xout)
                    a[i_clust].scale_box_by_overlap()
                    print "new box: %9.3f %9.3f %9.3f %9.3f"%(a[i_clust].box[0],a[i_clust].box[1],a[i_clust].box[2],a[i_clust].box[3] )

                    Xin=copy.deepcopy(a[i_clust].X[active_coords,:])
                    boxin=copy.deepcopy(a[i_clust].box)
                    Xout,boxout=rigid_body.optimize_packing(Xin,boxin, fin_factor*n_o,'pbc_rigid_box') #full_opt or pbc_rigid_box
                    a[i_clust].box=copy.deepcopy(boxout)
                    a[i_clust].X[active_coords,:] = copy.deepcopy(Xout)
                    print "before overlap box: %9.3f %9.3f %9.3f %9.3f"%(a[i_clust].box[0],a[i_clust].box[1],a[i_clust].box[2],a[i_clust].box[3] )
                    a[i_clust].scale_box_by_overlap()
                    print "after overlap box: %9.3f %9.3f %9.3f %9.3f"%(a[i_clust].box[0],a[i_clust].box[1],a[i_clust].box[2],a[i_clust].box[3] )

                    Xin=copy.deepcopy(a[i_clust].X[active_coords,:])
                    boxin=copy.deepcopy(a[i_clust].box)
                    Xout,boxout=rigid_body.optimize_packing(Xin,boxin, fin_factor*n_o,'condense_cluster')
                    a[i_clust].box=copy.deepcopy(boxout)
                    a[i_clust].X[active_coords,:] = copy.deepcopy(Xout)
                    print "before overlap box: %9.3f %9.3f %9.3f %9.3f"%(a[i_clust].box[0],a[i_clust].box[1],a[i_clust].box[2],a[i_clust].box[3] )
                    a[i_clust].scale_box_by_overlap()
                    print "after overlap box: %9.3f %9.3f %9.3f %9.3f"%(a[i_clust].box[0],a[i_clust].box[1],a[i_clust].box[2],a[i_clust].box[3] )

                filename+="-opt"
                print "writing "+ filename 
                print "===>fraction in cell = " + str( cluster.tetrahedron_volume*(a[i_clust].n_tet+1) /a[i_clust].box[3])
                print "\n"
                a[i_clust].print_pdb(filename,0) 
                a[i_clust].print_xtal_copies(filename) 
                a[i_clust].print_packing_fraction(filename)
                
                #resetting ALL openfaces after optimization
                a[i_clust].set_all_openfaces()

#JN DEBUG LINE ================================================
    import rlcompleter; import readline; readline.parse_and_bind("tab: complete")
    print "jndb\n"; import code; code.interact(local=locals())
#JN DEBUG LINE =================================================

if __name__=="__main__":
    main()
