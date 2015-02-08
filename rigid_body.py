import numpy as np
from cluster import edge,bond_length,base_to_c1_distance  #geometry
import copy

cutoff_sq = (0.45*edge)**2 # maximum distance for contact potentials.
k_spring=1.0
epsilon=1.0e-10           # numerical roundoff for distance comparisons.
rotational_scale=500.0    # bigger number= smaller stepsize (for all scales)
translational_scale=500.0 
volume_scale=10000.0

#gaussian potential parameters
sigma=0.005*edge
A0=0.25
alpha=1.0/2.0/sigma**2


def find_close_centroids(X,DXB):
    num_centroids=len(X)/5
    centroid_closelist=[ [] for i in range(27) ]; 
    # finding centroids that are close enough to have contact potentials
    for i in range(num_centroids):
        for k_cell in range(27):  
            for j in range( i+1,num_centroids):
                Xch = X[j*5,:] + DXB[k_cell][0,:] # creating symmetry copy
                dsq_check = np.sum( (X[i*5,:]-Xch)**2 )    
                if dsq_check < (2.*bond_length)**2 and dsq_check > epsilon :
                    centroid_closelist[k_cell].append( (i,j) )

    return centroid_closelist

def find_close_sites(X,DXB,centroid_closelist):
    site_closelist=[[] for i in range(len(centroid_closelist))]
    #centroid_closelist can be an array of length one for box_opt=False
    for j_cell in range(len(centroid_closelist)): 
        for k in centroid_closelist[j_cell]:
            i=k[0]; j=k[1]
            for l in range(1,5):
                for m in range(1,5): 
                    I=5*i+l; J=5*j+m
                    Xch=X[J,:]+DXB[j_cell][0,:] # creating symmetry copy
                    dsq_check = np.sum( (X[I,:]-Xch)**2 )          
                    if dsq_check<=cutoff_sq and dsq_check>=epsilon**2:
                        site_closelist[j_cell].append( ( I,J ) )
    return site_closelist


def quaternion_to_rotation_matrix(Q):
    q0 = Q[0]; q1 = Q[1]; q2 = Q[2]; q3 = Q[3]
    A = np.zeros((3,3))

    A[0,0] = ( q0**2 + q1**2 - q2**2 - q3**2 )/2.0
    A[0,1] = q1*q2 - q0*q3 
    A[0,2] = q1*q3 + q0*q2

    A[1,0] = q1*q2 + q0*q3
    A[1,1] = ( q0**2 - q1**2 + q2**2 - q3**2 )/2.0
    A[1,2] = q2*q3 - q0*q1 
 
    A[2,0] = q1*q3 - q0*q2
    A[2,1] = q0*q1 + q2*q3 
    A[2,2] = ( q0**2 - q1**2 - q2**2 + q3**2 )/2.0

    A=np.multiply(A,2.0)
    return A


def compute_forces_and_torques(X,T,F,P,DXB,closelist,box_opt):
    # searching those tetrahedra for vertices with contact potentials
    # full array need not be zeroed... next rev! 
       
    F[:]=0.0  # Force is a global array
    T[:]=0.0
    P[:]=0.0
    nonzero_forces = []
    for k_cell in range(len(closelist)):
        for i,j in closelist[k_cell]:
            fij=0.0
            Xch=X[j,:]+DXB[k_cell][0,:]; # creating symmetry copy
            dr = X[i,:]-Xch
            dsq_check = np.sum( dr**2 )
            fij=k_spring*dr 
            fij += A0*alpha*dr*np.exp(-alpha*dsq_check)

            I=i/5; J=j/5  #note the int trick works as floor
            F[ I ] +=  fij 
            F[ J ] += -fij 
            # each point has a lever arm about the tetrahedron center:
            r_i = X[i,:] - X[I*5,:]
            r_j = X[i,:] - X[J*5,:]
            
            if box_opt :#and :#not(k_cell==0):
                Xcom=np.divide(DXB[-1][0,:],-2.0) # center of prism obtained from last value of DXB
                for i in range(3):
                  if Xch[i] < Xcom[i]:  # force from minus side
                     P[0,i]+=fij[i]            
                  if Xch[i] > Xcom[i]:
                     P[0,i]-=fij[i]
                         
            # computing torques on each point
            tij    = np.cross(r_i,fij)
            tji    = np.cross(r_j,-fij)
            # accumulating torques for each centroid
            T[I]  += tij
            T[J]  += tji

    # returning a measure of the objective:
    return (np.linalg.norm(T)+np.linalg.norm(F), T,F,P )               

     
def propagate_rotational_dof(X,T,scale):
    A = np.zeros((3,3))
    for i in range(len(X)/5):
        t = -T[i,:] 
        norm=np.linalg.norm(t)
        if norm>=epsilon:
            s_th_2=norm/scale
            if (abs(s_th_2)>=1): 
                print "scale is not set correctly: s(th/2)= "+ str(s_th_2)
                s_th_2=0.1

            c_th_2=(1-s_th_2**2)**(0.5)

            q1,q2,q3=np.multiply(t, (s_th_2 /norm) )  #unit vector for last 3 components
            q0=c_th_2
            #q0=1.0 
            
            A[0,0] = ( q0**2 + q1**2 - q2**2 - q3**2 )/2.0
            A[0,1] = q1*q2 - q0*q3 
            A[0,2] = q1*q3 + q0*q2

            A[1,0] = q1*q2 + q0*q3
            A[1,1] = ( q0**2 - q1**2 + q2**2 - q3**2 )/2.0
            A[1,2] = q2*q3 - q0*q1 
         
            A[2,0] = q1*q3 - q0*q2
            A[2,1] = q0*q1 + q2*q3 
            A[2,2] = ( q0**2 - q1**2 - q2**2 + q3**2 )/2.0

            A=np.multiply(A,2.0)
        else:  #in case the prefilter doesn't catch it
            A[0,0]=1.0; A[1,1]=1.0; A[2,2]=1.0
        
        # updating global position in place
        Xc = X[5*i,:]  #location of center 
        Xr = X[range(5*i,5*i+5),:]-Xc #shifting rotation to center

        X[range(5*i,5*i+5),:]=np.dot(A,Xr.T).T + Xc

    return X
     
def propagate_translational_dof(X,F,scale): 
    for i in range(len(X)/5):
        dX=np.divide(-F[i],scale)
        for j in range(5):
            X[5*i+j,:]+=dX
    return X

def propagate_box_dof(X,P,box,scale):
# changing box dimensions
    db=np.zeros((1,3))
    box[3]=1.0 #volume vector updated
    for i in range(3):
        db[0,i]=-P[0,i]/scale
        box[i]+=db[0,i] 
        box[3]*=box[i] 

    DXB=update_symmetry_vectors(box)    
#moving coordinates to new center of mass    
    X=X+np.divide(db,2.0)
    if db: print "DXB :"+ str(DXB[-1])+ " in function: "+ str(X[-1,:])

    return X,box,DXB


def update_symmetry_vectors(box):
    #building cubic symmetries (27 reflections needed!) (same as cluster.init)     
    DXB = []; 
    for i in range(27): DXB.append( np.zeros((1,3)) )
    count=0; sft=[0,1.0,-1.0] 
   
    for i in range(3):
        for j in range(3):
            for k in range(3):
                DXB[count][0,:]=\
                  [ sft[i]*box[0],sft[j]*box[1],sft[k]*box[2] ]
                count+=1
    return DXB  

def optimize_packing(X,box,n_opt_steps,opt_flag):
    size = len(X);
    nt   = size/5
    F    = np.zeros( (nt,3 ) )  # resultant force
    T    = np.zeros( (nt,3 ) )  # torque on centroid
    P    = np.zeros( (1,3 ) )  # force on prism walls 

    DXB = update_symmetry_vectors(box)

    centroid_closelist = find_close_centroids(X,DXB) #one time centroid screen
    closelist=find_close_sites(X,DXB,centroid_closelist) #site screen
    print "optimizing " + opt_flag + " for "+str(n_opt_steps)
    X_orig=copy.deepcopy(X)
    box_opt=False
    if opt_flag=='full_opt' or opt_flag=='box_only':  box_opt=True
    S_orig,T,F,P = compute_forces_and_torques(X,T,F,P,DXB,closelist,box_opt)
    S=[]; S.append(S_orig)
    dX=[]; dX.append(0.0);

    if opt_flag=='condense_cluster':  #standard optimization (no pbc, no moving volume)
        box_opt=False
        for k_opt in range(n_opt_steps):
            X=propagate_rotational_dof(X,T,rotational_scale)
            X=propagate_translational_dof(X,F,translational_scale)
            closelist = find_close_sites(X,DXB,centroid_closelist)
            s,T,F,P   = compute_forces_and_torques(X,T,F,P,DXB,[closelist[0]],box_opt)
            dx = np.linalg.norm(X-X_orig) 
            S.append(s)
            dX.append(dx)

    if opt_flag=='box_only':
        box_opt=True
        for k_opt in range(n_opt_steps):
            X,box,DXB   = propagate_box_dof(X,P,box,volume_scale) 
            closelist = find_close_sites(X,DXB,centroid_closelist)
            s,T,F,P   = compute_forces_and_torques(X,T,F,P,DXB,closelist,box_opt)
            dx = np.linalg.norm(X-X_orig) 
            S.append(s)
            dX.append(dx)
            print "P: %9.3f %9.3f %9.3f |P| %9.3f"%( P[0,0],P[0,1],P[0,2],np.linalg.norm(P) )
            print "b: %9.3f %9.3f %9.3f  v %9.3f"%(box[0],box[1],box[2],box[3])


    if opt_flag=='pbc_rigid_box':
        box_opt=False
        for k_opt in range(n_opt_steps):
            X=propagate_rotational_dof(X,T,rotational_scale)
            X=propagate_translational_dof(X,F,translational_scale)
            closelist = find_close_sites(X,DXB,centroid_closelist)
            s,T,F,P   = compute_forces_and_torques(X,T,F,P,DXB,closelist,box_opt)
            dx = np.linalg.norm(X-X_orig) 
            S.append(s)
            dX.append(dx)


    if opt_flag=='full_opt':
        box_opt=False
        for k_opt in range(n_opt_steps):
            X=propagate_rotational_dof(X,T,rotational_scale)
            X=propagate_translational_dof(X,F,translational_scale)
            X,box,DXB   = propagate_box_dof(X,P,box,volume_scale) 
            closelist = find_close_sites(X,DXB,centroid_closelist)
            s,T,F,P   = compute_forces_and_torques(X,T,F,P,DXB,closelist,box_opt)
            dx = np.linalg.norm(X-X_orig) 
            S.append(s)
            dX.append(dx)



    print "k(fin): "+ str(k_opt) + " S: %7.5f dS %7.5f d|X-X0] %7.5f"%(s,s-S[-1], dx-dX[-1])
    return X, box

if __name__=="__main__":
      print "main called...nothing happening!"
#     optimize_packing(X,100)
