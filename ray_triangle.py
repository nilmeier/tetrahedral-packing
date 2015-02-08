# contains routines for computing ray triangle intersections following
#  http://en.wikipedia.org/wiki/Pl%C3%BCcker_coordinates
#  and 
#  http://www.researchgate.net/publication/241758122_Fast_Ray-Tetrahedron_Intersection_Using_Plucker_Coordinates




import numpy as np

db=False
#global tetrahedral indices
# face indexing
face_index=[]
#following indexing conventions of Platis & Theoharis
face_index.append([3,2,1])
face_index.append([2,3,0])
face_index.append([1,0,3])
face_index.append([0,1,2])
#we have to exclude the centroid coordinate at the 0 index
face_index=np.add(face_index,1)

# edge indexing
edge_index=[]
edge_index.append([1,2])
edge_index.append([2,0])
edge_index.append([0,1])


def plucker(x,y):
    """
    input:  2 numpy arrays with coordinates
    output: Plucker coordinates as an array of np arrays
    """
    return [(y-x),np.cross(x,y)]

def get_face(X,i):
    """
    uses face_index as a global
    """
    return X[face_index[i],:]

def get_edge(F,j):
    """
    following indexing convention in Platis & Theoharis
    """
    return plucker(F[ edge_index[j][0] ],F[edge_index[j][1] ] )


def get_ray(X,j_node,k_ray):
    """
    getting the ray that points from each of the 
    vertices of the face to the reference point.  I think we need only one, 
    but I'll return all 3 for now.
    """
    F=get_face(X,j_node)
    #remember that the first coordinate of X is the centroid
    x0=X[j_node+1,:]
    return plucker(x0,F[k_ray])

    
def permuted_inner_product(pr,ps):
    """
    permuted inner product
    """
    # pi_1 (.) :e pi_2 = u_r . v_s + u_s . v_r 
    return np.dot(pr[0],ps[1])+np.dot(ps[0],pr[1])

def evaluate_triangle_tests(t):
    case1 = False # case 1 is a ray entering a triangle
    case2 = False # case 2 is a ray exiting a triangle
    case3 = False # case 3 is a coplanar (intersecting) ray
    eps=1.0e-14 # numerical error 
    if t[0]>=eps and t[1]>=eps and t[2]>=eps:
            case1=True
            if abs(t[0])<=eps or abs(t[1])<=eps or abs(t[2])<=eps:
                case1=False           
    if t[0]<=0.0 and t[1]<=0.0 and t[2]<=0.0:
            case2=True
            if abs(t[0])<=eps or abs(t[1])<=eps or abs(t[2])<=eps:
                case2=False
    if abs(t[0])<=0 and abs(t[1])<=0 and abs(t[2])<=0:
        case3=True

    if db: print " c1: " + str(case1)+ " c2: " + str(case2) +  " c3: "+ str(case3)    

    if case1 or case2 or case3: 
        return True
    else:
        return False             


def check_overlap(X_new, new_node_index, X_ref,face_list):
    #use rays in X_new and triangles in X_ref
    # face_list is an array of faces to be checked (all faces=[0,1,2,3]
    # new_node_index is an integer 
    for i_face in face_list:        
        F_ref = get_face(X_ref,i_face)
        for j_point in [new_node_index]:   
            for k_ray in range(3): # may be restricted
                ray_new = get_ray(X_new,j_point,k_ray)
                if db: print "ray: "+ str(k_ray)+' = '+str(ray_new)

                # checking ray against all 3 sets (should fall out if overlap)
                if db: print 'edges in face '+str(i_face)
                testlist=[]
                for l_edge in range(3):
                    edge_ref = get_edge(F_ref,l_edge)
                    testlist.append(permuted_inner_product(ray_new,edge_ref))
                    if db: print "e"+str(l_edge) + "t="+str(testlist[:])
                
                intersect=evaluate_triangle_tests(testlist)
                if db: print "ray: "+str(k_ray) +" from node: "+ str(j_point) \
                      +" and triangle "+str(i_face)+ "intersect? "+str(intersect)  

                if intersect: return True      
    return False                      
                       
