from ray_triangle import *
import time

#plucker test    
x=np.array([2.0,3.0,7.0])
y=np.array([2.,1.0,0.0])    
print plucker(x,y)
# basic tetrahedroon for testing:
root2 = 2**(0.5)
X  = np.zeros((5,3))  # coordinates at origin
X[0,:]=[    0,    0,          0 ]
X[1,:]=[  0.5,    0, -0.5/root2 ]
X[2,:]=[ -0.5,    0, -0.5/root2 ]
X[3,:]=[  0.0,  0.5,  0.5/root2 ]
X[4,:]=[  0.0, -0.5,  0.5/root2 ]

# exact overlap
X_new=X
print "zero displacement"
for i_node in range(3):
    for j_face in range(3):
        print "checking node: "+ str(i_node) + " and face: "+ str(j_face) \
               +" "+str(check_overlap(X_new,i_node,X,[j_face]) )
        

# simple shift by 1/4 (guaranteed overlap)
d=0.25
X_new = np.add(X,d)
print "displacement ="+ str(d)
for i_node in range(3):
    for j_face in range(3):
        print "checking node: "+ str(i_node) + " and face: "+ str(j_face) \
               +" "+str(check_overlap(X_new,i_node,X,[j_face]) )
        



# simple shift by 5 (no overlap)
d=5.0
X_new = np.add(X,d)
print "displacement ="+ str(d)
for i_node in range(3):
    for j_face in range(3):
        print "checking node: "+ str(i_node) + " and face: "+ str(j_face) \
               +" "+str(check_overlap(X_new,i_node,X,[j_face]) )

