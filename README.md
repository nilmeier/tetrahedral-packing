This program will generate random tetrahedral packings in a cubic lattice.  It attempts
to find high density packings.  It is only in its first generation of development,
however, so it has some ways to go before it is complete.

To run with defaults, type python main.py

to see the options, add the -h flag:

Options:
  -h, --help            show this help message and exit
  -t N_T, --n_tet=N_T   number of trials per output cluster
  -c NC, --n_clusters=N_C
                        number of output clusters in set
  -b N_B, --n_buffer=N_B
                        number of trial (buffer) clusters
  -l L_T, --l_trial=L_T
                        number of tetrahedrons per trial
  -s RANDOM_SEED, --random_seed=RANDOM_SEED
                        random seed
  -o N_O, --n_opt_steps=N_O
                        # optimization steps

The output will go into a local directory that contains the 
options in the directory name.  Within the directory, you will 
find a list of pdb files.  These files contain methane like
molecules, with the hydrogens given as vertices.  To convert to 
the citrine format, there is a pdb_to_cit.py utility in this 
directory.  There are also a variety of scripts for generating
UCSF Chimera figures, which I will cover (later), once they are 
in a more stable form.

The following naming conventions apply to the output file:

clust-X-trial-Y-n-Z.XTAL_LOC.pdb

X:  the cluster number (goes from 0 to n_cluster-1)
Y:  the trial number (goes from 0 to n_tet-1)

Each trial consists of a series of l_trial attempts to place 
tetrahedra a given cluster.  As the packing becomes more dense, 
they are not always successful.

Z:  the number of tetrahedra in the pdb file 
    (up to Y*l_trial+1 maximum).

XTAL_LOC:  The label for the crystal copy.  

Since the symmetry is a simple rectangular prism displacement, 
each extension indicates the nature of the displacement relative 
to the central crystallographic unit with extension x0y0z0.  

The dimensions of the crystal are listed in the headers.  
For example xpymz0 means x is displaced by +DX y is displaced by
-DY, and z is unchanged, where DX,DY,DZ are the prism dimensions.



  

