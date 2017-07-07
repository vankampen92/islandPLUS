#include "HEADERS.h"
#include "da_IBT_Functions.h"

// #define VERBOSE

void UPGMA_CLUSTERING_PARTITION( double ** d, int No_of_SPECIES , 
				 int *** PARTITION, int * G, int ** K )
{
  /* This function is used to build a upgma based tree from a 
     distance matrix 'd'

     In addition, the function yields also the corresponding 
     'PARTITION' of this tree, which is the different ways in 
     which a set of 'No_of_SPECIES' objects can be subdivided 
     depending on a distance threshold defining the depth at which 
     tree branches are cut. 
 	
     So, the code is organized in two sequential calls:
     
     A. upgma_clustering(distance, No_of_SPECIES,
                   d, &N, Cluster, &No_of_NODES, Node_List );
     which yields the actual tree (stored in cluster). 
   
     B. upgma_cluster_to_partition (Cluster, No_of_SPECIES, 
                              PARTITION, G, K );
     which yields the actual partition set from the tree
     stored in Cluster. 
     
     The rest of code lines are sort of superfluous! 
  */
  int i,j; 
  int N, No_of_NODES; 

  N = No_of_SPECIES; 
  node ** Cluster = (node **)calloc( 2*No_of_SPECIES-1, sizeof(node *) );
  for(i = 0; i < (2*No_of_SPECIES-1); i++ )  {
    Cluster[i] = (node *)calloc(1, sizeof(node));
    Cluster[i]->left = Cluster[i]->right = NULL;
    Cluster[i]->No_of_SPECIES = 0;
    Cluster[i]->Species = (int *)calloc( No_of_SPECIES, sizeof(int) );
    Cluster[i]->d    = 0.0;
    if( i < No_of_SPECIES ) { 
      Cluster[i]->Sp = 1+i;
      Cluster[i]->No_of_SPECIES = 1; 
      Cluster[i]->Species[0]    = i;
    }
    else {
      Cluster[i]->Sp = 0; 
      Cluster[i]->No_of_SPECIES = 0; 
    }
  }
  // The first N nodes from Cluster correspond to the tree leaves 
  // with species labels in order {1, 2, 3, ..., No_of_SPECIES}.
  // Tree leaves point to NULL daughter branches. 
  // The leaf node Cluster[i] is set up with this label:
  // Cluster[i]->Sp = 1+i;
  // If Cluster[i]->Sp = 0, then this is an inner tree node with
  // two non-NULL daughter branches. 
  
  // The initial Node_List contains the tree leaves.
  int * Node_List = (int *)calloc( No_of_SPECIES, sizeof(int) );
  for(i = 0; i < No_of_SPECIES; i++ ) Node_List[i] = i;
  
  if (N != No_of_SPECIES) error(0,0,"Program aborted");
  No_of_NODES = No_of_SPECIES; //only leaves..

  double ** distance  = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  for(i=1; i<No_of_SPECIES; i++) { 
    distance[i] = (double *)calloc( i, sizeof(double) );
    for(j=0; j<i; j++) distance[i][j] = d[i][j]; 
  }

#if defined VERBOSE    
  Print_Triangular_Matrix( d,        No_of_SPECIES );
  Print_Triangular_Matrix( distance, No_of_SPECIES );
#endif

  Rprintf(" About to enter the upgma clustering algorithm...\n");
  
  upgma_clustering(distance, No_of_SPECIES,
		   d, &N, Cluster, &No_of_NODES, Node_List );
  
  Rprintf(" Just out from the upgma clustering algorithm...\n");   
  //getchar();

#if defined VERBOSE
  Rprintf("In Order Display\n");
  print_inorder(Cluster[No_of_NODES-1]);
  
  Rprintf("Pre Order Display\n");
  print_preorder(Cluster[No_of_NODES-1]);

  Rprintf("Post Order Display\n");
  print_postorder(Cluster[No_of_NODES-1]);

  Rprintf("Print out only all nodes\n");
  for(i=0; i<No_of_NODES; i++) {
    Rprintf("No of Species in cluster %d: %d ( Species: ", i, Cluster[i]->No_of_SPECIES);
    for(j=0; j<Cluster[i]->No_of_SPECIES; j++) Rprintf("%d ", Cluster[i]->Species[j]);
    Rprintf(")\n");
  }
#endif

  upgma_cluster_to_partition (Cluster, No_of_SPECIES, 
			      PARTITION, G, K );
	   
  //To delete the whole tree, we need to pass the tree root 
  //to this recurrent funcion. 
  deltree(Cluster[No_of_NODES-1]);
  free(Node_List);

  for(i=1; i<No_of_SPECIES; i++) free(distance[i]);
  free(distance);  
}

void upgma_cluster_to_partition ( node ** Cluster, int N, 
				  int *** PARTITION, int * G, int ** K )
{
  /* This function converts the information from a cluster tree into 
     the total set of partitions in which the whole set of N objects
     can be subdivided depending on a distance threshold defining 
     the depth at which tree branches are cut. 
  */
  int i, j, jj, k, n;
  int No_of_NODES;  /* Total number of nodes */

  node *** C_List = (node ***)calloc(N, sizeof(node **) );
  for(i = 0; i<N; i++) 
    C_List[i] = (node **)calloc(i+1, sizeof(node *) );

  No_of_NODES = 2 * N - 1;
  
  /* 0th Partition: one single set containing all objects */ 
  K[0][0] = N; /* No of elements of the set: all N */
  for(k=0; k<N; k++) PARTITION[0][0][k] = k;
  C_List[0][0] = Cluster[No_of_NODES-1]; 

  /* 1st Partition: by dividing the root node */
  n = No_of_NODES-1;
  K[1][0]      = Cluster[n]->right->No_of_SPECIES;
  C_List[1][0] = Cluster[n]->right; 
  for(k=0; k<Cluster[n]->right->No_of_SPECIES; k++)
    PARTITION[1][0][k] = Cluster[n]->right->Species[k];
  
  C_List[1][1] = Cluster[n]->left; 
  K[1][1] = Cluster[n]->left->No_of_SPECIES;
  for(k=0; k<Cluster[n]->left->No_of_SPECIES; k++)
    PARTITION[1][1][k] = Cluster[n]->left->Species[k];
  
  /* 2th to N-1 Partitions: */ 
  for (i=2; i<N; i++) { 
    n = No_of_NODES-i; /* Always choosing the deeper node 
			  to divide */
   
    K[i][0]      = Cluster[n]->right->No_of_SPECIES;
    C_List[i][0] = Cluster[n]->right; 
    for(k=0; k<Cluster[n]->right->No_of_SPECIES; k++)
      PARTITION[i][0][k] = Cluster[n]->right->Species[k];

    C_List[i][1] = Cluster[n]->left; 
    K[i][1] = Cluster[n]->left->No_of_SPECIES;
    for(k=0; k<Cluster[n]->left->No_of_SPECIES; k++)
      PARTITION[i][1][k] = Cluster[n]->left->Species[k];

    jj = 2;
    for(j=0; j<G[i-1]; j++) {
      if( Cluster[n] != C_List[i-1][j] ) {
	C_List[i][jj] = C_List[i-1][j];
	K[i][jj]      = C_List[i-1][j]->No_of_SPECIES;
	for(k=0; k<C_List[i-1][j]->No_of_SPECIES; k++)
	  PARTITION[i][jj][k] = C_List[i-1][j]->Species[k];
     
	jj++;
      }
    }
  }
    
  for(i = 0; i<N; i++) free(C_List[i]);  
  free(C_List);
}
