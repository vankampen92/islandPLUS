struct bin_tree {
  int Sp;                 /* Species Label : {1, 2, 3, ... } */
  double d;               
  int No_of_SPECIES;  
  int * Species;          /* List of species embraced by this node: { 0, 5, ...}*/
  struct bin_tree * right, * left;
};
typedef struct bin_tree node;

double Average_Node_Distance( node * N_1, node * N_2, 
			      double ** D, int S );  

void Print_Triangular_Matrix( double ** distance, int N );


void upgma_clustering(double ** D, int No_of_SPECIES,  
		      double ** distance, int * n, 
		      node ** T, int * No_of_NODES, int * Node_List ); 

void insert(node ** tree, double val); 

void print_preorder(node * tree); 

void print_inorder(node * tree); 

void print_postorder(node * tree);

void deltree(node * tree); 

node * search(node ** tree, double val);

void upgma_cluster_to_partition ( node ** Cluster, int N,
                                  int *** PARTITION, int * G, int ** K );

void UPGMA_CLUSTERING_PARTITION( double ** d, int No_of_SPECIES , 
				 int *** PARTITION, int * G, int ** K ); 
