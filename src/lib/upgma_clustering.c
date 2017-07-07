#include "HEADERS.h"
#include "da_IBT_Functions.h"

// #defined VERBOSE
double Average_Node_Distance( node * N_1, node * N_2, 
			      double ** D, int S )
{
  int i,j;
  int sp_i, sp_j; 

  double Ave_D;
  Ave_D = 0.0;
  for(i=0; i<N_1->No_of_SPECIES; i++)
    for(j=0; j<N_2->No_of_SPECIES; j++){
      sp_i = N_1->Species[i]; 
      sp_j = N_2->Species[j];

      if(sp_i >= S || sp_j >= S) error(0,0,"Program aborted");
      
      if ( sp_i > sp_j ) Ave_D += D[sp_i][ sp_j];
      else               Ave_D += D[sp_j][ sp_i];
    }
  
  Ave_D = Ave_D / (double)(N_1->No_of_SPECIES * N_2->No_of_SPECIES);

  return(Ave_D);
}

void Print_Triangular_Matrix( double ** distance, int N )
{
  int i,j;
  
  for(i=1; i<N; i++) {
    for(j=0; j<i; j++) 
       Rprintf(" [ d(%d, %d) = %g ]", i,j, distance[i][j] );
    Rprintf("\n");
  }
}

void upgma_clustering(double ** D, int No_of_SPECIES,  
		      double ** distance, int * n, 
		      node ** T, int * No_of_NODES, int * Node_List )
{
  /* Input:
     . D[][]: Initial distance matrix. Its entries remain unchanged
     through the whole algorithm.
     . S: Number of rows and columns of the initial distance matrix.
     . distance[][]: Distances between clusters as new clusters arise
     . (*n) Number of rows/columns of distance[][],  which is reduced by 
     one at each recurrent call
     . (* No_of_NODES): total number of nodes of the tree, which is 
     increased by a new inner node at each recurrent call.
     . Node_List[], list of nodes to compare in the "shrinking" input 
     distance[][] matrix, from first to last column/row. 
     
     Output:
     T: the tree countaning the cluster analysis. The root node is always
     node * root = T[(*No_of_NODES)-1]
  */
  int N, i, j, k, l, m, i_MIN, j_MIN;
  double x;

  N = (*n);
  if( N > 1) {
    double ** temp  = (double **)calloc( N, sizeof(double *) );
    for(i=1; i<N; i++) temp[i] = (double *)calloc( i, sizeof(double) );
    for(i=1; i<N; i++) 
      for(j=0; j<i; j++) temp[i][j] = distance[i][j]; 

    // Print_Triangular_Matrix( temp, N );
    // getchar();

    int * Index_List = (int *)calloc( N, sizeof(int) );
    int * List = (int *)calloc( N, sizeof(int) );
    for(i=0; i<N; i++) List[i] = Node_List[i]; 

    i_MIN = 1;
    j_MIN = 0;
    x     = temp[i_MIN][j_MIN]; 
    for(i=1; i<N; i++) 
      for(j=0; j<i; j++) { 
	if( temp[i][j] < x ) { 
	  i_MIN = i; 
	  j_MIN = j;
	}
	x = MIN( temp[i][j], x );
      }

    T[ * No_of_NODES ]->d     = temp[i_MIN][j_MIN]/2.0;
    T[ * No_of_NODES ]->left  = T[List[i_MIN]];
    T[ * No_of_NODES ]->right = T[List[j_MIN]];

    int S = 0;
    for(i=0; i<T[List[i_MIN]]->No_of_SPECIES; i++) T[ * No_of_NODES ]->Species[S++] = T[List[i_MIN]]->Species[i];
    for(i=0; i<T[List[j_MIN]]->No_of_SPECIES; i++) T[ * No_of_NODES ]->Species[S++] = T[List[j_MIN]]->Species[i];
    T[ * No_of_NODES ]->No_of_SPECIES =  S; 
    if ( S != (T[List[i_MIN]]->No_of_SPECIES + T[List[j_MIN]]->No_of_SPECIES) ) error(0,0,"Program aborted"); 
      
    Node_List[0] = (* No_of_NODES);
    
    k = 1; m = 0;
    for(i=0; i<N; i++) {
      if (i != i_MIN && i != j_MIN) {
	
	distance[k][0]  = Average_Node_Distance( T[ * No_of_NODES ], T[ List[ i ] ], 
						 D, No_of_SPECIES );
       
	Node_List[k]    =  List[i];
	Index_List[m++] =  i;
	k++;
      }
    }

    if (m != N-2) error("Program aborted");

    k = 1; 
    for(i=0; i<m; i++) {
      l = 1;
      for(j=i+1; j<m; j++) {
	if ( Index_List[i] > Index_List[j] ) x = temp[ Index_List[i] ][ Index_List[j] ];
	else                                 x = temp[ Index_List[j] ][ Index_List[i] ];
	
	distance[k+l][k] = x;
	l++;
      }
      k++;
    }

    * n = N-1;
    (* No_of_NODES) = (* No_of_NODES) + 1;

#if defined VERBOSE    
    Print_Triangular_Matrix( distance, * n );
    Rprintf("--------------------------------\n\n");
#endif
    upgma_clustering( D, No_of_SPECIES, 
		      distance, n, T, No_of_NODES, Node_List );
        
    for(i=1; i<N; i++) free(temp[i]);
    free(temp);
    free(List); free(Index_List);
  }
}

void insert(node ** tree, double val)
{
    node *temp = NULL;
    if(!(*tree))
    {
        temp = (node *)malloc(sizeof(node));
        temp->left = temp->right = NULL;
        temp->d    = val;
        *tree      = temp;
        return;
    }

    if(val < (*tree)->d)
    {
        insert(&(*tree)->left, val);
    }
    else if(val > (*tree)->d)
    {
        insert(&(*tree)->right, val);
    }

}

void print_preorder(node * tree)
{
  int i;

    if (tree)
    {
        Rprintf("d = %g\t",tree->d);
	for (i=0; i< tree->No_of_SPECIES; i++) Rprintf("%d ", tree->Species[i]);
	Rprintf("\n");
        print_preorder(tree->left);
        print_preorder(tree->right);
    }

}

void print_inorder(node * tree)
{
  int i;
    if (tree)
    {
        print_inorder(tree->left);
        Rprintf("d = %g\t",tree->d);
	for (i=0; i< tree->No_of_SPECIES; i++) Rprintf("%d ", tree->Species[i]);
	Rprintf("\n");
        print_inorder(tree->right);
    }
}

void print_postorder(node * tree)
{
  int i;

    if (tree)
    {
        print_postorder(tree->left);
        print_postorder(tree->right);
        Rprintf("d = %g\t",tree->d);
	for (i=0; i< tree->No_of_SPECIES; i++) Rprintf("%d ", tree->Species[i]);
	Rprintf("\n");
    }
}

void deltree(node * tree)
{
    if (tree)
    {
        deltree(tree->left);
        deltree(tree->right);
	free(tree->Species);
        free(tree);
    }
}

node * search(node ** tree, double val)
{
    if(!(*tree))
    {
        return NULL;
    }

    if(val < (*tree)->d)
    {
        search(&((*tree)->left), val);
    }
    else if(val > (*tree)->d)
    {
        search(&((*tree)->right), val);
    }
    else if(val == (*tree)->d)
    {
        return *tree;
    }
}

