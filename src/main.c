#include <stdio.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <gsl/gsl_rng.h>
#include <search.h>

#include "rgraph/src/tools.h"
#include "rgraph/src/graph.h"
#include "rgraph/src/datastruct.h"
#include "rgraph/src/modules.h"
#include "rgraph/src/bipartite.h"

/// HEADERS ///
// Functions that will be called from R
SEXP CNetcarto(SEXP nodes_in, SEXP nodes_out, SEXP weight,  SEXP r_coolingfac, SEXP r_seed,
			  SEXP r_iterfac, SEXP r_symmetric, SEXP r_auto_link, SEXP r_add_weight);
SEXP CBipartmod(SEXP nodes1, SEXP nodes2, SEXP weight, SEXP r_seed, SEXP r_iterfac,
			   SEXP r_coolingfac, SEXP r_degree_based, SEXP r_weighted, SEXP r_add_weight);

// Build networks from arrays
struct node_gra *ABuildNetwork(int E, int *nd_in, int *nd_out,
							   double *weight, int symmetric,
							   int auto_link_sw, int add_weight_sw);
struct binet *ABuildNetworkBipart(int E, int *nd_part1, int *nd_part2,
								  double *weight, int add_weight_sw);
////


/**
Netcarto function for R

Compute the modules, modularity and modularity roles of a network.

@param nodes_in,nodes_out Atomic integer vectors giving the id of the nodes.
@param weight Atomic numeric vector giving the weight of the edges.
@param r_seed Seed for the random number generator: Must be a positive
integer. 
@param r_iterfac At each temperature of the simulated annealing
(SA), the program performs fN^2 individual-node updates (involving the
movement of a single node from one module to another) and fN
collective updates (involving the merging of two modules and the split
of a module). The number "f" is the iteration factor.
@param r_symmetric If !=0 all edges a->b are copied b->a.
@param r_coolingfac Temperature cooling factor. 
@param r_auto_link If !=0 allows self looping edges a->a
@param r_add_weight If !=0 weights are summed if the edge already exist.
@return A list containing the 1) module, 2) within module z-score, 3)
        participation coefficient for each node and 4) modularity.
 */
SEXP CNetcarto(SEXP nodes_in, SEXP nodes_out, SEXP weight,  SEXP r_coolingfac, SEXP r_seed,
			  SEXP r_iterfac, SEXP r_symmetric, SEXP r_auto_link, SEXP r_add_weight){
 
  // Arguments
  int E = LENGTH(nodes_in); // Number of edges
  double iterfac = REAL(r_iterfac)[0];
  double coolingfac = REAL(r_coolingfac)[0];  
  long seed = INTEGER(r_seed)[0];
  int symmetric = INTEGER(r_symmetric)[0];
  int auto_link = INTEGER(r_auto_link)[0];
  int add_weight = INTEGER(r_add_weight)[0];

  struct node_gra *network = NULL;
  struct node_lis *p = NULL;
  struct group *part = NULL;
  struct group *g = NULL;
  gsl_rng *rand_gen;

  double Tf = 0.0;
  int N, mod_nb = 0;
  SEXP ans, module, z, P, modularity;

  
  //// RANDOM NUMBER GENERATOR INITIALIZATION
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);
  
  //// READ INPUT AND BUILD GRAPH 
 
  network = ABuildNetwork(E, INTEGER(nodes_in), INTEGER(nodes_out),
						  REAL(weight),
						  symmetric, auto_link, add_weight);
  N = CountNodes(network);
  
  //// BUILD OUTPUT STRUCTURE 
  PROTECT(ans = allocVector(VECSXP,4));
  PROTECT(module = allocVector(INTSXP,N));
  PROTECT(z = allocVector(REALSXP,N));
  PROTECT(P = allocVector(REALSXP,N));
  PROTECT(modularity = allocVector(REALSXP,1));
  SET_VECTOR_ELT(ans,0,module);
  SET_VECTOR_ELT(ans,1,z);
  SET_VECTOR_ELT(ans,2,P);
  SET_VECTOR_ELT(ans,3,modularity);

  //// COMPUTE RESULTS
 
  part = SACommunityIdent(network,
						  2.0 / (double)N, Tf, coolingfac,
						  iterfac, 0, 'o', 1, 'n', rand_gen);
 
  
  //// RETURN RESULTS
  // Get partition modularity. 
  REAL(modularity)[0] = Modularity(part);

  // Get the node metrics. 
  g = part;
  while ((g = g->next) != NULL) { // Modules
	mod_nb++;
    p = g->nodeList;
    while ((p = p->next) != NULL) { // Nodes in this module.
	  INTEGER(module)[p->node] = mod_nb;
      REAL(P)[p->node] = ParticipationCoefficient(p->ref);
      REAL(z)[p->node] = WithinModuleRelativeDegree(p->ref, g);
	}
  }
  
  /// Free Memory
  gsl_rng_free(rand_gen);
  RemovePartition(part);
  RemoveGraph(network);
  UNPROTECT(5);  
  return ans;
}

/**
Array Build Network
 
Build a network from a list of edges. Each edge i is defined by a
entry node, nd_in[i], a exit node nd_out[i] and a weight, weight[i].

The nodes id must be integer in [1, |Nodes|-1].

@param E number of edges
@param nd_in,nd_out list of size E, id of the nodes 
@param weight weight of the edges
@param symmetric if !=0 all edges a->b are copied b->a.
@param auto_link_sw allows auto link if !=0
@param add_weight if !=0, weights are summed if the edge already exist.
 */
struct node_gra *ABuildNetwork(int E,
							   int *nd_in,
							   int *nd_out,
							   double *weight,
							   int symmetric,
							   int auto_link_sw,
							   int add_weight_sw){
  struct node_gra *n1=NULL, *n2=NULL;
  struct node_gra *root=NULL, *last_add=NULL;
  void *node_dict=NULL;
  struct node_tree *n_tree=NULL, *ntree1=NULL, *ntree2=NULL;
  int N = 0;
  char buffer[255];
  
  // Create the header of the graph
  root = last_add = CreateHeaderGraph();

  // Get the number of nodes
  for (int i = 0; i < E; ++i)
	{
	  if(nd_in[i]>N)
		N = nd_in[i];
	  if(nd_out[i]>N)
		N = nd_out[i];
	}
  N ++; //Because node index start at 0.

  // Add the nodes 
  for (int i=0; i<N; ++i)
	{
	  sprintf (buffer,"nd_%d",i);
	  last_add = CreateNodeGraph(last_add, buffer);
	  //printf ("Create node %d, num:  %d, label %s\n", i, last_add->num,last_add->label);
	}

  // Add the edges 
  for (int i = 0; i < E; ++i)
	{
	  //printf ("### Add edge %d/%d between %d and %d \n",i,E,nd_in[i],nd_out[i]);
	  n1 = GetNode(nd_in[i],root);
	  n2 = GetNode(nd_out[i],root);
	  //printf ("Got nodes %d (%s) and %d (%s) \n",n1->num,n1->label,n2->num,n2->label);
	  AddAdjacency(n1, n2, auto_link_sw, add_weight_sw, weight[i], 0);
	  if (symmetric != 0 && n1 != n2)
		  AddAdjacency(n2, n1, auto_link_sw, add_weight_sw, weight[i], 0);
	  //printf ("OK.\n");
	}
  return root;
}


/**
Array Build Bipartite Network
Build a bipartite network from a list of edges. Each edge i is defined
by node in the first part 1, nd_part1[i], a node in part 2 nd_part2[i]
and a weight, weight[i].

The nodes id must be integer in [0, |Nodes|-1].

@param E number of edges
@param nd_part1,nd_part2 list of size E, id of the nodes 
@param weight weight of the edges
@param add_weight if !=0, weights are summed if the edge already exist.
*/
struct binet *ABuildNetworkBipart(int E,
								  int *nd_part1,
								  int *nd_part2,
								  double *weight,
								  int add_weight_sw)
{
  struct node_gra *root1=NULL, *root2=NULL;
  struct node_gra *last1=NULL, *last2=NULL;
  struct node_gra *n1=NULL, *n2=NULL;
  struct binet *net=NULL;
  int N1 = 0, N2 = 0;
  char buffer[255]; 
    
  // Initialize the networks
  net = CreateBipart();
  last1 = root1 = CreateHeaderGraph();
  last2 = root2 = CreateHeaderGraph();
  net->net1 = root1;
  net->net2 = root2;

  // Get the number of nodes
  for (int i = 0; i < E; ++i)
	{
	  //printf ("Edge %d, (%d,%d)\n",i,nd_part1[i],nd_part2[i]);
	  if(nd_part1[i]>N1)
		N1 = nd_part1[i];
	  if(nd_part2[i]>N2)
		N2 = nd_part2[i];
	}
  N1 ++; //Because node index start at 0.
  N2 ++; //Because node index start at 0.
  //printf ("N1 %d N2 %d\n",N1,N2 );

  // Add the nodes 
  for (int i=0; i<N1; ++i)
	{
	  sprintf (buffer,"nd_1_%d",i);
	  last1 = CreateNodeGraph(last1, buffer);
	  //printf ("Create node %d, num:  %d, label %s\n", i, last1->num,last1->label);
	}
  for (int i=0; i<N2; ++i)
	{
	  sprintf (buffer,"nd_2_%d",i);
	  last2 = CreateNodeGraph(last2, buffer);
	  //printf ("Create node %d, num:  %d, label %s\n", i, last2->num,last2->label);
	}

  // Add the edges
  for (int i = 0; i < E; ++i)
	{
	  n1 = GetNode(nd_part1[i],root1);
	  n2 = GetNode(nd_part2[i],root2); 
	  AddAdjacency(n1, n2, 0, add_weight_sw, weight[i], 0);
	  AddAdjacency(n2, n1, 0, add_weight_sw, weight[i], 0);
	}

  return net;
}

/**
Bipartmod function for R

Compute the modules, modularity and modularity roles of a bipartite network.

@param nodes1,nodes2 Atomic integer vectors giving the id of the nodes.
@param weight Atomic numeric vector giving the weight of the edges.
@param weighted_modularity if True the weighted formula for modularity is used.
@param r_seed Seed for the random number generator: Must be a positive
integer. 
@param r_iterfac At each temperature of the simulated annealing
(SA), the program performs fN^2 individual-node updates (involving the
movement of a single node from one module to another) and fN
collective updates (involving the merging of two modules and the split
of a module). The number "f" is the iteration factor.
@param r_coolingfac Temperature cooling factor. 
@param r_degree_based If true use the degree based formula.
@param r_weighted If true use the weighted version of the formula.
@return A list containing the 1) module, 2) within module z-score, 3)
        participation coefficient for each node and 4) modularity.
 */
SEXP CBipartmod(SEXP nodes1, SEXP nodes2, SEXP weight, SEXP r_seed,
				SEXP r_iterfac, SEXP r_coolingfac,
				SEXP r_degree_based, SEXP r_weighted, SEXP r_add_weight){

  SEXP ans, module, z, P, modularity;

  struct binet *network = NULL;
  struct node_lis *n = NULL;
  struct group *modules = NULL;
  struct group *g = NULL;
  struct node_gra *projected = NULL;
  gsl_rng *rand_gen;
  
  double Ti = 0.0, Tf = 0.0 ;
  int  N1,N2, mod_nb = 0;

  // Argument
  int E = LENGTH(nodes1); // Number of edges
  long seed = INTEGER(r_seed)[0];
  double iterfac = REAL(r_iterfac)[0];
  double Tsched = REAL(r_coolingfac)[0];
  int add_weight = INTEGER(r_add_weight)[0];
  int weighted = INTEGER(r_weighted)[0];
  int degree_based = INTEGER(r_degree_based)[0];
  
  //// RANDOM NUMBER GENERATOR INITIALIZATION
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);
  
  //// READ INPUT AND BUILD GRAPH 
  //printf("----- Building the network from %d edges.\n",E);
  network = ABuildNetworkBipart(E,INTEGER(nodes1),INTEGER(nodes2),
								REAL(weight),
								add_weight);
  N1 = CountNodes(network->net1);
  N2 = CountNodes(network->net2);
  //printf("----- The network has %d,%d nodes\n", N1,N2);
  
  //// BUILD OUTPUT STRUCTURE 
  PROTECT(ans = allocVector(VECSXP,4));
  PROTECT(module = allocVector(INTSXP,N1));
  PROTECT(z = allocVector(REALSXP,N1));
  PROTECT(P = allocVector(REALSXP,N1));
  PROTECT(modularity = allocVector(REALSXP,1));
  SET_VECTOR_ELT(ans,0,module);
  SET_VECTOR_ELT(ans,1,z);
  SET_VECTOR_ELT(ans,2,P);
  SET_VECTOR_ELT(ans,3,modularity);

  //// COMPUTE RESULTS
  //printf("----- Simulated Annealing \n");
  Ti = 1. / (double)N1;

  if (weighted == 1)
	modules = SACommunityIdentBipartWeighted(network,
										  Ti, Tf, Tsched, iterfac,
										  0, 'o', 1, 'n',
										  rand_gen);
  else
	modules = SACommunityIdentBipart(network,
								  Ti, Tf, Tsched, iterfac,
								  0, 'o', 1, 'n',
								  rand_gen);
  //printf("----- Done.\n");
  
  //// RETURN RESULTS
  // Get partition modularity. 
  REAL(modularity)[0] = ModularityBipart(network,modules);
  
  projected = ProjectBipart(network);
  MapPartToNetFast(modules, projected);

  // Get the node metrics. 
  g = modules;
  while ((g = g->next) != NULL) { // Modules
	//printf ("--Module %d\n",mod_nb);
	mod_nb++;
    n = g->nodeList;
    while ((n = n->next) != NULL) { // Nodes in this module.
	  //printf ("-----Node %d\n",n->node);
	  INTEGER(module)[n->node] = mod_nb;
	  if (degree_based == 1){
      REAL(P)[n->node] = ParticipationCoefficient(n->ref);
      REAL(z)[n->node] = WithinModuleRelativeDegree(n->ref, g);
	  }
	  else{
		REAL(P)[n->node] = WeightedParticipationCoefficient(n->ref,modules);
		REAL(z)[n->node] = WithinModuleRelativeStrength(n->ref, g);
	  }
	}
  }
  
  /// Free Memory
  gsl_rng_free(rand_gen);
  RemovePartition(modules);
  RemoveBipart(network);
  RemoveGraph(projected);
  UNPROTECT(5);  
  return ans;
}
