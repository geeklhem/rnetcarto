#include <stdio.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <gsl/gsl_rng.h>
#include <search.h>

#include "rgraph/src/sannealing.h"
#include "rgraph/src/io.h"
#include "rgraph/src/partition.h"
#include "rgraph/src/fillpartitions.h"

#define EPSILON 1.e-6

// Function that will be called from R
SEXP netcarto_binding(SEXP nodes_in, SEXP nodes_out, SEXP weight,
					  SEXP r_N, SEXP r_bipartite, SEXP r_clustering,
					  SEXP r_roles, SEXP r_diagonal_term,
					  SEXP r_coolingfac, SEXP r_seed, SEXP r_iterfac){

  // Arguments
  int E = LENGTH(nodes_in); // Number of edges
  double iterfac = REAL(r_iterfac)[0];
  double coolingfac = REAL(r_coolingfac)[0];
  long seed = INTEGER(r_seed)[0];
  int N = INTEGER(r_N)[0];
  int bipartite = INTEGER(r_bipartite)[0];
  int clustering = INTEGER(r_clustering)[0];
  int diagonal_term = INTEGER(r_diagonal_term)[0];
  int roles = INTEGER(r_roles)[0];
  Partition *part = NULL;
  AdjaArray *adj = NULL;
  gsl_rng *rand_gen;
  SEXP ans, module, z, P, modularity;
  unsigned int Ngroups=0;
  unsigned int i;
  int err = 0;
  //// RANDOM NUMBER GENERATOR INITIALIZATION
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  //// READ INPUT AND BUILD GRAPH
  if (!bipartite){
	if (Ngroups == 0) Ngroups = N;
	part = CreatePartition(N,Ngroups);
	adj = CreateAdjaArray(N,E);
	err = EdgeListToAdjaArray(INTEGER(nodes_in), INTEGER(nodes_out),
							  REAL(weight),	adj, part, 1);
	if (err != 0)
	  error("Initialisation error.\n");
  } else {
	ProjectBipartEdgeList(INTEGER(nodes_in), INTEGER(nodes_out), REAL(weight), E,
						  &part, &adj);
  }
  
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
  if (clustering){
		double Ti = 1. / (double)N;
		double Tf = 0.0;
		unsigned int nochange_limit=25;
		double proba_components = .5;
		AssignNodesToModules(part,rand_gen);
		err = GeneralSA(&part, adj, iterfac,
			  			Ti, Tf, coolingfac,
			  			proba_components, nochange_limit,
			  			rand_gen);
		if (err != 0)
		  error("Simulated annealing error.\n");
		CompressPartition(part);
		// Get partition modularity.
		REAL(modularity)[0] = PartitionModularity(part,adj,diagonal_term);
		for (i=0;i<part->N;i++)
			INTEGER(module)[i] = part->nodes[i]->module;
  }

  if(roles){
		double *connectivity, *participation;
		connectivity = (double*) calloc(part->N,sizeof(double));
		participation = (double*) calloc(part->N,sizeof(double));
		PartitionRolesMetrics(part, adj, connectivity, participation);
		for (i=0;i<part->N;i++){
				 REAL(z)[i] = connectivity[i];
				 REAL(P)[i] = participation[i];
		}
		free(connectivity);
		free(participation);
  }

  /// Free Memory
  gsl_rng_free(rand_gen);
	FreeAdjaArray(adj);
	FreePartition(part);
  UNPROTECT(5);
  return ans;
}
