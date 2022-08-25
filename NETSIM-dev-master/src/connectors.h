#ifndef _CONNECTORS_H_
#define _CONNECTORS_H_

#include "simulation.h"

double calculate_shortest( int i, int j, double dxs, double dxr, double L );
int calculate_closest( double distance, double position, double dxj, int N );

int gaussian_connect( struct network_parameters * net, struct simulation_parameters * sim,
		      struct synapse_parameters * syn, struct neuron * network, struct rng_state * rng );

int gaussian_connect_2D( struct network_parameters * net, struct simulation_parameters * sim,
			 struct synapse_parameters * syn, struct neuron * network, struct rng_state * rng );

int random_connect( struct network_parameters * net, struct neuron * network, struct rng_state * rng );

int test_connect( struct parameters * p, struct neuron * network );

int random_rewiring( struct parameters * p, struct neuron * network, struct rng_state * rng );


#endif
