#ifndef _NETSIM_H_
#define _NETSIM_H_

#include "simulation.h"

int initialize_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng );
int allocate_input_buffers( struct parameters * p, struct neuron * network );
int initialize_sustained_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng );
int initialize_2D_simulation( struct parameters * p, struct neuron * network );
int run_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng );

#endif
