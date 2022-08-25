#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "isaac64.h"
#include "gamma_dist_rng.h"
#include "rng.h"

#include "connectors.h"

double calculate_shortest( int i, int j, double dxs, double dxr, double L ) 
{
  // Calculates shortest distance bewteen two points on a ring with length L 
  // if distance bewteen each point is equally spaced, dx
  
  double opt1 = fabs(i*dxs - j*dxr);
  double opt2 = (L - opt1);
  return (( opt1 < opt2 ) ?  opt1 : opt2 );

}


/* function: calculate shortest distance between two nodes on a circle. */
/* warning: distance given has satisfy  distance/dxj > -N */
int calculate_closest( double distance, double position, double dxj, int N )
{

  double jdis;
  int j;

  jdis = position + distance;
  j = round ( jdis/dxj );
  
  j = (j+N)%N;

  return j;
  
}

/* gaussian connect */
int gaussian_connect( struct network_parameters * net, struct simulation_parameters * sim,
		      struct synapse_parameters * syn, struct neuron * network, struct rng_state * rng )
{

  /* init */
  /* special for calculating netsim $% */
  bool *temparr;

  /* same for calculating and saved netsim */
  int ii, jj;
  int Ne, Ni;
  int count, count_connections;
  double x_position, distance, dx_exc, dx_inh;

  double max_distance = 0;

  /* calculate */
  Ne = round( .8*net->N );
  Ni = round( .2*net->N );
  dx_exc = net->L / Ne;
  dx_inh = net->L / Ni;

  /* temporarily save connections in boolean size N, then translate to sparse */
  /* like this it is way faster to check for double connections */
  temparr = (bool *) malloc( net->N * sizeof(bool) ); /* $% */

  for ( ii = 0 ; ii < net->N ; ii++ )
  {

    x_position = (*(network+ii)).x_position;

    /* reset temparr */
    memset( temparr, (bool) 0, (size_t) net->N * sizeof(bool) );

    count = (int) round(.8*net->K);
    count_connections = 0;
    while ( count > 0 )
    {

      /* draw gaussian distance */
      distance = rng_gauss( rng ) * net->sigma_space;

      /* within range? */
      if ( fabs(distance) > (.5 * net->L) ) {  continue;  }

      jj = calculate_closest( distance, x_position, dx_exc , Ne );

      /* no self- or double connections */
      if ( ((ii==jj) && (*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) { max_distance = fabs(distance); }

      /* connect */
      *( ( ( *(network+ii) ).K ) + count_connections ) = jj;

      *(temparr+jj) = 1; /* $% */
      count--;
      count_connections++;

    }

    count = (int) round(.2*net->K);
    while ( count > 0 )
    {

      distance = rng_gauss( rng) * net->sigma_space;

      if ( fabs(distance) > (.5 * net->L) ) {  continue;  }

      jj = Ne + calculate_closest( distance, x_position, dx_inh, Ni );

      /* no self- or double connections */
      if ( ((ii==jj) && !(*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) { max_distance = fabs(distance); }
      
      /* connect */
      *( ( ( *(network+ii) ).K ) + count_connections ) = jj;

      *(temparr+jj) = 1; /* $% */
      count--;
      count_connections++;

    }

  }

  /* assign buffer length based on max distance */
  sim->buffer_length = (int) round( ( syn->synapse_delay + (max_distance / (syn->min_conduction_speed)) ) / sim->dt ) + 1; 

  /* cleanup */
  free( temparr );
  return NETSIM_NOERROR;

}

/* gaussian connect */
int gaussian_connect_2D( struct network_parameters * net, struct simulation_parameters * sim,
			 struct synapse_parameters * syn, struct neuron * network, struct rng_state * rng )
{

  /* init */
  /* special for calculating netsim $% */
  bool *temparr;

  /* same for calculating and saved netsim */
  int ii, jj, c1, c2;
  int Ne, Ni, NrowE, NrowI;
  int count, count_connections;
  double x_position, y_position, distance_x, distance_y, distance, dx_exc, dx_inh;

  double max_distance = 0;

  /* calculate */
  Ne = round( .8*net->N ); Ni = round( .2*net->N );
  NrowE = floor( sqrt(Ne) ); NrowI = floor( sqrt(Ni) );
  dx_exc = net->L / NrowE; dx_inh = net->L / NrowI;

  /* temporarily save connections in boolean size N, then translate to sparse */
  /* like this it is way faster to check for double connections */
  temparr = (bool *) malloc( net->N * sizeof(bool) ); /* $% */

  for ( ii = 0 ; ii < net->N ; ii++ )
  {

    x_position = (*(network+ii)).x_position; y_position = (*(network+ii)).y_position;

    /* reset temparr */
    memset( temparr, (bool) 0, (size_t) net->N * sizeof(bool) );

    count = (int) round(.8*net->K);
    count_connections = 0;
    while ( count > 0 )
    {

      /* draw gaussian distance */
      distance_x = rng_gauss( rng ) * net->sigma_space;
      distance_y = rng_gauss( rng ) * net->sigma_space;
      distance = sqrt( pow(distance_x,2) + pow(distance_y,2) );

      /* within range? */
      if ( fabs(distance) > (M_SQRT2 * .5 * net->L) ) {  continue;  }

      c1 = calculate_closest( distance_x, x_position, dx_exc, NrowE ); /* x-coordinate */
      c2 = calculate_closest( distance_y, y_position, dx_exc, NrowE ); /* y-coordinate */
      jj = c1 + ( c2 * NrowE ); 

      /* no self- or double connections */
      if ( ((ii==jj) && (*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) {  max_distance = fabs(distance); }

      /* connect */
      *( ( ( *(network+ii) ).K ) + count_connections ) = jj;

      *(temparr+jj) = 1; /* $% */
      count--; count_connections++;

    }

    count = (int) round(.2*net->K);
    while ( count > 0  )
    {

      /* draw gaussian distance */
      distance_x = rng_gauss( rng ) * net->sigma_space;
      distance_y = rng_gauss( rng ) * net->sigma_space;
      distance = sqrt( pow(distance_x,2) + pow(distance_y,2) );

      if ( fabs(distance) > (M_SQRT2 * .5 * net->L) ) {  continue;  }

      c1 = calculate_closest( distance_x, x_position, dx_inh, NrowI ); /* x-coordinate */
      c2 = calculate_closest( distance_y, y_position, dx_inh, NrowI ); /* y-coordinate */
      jj = Ne + c1 + ( c2 * NrowI );

      /* no self- or double connections */
      if ( ((ii==jj) && !(*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) { max_distance = fabs(distance); }
      
      /* connect */
      *( ( ( *(network+ii) ).K ) + count_connections ) = jj;

      *(temparr+jj) = 1; /* $% */
      count--;
      count_connections++;

    }

  }

  /* assign buffer length based on max distance */
  sim->buffer_length = (int) round( ( syn->synapse_delay + (max_distance / (syn->min_conduction_speed)) ) / sim->dt ) + 1; 

  /* cleanup */
  free( temparr );
  return NETSIM_NOERROR;

}

/* function: Erdos-Renyi random graph */
int random_connect( struct network_parameters * net, struct neuron * network, struct rng_state * rng )
{

  /* init */
  register int ii, jj;
  uint32_t target = 0;

  /* loop over network */
  for ( ii = 0 ; ii < net->N ; ii++ )
  {

    for ( jj = 0 ; jj < net->K ; jj++ )
    {

      /* connect */
      target = rng_uint32( rng ) % net->N;
      *( ( ( *(network+ii) ).K ) + jj ) = target;

    }

  }

  return NETSIM_NOERROR;

}

/* function: load test conections from file */
int test_connect( struct parameters * p, struct neuron * network )
{

  /* init */
  int ii, jj, id;
  FILE * file;
  file = fopen( p->simulation.optional_connection_path, "r" );
  int r = 0;
  
  /* loop over network */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    for ( jj = 0 ; jj < p->network.K ; jj++ ) 
    {  
      
      /* get connection from file */
      r = fscanf( file, "%d\n", &id );  
      if ( r < 0 ) {
	      //throw and error?
	      //fclose( file );
	      //return NETSIM_ERROR;
      }
      
      /* connect */
      *( ( ( *(network+ii) ).K ) + jj ) = id;

    }

  }

  /* clean up */
  fclose( file );
  return NETSIM_NOERROR;

}

/* function: random rewiring of EXCITATORY network connections */
int random_rewiring( struct parameters * p, struct neuron * network, struct rng_state * rng )
{

  /* init */
  int ii, jj, Ne;
  uint32_t target = 0;
  
  /* number of excitatory neurons */
  Ne = round( .8*p->network.N );

  /* loop over network */
  for ( ii = 0 ; ii < Ne ; ii++ )
  {

    for ( jj = 0 ; jj < p->network.K ; jj++ ) 
    {

      if ( rng_dbl64( rng ) < p->network.rewiring_probability )  
      {

        /* select random target */
        target = rng_uint32( rng ) % p->network.N;

        /* write random target */
        *( ( ( *(network+ii) ).K ) + jj ) = target;

      }
      
    }

  }

  /* assign buffer length as in random graph simulation */
  p->simulation.buffer_length = (int) ( ( p->synapse.synapse_delay + ((.5 * p->network.L) / (p->synapse.min_conduction_speed)) ) / p->simulation.dt ) + 1;     

  /* clean up */
  return NETSIM_NOERROR;

}
