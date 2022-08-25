#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include "NIDAQmx.h"
#include "isaac64.h"
#include "gamma_dist_rng.h"
#include "rng.h"

#include "connectors.h"
#include "netsim.h"


#include "interface_c.h"



/* define constants */
#define NETSIM_DIALOG_NOCOLOR__				"\x1b[0m"
#define NETSIM_DIALOG_INVCOLOR__			"\x1b[7m"

#define NETSIM_DIALOG_BG_COLOR__RED__		"\x1b[41;1m"

#define NETSIM_DIALOG_FG_COLOR__RED__		"\x1b[31;1m"
#define NETSIM_DIALOG_FG_COLOR__GREEN__		"\x1b[32;1m"
#define NETSIM_DIALOG_FG_COLOR__YELLOW__	"\x1b[33;1m"
#define NETSIM_DIALOG_FG_COLOR__BLUE__		"\x1b[34;1m"



double linspace_val( double start, double stop, int npts, int val_index )
{

	/* init */
	double val, interval;

	/* create output */
	interval = ( stop - start ) / ( npts - 1 );
	val = start + (double)(val_index)*interval;

	/* return */
	return val;

}

int sub2ind( unsigned int x1, unsigned int x2, unsigned int x3, unsigned int n1, unsigned int n2, unsigned int n3 )
{

	/* init */
	unsigned int index = 0;

	/* error checking */
	assert( (x1<=n1) && (x2<=n2) && (x3<=n3) );

	/* calculate linear index output */
	index = x1 + (x2-1)*n2 + (x3-1)*n1*n2;

	/* return */
	return index;

}

int ind2sub( unsigned int index, unsigned int dim, unsigned int n1, unsigned int n2, unsigned int n3 )
{

	/* init */
	unsigned int subscript_dim = 0;

	/* error checking */
	assert( index <= ((n1+1)*(n2+1)*(n3+1)) );	

	/* calculate subscript output */
	if ( dim == 1 ) {  subscript_dim = index % n1;  }
	else if ( dim == 2 ) {  subscript_dim = floor( index / n1 );  }
	else if ( dim == 3 ) {  subscript_dim = floor( index / (n1*n2) ); }

	/* return */
	return subscript_dim;

}


/* NETSIM startup dialog */
int startup_dialog()
{

	printf( "\n\n" );

	printf( "===============================\n" );
	printf( "============ %sNETSIM%s ===========\n", NETSIM_DIALOG_FG_COLOR__RED__, NETSIM_DIALOG_NOCOLOR__ );
	printf( "===============================\n" );

	printf( "\n\n" );

	return NETSIM_NOERROR;

}

int parse_input_file( struct parameters * p )
{

	/* init */
	FILE * file;
	char * line = NULL;
	char param[2048];
	char str[2048];
	double value, start, stop;
	int npts, val_index, loop_variable = 0, n1 = 0, n2 = 0, n3 = 0, n_loop_variables = 0;
	size_t len = 0;

	/* open file */
	file = fopen( p->simulation.parameter_file_path, "r" );
	if ( file == NULL ) {  printf( "Problem loading file\n\n" ); return 1;  }

	/* count number of loop variables in file */
	while ( getline( &line, &len, file ) != EOF )
	{
		if ( *line == '@' )
		{
			n_loop_variables++;
			sscanf( line+1, "%s = %lf:%lf:%d", param, &start, &stop, &npts );
			if ( n_loop_variables == 1 ) {  n1 = npts;  }
			else if ( n_loop_variables == 2 ) {  n2 = npts;  }
			else if ( n_loop_variables == 3 ) {  n3 = npts;  }
		}
	}
	assert( n_loop_variables <= 3 ); /* too many loop variables? */
	rewind( file );

	/* read file line by line */
	while ( getline( &line, &len, file ) != EOF )
	{

		if ( (*line != '#') && (*line != '$') && (*line != '\n') && (*line != ' ') )
		/* line comment character: #, string char: $, linspace char: @ */
		{

			/* parse values */
			if ( *line == '@' )
			{
				/* sscanf, linspace, assign value */
				loop_variable++;
				sscanf( line+1, "%s = %lf:%lf:%d", param, &start, &stop, &npts );
				assert( p->simulation.job_id != INT_MIN ); /* bad job ID? */
				val_index = ind2sub( p->simulation.job_id, loop_variable, n1, n2, n3 );
				value = linspace_val( start, stop, npts, val_index );
			}
			else
			{
				/* sscanf */
				sscanf( line, "%s = %lf", param, &value );
			}		

			/* translation layer */
			if ( strcmp( param, "N" ) == 0 ) {  p->network.N = (int)round(value);  }
			else if ( strcmp( param, "K" ) == 0 ) {  p->network.K = (int)round(value);  }
			else if ( strcmp( param, "L" ) == 0 ) {  p->network.L = value;  }
			else if ( strcmp( param, "sigma_space" ) == 0 ) {  p->network.sigma_space = value;  }
			else if ( strcmp( param, "poisson_rate" ) == 0 ) {  p->network.poisson_rate = value;  }
			else if ( strcmp( param, "poisson_start_time" ) == 0 ) {  p->network.poisson_start_time = value;  }
			else if ( strcmp( param, "poisson_stop_time" ) == 0 ) {  p->network.poisson_stop_time = value;  }
			else if ( strcmp( param, "rewiring_probability" ) == 0 ) {  p->network.rewiring_probability = value;  }

			else if ( strcmp( param, "T" ) == 0 ) {  p->simulation.T = value;  }
			else if ( strcmp( param, "dt" ) == 0 ) {  p->simulation.dt = value;  }
			else if ( strcmp( param, "vm_sigma" ) == 0 ) {  p->simulation.vm_sigma = value;  }
			else if ( strcmp( param, "vm_mean" ) == 0 ) {  p->simulation.vm_mean = value;  }
			else if ( strcmp( param, "ge_sigma" ) == 0 ) {  p->simulation.ge_sigma = value;  }
			else if ( strcmp( param, "ge_mean" ) == 0 ) {  p->simulation.ge_mean = value;  }
			else if ( strcmp( param, "gi_sigma" ) == 0 ) {  p->simulation.gi_sigma = value;  }
			else if ( strcmp( param, "gi_mean" ) == 0 ) {  p->simulation.gi_mean = value;  }
			else if ( strcmp( param, "bin_size" ) == 0 ) {  p->simulation.bin_size = (int)round(value);  }
			else if ( strcmp( param, "start_record_time" ) == 0 ) {  p->simulation.start_record_time = value;  }
			else if ( strcmp( param, "stop_record_time" ) == 0 ) {  p->simulation.stop_record_time = value;  }
			else if ( strcmp( param, "record_downsample_factor" ) == 0 ) {  p->simulation.record_downsample_factor = (int)round(value);  }
			else if ( strcmp( param, "connector" ) == 0 ) {  p->simulation.connector = (int)round(value);  }
			else if ( strcmp( param, "initiator" ) == 0 ) { p->simulation.initiator = (int)round(value);  }
			else if ( strcmp( param, "output_file_code" ) == 0 ) {  p->simulation.output_file_code = (int)round(value);  }
			else if ( strcmp( param, "report_minutes" ) == 0 ) { p->simulation.report_minutes = value;  }
			else if ( strcmp( param, "cutoff_frequency" ) == 0 ) { p->simulation.cutoff_frequency = value;  }
			else if ( strcmp( param, "conduction_speed_code" ) == 0 ) { p->simulation.conduction_speed_code = value;  }
			else if ( strcmp( param, "gamma_parameter_shape" ) == 0 ) { p->simulation.gamma_parameter_shape = value;  }
			else if ( strcmp( param, "gamma_parameter_scale" ) == 0 ) { p->simulation.gamma_parameter_scale = value;  }
			else if ( strcmp( param, "save_connectivity" ) == 0 ) { p->simulation.save_connectivity = value;  }

			else if ( strcmp( param, "ge" ) == 0 ) {  p->synapse.ge = value;  }
			else if ( strcmp( param, "gi" ) == 0 ) {  p->synapse.gi = value;  }
			else if ( strcmp( param, "Ee" ) == 0 ) {  p->synapse.Ee = value;  }
			else if ( strcmp( param, "Ei" ) == 0 ) {  p->synapse.Ei = value;  }
			else if ( strcmp( param, "synapse_delay" ) == 0 ) {  p->synapse.synapse_delay = value;  }
			else if ( strcmp( param, "min_conduction_speed" ) == 0 ) {  p->synapse.min_conduction_speed = value;  }
			else if ( strcmp( param, "max_conduction_speed" ) == 0 ) {  p->synapse.max_conduction_speed = value;  }
			else if ( strcmp( param, "taue" ) == 0 ) {  p->synapse.taue = value;  }
			else if ( strcmp( param, "taui" ) == 0 ) {  p->synapse.taui = value;  }
			else if ( strcmp( param, "min_uniform_delay" ) == 0 ) {  p->synapse.min_uniform_delay = value;  }
			else if ( strcmp( param, "max_uniform_delay" ) == 0 ) {  p->synapse.max_uniform_delay = value;  }
			else if ( strcmp( param, "p_release" ) == 0 ) {  p->synapse.p_release = value;  }

			else if ( strcmp( param, "taum" ) == 0 ) {  p->neuron.taum = value; }
			else if ( strcmp( param, "vr" ) == 0 ) {  p->neuron.vr = value;  }
			else if ( strcmp( param, "vreset" ) == 0 ) {  p->neuron.vreset = value;  }
			else if ( strcmp( param, "vth" ) == 0 ) {  p->neuron.vth = value; }
			else if ( strcmp( param, "taur" ) == 0 ) {  p->neuron.taur = value;  }
			else if ( strcmp( param, "El" ) == 0 ) {  p->neuron.El = value;  }
			else if ( strcmp( param, "Ie" ) == 0 ) {  p->neuron.Ie = value;  }
			else if ( strcmp( param, "Cm" ) == 0 ) {  p->neuron.Cm = value;  }

			else {  printf( "Improper parameter read for %s.\n", param ); return NETSIM_BAD_PARAMETER;  }

		}
		else if ( *line == '$' )
		{

			/* sscanf */
			sscanf( line+1, "%s = %s", param, str );
			printf( "Reading %s: %s \n", param, str );

			/* translation layer */
			if ( strcmp( param, "optional_connection_path" ) == 0 ) {  strcpy( p->simulation.optional_connection_path, str );  }
			else if ( strcmp( param, "output_path" ) == 0 ) {
				if ( p->simulation.output_path_override == 0 ) { strcpy( p->simulation.output_path, str ); }
			}
			else {  printf( "Improper parameter read (string).\n" ); return NETSIM_BAD_PARAMETER;  }

		}

	}

	/* if part of a PARAMETER SCAN, modify OUPUT FILE CODE with JOB ID */
	if ( n_loop_variables > 0 ) {  p->simulation.output_file_code = p->simulation.job_id;  }

	/* cleanup */
	fclose( file );
	free( line );

	/* return */
	return NETSIM_NOERROR;

}

int parse_inputs( int argc, char ** argv, struct parameters * p )
{

	/* init options */
	int option = 0;

	/* parse inputs */
	while ( ( option = getopt( argc, argv, "f:j:o:" ) ) != -1 ) {
	    switch (option) {
	        case 'f'  : strcpy( p->simulation.parameter_file_path, optarg );
	                    break;
	        case 'j'  : p->simulation.job_id = (unsigned int) atoi(optarg);
			    break;
	        case 'o'  : strcpy( p->simulation.output_path, optarg );
	        	    p->simulation.output_path_override = 1;
			    printf("Command line output_path: %s overriding parameter file\n", optarg);
			    break;
	        default   : printf( "USAGE: \"netsim -j job_id -f parameter_file -o output_directory\"\n" );
	                    return NETSIM_BAD_PARAMETER;
	    }
	}

	/* return */
	return NETSIM_NOERROR;

}

int load_array_from_file( double * array, int length, char * file_path )
{

	/* init */
	int nbytes;
	FILE * file;

	/* open file */
	file = fopen( file_path, "rb" );
	if ( file == NULL ) {  printf( "Problem loading file\n\n" ); return NETSIM_BAD_PARAMETER;  }

	/* read file */
	nbytes = fread( array, sizeof(double), length, file );
	if ( nbytes != length ) {  printf( "Problem reading file\n\n" ); return NETSIM_BAD_PARAMETER;  }

	/* clean up */
	fclose( file );
	return NETSIM_NOERROR;

}

/* function: load file connections (binary) */
int file_connect( struct parameters * p, struct neuron * network )
{

  /* init */
  int ii, K;
  FILE * file;

  /* prepare */
  K = p->network.K;
  file = fopen( p->simulation.optional_connection_path, "rb" );

  /* loop over network */
  for ( ii = 0; ii < p->network.N ; ii++ )
  {

    if ( fread( ( *(network+ii) ).K, sizeof(uint32_t), K, file ) != 0 )
    {

      printf("Error in reading synaptic connections (ii=%d)\n", ii);
      return NETSIM_ERROR;

    }

  }

  /* clean up */
  fclose( file );
  return NETSIM_NOERROR;

}

/* function: intialize for test simulations */
int initialize_test_simulation( struct parameters * p, struct neuron * network )
{

  int ii;

  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    ( *(network+ii) ).excitatory = 1;
    if ( ii == 0 ) {  ( *(network+ii) ).v = -49e-3;  }
    else {  ( *(network+ii) ).v = p->neuron.vreset;  }
    ( *(network+ii) ).x_position = ii * ( p->network.L / p->network.N );

  }

  return NETSIM_NOERROR;

}

int update_sdi( char * sdi, int up )
{

	/* declare */
	int len;

	/* init */
	len = 0;

	while ( *(sdi + len) != '\0' ) {  printf("%d", len); len++;  }

	assert( (len + up) > 0 );

	sdi[len + up - 1 ] = '-';
	sdi[len + up] = '\0';

	return NETSIM_NOERROR;

}

int save_connectivity( struct parameters * p, struct neuron * network )
{
  
  printf("Saving connectivity...\n");

  /* init */
  int ii, jj;
  char str[2048];
  FILE *s, *r;

  /* bitfield variables */
  uint32_t * send, * receive;
  send = (uint32_t *) malloc( sizeof(uint32_t) * p->network.K );
  receive = (uint32_t *) malloc( sizeof(uint32_t) * p->network.K );

  sprintf( str, "%s/%08dii.bin", p->simulation.output_path, p->simulation.output_file_code );
  
  printf( "Saving sending connections in: %s\n", str );
  s = fopen( str, "wb" );
  sprintf( str, "%s/%08djj.bin", p->simulation.output_path, p->simulation.output_file_code );

  printf( "Saving receiving connections in: %s\n", str );
  r = fopen( str, "wb" );

  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {
    
    // write array of repeating ii of same size as receiving neurons
    for ( jj = 0; jj < p->network.K ; jj++ ) 
    { 
      *( send + jj ) = ii;      
      *( receive + jj ) = *( ( ( *(network+ii) ).K ) + jj );
    }

    fwrite( send, sizeof(uint32_t), p->network.K, s );
    fwrite( receive, sizeof(uint32_t), p->network.K, r );  

  }

  free( send ); free( receive );
  fclose( s );
  fclose( r );

  return NETSIM_NOERROR;

}

/* save incoming connections */
int save_incoming_connections( struct parameters * p, struct neuron * network )
{
  
  /* bitfield variables */
  uint32_t target = 0;

  /* init */
  int * count; int ii, jj; 
  FILE * c; char str[2048];
  count = (int *) malloc( sizeof(int) * p->network.N );

  /* set up output file */
  sprintf( str, "%s/%08dconn.bin", p->simulation.output_path, p->simulation.output_file_code );
  c = fopen( str, "wb" );

  /* count incoming connections */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {
    
    for ( jj = 0; jj < p->network.K ; jj++ ) 
    {
      target = *( ( ( *(network+ii) ).K ) + jj );
      ( *( count + target ) )++;
    }

  }

  /* write output file */
  fwrite( count, sizeof(int), p->network.N, c );
  free( count );
  fclose( c );

  /* return */
  return NETSIM_NOERROR;

}


/* main function */
int main( int argc, char ** argv )
{

  /* print startup dialog */
  assert( startup_dialog() >= NETSIM_NOERROR );

  /* init */
  int ii;
  char * file_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
  clock_t begin, end;
  double buildtime, runtime, totaltime;
  struct parameters p = {{0}};

  /* allocate space for output path before parsing args */
  p.simulation.output_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
    
  /* parse input file path */
  p.simulation.parameter_file_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
  parse_inputs( argc, argv, &p );

  /* parse input file */
  p.simulation.optional_connection_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
  assert( parse_input_file( &p ) >= NETSIM_NOERROR );

  /* set defaults (not network parameters!) */
  if ( p.simulation.report_minutes == 0 ) {  p.simulation.report_minutes = 5;  }
  if ( p.neuron.vr != 0 && p.neuron.vreset == 0 ) {  p.neuron.vreset = p.neuron.vr;  }
  if ( p.network.poisson_rate>0 && p.network.poisson_stop_time == 0 ) {  p.network.poisson_stop_time = p.simulation.T;  }
  if ( p.network.poisson_rate == 0 ) {  p.network.poisson_stop_time = 0;  }
  if ( p.network.poisEIratio < .001 ) {  p.network.poisEIratio = .8;  }
  if ( p.network.spatial_dimensions == 0 ) {  p.network.spatial_dimensions = 1;  } /* 1D by default */
  if ( p.simulation.bin_size == 0 ) {  p.simulation.bin_size = 100;  }
  if ( p.simulation.cutoff_frequency == 0 ) {  p.simulation.cutoff_frequency = 160;  }
  if ( p.simulation.job_id == 0 ) {  p.simulation.job_id = INT_MIN;  }
  assert( (p.synapse.min_uniform_delay == 0) && (p.synapse.max_uniform_delay == 0) ); /* wrong simulator ^& */

  /* init network */
  struct neuron * network;
  network = malloc( sizeof(struct neuron) * p.network.N );

  /* init RNG */
  unsigned int seed = p.simulation.output_file_code;
  struct rng_state * rng;
  rng = malloc( sizeof(struct rng_state) );
  rng_init( rng, seed );

  /* init network */
  struct output o;
  o.spikecount = 0;
  assert( initialize_simulation( &p, network, &o, rng ) >= NETSIM_NOERROR );

  /* make connections */
  printf( "     Connecting network...\n" ); fflush( stdout ); begin = clock();

  if ( p.simulation.connector == NETSIM_RANDOM_CONNECT )  {  random_connect( &(p.network), network, rng );  }
  else if ( p.simulation.connector == NETSIM_GAUSSIAN_CONNECT )  {  gaussian_connect( &(p.network), &(p.simulation), &(p.synapse), network, rng );  }
  else if ( p.simulation.connector == NETSIM_TEST_CONNECT ) {  test_connect( &p, network );  }
  else if ( p.simulation.connector == NETSIM_GAUSSIAN_2D_CONNECT ) 
  {  
  	initialize_2D_simulation( &p, network ); 
    gaussian_connect_2D( &(p.network), &(p.simulation), &(p.synapse), network, rng );
    p.network.spatial_dimensions = 2;
  }  
  else {  printf("Given connector not specified\n"); return NETSIM_BAD_PARAMETER;  }

  /* random rewiring if specified */
  if ( p.network.rewiring_probability > 0 ) 
  {   
    printf( "     Rewiring network with probability %2.2f...\n", p.network.rewiring_probability ); fflush( stdout ); 
    random_rewiring( &p, network, rng ); 
  }

  /* save connectivity */
  if ( p.simulation.save_connectivity == 1 ) 
  	{  save_connectivity( &p, network ); save_incoming_connections( &p, network );  }

  end = clock(); buildtime = (double) (end - begin) / CLOCKS_PER_SEC;

  /* search for the right conductances to get desired rate */
  begin = clock();

  /* initialize input buffers */
  allocate_input_buffers( &p, network );

  /* is it a self-sustained simulation */
  if ( p.simulation.vm_mean != 0 || p.simulation.initiator == NETSIM_INIT_SUSTAINED )
  {  assert( initialize_sustained_simulation( &p, network, &o, rng ) <= NETSIM_NOERROR );  }

  /* start clock */
  printf( "Running network with ge %.10f nS and gi %.10f nS\n", p.synapse.ge*1e9, p.synapse.gi*1e9 );
  fflush( stdout );
  begin = clock();

  /* run simulation */
  assert( run_simulation( &p, network, &o, rng ) <= NETSIM_NOERROR );

  /* calculate simulation time */
  end = clock(); runtime = (double)( end - begin ) / CLOCKS_PER_SEC;
  totaltime = buildtime + runtime;

  /* print summary values to screen */
  printf( "Endtime: %g\nRate (Hz): %g\nBuildtime (s): %g, Runtime (s): %g, Totaltime (s): %g\n",
    o.endtime, ((double)o.spikecount)/p.network.N/o.endtime, buildtime, runtime, totaltime );
  printf( "Total number of spikes: %d\n", o.spikecount );
  printf( "Mean membrane potential: %f (mV)\n", (o.vm_mean_all)*1000 );

  /* clean up */
  fclose( o.spikefile );
  fclose( o.vfile );
  fclose( o.gefile );
  fclose( o.gifile );
  fclose( o.individual_traces_vms );
  fclose( o.individual_traces_ge );
  fclose( o.individual_traces_gi );

  for ( ii = 0 ; ii < p.network.N ; ii++ ) {  free( (*(network+ii)).K );  }
  free( network ); free( rng );
  free( file_path );

  free( p.simulation.optional_connection_path );
  free( p.simulation.output_path );
  free( p.simulation.parameter_file_path );

  return NETSIM_NOERROR;

}
