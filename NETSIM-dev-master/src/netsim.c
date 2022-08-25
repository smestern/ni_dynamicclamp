/********************************************************************************
 *                                                                              *
 *  FILE:   netsim.c                                                            *
 *  VERSION:  0.1                                                               *
 *  PROGRAM:  NETSIM                                                            *
 *                                                                              *
 *  PURPOSE:  A fast, large-scale simulator for topographic spiking networks    *
 *                                                                              *
 *  Copyright (C) 2016-2020 Lyle Muller                                         *
 *  http://mullerlab.ca                                                         *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 *  DEVELOPMENT: Lyle Muller, Charlee Fletterman, Theo Desbordes, Gabriel       *
 *  Benigno, Christopher Steward                                                *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 * This file is part of NETSIM.                                                 *
 *                                                                              *
 *     NETSIM is free software: you can redistribute it and/or modify           *
 *     it under the terms of the GNU General Public License as published by     *
 *     the Free Software Foundation, either version 3 of the License, or        *
 *     (at your option) any later version.                                      *
 *                                                                              *
 *     NETSIM is distributed in the hope that it will be useful,                *
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *     GNU General Public License for more details.                             *
 *                                                                              *
 *     You should have received a copy of the GNU General Public License        *
 *     along with NETSIM.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                              *
 ********************************************************************************/


/* global includes */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* local includes */
#include "isaac64.h"
#include "gamma_dist_rng.h"
#include "rng.h"
#include "interface_c.h"
#include "netsim.h"

/* function: init to zero */
int initialize_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng )
{

  /* init */
  register int ii;
  int Ne = (int) round( .8 * p->network.N ) ;
  int buffer_length;
  double dx_exc, dx_inh;
  double tmp_conduction_speed;
  char str[250];

  /* calculate */
  dx_exc = p->network.L / Ne;
  dx_inh = p->network.L / ( .2 * p->network.N );

  /* calculate needed buffer size */
  if ( p->synapse.max_uniform_delay > 0 ) {  buffer_length = (int) ( p->synapse.max_uniform_delay / p->simulation.dt ) + 1;  }
  else {  buffer_length = (int) ( ( p->synapse.synapse_delay + ((.5 * p->network.L) / (p->synapse.min_conduction_speed)) ) / p->simulation.dt ) + 1;  }
  p->simulation.buffer_length = buffer_length;

  /* initialize network */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    ( *(network+ii) ).excitatory    = ( ii < Ne );
    ( *(network+ii) ).v             = p->neuron.El;
    ( *(network+ii) ).ge            = 0;
    ( *(network+ii) ).gi            = 0;
    ( *(network+ii) ).K             = (uint32_t *) malloc( sizeof(uint32_t) * p->network.K );
    ( *(network+ii) ).lastSpike     = -1;
    ( *(network+ii) ).lastPoisE     = -1;
    ( *(network+ii) ).lastPoisI     = -1;

    /* calculate neuron's x_position in the ring */
    if ( ( *(network+ii) ).excitatory == 1 ) {  (*(network+ii) ).x_position = ii * dx_exc;  }
    else {  ( *(network+ii) ).x_position = ( ii - Ne ) * dx_inh;  }
    ( *(network+ii) ).y_position = 0;

    /* calculate conduction velocity: code 0 -> constant speed, 1 -> uniform, 2 -> gamma distribution  */
    if ( p->simulation.conduction_speed_code == 0 ) {  ( *(network+ii) ).conduction_speed = p->synapse.min_conduction_speed;  }
    else if ( p->simulation.conduction_speed_code == 1 ) {  ( *(network+ii) ).conduction_speed =
      ( p->synapse.max_conduction_speed - p->synapse.min_conduction_speed ) * rng_dbl64(rng) + p->synapse.min_conduction_speed;  }
    else if ( p->simulation.conduction_speed_code == 2 ) {  

      do { tmp_conduction_speed = gamma_sampling( p->simulation.gamma_parameter_shape, p->simulation.gamma_parameter_scale, rng ); }
      while ( (tmp_conduction_speed < p->synapse.min_conduction_speed) || (tmp_conduction_speed > p->synapse.max_conduction_speed) );

      ( *(network+ii) ).conduction_speed = tmp_conduction_speed;

    }

    assert( ( *(network+ii) ).K!=NULL );

  }
  
  /* initialize output files */
  if ( p->simulation.stop_record_time > 0 )
  {

    sprintf( str,  "%s/%08dvms.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->vfile = fopen( str, "w" );
    sprintf( str, "%s/%08dge.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->gefile = fopen( str, "w" );
    sprintf( str, "%s/%08dgi.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->gifile = fopen( str, "w" );

  }
  /* if stop_record_time == 0 nothing will be recorded so no point in making the file */
  else
  {

    o->vfile = fopen( "/dev/null", "w" );
    o->gefile = fopen( "/dev/null", "w" );
    o->gifile = fopen( "/dev/null", "w" );

  }

  sprintf( str,  "%s/%08dspk.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->spikefile = fopen( str, "wb" );
  sprintf( str,  "%s/%08dindividualvms.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->individual_traces_vms = fopen( str, "wb" );
  sprintf( str,  "%s/%08dindividualge.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->individual_traces_ge = fopen( str, "wb" );
  sprintf( str,  "%s/%08dindividualgi.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->individual_traces_gi = fopen( str, "wb" );

  assert( o->vfile!=NULL && o->spikefile!=NULL && o->individual_traces_vms!=NULL && o->individual_traces_ge!=NULL && o->individual_traces_gi!=NULL );

  return 0;

}

/* function: allocate input buffers */
int allocate_input_buffers( struct parameters * p, struct neuron * network )
{

  /* init */
  int ii, buffer_length; 
  buffer_length = p->simulation.buffer_length;

  /* initialize buffers in network */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    ( *(network+ii) ).buff_input_E  = (double *) calloc( buffer_length, sizeof(double) );
    ( *(network+ii) ).buff_input_I  = (double *) calloc( buffer_length, sizeof(double) );

    /* check memory allocation */
    assert( ( *(network+ii) ).buff_input_E!=NULL && ( *(network+ii) ).buff_input_I!=NULL );

  }

  return NETSIM_NOERROR;

}

/* function: intialize for self-sustained activity */
int initialize_sustained_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng )
{

  printf( "> Initializing self-sustained activity " );

  if ( p->simulation.initiator == 0 )
  {

    printf( "with gaussian vm, ge, and gi\n" );

    /* init */

    register int ii;

    for ( ii = 0 ; ii < p->network.N ; ii++ )
    {
      /* init v, ge and gi with normal distribution for sustained (emperical) */
      ( *(network+ii) ).v             = rng_gauss( rng ) * p->simulation.vm_sigma + p->simulation.vm_mean;
      ( *(network+ii) ).ge            = rng_gauss( rng ) * p->simulation.ge_sigma + p->simulation.ge_mean;
      ( *(network+ii) ).gi            = rng_gauss( rng ) * p->simulation.gi_sigma + p->simulation.gi_mean;
    }

  }
  else if ( p->simulation.initiator == 1 )
  {

    printf( "with poisson input\n" );

    /* init */
    int ii, jj;
    double hold_poisEIratio, hold_poisson_rate;
    double hold_poisson_start_time, hold_poisson_stop_time, hold_T;
    FILE * hold_spikefile;
    FILE * hold_vfile;
    FILE * hold_gefile;
    FILE * hold_gifile;
    FILE * hold_individual_traces_vms;
    FILE * hold_individual_traces_ge;
    FILE * hold_individual_traces_gi;
    int hold_spikecount;
    int buffer_length, buffer_position;
    char str[2048];

    /* calculate needed buffer size */
    buffer_length = p->simulation.buffer_length;
    double * temparrE = calloc( buffer_length, sizeof(double) );
    double * temparrI = calloc( buffer_length, sizeof(double) );

    /* save for later to conserve original settings */
    hold_poisson_rate = p->network.poisson_rate;
    hold_poisEIratio = p->network.poisEIratio;
    hold_poisson_start_time = p->network.poisson_start_time;
    hold_poisson_stop_time = p->network.poisson_stop_time;
    hold_spikefile = o->spikefile;
    hold_vfile = o->vfile;
    hold_gefile = o->gefile;
    hold_gifile = o->gifile;
    hold_individual_traces_vms = o->individual_traces_vms;
    hold_individual_traces_ge = o->individual_traces_ge;
    hold_individual_traces_gi = o->individual_traces_gi;
    hold_spikecount = o->spikecount;
    hold_T = p->simulation.T;

    /* set */
    p->network.poisson_rate = 2000 / p->network.K;
    p->network.poisEIratio = 1;
    p->network.poisson_start_time = 0;
    p->network.poisson_stop_time = 200e-3;
    sprintf( str,  "%s/%08dspk_init.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->spikefile = fopen( str, "wb" ); /* reset this file every time initializer runs */
    p->simulation.T = 300e-3;
    o->vfile = fopen( "/dev/null", "rb" );
    o->gefile = fopen( "/dev/null", "rb" );
    o->gifile = fopen( "/dev/null", "rb" );
    o->individual_traces_vms = fopen( "/dev/null", "rb" );
    o->individual_traces_ge = fopen( "/dev/null", "rb" );
    o->individual_traces_gi = fopen( "/dev/null", "rb" );

    /* consistency check */
    if ( p->synapse.p_release > 0 && p->synapse.p_release < 1)
    {  printf( "Currently no external_input with p_release implemented. Exiting..\n" ); return -1;  }

    fclose( o->spikefile );

    buffer_position = (int) round( o->endtime / p->simulation.dt ) % buffer_length;

    for ( ii=0; ii<p->network.N; ii++ )
    {

      (*(network+ii)).lastSpike = (*(network+ii)).lastSpike - p->simulation.T;
      (*(network+ii)).lastPoisE = (*(network+ii)).lastPoisE - (int) round( p->simulation.T/p->simulation.dt );
      (*(network+ii)).lastPoisI = (*(network+ii)).lastPoisI - (int) round( p->simulation.T/p->simulation.dt );

      for ( jj=0; jj<buffer_length; jj++ )
      {
        *(temparrE+jj) = *( (*(network+ii)).buff_input_E + ( buffer_position + jj ) % buffer_length );
        *(temparrI+jj) = *( (*(network+ii)).buff_input_I + ( buffer_position + jj ) % buffer_length );
      }

      for ( jj=0; jj<buffer_length; jj++ )
      {
        *( (*(network+ii)).buff_input_E + jj) = *(temparrE+jj);
        *( (*(network+ii)).buff_input_I + jj) = *(temparrI+jj);
      }

    }

    p->network.poisson_rate = hold_poisson_rate;
    p->network.poisEIratio = hold_poisEIratio;
    p->network.poisson_start_time = hold_poisson_start_time;
    p->network.poisson_stop_time = hold_poisson_stop_time;
    p->simulation.T = hold_T;
    o->spikefile = hold_spikefile;
    o->spikecount = hold_spikecount;
    o->vfile = hold_vfile;
    o->gefile = hold_gefile;
    o->gifile = hold_gifile;
    o->individual_traces_vms = hold_individual_traces_vms;
    o->individual_traces_ge = hold_individual_traces_ge;
    o->individual_traces_gi = hold_individual_traces_gi;

  }

  return 0;

}

/* function initialize 2D simulation */
int initialize_2D_simulation( struct parameters * p, struct neuron * network )
{

  /* init */
  int Ne = (int) round( .8 * p->network.N ); 
  int Ni = p->network.N - Ne;
  register int ii, jj;
  int NrowE, NrowI, index; 

  double dx_exc, dx_inh;

  /* calculate */
  NrowE = floor( sqrt(Ne) ); NrowI = floor( sqrt(Ni) );
  dx_exc = p->network.L / NrowE; dx_inh = p->network.L / NrowI;

  /* loop over neurons and calculate lattice indices */
  for ( ii = 0 ; ii < NrowE ; ii++ )
  {
    for ( jj = 0 ; jj < NrowE ; jj++ )
    {
      /* calculate index */
      index = jj + ( ii * NrowE );

      ( *(network+index) ).x_position = jj * dx_exc; 
      ( *(network+index) ).y_position = ii * dx_exc;
    }
  }

  /* loop over neurons and calculate lattice indices */
  for ( ii = 0 ; ii < NrowI ; ii++ )
  {
    for ( jj = 0 ; jj < NrowI ; jj++ )
    {
      /* calculate index */
      index = jj + ( ii * NrowI ) + Ne;

      ( *(network+index) ).x_position = jj * dx_inh; 
      ( *(network+index) ).y_position = ii * dx_inh;
    }
  }

  return NETSIM_NOERROR;

}

/* function: run_simulation */
int run_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng )
{

  /* make local structure copies */
  struct network_parameters net = p->network;
  struct neuron_parameters neu = p->neuron;
  struct synapse_parameters syn = p->synapse;
  struct simulation_parameters sim = p->simulation;

  /* init local variables */
  register int ii, jj, kk;
  int tt, N, K, Nsteps;
  int Ne;
  int buffer_length, buffer_position, insert_position;
  int number_of_seconds_between_progress_report;
  int spikecount_first_timecheck;
  double t, dt, time_first_timecheck, rate;
  double taue, taui, taur, taum, Ee, Ei, vreset, vth, ge, gi, Gl, El, Ie, Cm;
  clock_t begin, intermediate;
  uint32_t target = 0;
  float64 i_in = 0;

  /* init special for recording ^& */
  double * record_1_timestep_vms, * record_1_timestep_binned_vms;
  double * record_1_timestep_ge, * record_1_timestep_binned_ge;
  double * record_1_timestep_gi, * record_1_timestep_binned_gi;
  double totalsum_v, sum_v, sum_ge, sum_gi;
  char * spike_array;
  char * start_spike_array;
  int spike_index, vsteps, nvsteps;
  uint32_t origin;
  int number_of_recorded_individuals, record_individuals_jump;
  int bin_size, record_downsample_factor, record_start_bin, record_stop_bin;
  int spikes_in_this_round;
  unsigned long long totalcount_v;
  int NrowE;

  /* special init for p relase */
  //double p_release;
  //unsigned int spike_transmitted, spike_not_transmitted;
  //spike_transmitted = spike_not_transmitted = 0;

  /* special init for recalculating delay */
  int delay_in_bins;
  double distance_x, distance_y, distance, L, synapse_delay;

  /* make local variable copies */
  N = net.N; K = net.K; L = net.L;
  Ne = (int) round( .8 * N ) ;
  taur = neu.taur; taum = neu.taum; vreset = neu.vreset; vth = neu.vth; El = neu.El; Ie = neu.Ie;
  taue = syn.taue; taui = syn.taui; Ee = syn.Ee; Ei = syn.Ei; ge = syn.ge; gi = syn.gi; synapse_delay = syn.synapse_delay; //p_release = syn.p_release;
  dt = sim.dt; bin_size = sim.bin_size; record_downsample_factor = sim.record_downsample_factor;

  #ifdef EXTERNAL_INPUT
  /* start POISSON */
  int N_poisson_E, N_poisson_I, delay_in_bins_poisson;
  int stop_poisson_bin, start_poisson_bin;
  int count_poisE, count_poisI;
  double poisson_rate, poisson_ISI;
  poisson_rate = net.poisson_rate;
  start_poisson_bin = (int) round( net.poisson_start_time / sim.dt );
  stop_poisson_bin = (int) round( net.poisson_stop_time / sim.dt );
  N_poisson_E = (int) round( net.poisEIratio * net.K );
  N_poisson_I = net.K - N_poisson_E;
  count_poisE = 0; count_poisI = 0;
  /* end POISSON */
  #endif

  /* init recording */
  number_of_recorded_individuals = 10;
  record_individuals_jump = (int) floor( Ne / number_of_recorded_individuals );

  if ( p->network.spatial_dimensions == 1 ) {  vsteps = (int) floor( N / bin_size ); }
  if ( p->network.spatial_dimensions == 2 ) 
  	{  vsteps = (int) floor( floor(0.8*N) / (bin_size*bin_size) ); nvsteps = (int) floor( sqrt(vsteps) );  }
  NrowE = floor( sqrt(floor(0.8*N)) );

  record_start_bin = (int) floor( sim.start_record_time / dt );
  record_stop_bin = (int) floor( sim.stop_record_time / dt );
  record_1_timestep_vms = (double *) malloc ( sizeof(double) * N );
  record_1_timestep_binned_vms = (double *) malloc ( sizeof(double) * (int) vsteps );
  record_1_timestep_ge = (double *) malloc ( sizeof(double) * N );
  record_1_timestep_binned_ge = (double *) malloc ( sizeof(double) * (int) vsteps );
  record_1_timestep_gi = (double *) malloc ( sizeof(double) * N );
  record_1_timestep_binned_gi = (double *) malloc ( sizeof(double) * (int) vsteps );
  spike_array = (char *)malloc( ( sizeof(int) + sizeof(double) ) * N );

  start_spike_array = spike_array;
  fwrite( &number_of_recorded_individuals, sizeof(int), 1, o->individual_traces_vms );
  fwrite( &number_of_recorded_individuals, sizeof(int), 1, o->individual_traces_ge );
  fwrite( &number_of_recorded_individuals, sizeof(int), 1, o->individual_traces_gi );

  /* calculate parameters */
  Cm = neu.Cm; Gl = Cm / taum;
  Nsteps = round( sim.T / dt );
  number_of_seconds_between_progress_report = sim.report_minutes * 60;
  buffer_length = sim.buffer_length;

  /* initialize variables and arrays */
  begin = clock(); intermediate = clock();
  spikes_in_this_round = 0; totalsum_v = 0; totalcount_v = 0;
  spikecount_first_timecheck = 0; time_first_timecheck = 0;
  init_ni(dt, 1,1);
  /* MAIN SIMULATION LOOP */
  for ( tt = 0 ; tt < Nsteps ; tt++ )
  {

    t = dt * tt;
    buffer_position = tt % buffer_length;

    spike_array = start_spike_array; /* ^& */
    spikes_in_this_round = 0; /* ^& */

    if ( ( time_first_timecheck == 0 ) && ( t > 5e-3 ) )
    {
      spikecount_first_timecheck = o->spikecount; time_first_timecheck = t;
    }

    /* near to no slow-down (between 4e-5 and 5e-4 sec per timestep) */
    if ( ( ((int) ( clock() - intermediate ) / CLOCKS_PER_SEC ) > number_of_seconds_between_progress_report ) && ( t > 10e-3 ) )
    {

      rate = (double) ( o->spikecount - spikecount_first_timecheck ) / N / ( t - time_first_timecheck );
      printf("%f out of %f seconds simulated in %f (number of spikes so far: %d -> %f (Hz) )\n",
        t, sim.T, (double)( clock() - begin)/ CLOCKS_PER_SEC, o->spikecount, rate );
      fflush( stdout );
      intermediate = clock();

      /* if the rate is too high, stop the simulation */
      if ( rate > sim.cutoff_frequency )
      {
          printf( "Early exit because firing rate > %d Hz \n", sim.cutoff_frequency);
          o->endtime = t;
          o->vm_mean_all = totalsum_v / totalcount_v;

          return NETSIM_EARLY_EXIT;
      }

    }

    for ( ii = 0 ; ii < N ; ii++ )
    {

      #ifdef EXTERNAL_INPUT
      /* start POISSON */
      if ( ( (tt >= start_poisson_bin) && (tt < stop_poisson_bin) ) )
      {

        if ( ( N_poisson_E==0) || (*(network+ii)).lastPoisE > tt ) {  ;  }
        else
        {

          while ( (*(network+ii)).lastPoisE <= tt )
          {

            if ( (*(network+ii)).lastPoisE == tt ) /* need to do this for edge case of start, not too happy with it */
            {
              *(( *(network + ii)).buff_input_E + buffer_position ) += ge;
              count_poisE++;
            }

            poisson_ISI = ( -log( rng_dbl64( rng ) ) / ( poisson_rate * N_poisson_E ) );
            delay_in_bins_poisson = (int) round( poisson_ISI / dt );
            ( *(network+ii) ).lastPoisE = tt + delay_in_bins_poisson;

          }

        }

        if ( ( N_poisson_I==0) || (*(network+ii)).lastPoisI > tt ) {  ;  }
        else
        {

          while ( ( (*(network+ii)).lastPoisI <= tt  ) )
          {

            if ( (*(network+ii)).lastPoisI == tt ) /* need to do this for edge case of start, not too happy with it */
            {
              *(( *(network + ii)).buff_input_I + buffer_position ) += gi;
              count_poisI++;
            }

            poisson_ISI =  ( -log( rng_dbl64( rng ) ) / ( poisson_rate * N_poisson_I ) );
            delay_in_bins_poisson = (int) round( poisson_ISI / dt );
            ( *(network+ii) ).lastPoisI = tt + delay_in_bins_poisson;

          }

        }

      }
      /* end POISSON */
      #endif

      /* integrate synaptic input */

      ( *(network+ii) ). ge +=  -dt*( ( *(network+ii) ).ge )/taue +
                                *(( *(network+ii) ).buff_input_E + buffer_position);
      ( *(network+ii) ). gi +=  -dt*( ( *(network+ii) ).gi )/taui +
                                *(( *(network+ii) ).buff_input_I + buffer_position);


      /* reset buffer */
      *( ( *(network+ii) ).buff_input_E + buffer_position ) = 0;
      *( ( *(network+ii) ).buff_input_I + buffer_position ) = 0;

      if ( t > ( (*(network+ii)).lastSpike + taur ) )
      {
        if (ii==0)
        {
          i_in= (( (*(network+ii)).ge*(Ee-(*(network+ii)).v) +  /* excitatory synapse */
                      (*(network+ii)).gi*(Ei-(*(network+ii)).v) +             /* leak conductance */
                      Ie ) ); 
          //printf i_in for debugging
          //printf("%f\n", i_in*1e9);
          ( *(network+ii) ).v = step_clamp(t, i_in);
          printf("%f\n", (*(network+ii)).v);
        }
        else {
        ( *(network+ii) ).v +=
                dt * (
                  ( (*(network+ii)).ge*(Ee-(*(network+ii)).v) +  /* excitatory synapse */
                    (*(network+ii)).gi*(Ei-(*(network+ii)).v) +  /* inhibitory synapse */
                    Gl * (El - (*(network+ii)).v ) +             /* leak conductance */
                    Ie ) / Cm ) ;                                /* current injection */
        
        }
        /* save membrane potential into temporary array ^& */
        *( record_1_timestep_vms + ii ) = (*(network+ii)).v;
        *( record_1_timestep_ge + ii ) = (*(network+ii)).ge;
        *( record_1_timestep_gi + ii ) = (*(network+ii)).gi;

        if ( ( *(network+ii) ).v >= vth )
        {

          /* increment spike counter ^& */
          o->spikecount++;

          /* complicated way to save the spike index and spike time in a char array ^& */
          spike_index = ii;
          memcpy( (void *)spike_array, (void *)(&spike_index), sizeof(int) ); spike_array += sizeof(int);
          memcpy( (void *)spike_array, (void *)(&t), sizeof(double) ); spike_array += sizeof(double);
          spikes_in_this_round++;
          /* end ^& */

          ( *(network+ii) ).v = vreset;
          ( *(network+ii) ).lastSpike = t;

          if ( ( *(network+ii) ).excitatory )
          {

            for ( jj = 0 ; jj < K ; jj++ )
            {

              #include "templates/release_probability_template_start.h"
              target = *( ( ( *(network+ii) ).K ) + jj );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_E + insert_position ) += ge;
              #include "templates/release_probability_template_stop.h"

            }            

          }
          else
          {

            for ( jj = 0 ; jj < K ; jj++ )
            {

              #include "templates/release_probability_template_start.h"
              target = *( ( ( *(network+ii) ).K ) + jj );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_I + insert_position ) += gi;
              #include "templates/release_probability_template_stop.h"

            }

          }

        }

      }
      else
      {

        /* special for recording ^& */
        *( record_1_timestep_vms + ii ) = (*(network+ii)).v;
        *( record_1_timestep_ge + ii ) = (*(network+ii)).ge;
        *( record_1_timestep_gi + ii ) = (*(network+ii)).gi;

      }
      

    }

    /* save spikes outside of record times too ^& */
    spike_array = start_spike_array;
    fwrite( spike_array, sizeof(int) + sizeof(double), spikes_in_this_round, o->spikefile );

    /* write values to file */
    /* special for recording ^& */
    if (  (tt >= record_start_bin) && (tt < record_stop_bin) &&
          (( (record_downsample_factor) == 0) || ((tt % record_downsample_factor) == 0 ))  )
    {

      /* write binned state variables to file */
  		if ( p->network.spatial_dimensions == 1 )
      {

    		for ( ii = 0 ; ii < vsteps ; ii++ )
    		{

     			sum_v = sum_ge = sum_gi = 0;
      		for ( jj = 0 ; jj < bin_size ; jj++ )
      		{
        		sum_v += *(record_1_timestep_vms+(ii*bin_size)+jj);
        		sum_ge += *(record_1_timestep_ge+(ii*bin_size)+jj);
        		sum_gi += *(record_1_timestep_gi+(ii*bin_size)+jj);
      		}

        	*(record_1_timestep_binned_vms+ii) = sum_v / bin_size;
    			*(record_1_timestep_binned_ge+ii) = sum_ge / bin_size;
    			*(record_1_timestep_binned_gi+ii) = sum_gi / bin_size;

    			totalsum_v += sum_v;
    			totalcount_v += bin_size;

        }

    	}
    	else if ( p->network.spatial_dimensions == 2 )
    	{

    		for ( ii = 0 ; ii < vsteps ; ii++ )
    		{

     			sum_v = sum_ge = sum_gi = 0;
      		origin = floor( ii / nvsteps ) * ( bin_size * NrowE ) + ( ii % nvsteps ) * ( bin_size );            
      		
      		for ( jj = 0 ; jj < bin_size ; jj++ )
          {
     	  		for ( kk = 0 ; kk < bin_size ; kk++ )
     	  		{
	        		sum_v += *(record_1_timestep_vms+origin+(jj*NrowE)+kk);
	        		sum_ge += *(record_1_timestep_ge+origin+(jj*NrowE)+kk);
	        		sum_gi += *(record_1_timestep_gi+origin+(jj*NrowE)+kk);
    	  		}
      		}

      		*(record_1_timestep_binned_vms+ii) = sum_v / (bin_size*bin_size);
          *(record_1_timestep_binned_ge+ii) = sum_ge / (bin_size*bin_size);
          *(record_1_timestep_binned_gi+ii) = sum_gi / (bin_size*bin_size);

          totalsum_v += sum_v;
          totalcount_v += (bin_size*bin_size);

      	}

    	}

      fwrite( record_1_timestep_binned_vms, sizeof(double), vsteps, o->vfile );
      fwrite( record_1_timestep_binned_ge, sizeof(double), vsteps, o->gefile );
      fwrite( record_1_timestep_binned_gi, sizeof(double), vsteps, o->gifile );

    } 

    /* keep track of a few unbinned neurons to verify that everything looks normal */
    for ( jj = 0 ; jj < Ne ; jj+=record_individuals_jump )
    {
      fwrite( ( record_1_timestep_vms+jj ), sizeof(double), 1, o->individual_traces_vms );
      fwrite( ( record_1_timestep_ge+jj ), sizeof(double), 1, o->individual_traces_ge );
      fwrite( ( record_1_timestep_gi+jj ), sizeof(double), 1, o->individual_traces_gi );
    }
    /* end ^& */

  }

  /* write last values to the output struct */
  o->endtime = sim.T;
  o->vm_mean_all = totalsum_v / totalcount_v;
  
  //clean up the NI
  clean_up();
  /* return value */
  return NETSIM_NOERROR;

}


