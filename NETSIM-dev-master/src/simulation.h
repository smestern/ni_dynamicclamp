#ifndef _SIMULATION_H_
#define _SIMULATION_H_

/* defined constants */
#define NETSIM_NOERROR            0
#define NETSIM_ERROR              1
#define NETSIM_BAD_PARAMETER      2
#define NETSIM_EARLY_EXIT         -1

#define NETSIM_RANDOM_CONNECT     1
#define NETSIM_GAUSSIAN_CONNECT   2
#define NETSIM_TEST_CONNECT       3
#define NETSIM_GAUSSIAN_2D_CONNECT   4

#define NETSIM_INIT_STANDARD      0
#define NETSIM_INIT_SUSTAINED     1

#define FILE_BUFFER_SIZE          2048

/* individual parameter structures */
struct neuron_parameters {
  double taum, taur, tauref, vr, vth, vreset, Gl, El, Ie, Cm;
};

struct synapse_parameters {
  double ge, gi, taue, taui, Ee, Ei, synapse_delay, min_conduction_speed, max_conduction_speed, p_release, min_uniform_delay, max_uniform_delay;
};

struct network_parameters {
  int N, K;
  int spatial_dimensions;
  double L, sigma_space, poisson_rate, poisson_start_time, poisson_stop_time, poisEIratio, rewiring_probability;
};

struct simulation_parameters {

  /* simulation control parameters */
  int connector, initiator, output_file_code, bin_size, record_downsample_factor, buffer_length, cutoff_frequency;
  int conduction_speed_code, gamma_parameter_shape, output_path_override;
  double T, dt, vm_mean, vm_sigma, ge_mean, ge_sigma, gi_mean, gi_sigma, gamma_parameter_scale;
  char * parameter_file_path, * optional_connection_path, * output_path;

  /* recording parameters */
  int spikecount, save_connectivity;
  double start_record_time, stop_record_time, report_minutes;
  FILE * vfile, * spikefile, * gefile, * gifile;

  /* scan parameters */
  int job_id;

};

/* umbrella parameter structure */
struct parameters {
  struct neuron_parameters neuron;
  struct synapse_parameters synapse;
  struct network_parameters network;
  struct simulation_parameters simulation;
};

// /* synapse bitfield structure */
// struct trip {
//   unsigned long offset : 64;
// };

/* neuron structure */
struct neuron {
  bool excitatory;
  double v;
  double ge;
  double gi;
  uint32_t * K;
  double conduction_speed;
  double x_position;
  double y_position;
  double lastSpike;
  int lastPoisE;
  int lastPoisI;
  double * buff_input_E;
  double * buff_input_I;
};

/* output structure */
struct output {
  int spikecount, spikecount_init;
  double endtime, vm_mean_all, p_release_out;
  double vmmean, vmsigma, gemean, gesigma, gimean, gisigma; /* end-state of the network */
  int NsideE; /* for 2d lattice (excitatory) */
  int NsideI; /* for 2d lattice (inhibitory) */
  FILE * vfile;
  FILE * spikefile;
  FILE * individual_traces_vms;
  FILE * individual_traces_ge;
  FILE * individual_traces_gi;
  FILE * gefile;
  FILE * gifile;
};

#endif
