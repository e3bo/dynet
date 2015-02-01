/* igraph library headers */

#include <igraph/igraph.h>

/* gsl headers */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

/*headers for hashing */
#include "uthash.h"		/* hash table macros */
#include <stddef.h>		/* offset of */

/* headers for argument parsing */
#include <stdlib.h>
#include <error.h>
#include <argp.h>

/* errno.h has declaration of errno variable */
#include <errno.h>

/** more stuff for argument parsing **/

const char *argp_program_version = "dynet 0.2";
const char *argp_program_bug_address = "<ebo@mail.utexas.edu>";

static char doc[] = "DyNet -- simulator for disease spreading\
 on dynamic networks\
\vThe simulator allows the network to relax from an initial\
 degree distribuion toward a target degree distribution. The\
 condition that degree correlations are close to zero is used\
 to determine the types of edges that are added and deleted.\
 No further documentation is available yet.";

/* A description of the arguments we accept. */
static char args_doc[] = " ";

/* Keys for options without short-options */
#define OPT_ABORT 1

/* macros and functions for network generation*/
#define DYNET_EXP (int) 1
#define DYNET_REG (int) 2
#define DYNET_PARETO (int) 3
#define DYNET_ER (int) 4

/* The options we understand. */
static struct argp_option options[] = {
  {0, 0, 0, 0, "General Options"},
  {"verbose", 'v', 0, 0, "Produce verbose output"},
  {"quiet", 'q', 0, 0, "Don't produce any outpu"},
  {"silent", 's', 0, OPTION_ALIAS},
  {"output", 'o', "FILE", 0,
   "Output to FILE instead of standard output"},
  {0, 0, 0, 0, "Network Model Options"},
  {"network_rate", 'n', "RATE", 0,
   "Network degree distribution approaches equilibirium "
   "at rate RATE (default 0.4)"},
  {"network_tension", 'k', "RATE", 0,
   "Network degree correlations are removed in proportion to "
   "thier size at rate RATE (default 50)"},
  {"network_size", 'N', "COUNT", 0,
   "the network is made to have about COUNT (default 10000) nodes "},
  {"network_type", 'y', "INTEGER_CODE", 0,
   "INTEGER_CODE (default 2) determines the limiting degree distribution:"
    " 1==exponential, 2==regular, 3==Pareto, 4==Poisson"},
  {"par1", 'm', "RATIONAL", 0,
   "If network_type is in {1, 2, 4}, RATIONAL (default 4) is "
    " the mean of the degree distribution. If network_type is 3,"
   " RATIONAL is A in definition of Pareto pdf. in GSL "
   "documentation"},
  {"par2", 'b', "RATIONAL", 0,
   "If network_type is PARETO, RATIONAL is B in the definition of"
   " the Pareto pdf in the GSL documentation"},
  {0, 0, 0, 0, "Disease Model Options"},
  {"trans_rate", 't', "RATE", 0,
   "Disease moves at rate RATE (default 2) across an edge"},
  {"recov_rate", 'r', "RATE", 0,
   "Infected nodes reacover at rate RATE (default 1)"},
  {"interval", 'i', "LENGTH", 0,
   "State variables are printed to output at every LENGTH\
 (default 0.05) time units"},
  {"epsilon", 'e', "FRACTION", 0, "Initial fraction FRACTION "
   "(default 0.1) of edges from suscetible nodes pointing "
   "to infected nodes"},
  {"finish_time", 'f', "LENGTH", 0, "Simulation time will start "
   "at zero and run LENGTH (default 5) time units"},
  {0, 0, 0, 0, "Evolutionary Model Options"},
  {"sample_fraction", 'p', "FRACTION", 0,
    "FRACTION (default 0.001) nodes in virus gene tree are selected"
      "for inclusion in newick output"},
  {0}
};

/* maximum rate allowed in arguments */
#define MAX_RATE 10000.

/* Used by `main' to communicate with `parse_opt'. */
struct arguments
{
  int silent, verbose;		/* `-s', `-v' */
  char *output_file;		/* FILE arg to `--output' */
  double trans_rate;		/* RATE arg to `--trans_rate' */
  double recov_rate;		/* RATE arg to `--recov_rate' */
  double interval;		/* LENGTH arg to `--interval' */
  double epsilon;		/* FRACTION arg to `--epsilon' */
  double finish_time;		/* LENGTH arg to `--finish_time' */
  double network_rate;		/* RATE arg to `--network_rate' */
  double network_tension;	/* RATE arg to `--network_tension' */
  int network_size;		/* COUNT arg to `--network_size' */
  int network_type;             /* INTEGER_CODE arg to `--network_type' */
  double par1;                  /* RATIONAL arg to `--par1' */
  double par2;                  /* RATIONAL arg to `--par2' */
  double sample_fraction;	/* FRACTION arg to `--sample_fraction' */


};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the INPUT argument from `argp_parse', which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = (struct arguments *) state->input;

  switch (key)
    {
    case 'q':
    case 's':
      arguments->silent = 1;
      break;
    case 'v':
      arguments->verbose = 1;
      break;
    case 'o':
      arguments->output_file = arg;
      break;
    case 't':
      arguments->trans_rate = strtod (arg, NULL);
      if (arguments->trans_rate < 0 || arguments->trans_rate > MAX_RATE)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "trans_rate=%g, should be in [0, %g]",
			 arguments->trans_rate, MAX_RATE);

	}
      break;
    case 'r':
      arguments->recov_rate = strtod (arg, NULL);
      if (arguments->recov_rate < 0 || arguments->recov_rate > MAX_RATE)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "recov_rate=%g, should be in [0, %g]",
			 arguments->recov_rate, MAX_RATE);

	}
      break;
    case 'i':
      arguments->interval = strtod (arg, NULL);
      if (arguments->interval <= 0 || arguments->interval > 1e2)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "interval=%g, should be in [0, %g]",
			 arguments->recov_rate, 1e2);

	}
      break;
    case 'e':
      arguments->epsilon = strtod (arg, NULL);
      if (arguments->epsilon < 0 || arguments->epsilon > 1)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "epsilon=%g, should be in [0, 1]",
			 arguments->epsilon);

	}
      break;
    case 'f':
      arguments->finish_time = strtod (arg, NULL);
      if (arguments->finish_time < 0 || arguments->finish_time > 100)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "finish_time=%g, should be in [0, 100]",
			 arguments->finish_time);

	}
      break;
    case 'n':
      arguments->network_rate = strtod (arg, NULL);
      if (arguments->network_rate < 0 || arguments->network_rate > MAX_RATE)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "network_rate=%g, should be in [0, %g]",
			 arguments->network_rate, MAX_RATE);

	}
      break;
    case 'k':
      arguments->network_tension = strtod (arg, NULL);
      if (arguments->network_tension < 0
	  || arguments->network_tension > MAX_RATE)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "network_tension=%g, should be in [0, %g]",
			 arguments->network_tension, MAX_RATE);

	}
      break;
    case 'N':
      arguments->network_size = atoi (arg);
      if (arguments->network_size < 1 || arguments->network_size > 1e6)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "network_size=%d, should be in [1, %d]",
			 arguments->network_size, (int) 1e6);

	}
      break;
    case 'p':
      arguments->sample_fraction = strtod (arg, NULL);
      if (arguments->sample_fraction < 0 || arguments->sample_fraction > 1)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "sample_fraction =%g, should be in [0, 1]",
			 arguments->sample_fraction);

	}
      break;
    case 'y':
      arguments->network_type =  atoi (arg);
      if (arguments->network_type != DYNET_REG 
          && arguments->network_type != DYNET_ER
          && arguments->network_type != DYNET_PARETO
          && arguments->network_type != DYNET_EXP)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "network type =%d, should be in {%d, %d, %d, %d}",
			 arguments->network_type, DYNET_EXP, DYNET_REG,
                         DYNET_PARETO, DYNET_ER);

	}
      break;
    case 'm':
      arguments->par1 = strtod (arg, NULL);
      if (arguments->par1< 0 || arguments->par1 > 100)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "par1 =%g, should be in [0, 100]",
			 arguments->par1);

	}
      break;
    case 'b':
      arguments->par2 = strtod (arg, NULL);
      if (arguments->par2< 1 || arguments->par2 > 100)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "par2 =%g, should be in [1, 100]",
			 arguments->par2);

	}
      break;
    case ARGP_KEY_ARG:
      /* Too many arguments. */
      argp_usage (state);
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/*structure for events in simulation */

#define MUTATE (char) 'm'
#define INFECT (char) 'i'
#define RECOVER (char) 'r'
#define EDGE_ADD (char) 'a'
#define EDGE_DEL (char) 'd'


struct event
{
  /* key is an aggregate of event code, ego_id, and alter_id */
  /* TODO use separate structures for different event_codes */
  char event_code;
  int ego_id;
  /* The type of the last part of the key must be consistent 
   * with the keylen calculation below */
  int alter_id;

  int i_deg;
  int j_deg;
  int N_ij;			/* number of edges between degree i and j nodes */
  double rate;
  int strain_id;
  int phylo_id;

  /* UT_hash_handle member makes structure a hashtable */
  UT_hash_handle hh;
};

/** genealogy struct **/
struct node
{
  int abv;
  int ndes;
  float time;
  int flag;
};


/*hash table */
struct event *event_table = NULL;

/*keylen holds the size of the hash key */
unsigned keylen;

/* sum of rate of all events in hash table */
double rate_sum;

/* array to track state of the each node */
char *node_states;

/* arrays to track states of nodes of each degree class */
int *S_k;
int *I_k;

/* maximum number of stochastic sim steps */
#define STEPMAX 50000

/*maximum degree correlations allowed */
#define MAXREDGE 0.05

/* SIR model parameters */
double trans_rate;
double recov_rate;

/* OUTSIDE is the ID of any host outside of population that transmit 
 * the initial infections */
#define OUTSIDE -99

/* rounding doubles to nearest integer */
/* source http://www.cs.tut.fi/~jkorpela/round.html */
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/** function prototypes **/

int infect (int infector_id, int infectee_id, igraph_t *g, 
            igraph_vector_t *degree);
int recover (int recoverer_id, igraph_t *g, igraph_vector_t *degree);
void delete_all_events ();
void delete_event (struct event *ev);
void delete_event_by_key (char event_code, int ego_id, int alter_id);
void print_rates (struct event *event_table);
int check_event_by_key (char event_code, int ego_id, int alter_id);
int check_infector_events (int infector_id, igraph_t *g);
double get_rate_sum ();
double get_infection_rate_sum ();
int neighbor_state_count (int susceptible_id, int *num_S_nb, int *num_I_nb,
                          igraph_t *g);
int get_num_infected (int max_degree);
int get_approx_case_prod_variance (double trans, double gen_time, int num_infected, 
                               double *icpm, double *icpv, igraph_t *g);
double get_approx_case_prod_mean (double gen_time, int num_infected);
FILE *open_output (const char *filename);
int close_file (FILE *f);
void prtree3 (struct node *ptree, int num_nodes, 
              FILE** f_tree, FILE** f_coal_times, FILE** f_samp_times);
int sequence_to_graph(int *sequence, int num_nodes, igraph_t *g);
int get_degree_seq (int *resulting_degree_seq, int num_nodes, 
                    double  par1, double par2, int net_type);

     /** functions **/

void
delete_all_events ()
{
  struct event *current_event;

  while (event_table)
    {
      current_event = event_table;
      HASH_DEL (event_table, current_event);
      free (current_event);
    }
}

void
print_rates (struct event *event_table)
{
  struct event *ev;
  unsigned z = 0;

  for (ev = event_table; ev != NULL; ev = (struct event *) ev->hh.next)
    {
      printf ("rate is %g\n", ev->rate);
    }
}


void
delete_event (struct event *ev)
{
  HASH_DEL (event_table, ev);
  free (ev);
}

int
check_event_by_key (char event_code, int ego_id, int alter_id)
{

  struct event *ev1, ev2;

  memset (&ev2, 0, sizeof (struct event));
  ev2.event_code = event_code;
  ev2.ego_id = ego_id;
  ev2.alter_id = alter_id;
  HASH_FIND (hh, event_table, &ev2.event_code, keylen, ev1);
  if (!ev1)
    {
      error_at_line (0, errno, __FILE__, __LINE__,
		     "failed to find hash with key %c_%d_%d",
		     event_code, ego_id, alter_id);
    }
  else if (event_code == INFECT && ev1->rate != trans_rate)
    {
      error_at_line (0, 0, __FILE__, __LINE__,
		     "infection event with rate %g, should be %g",
		     ev1->rate, trans_rate);
    }
  else if (event_code == RECOVER && ev1->rate != recov_rate)
    {
      error_at_line (0, 0, __FILE__, __LINE__,
		     "recovery event with rate %g, should be %g",
		     ev1->rate, recov_rate);
    }
  return 0;
}

void
delete_event_by_key (char event_code, int ego_id, int alter_id)
{

  struct event *ev1, ev2;

  memset (&ev2, 0, sizeof (struct event));
  ev2.event_code = event_code;
  ev2.ego_id = ego_id;
  ev2.alter_id = alter_id;
  HASH_FIND (hh, event_table, &ev2.event_code, keylen, ev1);
  if (!ev1)
    {
      fprintf (stderr, "error: %s: %d: failed to find hash key \n",
	       __FILE__, __LINE__);
      exit (1);
    }
  HASH_DEL (event_table, ev1);
  free (ev1);
}

int
main (int argc, char *argv[])
{

  /* index for loops */
  size_t i;

  /* value returned by functions called */
  int ret = 0;

  /* contact network */
  igraph_t g;

  /* number of hosts of degree i*/
  int N_i;

  /* vector with degree of each host  */
  igraph_vector_t degree;

  /* maximum degree of nodes in contact network */
  int max_degree;

  /* initial degree sequence of nodes in contact network */
  int *initial_degree_seq;

  /* element I of degree_dist contains the number of nodes 
   * with degree I */
  gsl_histogram *degree_dist;

  /* counter of number of errors */
  errno = 0;

  /* epsilon is the initial fraction of infected nodes */
  double epsilon;

  /* I_initial holds the number of initially infected hosts in 
   * based on epsilon */
  int I_initial;

  /* array holding IDS of initially infected hosts */
  int *I_initial_IDs;
  int *all_nodes_IDs;

  /* time in simulation */
  double time = 0;

  /* count of steps in simulation */
  int step_count = 0;

  double rand;
  double cum_density, p_I, p_S;

  struct arguments arguments;

/* Default values for variable set by arguments */
  arguments.silent = 0;
  arguments.verbose = 0;
  arguments.output_file = (char *) "-";
  arguments.trans_rate = 2;
  arguments.recov_rate = 1;
  arguments.interval = 0.05;
  arguments.epsilon = 0.1;
  arguments.finish_time = 5;
  arguments.network_rate = 0.4;
  arguments.network_tension = 50;
  arguments.network_size = 10000;
  arguments.network_type = DYNET_REG;
  arguments.par1 = 4;
  arguments.par2 = 1;
  arguments.sample_fraction = 0.001;

  /*TODO add option to use high tension to eliminate degree correlations
   * in initial network before beginning simulation */

  struct event *ev, *ev1, ev2;
  unsigned z;

  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  fprintf (stdout, "trans_rate = %g\nrecov_rate = %g\ninterval = %g\n"
	   "epsilon = %g\nfinish_time = %g\nnetwork_rate = %g\n"
	   "network_tension = %g\nnetwork_size = %d\nnetwork_type = %d\n"
           "par1 = %g\npar2 = %g\n"
           "sample_fraction= %g\nOUTPUT_FILE = %s\nVERBOSE = %s\nSILENT = %s\n",
	   arguments.trans_rate,
	   arguments.recov_rate,
	   arguments.interval,
	   arguments.epsilon,
	   arguments.finish_time,
	   arguments.network_rate,
	   arguments.network_tension,
	   arguments.network_size,
           arguments.network_type,
           arguments.par1,
           arguments.par2,
           arguments.sample_fraction,
	   arguments.output_file,
	   arguments.verbose ? "yes" : "no", arguments.silent ? "yes" : "no");
#ifdef    DYNET_DEBUG
  fprintf (stdout, "dynet debugging level = %d\n", DYNET_DEBUG);
#endif

  /* node_hoste[I] contains the id of the genealogy node in host  with id I having the greatest time */
  int *node_host;		

  struct node *ptree;

  ptree =
    (struct node *) calloc ((unsigned) (2 * arguments.network_size),
			    sizeof (struct node));
  if (!ptree)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  node_host =
    (int *) calloc ((unsigned) (arguments.network_size), sizeof (int));

  if (!node_host)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  int gen_ind = 0;

  trans_rate = arguments.trans_rate;
  recov_rate = arguments.recov_rate;
  epsilon = arguments.epsilon;

  keylen =
    offsetof (struct event, alter_id) + sizeof (int) - offsetof (struct event,
								 event_code);

  initial_degree_seq = (int *) malloc ((arguments.network_size) * sizeof (int));

  if (!initial_degree_seq)
    {
      fprintf (stderr, "Error: %s: %d: Malloc failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  if (arguments.network_type != DYNET_ER)
    {
  get_degree_seq (initial_degree_seq, arguments.network_size,
                 arguments.par1, arguments.par2, arguments.network_type); 
  sequence_to_graph (initial_degree_seq, arguments.network_size, &g);
    }
  else
    {
      igraph_erdos_renyi_game (&g, IGRAPH_ERDOS_RENYI_GNM,
                               arguments.network_size,
                               round (arguments.network_size * arguments.par1 /2),
                               0, 0);
    }



  /* get the degree distribution */
  igraph_vector_init (&degree, 0);
  igraph_degree (&g, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
  max_degree = (int) igraph_vector_max (&degree);
  degree_dist = gsl_histogram_alloc (max_degree + 1);
  gsl_histogram_set_ranges_uniform (degree_dist, -0.5, max_degree + 0.5);
  for (i = 0; i < igraph_vector_size (&degree); i++)
    {
      gsl_histogram_increment (degree_dist, VECTOR (degree)[i]);
    }

  // open filestreams and print headers
  FILE *fp, *of;

  if (arguments.output_file == "-")
    {
      of = stdout;
    }
  else
    {
      of = open_output (arguments.output_file);
    }
  fprintf (of, "        time   ");
  for (i = 0; i <= max_degree; i++)
    {
      N_i = round(gsl_histogram_get (degree_dist, i));
      if (N_i)
        {
          fprintf (of, "N_%05u S_%05u I_%05u   ", i + 1, i + 1, i + 1);
        }
    }
fprintf (of, "   incidence  num_infected      gen_time          icpv          icpm      vm_ratio\n"); 

/*  fprintf (of, "r_edge p_I p_S\n"); */

  const gsl_rng_type *T;
  gsl_rng *rng;

  /* create a generator by the environment variable 
   * GSL_RNG_TYPE */

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);

  double write_point = 0;

  /* put intial infection and recovery event into table */

  node_states = (char *) malloc ((arguments.network_size + 1) * sizeof (char));
  if (!node_states)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  memset (node_states, 's', arguments.network_size * sizeof (char));
  node_states[arguments.network_size * sizeof (char)] = '\0';

  S_k = (int *) malloc ((max_degree + 1) * sizeof (int));
  if (!S_k)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  for (i = 0; i <= max_degree; i++)
    {
      S_k[i] = gsl_histogram_get (degree_dist, i);
    }

  I_k = (int *) calloc ((max_degree + 1), sizeof (int));
  if (!I_k)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  /* TODO select host at random from network */
  infect (OUTSIDE, 0, &g, &degree);
  ptree->time = 0;
  node_host[0] = 0;
  gen_ind++;

  while (time < arguments.finish_time && rate_sum > 0)
    {
      step_count++;

#if DYNET_DEBUG > 0
      /* check that rate sum is being updated correctly */

      if (fabs (get_rate_sum () - rate_sum) > 1e-6)
	{
	  error_at_line (EXIT_FAILURE, 0, __FILE__, __LINE__,
			 "rate_sum - get_rate_sum ()= %g\n",
			 rate_sum - get_rate_sum ());

	}

      /* check that counts of susceptibles and infecteds are 
       * correct when summed over all degree classes */

      int Scount, Icount, Scount2, Icount2;
      Scount = Icount = Scount2 = Icount2 = 0;
      for (i = 0; node_states[i] != 0; i++)
	{
	  if (node_states[i] == 's')
	    Scount++;
	  else if (node_states[i] == 'i')
	    Icount++;
	}
      for (i = 0; i <= max_degree; i++)
	{
	  Scount2 += S_k[i];
	  Icount2 += I_k[i];
	}
      if (Scount != Scount2 || Icount != Icount2)
	{
	  fprintf (stderr, "Debug: %s: %d: counts don't match\n",
		   __FILE__, __LINE__);
	  fprintf (stderr, "\tSc=%d Sc2=%d, Ic=%d, Ic2=%d\n", Scount, Scount2,
		   Icount, Icount2);
	}

      /* Check that counts of Susceptiable and infecteds by degree
       * class are being updated correctly. */

      int deg, *I_k_check, *S_k_check;
      I_k_check = (int *) calloc ((max_degree + 1), sizeof (int));
      if (!I_k_check)
	{
	  error_at_line (EXIT_FAILURE, errno, __FILE__, __LINE__,
			 "I_k_check");
	}
      S_k_check = (int *) calloc ((max_degree + 1), sizeof (int));
      if (!S_k_check)
	{
	  error_at_line (EXIT_FAILURE, errno, __FILE__, __LINE__,
			 "S_k_check");
	}
      for (i = 0; i < arguments.network_size; i++)
	{
	  if (node_states[i] == 's')
	    {
	      deg = VECTOR (degree)[i];
	      S_k_check[deg]++;
	    }
	  if (node_states[i] == 'i')
	    {
	      deg = VECTOR (degree)[i];
	      I_k_check[deg]++;
	    }
	}
      error_message_count = 0;
      for (i = 0; i <= max_degree; i++)
	{
	  if (S_k_check[i] != S_k[i])
	    {
	      error_at_line (0, 0, __FILE__, __LINE__,
			     "S_k[%d] = %d, should be %d",
			     i, S_k[i], S_k_check[i]);
	    }
	  if (I_k_check[i] != I_k[i])
	    {
	      error_at_line (0, 0, __FILE__, __LINE__,
			     "I_k[%d] = %d, should be %d",
			     i, I_k[i], I_k_check[i]);
	    }
	}
      if (error_message_count != 0)
	{
	  error (EXIT_FAILURE, 0, "%u errors found", error_message_count);
	}
      free (I_k_check);
      free (S_k_check);

      /* calculate p.I, i.e. fraction of arcs from susctpibles
       * going to infecteds */
      int num_S_nb, num_I_nb, sum_S_nb, sum_I_nb, sum_nb, susceptible_id;
      sum_nb = sum_S_nb = sum_I_nb = 0;
      for (i = 0; i < arguments.network_size; i++)
	{
	  if (node_states[i] == 's')
	    {
	      num_S_nb = num_I_nb = 0;
	      sum_nb += VECTOR (degree)[i];
	      susceptible_id = i;
	      neighbor_state_count (susceptible_id, &num_S_nb, &num_I_nb, &g);
	      sum_S_nb += num_S_nb;
	      sum_I_nb += num_I_nb;
	    }
	}
      p_S = (double) sum_S_nb / sum_nb;
      p_I = (double) sum_I_nb / sum_nb;
#endif
#if DYNET_DEBUG > 1

      /* Check for all events in the hash table, making sure the 
       * hash table is consistent with the current state of the network. */
      error_message_count = 0;
      int infector_id;
      int event_count = 0;

      /*first check for and count one edge addition or deletion event that 
       * can always happen */

      memset (&ev2, 0, sizeof (struct event));
      ev2.event_code = EDGE_ADD;
      HASH_FIND (hh, event_table, &ev2.event_code, keylen, ev1);
      if (ev1)
	{
	  event_count++;
	}
      else
	{
	  ev2.event_code = EDGE_DEL;
	  HASH_FIND (hh, event_table, &ev2.event_code, keylen, ev1);
	  if (ev1)
	    {
	      event_count++;
	    }
	}

      for (i = 0; i < arguments.network_size; i++)
	{
	  if (node_states[i] == 'i')
	    {
	      infector_id = i;
	      event_count += check_infector_events (infector_id, &g);
	    }
	}
      int hash_count = HASH_COUNT (event_table);
      if (error_message_count != 0)
	{
	  error (EXIT_FAILURE, 0, "%u events are missing in hash table "
		 "or are present but occuring at the wrong rate",
		 error_message_count);
	}
      if (hash_count != event_count)
	{
	  error (EXIT_FAILURE, 0, "%d extra events in hash table\n",
		 hash_count - event_count);
	}

#endif
#if DYNET_DEBUG_TEST == 1
      if (step_count == 10)
	{
	  /* set the rate of the first event in the table to 0 */
	  fprintf (stderr,
		   "Changing the rate of an event to test debugging tests\n");
	  event_table->rate = 0;
	}
#endif
#if DYNET_DEBUG_TEST == 2
      if (step_count == 70)
	{
	  /* add an extra event to the table */
	  fprintf (stderr, "Adding an event to test debugging tests\n");
	  ev1 = (struct event *) malloc (sizeof (struct event));
	  if (!ev1)
	    {
	      fprintf (stderr,
		       "Error: %s: %d: malloc failed\n", __FILE__, __LINE__);
	      return (1);
	    }
	  memset (ev1, 0, sizeof (struct event));
	  ev1->event_code = INFECT;
	  ev1->rate = 0;
	  HASH_ADD (hh, event_table, event_code, keylen, ev1);
	}
#endif
#if DYNET_DEBUG_TEST == 3
      if (step_count == 50)
	{
	  /* remove an event at beginning of the table */
	  fprintf (stderr, "Removing an event to test debugging tests\n");
	  delete_event (event_table);
	}
#endif
#if DYNET_DEBUG_TEST == 4
      if (step_count == 50)
	{
	  /* give two events rates with compensating erros */
	  fprintf (stderr, "Swapping event rates to test debugging tests\n");
	  ev = event_table;
	  ev->rate += 0.001;
	  ev = (struct event *) ev->hh.next;
	  ev->rate -= 0.001;
	}
#endif
#if DYNET_DEBUG_TEST == 5
      if (step_count == 50)
	{
	  /* change event code  */
	  fprintf (stderr,
		   "Changing event code of an event to test debugging tests\n");
	  ev = event_table;
	  ev->event_code = 'Z';
	}
#endif


      if (step_count > STEPMAX)
	{
	  fprintf (stderr, "Warning: %s: %d: reached STEPMAX, breaking\n",
		   __FILE__, __LINE__);
	  break;
	}

      /* determine time of next event */
      time += gsl_ran_exponential (rng, 1 / rate_sum);
      if (time > write_point)
	{
	  fprintf (of, "%9e   ", write_point);
	  write_point += arguments.interval;
	  for (i = 0; i <= max_degree; i++)
	    {
              int N_i = round(gsl_histogram_get (degree_dist, i));
              if (N_i)
                {
                  fprintf (of, "%7ld %7d %7d   ",  N_i, S_k[i], I_k[i]);
                }
	    }
          double incidence = get_infection_rate_sum();
          int num_infected = get_num_infected(max_degree);
          double gen_time = (double) num_infected / incidence;
/*          double icpm = get_approx_case_prod_mean (gen_time, num_infected);*/
          double icpm;
          double icpv;
          get_approx_case_prod_variance (arguments.trans_rate,
                                         gen_time,
                                         num_infected, 
                                         &icpm, &icpv, &g);
          double vm_ratio = icpv / icpm;

          fprintf (of, "%9e  %12d  %9e  %9e  %9e  %9e\n", 
                   incidence, num_infected, gen_time, icpv, icpm, vm_ratio );
/*	  fprintf (of, " %g %g\n", p_I, p_S); */
	}

      /* select one event to be next proportionally to it's rate */

      rand = rate_sum * gsl_rng_uniform (rng);
      cum_density = 0;
      for (ev = event_table; ev != NULL; ev = (struct event *) ev->hh.next)
	{
	  cum_density += ev->rate;
	  if (cum_density > rand)
	    {
	      break;
	    }
	}

      switch (ev->event_code)
	{
	case RECOVER:
	  (ptree + gen_ind)->time = time;
	  (ptree + gen_ind)->abv = node_host[ev->ego_id];
	  node_host[ev->ego_id] = gen_ind;
	  gen_ind++;
	  recover (ev->ego_id, &g, &degree);
	  break;
	case INFECT:
	  (ptree + gen_ind)->time = time;
	  if (gen_ind > 0)
	    {
	      (ptree + gen_ind)->abv = node_host[ev->ego_id];
	    }
	  node_host[ev->ego_id] = gen_ind;
	  node_host[ev->alter_id] = gen_ind;
	  gen_ind++;
	  infect (ev->ego_id, ev->alter_id, &g, &degree);
	  break;
	default:
	  fprintf (stderr, "Error: %s: %d: Illegal event code\n",
		   __FILE__, __LINE__);
	  ret = 1;
	}
      if (ret == 1)
	{
	  break;
	}
    }				/* end while() */
  delete_all_events ();
  close_file (of);
  free (node_states);
  free (S_k);
  free (I_k);

  FILE *f_tree = open_output ("tree");
  FILE *f_coal_times = open_output ("coal_times");
  FILE *f_samp_times = open_output ("samp_times");

  int *sampled_pnode_IDs, *all_pnode_IDs;
  int sample_size = round (arguments.sample_fraction * gen_ind);
  if (sample_size == 0)
    {
      sample_size = 1;
    }
  sampled_pnode_IDs = (int *) malloc (sample_size * sizeof (int));
  all_pnode_IDs = (int *) malloc (gen_ind * sizeof (int));
  if (!sampled_pnode_IDs || !all_pnode_IDs)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  for (i = 0; i < gen_ind; i++)
    {
      all_pnode_IDs[i] = i;
      (ptree + i)->flag = 0;
    }
  gsl_ran_choose (rng, sampled_pnode_IDs, sample_size, all_pnode_IDs, gen_ind,
		  sizeof (int));
  for (i = 0; i < sample_size; i++)
    {
      int w = sampled_pnode_IDs[i];
      (ptree + w)->flag = 1;
      while (w != 0)
	{
	  w = (ptree + w)->abv;
	  (ptree + w)->flag = 1;
	}
    }
  prtree3 (ptree, gen_ind, &f_tree, &f_coal_times, &f_samp_times);
  gsl_rng_free (rng);
  free (all_pnode_IDs);
  free (sampled_pnode_IDs);


  free (ptree);
  free (initial_degree_seq);
  gsl_histogram_free (degree_dist);
  igraph_destroy (&g);
  close_file (f_tree);
  close_file (f_coal_times);
  close_file (f_samp_times);
  return 0;
}


int
infect (int infector_id, int infectee_id, igraph_t *g, 
        igraph_vector_t *degree)
{
  igraph_vector_t neighbors;
  int neighbor_id;
  int infectee_id_deg;
  struct event *ev1;
  size_t i;

  /* update system dynamic variable */
  node_states[infectee_id] = 'i';
  infectee_id_deg = round (igraph_vector_e (degree, infectee_id));
  S_k[infectee_id_deg]--;
  I_k[infectee_id_deg]++;

  /* add event of infectee recovering from the table */
  ev1 = (struct event *) malloc (sizeof (struct event));
  if (!ev1)
    {
      fprintf (stderr, "Error: %s: %d: malloc failed\n", __FILE__, __LINE__);
      return (1);
    }
  memset (ev1, 0, sizeof (struct event));
  ev1->ego_id = infectee_id;
  ev1->rate = recov_rate;
  ev1->event_code = RECOVER;
  HASH_ADD (hh, event_table, event_code, keylen, ev1);
  rate_sum += recov_rate;

  /* add events of infectee infecting susceptible neighbors and
   * remove events of infected host getting infeced by infectious
   * neighbors */

  igraph_vector_init (&neighbors, 0);
  igraph_neighbors (g, &neighbors, (igraph_integer_t) infectee_id, IGRAPH_ALL);
  for (i = 0; i < igraph_vector_size (&neighbors); i++)
    {
      neighbor_id = VECTOR (neighbors)[i];
      switch (node_states[neighbor_id])
	{
	case 's':
	  ev1 = (struct event *) malloc (sizeof (struct event));
	  if (!ev1)
	    {
	      fprintf (stderr, "Error: %s: %d: malloc failed\n", __FILE__,
		       __LINE__);
	      return (1);
	    }
	  memset (ev1, 0, sizeof (struct event));
	  ev1->ego_id = infectee_id;
	  ev1->alter_id = neighbor_id;
	  ev1->rate = trans_rate;
	  ev1->event_code = INFECT;
	  HASH_ADD (hh, event_table, event_code, keylen, ev1);
	  rate_sum += trans_rate;
	  break;
	case 'i':
	  delete_event_by_key (INFECT, neighbor_id, infectee_id);
	  rate_sum -= trans_rate;
	  break;
	case 'r':
	  break;
	default:
	  fprintf (stderr,
		   "Error: %s: %d: invalid state for node\n",
		   __FILE__, __LINE__);
	  exit (1);
	}
    }
  igraph_vector_destroy (&neighbors);
  return 0;
}

int
neighbor_state_count (int susceptible_id, int *num_S_nb, int *num_I_nb,
                      igraph_t *g)
{

  int neighbor_id;
  igraph_vector_t neighbors;
  size_t i;

  *num_S_nb = 0;
  *num_I_nb = 0;


  /* check for event of infector infecting neighbors */

  igraph_vector_init (&neighbors, 0);
  igraph_neighbors (g, &neighbors, susceptible_id, IGRAPH_ALL);

  for (i = 0; i < igraph_vector_size (&neighbors); i++)
    {
      neighbor_id = VECTOR (neighbors)[i];
      switch (node_states[neighbor_id])
	{
	case 's':
	  ++*num_S_nb;
	  break;
	case 'i':
	  ++*num_I_nb;
	  break;
	case 'r':
	  break;
	default:
	  fprintf (stderr,
		   "Error: %s: %d: invalid state for node\n",
		   __FILE__, __LINE__);
	  exit (1);
	}
    }
  igraph_vector_destroy (&neighbors);
  return 0;

}

int
check_infector_events (int infector_id, igraph_t *g)
{

  int neighbor_id, event_count = 0;
  struct event *ev1;
  igraph_vector_t neighbors;
  size_t i;

  /* check for recovery event */

  check_event_by_key (RECOVER, infector_id, 0);
  event_count++;

  /* check for event of infector infecting neighbors */

  igraph_vector_init (&neighbors, 0);
  igraph_neighbors (g, &neighbors, infector_id, IGRAPH_ALL);
  
  for (i = 0; i < igraph_vector_size (&neighbors); i++)
    {
      neighbor_id = VECTOR (neighbors)[i];
      switch (node_states[neighbor_id])
	{
	case 's':
	  check_event_by_key (INFECT, infector_id, neighbor_id);
	  event_count++;
	  break;
	case 'i':
	  break;
	case 'r':
	  break;
	default:
	  fprintf (stderr,
		   "Error: %s: %d: invalid state for node\n",
		   __FILE__, __LINE__);
	  exit (1);
	}
    }
  igraph_vector_destroy (&neighbors);
  return event_count;
}


int
recover (int recoverer_id, igraph_t *g, igraph_vector_t *degree)
{

  int recoverer_id_deg, neighbor_id;
  struct event *ev1;
  igraph_vector_t neighbors;
  size_t i;

  /* update system dynamic variable */

  node_states[recoverer_id] = 'r';
  recoverer_id_deg = round (igraph_vector_e (degree, recoverer_id));
  I_k[recoverer_id_deg]--;

  /* remove event of recoverer recovering again */

  delete_event_by_key (RECOVER, recoverer_id, 0);
  rate_sum -= recov_rate;

  /* remove event of recoverer infecting neighbors */

  igraph_vector_init (&neighbors, 0);
  igraph_neighbors (g, &neighbors, recoverer_id, IGRAPH_ALL);
  for (i = 0; i < igraph_vector_size (&neighbors); i++)
    {
      neighbor_id = VECTOR (neighbors)[i];
      switch (node_states[neighbor_id])
	{
	case 's':
	  delete_event_by_key (INFECT, recoverer_id, neighbor_id);
	  rate_sum -= trans_rate;
	  break;
	case 'i':
	  break;
	case 'r':
	  break;
	default:
	  fprintf (stderr,
		   "Error: %s: %d: invalid state for node\n",
		   __FILE__, __LINE__);
	  exit (1);
	}
    }
  igraph_vector_destroy (&neighbors);
  return 0;
}

double
get_infection_rate_sum ()
{
  struct event *ev;
  double infection_rate_sum = 0;

  for (ev = event_table; ev != NULL; ev = (struct event *) ev->hh.next)
    {
      if ( ev->event_code == INFECT)
        {
          infection_rate_sum += ev->rate;
        }
    }
  return infection_rate_sum;
}

int
get_num_infected (int max_degree)
{
  int i;
  int num_infected = 0;

  for (i = 0; i <= max_degree; i++)
    {
      num_infected += I_k[i];
    }
  return num_infected;
}

double
get_approx_case_prod_mean (double gen_time, int num_infected)
{
  struct event *ev;
  double mean_sum = 0;

  /* probability of transmission over edge within GEN_TIME */
  double p;
  
  /* individual case production mean */
  double icpm;

  for (ev = event_table; ev != NULL; ev = (struct event *) ev->hh.next)
    {
      if ( ev->event_code == INFECT)
        {
          p = 1 - exp (-ev->rate * gen_time);
          mean_sum += p;
        }
    }
  icpm = mean_sum / num_infected;
  return icpm;
}

int
get_approx_case_prod_variance (double trans, double gen_time, int num_infected, 
                               double *icpm, double *icpv, igraph_t *g)
{
  int num_S_nb;
  int num_I_nb;
  int i;

  double p;
  double np;

  double E_y = 0;
  double E_y_sq = 0;

  p = 1 - exp (-trans * gen_time ); 

  for (i = 0; node_states[i] != 0; i++)
    {
      if (node_states[i] == 'i')
        {
          neighbor_state_count (i, &num_S_nb, &num_I_nb, g);
          np = num_S_nb * p;
          E_y += np;
          E_y_sq += np * (np - p + 1);
        }
    }
  *icpm = (double) E_y / num_infected;
  *icpv = (double) E_y_sq / num_infected  - (*icpm) * (*icpm);
  return 0;
}

double
get_rate_sum ()
{
  struct event *ev;
  double rate_sum_check = 0;

  for (ev = event_table; ev != NULL; ev = (struct event *) ev->hh.next)
    {
      rate_sum_check += ev->rate;
    }
  return rate_sum_check;
}

void
prtree3 (struct node *ptree, int num_nodes, 
         FILE** f_tree, FILE** f_coal_times, FILE** f_samp_times)
  /* print genealogy tree in Newick format */
  /* after Hudson's ms code */
{
  double t;
  int i, *descl, *descr, w, tmp;
  void parens (struct node *ptree, int *descl, int *descr, int noden, 
               FILE** f_tree, FILE** f_coal_times, FILE** f_samp_times);

  descl = (int *) malloc ((unsigned) (num_nodes) * sizeof (int));
  descr = (int *) malloc ((unsigned) (num_nodes) * sizeof (int));
  for (i = 0; i < num_nodes; i++)
    descl[i] = descr[i] = -1;
  for (i = 1; i < num_nodes; i++)
    {
      if (descl[(ptree + i)->abv] == -1 && (ptree + i)->flag == 1)
	descl[(ptree + i)->abv] = i;
      else if ((ptree + i)->flag == 1)
	descr[(ptree + i)->abv] = i;
    }
  for (i = 1; i < num_nodes; i++)
    {
      if (descl[i] != -1 && descr[i] == -1)
	{
	  w = descl[i];
	  while (descr[w] == -1)
	    {
	      tmp = descl[w];
	      descl[w] = -1;
	      if (tmp == -1)
		{
		  break;
		}
	      w = tmp;
	    }
	  if (tmp != -1)
	    {
	      descl[i] = descl[w];
	      descr[i] = descr[w];
              (ptree + descl[i])->abv = i;
              (ptree + descr[i])->abv = i;
              (ptree + i)->time = (ptree + w)->time;
	      descl[w] = -1;
	      descr[w] = -1;
	    }
	  if (tmp == -1)
	    {
	      descl[i] = -1;
	      tmp = 0;
	      while (descl[tmp] != i && descr[tmp] != i)
		{
		  tmp++;
		}
	      if (descl[tmp] == i)
		{
		  descl[tmp] = w;
                  (ptree + w)->abv = tmp;
		}
	      if (descr[tmp] == i)
		{
		  descr[tmp] = w;
                  (ptree + w)->abv = tmp;
		}
	    }
	}
    }
  parens (ptree, descl, descr, 1, f_tree, f_coal_times, f_samp_times);
  fprintf (*f_coal_times, "\n");
  fprintf (*f_samp_times, "\n");

  free (descl);
  free (descr);
}

void
parens (struct node *ptree, int *descl, int *descr, int noden, 
        FILE** f_tree, FILE** f_coal_times, FILE** f_samp_times)
  /* print genealogy tree in Newick format */
  /* after Hudson's ms code */
{
  double time;

  if (descl[noden] == -1)
    {
      time = (ptree + (ptree + noden)->abv)->time - (ptree + noden)->time;
      fprintf (*f_tree, "%d_%g:%g", noden + 1, (ptree + noden)->time, -1*time);
      fprintf (*f_samp_times, "%g ", (ptree + noden)->time);
    }
  else
    {
      fprintf (*f_coal_times, "%g ", (ptree + noden)->time);
      fprintf (*f_tree, "(");
      parens (ptree, descl, descr, descl[noden], f_tree, 
              f_coal_times, f_samp_times);
      fprintf (*f_tree, ",");
      parens (ptree, descl, descr, descr[noden], f_tree, 
              f_coal_times, f_samp_times);
      if ((ptree + noden)->abv == 0)
	{
	  fprintf (*f_tree, ");\n");
	}
      else
	{
	  time = (ptree + (ptree + noden)->abv)->time - (ptree + noden)->time;
	  fprintf (*f_tree, "):%g", -1*time);
	}
    }
}

int
get_degree_seq (int *resulting_degree_seq, int num_nodes, 
                    double  par1, double par2, int net_type)
{
  int i, degree_sum;
  double mean_degree;
  const gsl_rng_type * T;
  gsl_rng *r;

  if (net_type == DYNET_PARETO)
    {
      gsl_rng_env_setup ();
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);
      degree_sum = 0;
      for (i = 0; i < num_nodes; i++)
      {
        resulting_degree_seq[i] = gsl_ran_pareto (r, par1, par2);
        degree_sum += resulting_degree_seq[i];
      }
      if (degree_sum % 2 != 0)
        {
          fprintf (stderr,
                   "Warning: adjusting degree sequence so that it's even\n");
          resulting_degree_seq[i-1]++;
        }
      gsl_rng_free (r);
    }
  if (net_type == DYNET_EXP)
    {
      mean_degree = par1;
      gsl_rng_env_setup ();
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);
      degree_sum = 0;
      for (i = 0; i < num_nodes; i++)
      {
        resulting_degree_seq[i] = round(gsl_ran_exponential (r, mean_degree));
        if (resulting_degree_seq[i] == 0)
          {
            resulting_degree_seq[i]++;
          }
        degree_sum += resulting_degree_seq[i];
      }
      if (degree_sum % 2 != 0)
        {
          fprintf (stderr,
                   "Warning: adjusting degree sequence so that it's even\n");
          resulting_degree_seq[i-1]++;
        }
      gsl_rng_free (r);
    }
  if (net_type == DYNET_REG)
    {
      mean_degree = par1;
      if (mean_degree - (int) mean_degree > 0 || (int) mean_degree *num_nodes % 2 !=0)
        {
          fprintf (stderr, "Error: %s: %d:\
 mean degree impossible for regular net\n", __FILE__, __LINE__);
          return 1;
        }
      for (i = 0; i <num_nodes; i++)
        {
          resulting_degree_seq[i] = (int) mean_degree;
        }
    }
  return 0;
}

int
sequence_to_graph(int *sequence, int num_nodes, igraph_t* g) 
{

/* construct a graph from a degree sequence*/

  igraph_vector_t v;
  size_t i;

  igraph_vector_init(&v, 0);
  igraph_vector_resize(&v, num_nodes);
  
  for (i=0; i<num_nodes; i++)
    {
      VECTOR(v)[i]=sequence[i];
    }
  igraph_degree_sequence_game(g, &v, NULL, IGRAPH_DEGSEQ_VL);
/* TODO 
 * figure out how to use these function so that initial graph can have multiple components and isolates 
  igraph_degree_sequence_game(g, &v, NULL, IGRAPH_DEGSEQ_SIMPLE);
  igraph_simplify (g, 1, 1);
  */
  igraph_vector_destroy(&v);
  
  return 0;
}

FILE *
open_output (const char *filename)
     /* open filename for output; return NULL if problem */
     /* after p. 367, C A Reference Manual 5e */
{
  FILE *f;
  errno = 0;
  /* Functions below might choke on a NULL filename. */
  if (filename == NULL) filename = "";
  f = fopen (filename, "w");
  if (f == NULL)
    {
      error_at_line (EXIT_FAILURE, errno, __FILE__, __LINE__,
                     "Failed to open file %s", filename);
    }
  return f;
}

int 
close_file (FILE *f)
     /* Close file f */
     /* after p. 367, C A Reference Manual 5e */
{
  int s = 0;
  if (f == NULL)
    {
      return 0;
    }
  errno = 0;
  s = fclose (f);
  if (s == EOF) perror ("Close failed");
  return s;
}

