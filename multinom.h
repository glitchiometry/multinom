#include "basics.h"
#include "stdio.h"

#define DTYPE_ARRAY_VOIDSTAR 4
#define DTYPE_CONST_INT 1
#define DTYPE_CONST_FLOAT 2
#define DTYPE_INDETERMINATE_CONST 3
#define DTYPE_CONST_ZERO 0
#define M_EV_ROOT_TYPE 0
#define M_EV_COMP_TYPE 1
#define MULTINOM_ADD 0
#define MULTINOM_SUB 1
#define MULTINOM_MULT 2

/*
 * NOTES: 
 * 	The goal of this (WIP) library is to implement fundamental routines for performing calculations with polynomials in several variables.
 * 	Eventually, the hope is to accommodate polynomials over arbitrary (commutative) rings (but maybe even non-commutative rings.)
 * 	One important motive for this library is to facilitate the discovery of efficient algorithms for evaluating multivariable polynomials.
 * 	The na'ive approach (generalizing Horner's rule for dense single variable polynomials) may not be generically optimal in the multivariate
 * 	case (possibly even for single variable polynomials.) There does not appear to be an obvious structure to optimal 
 * 	polynomial evaluation schemes (at least as far as this author can tell), but this does not imply that such structure does not exist.
 * 	The hope is that by facilitating discovery of optimal or near-optimal evaluation algorithms, any (possibly disordered or semi-random)
 * 	structure can be examined statistically.
 * 	Toward this end, in addition to performing basic arithmetic operations on multivariable polynomials (addition, subtraction, multiplication)
 * 	it is necessary to represent the general form of a (multivariable) polynomial with unknown/variable coefficients, as well as the general forms
 * 	of polynomials that can be obtained by combining such general polynomials.
 * 	Perhaps the simplest approach to achieving this goal is to express the unknown/variable coefficients as additional variables (expanding the 
 * 	'dimensionality' of the multivariable polynomial) and performing the usual arithmetic operations. The challenge then becomes how to decompose
 * 	the set of polynomials in the last 'm' dimensions say into elementary (algebraically independent) factors/components.  
 * 	To be continued...
 */

typedef struct
{
	int sym_id;
	void *coeff;
} indeterminate;

typedef struct
{
	int dim;
	int deg;
	// An array of 'coefficients' consisting of multinomials of lower dimensionality (in dim - 1 variables)
	void *data;
	int dtype;
	//array_voidstar coeff;
	void *src;
} multinomial_exp;

void multinomial_exp_init(multinomial_exp *mne, int dim, int deg, void *data, int dtype);
void free_multinomial_exp(void *mne);
void negative_multinomial_exp(multinomial_exp *mne);
void add_sub_multinomial_data(multinomial_exp *mne, multinomial_exp *sub_mne, char mode);
void add_sub_multinomial_exp(multinomial_exp *mne, multinomial_exp *sub_mne, char mode);
void add2multinomial_exp(multinomial_exp *mne, multinomial_exp *sub_mne);
void multiply_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3);
void compose_multinomial_exp(multinomial_exp **mne1, int dim1, multinomial_exp **mne2, int dim2, multinomial_exp **mne3);
void transcribe_multinomial_exp(multinomial_exp *src, multinomial_exp *dest);
void display_multinomial_exp(multinomial_exp *mne);
void multinomial_exp2str(multinomial_exp *mne, array_char *mne_str);
void pow_mutinomial_exp(multinomial_exp *mne, multinomial_exp *mne_p, int p);
char nonzero_test(multinomial_exp *mne);
void polynomial_int2multinomial_exp(array_int *p, multinomial_exp *mne, int dim);
void polynomial_double2multinomial_exp(array_double *p, multinomial_exp *mne, int dim);
void add_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3);
int multinomial_exp_eval_int(multinomial_exp *mne, int *vals);
double multinomial_exp_eval_double(multinomial_exp *mne, double *vals);
int count_indeterminates(multinomial_exp *mne);
char multinomial_exp_equals(multinomial_exp *mne1, multinomial_exp *mne2);
void truncate_multinomial_exp(multinomial_exp *mne, multinomial_exp *tmne, int dim);
void indetermify_multinomial_exp(multinomial_exp *mne, multinomial_exp *imne);

int multinomial_exp_complexity(multinomial_exp *mne);

typedef struct
{
	multinomial_exp *expr;
	void *expr_1;
	void *expr_2;
	char op;
	char type;
} multinomial_evaluator;

void multinomial_evaluator_init(multinomial_evaluator *mnev, multinomial_exp *expr);
void free_multinomial_evaluator(multinomial_evaluator *mnev);
void multinomial_evaluator_combine(multinomial_evaluator *mnev1, multinomial_evaluator *mnev2, char op, multinomial_evaluator *comp);
void multinomial_evaluator_combine_exp(multinomial_evaluator *mnev1, multinomial_evaluator *mnev2, char op, multinomial_evaluator *comp);
int multinomial_evaluator_complexity(multinomial_evaluator *mnev);
void multinomial_evaluator_eval_double(multinomial_evaluator *mnev, double *vals, double *cval);
void multinomial_evaluator_eval_int(multinomial_evaluator *mnev, int *vals, int *cval);
