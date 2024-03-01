#include "multinom.h"

char recognized_type(int t)
{
	if (t == DTYPE_MULTINOMIAL ||
		t == DTYPE_CONST_INT ||
		t == DTYPE_CONST_FLOAT ||
		t == DTYPE_CONST_COMPLEX ||
		t == DTYPE_INDETERMINATE_CONST ||
		t == DTYPE_CONST_ZERO ||
		t == DTYPE_CONST) return 1;
	else return 0;
}

char check_types(multinomial_exp *mne)
{
	if (mne != NULL)
	{
		if (recognized_type((*mne).dtype))
		{
			if ((*mne).dtype == DTYPE_MULTINOMIAL)
			{
				array_voidstar *mne_data = (array_voidstar *) (*mne).data;
				for (int i = 0; i < (*mne_data).len; i++)
				{
					multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[i];
					if (check_types(aux)) {}
					else return 0;
				}
			}
		}
		else return 0;
	}
	return 1;
}

void multinomial_exp_init_zero(multinomial_exp *mne)
{
	(*mne).dtype = DTYPE_CONST_ZERO;
	(*mne).deg = 0;
	(*mne).dim = -1;
	(*mne).data = NULL;
}

void multinomial_exp_init_int(multinomial_exp *mne, int val)
{
	multinomial_exp_init(mne, -1, 0, &val, DTYPE_CONST_INT);
}

void multinomial_exp_init_double(multinomial_exp *mne, double val)
{
	multinomial_exp_init(mne, -1, 0, &val, DTYPE_CONST_FLOAT);
}

void multinomial_exp_init_complex(multinomial_exp *mne, complex val)
{
	multinomial_exp_init(mne, -1, 0, &val, DTYPE_CONST_COMPLEX);
}


void multinomial_exp_init(multinomial_exp *mne, int dim, int deg, void *data, int dtype)
{
	(*mne).dtype = dtype;
	(*mne).deg = deg;
	(*mne).dim = dim;
	if ((dtype == DTYPE_MULTINOMIAL) && (data != NULL))
	{
		int data_len = deg + 1;
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *in_data = (array_voidstar *) data;
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		if ((*in_data).len < data_len)
		{
			printf("Warning (multinomial_exp_init): Mismatch between degree + 1 = %d and data length = %d\n", deg + 1, (*in_data).len); 
		}
		array_voidstar_init(mne_data, data_len);
		(*mne_data).len = data_len;
		for (int i = 0; i < data_len; i++)
		{
			if ((*in_data).e[i] != NULL) 
			{
				multinomial_exp *in_m_exp = (multinomial_exp *) (*in_data).e[i];
				(*mne_data).e[i] = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				transcribe_multinomial_exp(in_m_exp, (multinomial_exp *) (*mne_data).e[i]);
			}
			else
			{
				(*mne_data).e[i] = NULL;
			}
		}
		for (int i = (*in_data).len; i <= deg; i++)
		{
			(*mne_data).e[i] = NULL;
		}
	}
	else 
	{
		if (dtype == DTYPE_MULTINOMIAL)
		{
			(*mne).data = NULL;
		}
		if (dtype == DTYPE_INDETERMINATE_CONST)
		{
			if (data != NULL)
			{
				(*mne).data = (indeterminate *) calloc(1, sizeof(indeterminate));
				*((indeterminate *) (*mne).data) = *((indeterminate *) data);
			}
			else (*mne).data = NULL;
			//(*mne).data = NULL;
		}
		if (dtype == DTYPE_CONST_INT)
		{
			(*mne).data = (int *) calloc(1, sizeof(int));
			*((int *) (*mne).data) = *((int *) data);
		}
		if (dtype == DTYPE_CONST_FLOAT)
		{
			(*mne).data = (double *) calloc(1, sizeof(double));
			*((double *) (*mne).data) = *((double *) data);
		}
		if (dtype == DTYPE_CONST_COMPLEX)
		{
			(*mne).data = (complex *) calloc(1, sizeof(complex));
			*((complex *) (*mne).data) = *((complex *) data);
		}
	}
}

void permute_multinomial_exp_shallow(multinomial_exp *mne, int (*perm)(int))
{
	if ((*mne).dim > -1)
	{
		(*mne).dim = perm((*mne).dim);
		if ((*mne).dtype == DTYPE_MULTINOMIAL)
		{
			if ((*mne).data != NULL) {}
			else
			{
				printf("Something weird happened in permute_multinomial_exp_shallow: multinomials with DTYPE_MULTINOMIAL shouldn't have NULL data fields\n");
				exit(EXIT_FAILURE);
			}
			array_voidstar *mne_data = (array_voidstar *) (*mne).data;
			for (int i = 0; i < (*mne_data).len; i++)
			{
				if ((*mne_data).e[i] != NULL) permute_multinomial_exp_shallow((multinomial_exp *) (*mne_data).e[i], perm);
			}
		}
	}
}


void aux_func_simplify_multinomial_exp(multinomial_exp *mne)
{
	if ((*mne).dtype == DTYPE_MULTINOMIAL)
	{
		int deg = -1;
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		for (int i = 0; i < (*mne_data).len; i++)
		{
			if (nonzero_test((multinomial_exp *) (*mne_data).e[i]))
			{
				deg = i;
			}
		}
		if (deg > -1)
		{
			for (int i = (*mne).deg; i > deg; i--)
			{
				remove_array_voidstar(mne_data, i, free_multinomial_exp);
			}
			(*mne).deg = deg;
			if ((*mne).deg == 0)
			{
				// Replace mne by its (simplified) 0th degree term
				multinomial_exp *aux_mne = (multinomial_exp *) (*mne_data).e[0];
				aux_func_simplify_multinomial_exp(aux_mne);
				(*mne) = (*aux_mne);
				free_array_voidstar(mne_data, free_multinomial_exp);
				free(aux_mne); // RESUME: check this! (Freeing this pointer may be unnecessary)
			}
		}
		else if (deg == -1)
		{
			free_multinomial_exp(mne);
			multinomial_exp_init_zero(mne);
		}
	}
	else
	{
		if (nonzero_test(mne)) {}
		else if ((*mne).dtype != DTYPE_CONST_ZERO)
		{
			free_multinomial_exp(mne);
			multinomial_exp_init_zero(mne);
		}
	}
}

// NOTE: this could probably be made more efficient by passing the multinomial_exp as a pointer-of-pointer
// 	(i.e. multinomial_exp **mne); for most applications, however, (e.g. 'sorting' multinomials after permuting
// 	their variables) the gain in performance would probably be negligible.
void add2multinomial_exp_monomial_rec(multinomial_exp *mne, int *vars, int *pows, int len, multinomial_exp *coeff)
{
	if (mne != NULL) {}
	else
	{
		printf("Error (add2multinomial_exp_monomial_rec): specified null multinomial_exp pointer\n");
		exit(EXIT_FAILURE);
	}
	if (len > 0)
	{
		if (has_constant_type(mne) || ((*mne).dim > vars[0]))
		{
			// Replace mne with a multinomial whose constant term points to mne
			multinomial_exp aux;
			array_voidstar aux_data;
			array_voidstar_init(&aux_data, pows[0] + 1);
			aux_data.len = pows[0] + 1;
			// NOTE: consider replacing array_voidstar entries for each power with occupancy lists
			// 	for sparse multinomials.
			if (nonzero_test(mne))
			{
				multinomial_exp *aux0 = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				transcribe_multinomial_exp(mne, aux0);
				aux_data.e[0] = (void *) aux0;
			}
			free_multinomial_exp(mne);
			multinomial_exp_init(&aux, vars[0], pows[0], &aux_data, DTYPE_MULTINOMIAL);
			free_array_voidstar(&aux_data, free_multinomial_exp);
			(*mne) = aux;
		}
		if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).dim < vars[0])
		{
			array_voidstar *mne_data = (array_voidstar *) (*mne).data;
			if ((*mne_data).len > 0) {}
			else add2array_voidstar(mne_data, NULL);
			if ((*mne_data).e[0] != NULL) {}
			else
			{
				multinomial_exp *aux = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init_zero(aux);
				(*mne_data).e[0] = aux;
			}
			add2multinomial_exp_monomial_rec((multinomial_exp *) (*mne_data).e[0], vars, pows, len, coeff);
			if (nonzero_test((multinomial_exp *) (*mne_data).e[0])) {}
			else
			{
				free_multinomial_exp((multinomial_exp *) (*mne_data).e[0]);
				(*mne_data).e[0] = NULL;
			}
		}
		else if ((*mne).dim == vars[0])
		{
			array_voidstar *mne_data = (array_voidstar *) (*mne).data;
			if ((*mne_data).len <= pows[0])
			{
				while ((*mne_data).len <= pows[0]) add2array_voidstar(mne_data, NULL);
			}
			if ((*mne_data).e[pows[0]] != NULL) {}
			else
			{
				multinomial_exp *aux1 = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init_zero(aux1);
				(*mne_data).e[pows[0]] = (void *) aux1;
			}
			add2multinomial_exp_monomial_rec((multinomial_exp *) (*mne_data).e[pows[0]], &(vars[1]), &(pows[1]), len - 1, coeff);
			// NOTE: until nonzero_test is improved, this step will be somewhat slow for large multinomials
			if (nonzero_test((multinomial_exp *) (*mne_data).e[pows[0]])) 
			{
				(*mne).deg = (*mne).deg > pows[0] ? (*mne).deg : pows[0];
			} 
			else 
			{
				free_multinomial_exp((multinomial_exp *) (*mne_data).e[pows[0]]);
				(*mne_data).e[pows[0]] = NULL;
				if ((*mne).deg == pows[0]) 
				{
					while ((*mne).deg > -1 && (*mne_data).e[(*mne).deg] == NULL) (*mne).deg -= 1;
					if ((*mne).deg > 0) {}
					else
					{
						aux_func_simplify_multinomial_exp(mne); // RESUME: test this!
					}
				}
			}
		}
		else
		{
			printf("Error (add2multinomial_exp_monomial_rec): unhandled case\n");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		// Add coeff to the constant term of mne
		add2multinomial_exp(mne, coeff);
	}
}

void add2multinomial_exp_monomial(multinomial_exp *mne, int *vars, int *pows, int len, multinomial_exp *coeff)
{
	array_int pi;
	array_int_init(&pi, len);
	for (int i = 0; i < len; i++) pi.e[i] = i;
	pi.len = len;
	array_int avars;
	avars.e = vars;
	avars.len = len;
	sort_array_int_permutation(&avars, &pi);
	int aux_v[len], aux_p[len];
	for (int i = 0; i < (avars).len; i++) aux_v[i] = (avars).e[pi.e[i]];
	for (int i = 0; i < (avars).len; i++) aux_p[i] = pows[pi.e[i]];
	add2multinomial_exp_monomial_rec(mne, &(aux_v[0]), &(aux_p[0]), len, coeff);
	free_array_int(&pi);
}

void sort_multinomial_exp_rec(multinomial_exp *mne, array_int *vars, array_int *pows, multinomial_exp *smne)
{
	if ((*mne).dtype == DTYPE_MULTINOMIAL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		add2array_int(vars, (*mne).dim);
		for (int i = 0; i < (*mne_data).len; i++)
		{
			if ((*mne_data).e[i] != NULL)
			{
				int pows_len = (*pows).len;
				add2array_int(pows, i);
				multinomial_exp *term = (multinomial_exp *) (*mne_data).e[i];
				sort_multinomial_exp_rec(term, vars, pows, smne);
				remove_array_int(pows, pows_len);
			}
		}
		remove_array_int(vars, (*vars).len - 1);
	}
	else
	{
		// Define a monomial and add it to (*smne)
		if (nonzero_test(mne))
		{
			add2multinomial_exp_monomial(smne, (*vars).e, (*pows).e, (*vars).len, mne); 
		}
	}
}

void sort_multinomial_exp(multinomial_exp *mne, multinomial_exp *smne)
{
	if ((*mne).dtype == DTYPE_MULTINOMIAL)
	{
		// RESUME
		// Collect terms multiplying each power of the lowest dimension (or variable of smallest index.)
		// Simple (but possibly slow) approach: 
		// 	Add multinomials to a graded multinomial
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		array_int vars;
		array_int pows;
		array_int_init(&vars, 1);
		array_int_init(&pows, 1);
		multinomial_exp_init_zero(smne);
		sort_multinomial_exp_rec(mne, &vars, &pows, smne);
	}
}

void permute_multinomial_exp(multinomial_exp *mne, int (*perm)(int), multinomial_exp *pmne)
{
	if (mne != NULL) {}
	else return;
	if (perm != NULL && (*mne).dtype == DTYPE_MULTINOMIAL) {}
	else
	{
		printf("Permute: Edge case ");
		display_multinomial_exp(mne);
		transcribe_multinomial_exp(mne, pmne);
		return;
	}
	multinomial_exp spmne;
	transcribe_multinomial_exp(mne, &spmne);
	permute_multinomial_exp_shallow(&spmne, perm);
	display_multinomial_exp(&spmne);
	sort_multinomial_exp(&spmne, pmne);
	free_multinomial_exp(&spmne);
}

void free_multinomial_exp(void *mne)
{
	if (mne != NULL)
	{
		multinomial_exp *aux = (multinomial_exp *) mne;
		if (nonzero_test(aux))
		{
			if ((*aux).dtype == DTYPE_MULTINOMIAL && (*aux).data != NULL)
			{
				free_array_voidstar((array_voidstar *) (*aux).data, free_multinomial_exp);
				free((*aux).data);
			}
			else 
			{
				// Address other nontrivial cases here
				free((*aux).data);
			}
		}
		else if ((*aux).dtype != DTYPE_CONST_ZERO)
		{
			printf("Warning: freeing nontrivial zero multinomial\n");
			if ((*aux).dtype == DTYPE_MULTINOMIAL && (*aux).data != NULL)
			{
				array_voidstar *aux_data = (array_voidstar *) (*aux).data;
				free_array_voidstar(aux_data, free_multinomial_exp);
				free(aux_data);
			}
			else if ((*aux).data != NULL)
			{
				free((*aux).data);
			}
		}
	}
}

void transcribe_multinomial_exp(multinomial_exp *src, multinomial_exp *dest)
{
	multinomial_exp_init(dest, (*src).dim, (*src).deg, (*src).data, (*src).dtype);
}

char indeterminate_const_test(multinomial_exp *mne)
{
	return (*mne).dtype == DTYPE_INDETERMINATE_CONST;
}

char valid_numeric_test(multinomial_exp *mne)
{
	return ((*mne).dtype == DTYPE_CONST_INT || (*mne).dtype == DTYPE_CONST_FLOAT || (*mne).dtype == DTYPE_CONST_COMPLEX) && (*mne).data != NULL;
}

void add2multinomial_data(multinomial_exp *mne, multinomial_exp *sub_mne)
{
	//printf("Adding constants of type %d, %d\n", (*mne).dtype, (*sub_mne).dtype);
	if ((*mne).deg == 0 && (*sub_mne).deg == 0) {}
	else printf("Something weird happened!\n");
	if (indeterminate_const_test(mne)) return;
	else if (indeterminate_const_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, -1, 0, NULL, DTYPE_INDETERMINATE_CONST);
		return;
	}
	else if (valid_numeric_test(mne) && valid_numeric_test(sub_mne))
	{
		if ((*mne).dtype == DTYPE_CONST_INT)
		{
			int *mne_data = (int *) (*mne).data;
			if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				double *sub_mne_data = (double *) (*sub_mne).data;
				double *new_data = (double *) calloc(1, sizeof(double));
				(*new_data) = (*sub_mne_data) + (*mne_data);
				free(mne_data);
				(*mne).data = (void *) new_data;
				(*mne).dtype = DTYPE_CONST_FLOAT;
			}
			if ((*sub_mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex *sub_mne_data = (complex *) (*sub_mne).data;
				complex *new_data = (complex *) calloc(1, sizeof(complex));
				(*new_data) = (*sub_mne_data) + (*mne_data);
				free(mne_data);
				(*mne).data = (void *) new_data;
				(*mne).dtype = DTYPE_CONST_COMPLEX;
			}
		}
		else if ((*mne).dtype == DTYPE_CONST_FLOAT)
		{
			double *mne_data = (double *) (*mne).data;
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				double *sub_mne_data = (double *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex *sub_mne_data = (complex *) (*sub_mne).data;
				complex *new_data = (complex *) calloc(1, sizeof(complex));
				(*new_data) = (*sub_mne_data) + (*mne_data);
				free(mne_data);
				(*mne).data = (void *) new_data;
				(*mne).dtype = DTYPE_CONST_COMPLEX;
			}
		}
		else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
		{
			complex *mne_data = (complex *) (*mne).data;
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				double *sub_mne_data = (double *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex *sub_mne_data = (complex *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
		}
	}
	else if (!nonzero_test(mne) && nonzero_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, (*sub_mne).dim, (*sub_mne).deg, (*sub_mne).data, (*sub_mne).dtype);
	}
}

void sub_from_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2)
{
	add_sub_multinomial_exp(mne1, mne2, MULTINOM_SUB);
}

void sub_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3)
{
	transcribe_multinomial_exp(mne1, mne3);
	sub_from_multinomial_exp(mne3, mne2);
}

void add_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3)
{
	transcribe_multinomial_exp(mne1, mne3);
	add2multinomial_exp(mne3, mne2);
}

void negative_multinomial_exp(multinomial_exp *mne)
{
	// TBD
	if (mne != NULL) {}
	else return;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		for (int i = 0; i < (*mne_data).len; i++)
		{
			if ((*mne_data).e[i] != NULL) negative_multinomial_exp((multinomial_exp *) (*mne_data).e[i]);
		}
	}
	else if ((*mne).dtype == DTYPE_CONST_INT && (*mne).data != NULL)
	{
		int *val = (int *) (*mne).data;
		(*val) = -(*val);
	}
	else if ((*mne).dtype == DTYPE_CONST_FLOAT && (*mne).data != NULL)
	{
		double *val = (double *) (*mne).data;
		(*val) = -(*val);
	}
	else if ((*mne).dtype == DTYPE_CONST_COMPLEX && (*mne).data != NULL)
	{
		complex *val = (complex *) (*mne).data;
		(*val) = -(*val);
	}
	// RESUME: address indeterminate constant case to allow symbolic arithmetic (although this may be easier to accomplish simply by expanding the space of variables)
}

void add_sub_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, char mode)
{
	if (nonzero_test(mne2)) {}
	else return;
	if (check_types(mne1)) {}
	else
	{
		printf("Failed type check for mne1 at beginning of add_sub_multinomial_exp\n");
		exit(EXIT_FAILURE);
	}
	if (check_types(mne2)) {}
	else
	{
		printf("Failed type check for mne2 at beginning of add_sub_multinomial_exp\n");
		exit(EXIT_FAILURE);
	}
	if ((((*mne1).dim < (*mne2).dim) || ((*mne2).dim < 0)) && (*mne1).dtype == DTYPE_MULTINOMIAL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne1).data;
		if ((*mne_data).len > 0) {}
		else add2array_voidstar(mne_data, NULL);
		if ((*mne_data).e[0] != NULL)
		{
			multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[0];
			//add2multinomial_exp(aux, sub_mne);
			add_sub_multinomial_exp(aux, mne2, mode);
		}
		else
		{
			(*mne_data).e[0] = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
			transcribe_multinomial_exp(mne2, (multinomial_exp *) (*mne_data).e[0]);
			if (mode == MULTINOM_ADD) {}
			else if (mode == MULTINOM_SUB) 
			{
				negative_multinomial_exp((multinomial_exp *) (*mne_data).e[0]); 
			}
		}
	}
	else if ((((*mne1).dim > (*mne2).dim) || ((*mne1).dim == -1)) && ((*mne2).dtype == DTYPE_MULTINOMIAL))
	{
		multinomial_exp prec;
	       	transcribe_multinomial_exp(mne2, &prec);
		if (mode == MULTINOM_SUB) negative_multinomial_exp(&prec);
		add_sub_multinomial_exp(&prec, mne1, MULTINOM_ADD);
		free_multinomial_exp(mne1);
		(*mne1) = prec;
		if (check_types(mne1)) {}
		else
		{
			printf("Failed type check after swap\n");
			exit(EXIT_FAILURE);
		}
	}
	else if ((*mne1).dim == (*mne2).dim)
	{
		if ((*mne1).dim > -1)
		{
			char nonzero_status = 0; // RESUME: update this at each step of the calculation where mne_data is changed, and set mne1 to zero if non_zero_status is still zero
			int min_deg;
			array_voidstar *mne_data = (array_voidstar *) (*mne1).data;
			array_voidstar *sub_mne_data = (array_voidstar *) (*mne2).data;
			if ((*mne1).deg < (*mne2).deg)
			{
				min_deg = (*mne1).deg;
				add_mem_array_voidstar_until(mne_data, (*mne2).deg + 1);
				(*mne_data).len = (*mne2).deg + 1;
				for (int i = (*mne1).deg + 1; i <= (*mne2).deg; i++)
				{
					if (i < (*sub_mne_data).len && (*sub_mne_data).e[i] != NULL) 
					{
						multinomial_exp *sub_mne_term = (multinomial_exp *) (*sub_mne_data).e[i];
						if (check_types(sub_mne_term)) {}
						else
						{
							printf("Error: sub_mne_term has unrecognized type\n");
							exit(EXIT_FAILURE);
						}
						if (nonzero_test(sub_mne_term))
						{
							nonzero_status = 1;
							(*mne_data).e[i] = (void *) ((multinomial_exp *) calloc(1, sizeof(multinomial_exp)));
							transcribe_multinomial_exp((multinomial_exp *) (*sub_mne_data).e[i], (multinomial_exp *) (*mne_data).e[i]);
							if (mode == MULTINOM_ADD) {}
							else if (mode == MULTINOM_SUB)
							{
								negative_multinomial_exp((multinomial_exp *) (*mne_data).e[i]);
							}
						}
					}
					else (*mne_data).e[i] = NULL;
				}
				(*mne1).deg = (*mne2).deg;
			}
			else min_deg = (*mne2).deg;
			for (int i = 0; i <= min_deg; i++)
			{
				multinomial_exp *mne_desc = (multinomial_exp *) (*mne_data).e[i];
				multinomial_exp *sub_mne_desc = (multinomial_exp *) (*sub_mne_data).e[i];
				if (mne_desc != NULL && sub_mne_desc != NULL) 
				{
					add_sub_multinomial_exp(mne_desc, sub_mne_desc, mode);
					if (nonzero_test(mne_desc)) 
					{
						nonzero_status = 1;
					}
					else
					{
						free_multinomial_exp(mne_desc);
						(*mne_data).e[i] = NULL;
					}
				}
				else
				{
					if (mne_desc == NULL && sub_mne_desc != NULL)
					{
						if (nonzero_test(sub_mne_desc))
						{
							nonzero_status = 1;
							(*mne_data).e[i] = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
							mne_desc = (multinomial_exp *) (*mne_data).e[i];
							transcribe_multinomial_exp(sub_mne_desc, mne_desc);
							if (mode == MULTINOM_ADD) {}
							else
							{
								negative_multinomial_exp(mne_desc);
							}
						}
					}
				}
			}
			if (nonzero_status) {}
			else
			{
				free_array_voidstar(mne_data, free_multinomial_exp);
				(*mne1).data = NULL;
				(*mne1).dtype = DTYPE_CONST_ZERO;
			}
		}
		else
		{
			add_sub_multinomial_data(mne1, mne2, mode);
			// add2multinomial_data(mne, sub_mne); 
		}
	}
	else
	{
		printf("Error (add_sub_multinomial_exp): unhandled case dim1 = %d, dim2 = %d\n", (*mne1).dim, (*mne2).dim);
		exit(EXIT_FAILURE);
	}
	if (check_types(mne1)) {}
	else
	{
		printf("Failed type check in add_sub_multinomial_exp\n");
		exit(EXIT_FAILURE);
	}
} // END add_sub_multinomial_exp

void add_sub_multinomial_data(multinomial_exp *mne, multinomial_exp *sub_mne, char mode)
{
	if ((*mne).deg == 0 && (*sub_mne).deg == 0) {}
	else printf("Something weird happened! (add_sub_multinomial_data)\n");
	if (indeterminate_const_test(mne)) return;
	else if (indeterminate_const_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, -1, 0, NULL, DTYPE_INDETERMINATE_CONST);
		return;
	}
	else if (valid_numeric_test(mne) && valid_numeric_test(sub_mne))
	{
		if ((*mne).dtype == DTYPE_CONST_INT)
		{
			int *mne_data = (int *) (*mne).data;
			if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
				if ((*mne_data) == 0)
				{
					free(mne_data);
					(*mne).data = NULL;
					(*mne).dtype = DTYPE_CONST_ZERO;
				}
			}
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				double *sub_mne_data = (double *) (*sub_mne).data;
				double *new_data = (double *) calloc(1, sizeof(double));
				if (mode == MULTINOM_ADD) (*new_data) = (*mne_data) + (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*new_data) = (*mne_data) - (*sub_mne_data);
				else 
				{
					printf("Mode must correspond to either addition or subtraction\n");
					exit(EXIT_FAILURE);
				}
				free(mne_data);
				if ((*new_data) != 0)
				{
					(*mne).data = (void *) new_data;
					(*mne).dtype = DTYPE_CONST_FLOAT;
				}
				else
				{
					free(new_data);
					(*mne).data = NULL;
					(*mne).dtype = DTYPE_CONST_ZERO;
				}
			}
			if ((*sub_mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex *sub_mne_data = (complex *) (*sub_mne).data;
				complex *new_data = (complex *) calloc(1, sizeof(complex));
				if (mode == MULTINOM_ADD) (*new_data) = (*mne_data) + (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*new_data) = (*mne_data) - (*sub_mne_data);
				else 
				{
					printf("Mode must correspond to either addition or subtraction\n");
					exit(EXIT_FAILURE);
				}
				free(mne_data);
				if ((*new_data) != 0)
				{
					(*mne).data = (void *) new_data;
					(*mne).dtype = DTYPE_CONST_COMPLEX;
				}
				else
				{
					free(new_data);
					(*mne).data = NULL;
					(*mne).dtype = DTYPE_CONST_ZERO;
				}
			}
		}
		else if ((*mne).dtype == DTYPE_CONST_FLOAT)
		{
			double *mne_data = (double *) (*mne).data;
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				double *sub_mne_data = (double *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
				if ((*mne_data) != 0) {}
				else
				{
					free(mne_data);
					(*mne).data = NULL;
					(*mne).dtype = DTYPE_CONST_ZERO;
				}
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
				if ((*mne_data) != 0) {}
				else
				{
					free(mne_data);
					(*mne).data = NULL;
					(*mne).dtype = DTYPE_CONST_ZERO;
				}
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex *sub_mne_data = (complex *) (*sub_mne).data;
				complex *new_data = (complex *) calloc(1, sizeof(complex));
				(*new_data) = (*mne_data);
				if (mode == MULTINOM_ADD) (*new_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*new_data) -= (*sub_mne_data);
				free(mne_data);
				if ((*new_data) != 0) 
				{
					(*mne).data = (void *) new_data;
					(*mne).dtype = DTYPE_CONST_COMPLEX;
				}
				else
				{
					free(new_data);
					(*mne).data = NULL;
					(*mne).dtype = DTYPE_CONST_ZERO;
				}
			}
		}
		else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
		{
			complex *mne_data = (complex *) (*mne).data;
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				double *sub_mne_data = (double *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex *sub_mne_data = (complex *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
			}
			if ((*mne_data) != 0) {}
			else
			{
				free(mne_data);
				(*mne).data = NULL;
				(*mne).dtype = DTYPE_CONST_ZERO;
			}
		}
		else
		{
			printf("Something weird happened\n");
			exit(0);
		}
	}
	else if (!nonzero_test(mne) && nonzero_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, (*sub_mne).dim,(*sub_mne).deg, (*sub_mne).data, (*sub_mne).dtype);
		if (mode == MULTINOM_SUB) negative_multinomial_exp(mne);
	}
}

void add2multinomial_exp(multinomial_exp *mne, multinomial_exp *sub_mne)
{
	add_sub_multinomial_exp(mne, sub_mne, MULTINOM_ADD);
	//printf("\t\tResult: ");
	//display_multinomial_exp(mne);
}

void multinomial_exp2str(multinomial_exp *mne, array_char *mne_str)
{
	if (nonzero_test(mne)) {}
	else
	{
		add2array_char(mne_str, '0');
		return;
	}
	if ((*mne).deg > 0 && (*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		int k = (*mne).dim;
		int deg_ulim = (*mne).deg;
		while ((*mne_data).e[deg_ulim] == NULL && deg_ulim > 0)
		{
			deg_ulim -= 1;
		}
		add2array_char(mne_str, '(');
		for (int i = 0; i < deg_ulim; i++)
		{
			multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[i];
			if (aux != NULL && nonzero_test(aux))
			{
				multinomial_exp2str(aux, mne_str);
				if (i > 0)
				{
					char buf[256];
					sprintf(buf, ".x_%d^%d + ", k, i);
					append_array_char(mne_str, buf, strlen(buf));
				}
				else
				{
					char buf[256];
					sprintf(buf, " + ");
					append_array_char(mne_str, buf, strlen(buf));
				}
			}
		}
		multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[deg_ulim];
		if (aux != NULL && nonzero_test(aux))
		{
			multinomial_exp2str(aux, mne_str);
			char buf[256];
			sprintf(buf, ".x_%d^%d)", k, deg_ulim);
			append_array_char(mne_str, buf, strlen(buf));
		}
	}
	else
	{
		if ((*mne).dtype == DTYPE_INDETERMINATE_CONST)
		{
			if ((*mne).data == NULL) add2array_char(mne_str, 'C');
			else 
			{
				indeterminate c = *((indeterminate *) (*mne).data);
				char buf[256];
				sprintf(buf, "C_%d", c.sym_id);
				append_array_char(mne_str, buf, strlen(buf));
			}
		}
		else if (nonzero_test(mne))
		{
			char buf[256];
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				sprintf(buf, "%d", *((int *) (*mne).data));
			}
			else if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				sprintf(buf, "%g", *((double *) (*mne).data));
			}
			else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex val = *((complex *) (*mne).data);
				sprintf(buf, "%g + %gj", creal(val), cimag(val));
			}
			append_array_char(mne_str, buf, strlen(buf));
		}
	}
}

void display_multinomial_exp(multinomial_exp *mne)
{
	array_char buf;
	array_char_init(&buf, 1);
	multinomial_exp2str(mne, &buf);
	for (int i = 0; i < buf.len; i++)
	{
		fputc(buf.e[i], stdout);
	}
	printf("\n");
	free_array_char(&buf);
}

void multiply_by_monomial_multinomial_exp(multinomial_exp *mne, int dim, int p)
{
	if (nonzero_test(mne) && mne != NULL && p > 0) {}
	else return;
	int min_dim = (*mne).dim;
	if (dim > min_dim && (*mne).dtype == DTYPE_MULTINOMIAL)
	{
		printf("multiply_by_monomial_multinomial_exp: dim = %d > %d\n", dim, min_dim);
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		for (int i = 0; i < (*mne_data).len; i++)
		{
			multiply_by_monomial_multinomial_exp((multinomial_exp *) (*mne_data).e[i], dim, p);
		}
	}
	else if (dim == min_dim)
	{
		printf("multiply_by_monomial_multinomial_exp: equal case\n");
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		int i1 = p + (*mne_data).len;
		(*mne).deg = i1 - 1;
		int i0 = (*mne_data).len;
		add_mem_array_voidstar_until(mne_data, i1);
		(*mne_data).len = i1;
		while (i0 > 0)
		{
			i0 -= 1;
			i1 -= 1;
			(*mne_data).e[i1] = (*mne_data).e[i0];
		}
		while (i1 > 0)
		{
			i1 -= 1;
			(*mne_data).e[i1] = NULL;
		}
	}
	else
	{
		printf("Creating an auxiliary multinomial of dimension %d (c.f. dimension %d)\n", dim, (*mne).dim);
		multinomial_exp *aux = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
		(*aux) = (*mne);
		(*mne).dim = dim;
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		int pp1 = p + 1;
		array_voidstar_init(mne_data, pp1);
		(*mne_data).len = pp1;
		(*mne_data).e[p] = (void *) aux;
		(*mne).dtype = DTYPE_MULTINOMIAL;
	}
}

void multiply_multinomial_exp_data(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3)
{
	//printf("\tMultiplying const data\n");
	int mne3_dtype;
	void *mne3_data;
	if (indeterminate_const_test(mne1) || indeterminate_const_test(mne2))
	{
		mne3_dtype = DTYPE_INDETERMINATE_CONST;
		mne3_data = NULL;
	}
	else if ((*mne1).dtype == DTYPE_CONST_INT)
	{
		//printf("\tInteger case\n");
		int mne1_val = *((int *) (*mne1).data);
		//printf("\tmne1_val = %d\n", mne1_val);
		if ((*mne2).dtype == DTYPE_CONST_INT)
		{
			mne3_dtype = DTYPE_CONST_INT;
			int mne2_val = *((int *) (*mne2).data);
			//printf("\tmne2_val = %d\n", mne2_val);
			int *aux = (int *) calloc(1, sizeof(int));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
			//printf("\tmne3_val = %d\n", (*aux));
		}
		else if ((*mne2).dtype == DTYPE_CONST_FLOAT)
		{
			mne3_dtype = DTYPE_CONST_FLOAT;
			double mne2_val = *((double *) (*mne2).data);
			double *aux = (double *) calloc(1, sizeof(int));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
		}
		else if ((*mne2).dtype == DTYPE_CONST_COMPLEX)
		{
			mne3_dtype = DTYPE_CONST_COMPLEX;
			complex mne2_val = *((complex *) (*mne2).data);
			complex *aux = (complex *) calloc(1, sizeof(complex));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
		}
		else
		{
			printf("Error (multiply_multinomial_exp_data): unhandled case (mne2_type) = %d\n", (*mne2).dtype);
			exit(EXIT_FAILURE);
		}
	}
	else if ((*mne1).dtype == DTYPE_CONST_FLOAT)
	{
		mne3_dtype = DTYPE_CONST_FLOAT;
		double mne1_val = *((double *) (*mne1).data);
		if ((*mne2).dtype == DTYPE_CONST_INT)
		{
			int mne2_val = *((int *) (*mne2).data);
			double *aux = (double *) calloc(1, sizeof(double));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
		}
		else if ((*mne2).dtype == DTYPE_CONST_FLOAT)
		{
			double mne2_val = *((double *) (*mne2).data);
			double *aux = (double *) calloc(1, sizeof(double));
			(*aux) = mne2_val * mne1_val;
			mne3_data = (void *) aux;
		}
		else if ((*mne2).dtype == DTYPE_CONST_COMPLEX)
		{
			mne3_dtype = DTYPE_CONST_COMPLEX;
			complex mne2_val = *((complex *) (*mne2).data);
			complex *aux = (complex *) calloc(1, sizeof(complex));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
		}
		else printf("Something weird happened!\n");
	}
	else if ((*mne1).dtype == DTYPE_CONST_COMPLEX)
	{
		mne3_dtype = DTYPE_CONST_COMPLEX;
		complex mne1_val = *((complex *) (*mne1).data);
		if ((*mne2).dtype == DTYPE_CONST_INT)
		{
			int mne2_val = *((int *) (*mne2).data);
			complex *aux = (complex *) calloc(1, sizeof(complex));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
		}
		else if ((*mne2).dtype == DTYPE_CONST_FLOAT)
		{
			double mne2_val = *((double *) (*mne2).data);
			complex *aux = (complex *) calloc(1, sizeof(complex));
			(*aux) = mne2_val * mne1_val;
			mne3_data = (void *) aux;
		}
		else if ((*mne2).dtype == DTYPE_CONST_COMPLEX)
		{
			mne3_dtype = DTYPE_CONST_COMPLEX;
			complex mne2_val = *((complex *) (*mne2).data);
			complex *aux = (complex *) calloc(1, sizeof(complex));
			(*aux) = mne1_val * mne2_val;
			mne3_data = (void *) aux;
		}
		else printf("Something weird happened!\n");
	}
	else
	{
		printf("Error (multiply_multinomial_exp_data): unhandled case\n");
		exit(EXIT_FAILURE);
	}
	multinomial_exp_init(mne3, -1, 0, mne3_data, mne3_dtype);
	//printf("\tAttempting to free mne3_data\n");
	free(mne3_data);
	//printf("\t(done)\n");
}

void check_null(array_voidstar *mne3_incr_data)
{
	for (int i = 0; i < (*mne3_incr_data).mem; i++)
	{
		if ((*mne3_incr_data).e[i] == NULL) {}
		else
		{
			printf("Error: data array contains non-null entries\n");
			exit(EXIT_FAILURE);
		}
	}
}

void multiply_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3)
{
	//printf("multiply_multinomial_exp:\n");
	if (nonzero_test(mne1) && nonzero_test(mne2)) {}
	else
	{
		multinomial_exp_init(mne3, -1, 0, NULL, DTYPE_CONST_ZERO);
		return;
	}
	if (((*mne1).dtype == DTYPE_MULTINOMIAL) && ((*mne1).dim <= (*mne2).dim || has_constant_type(mne2)))
	{
		//printf("\tmultiplying multinomials of dimensions %d, %d\n", (*mne1).dim, (*mne2).dim);
		int mne3_deg;
		char status = ((*mne1).dim == (*mne2).dim);
		multinomial_exp_init(mne3, -1, 0, NULL, DTYPE_CONST_ZERO);
		array_voidstar *mne1_data = (array_voidstar *) (*mne1).data;
		multinomial_exp mne3_incr;
		multinomial_exp_init(&mne3_incr, (*mne1).dim, 0, NULL, DTYPE_MULTINOMIAL);
		mne3_incr.data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne3_incr_data = (array_voidstar *) mne3_incr.data;
		mne3_deg = !status ? (*mne1).deg : (*mne1).deg + (*mne2).deg;
		array_voidstar_init(mne3_incr_data, mne3_deg + 1);
		(*mne3_incr_data).len = mne3_deg + 1;
		mne3_incr.deg = mne3_deg;
		mne3_incr.dim = (*mne1).dim;
		for (int i = 0; i <= (*mne1).deg; i++)
		{
			if ((*mne1_data).e[i] != NULL)
			{
				// RESUME: Fix this! Address the case where (*mne1).dim == (*mne2).dim,
				// 	and make sure that the degree of mne3_incr is adjusted accordingly.
				// UPDATE: A first attempt has been made, but is yet to be tested.
				multinomial_exp aux;
				multinomial_exp *mne1_desc = (multinomial_exp *) (*mne1_data).e[i];
				multiply_multinomial_exp(mne1_desc, mne2, &aux);
				//display_multinomial_exp(&aux);
				if (!status || has_constant_type(&aux))
				{
					if (nonzero_test(&aux))
					{
						//mne3_incr.deg = i;
						//(*mne3_incr_data).len = i + 1;
						(*mne3_incr_data).e[i] = (void *) &aux;
						mne3_incr.deg = i;
						(*mne3_incr_data).len = i + 1;
						//mne3_incr.deg = i;
						add2multinomial_exp(mne3, &mne3_incr);
						//display_multinomial_exp(mne3);
						(*mne3_incr_data).e[i] = NULL;
						check_null(mne3_incr_data);
					}
				}
				else
				{
					// shift aux_data over by i steps
					if (nonzero_test(&aux))
					{
						array_voidstar *aux_data = (array_voidstar *) aux.data;
						void **incr_data = &((*mne3_incr_data).e[i]);
						int ii_max = -1;
						for (int ii = 0; ii <= aux.deg; ii++)
						{
							multinomial_exp *aux_di = (multinomial_exp *) (*aux_data).e[ii];
							if (nonzero_test(aux_di)) ii_max = ii;
							incr_data[ii] = (*aux_data).e[ii];
						}
						if (ii_max > -1)
						{
							mne3_incr.deg = ii_max + i;
							(*mne3_incr_data).len = mne3_incr.deg + 1;
							add2multinomial_exp(mne3, &mne3_incr);
						}
						//(*mne3_incr_data).len = (*aux_data).len + i;
						//mne3_incr.deg = aux.deg + i;
						for (int ii = 0; ii <= aux.deg; ii++) incr_data[ii] = NULL;
						check_null(mne3_incr_data);
					}
				}
				free_multinomial_exp(&aux);
			}
		}
		free_array_voidstar(mne3_incr_data, NULL);
		free(mne3_incr.data);
		//free_multinomial_exp(&mne3_incr);
	}
	else if ((*mne2).dtype == DTYPE_MULTINOMIAL && ((*mne2).dim < (*mne1).dim || has_constant_type(mne1)))
	{
		multiply_multinomial_exp(mne2, mne1, mne3);
	}
	else if (has_constant_type(mne1) && has_constant_type(mne2))
	{
		multiply_multinomial_exp_data(mne1, mne2, mne3);
	}
	else
	{
		printf("Something weird happened!\n"); 
	}
	if (check_types(mne3)) {}
	else
	{
		printf("Failed type check in multiply_multinomial_exp\n");
		exit(EXIT_FAILURE);
	}
	//printf("(done) (multiply_multinomial_exp)\n");
}

void pow_mutinomial_exp(multinomial_exp *mne, multinomial_exp *mne_p, int p)
{
	if ((*mne).dtype == DTYPE_MULTINOMIAL)
	{
		transcribe_multinomial_exp(mne, mne_p);
		multinomial_exp *aux = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
		for (int i = 0; i < p; i++)
		{
			multiply_multinomial_exp(mne, mne_p, aux);
			multinomial_exp *temp = mne_p;
			mne_p = aux;
			aux = temp;
			free_multinomial_exp((void *) aux);
		}
		free(aux);
	}
	else
	{
		if ((*mne).dtype == DTYPE_INDETERMINATE_CONST)
		{
			multinomial_exp_init(mne_p, -1, 0, NULL, DTYPE_INDETERMINATE_CONST);
		}
		else
		{
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				int mne_val = *((int *) (*mne).data);
				int mne_p_val = 1;
				if (p > 0 && p < 11)
				{
					mne_p_val = mne_val;
					// Consider replacing this with an 'optimized' evaluator for large 'p'
					for (int i = 1; i < p; i++)
					{
						mne_p_val *= mne_val;
					}
				}
				else if (p > 10)
				{
					double pre = pow((double) mne_val, p);
					mne_p_val = pre > 0 ? (int) (pre + 0.01) : (int) (pre - 0.01);
				}
				multinomial_exp_init(mne_p, -1, 0, &mne_p_val, DTYPE_CONST_INT);
			}
			if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				double mne_val = *((double *) (*mne).data);
				double mne_p_val = 1;
				if (p > 0 && p < 11)
				{
					mne_p_val = mne_val;
					for (int i = 1; i < p; i++) mne_p_val *= mne_val;
				}
				else if (p > 10)
				{
					mne_p_val = pow(mne_val, p);
				}
				multinomial_exp_init(mne_p, -1, 0, &mne_p_val, DTYPE_CONST_FLOAT);
			}
			if ((*mne).dtype == DTYPE_CONST_COMPLEX)
			{
				complex mne_val = *((complex *) (*mne).data);
				complex mne_p_val = 1;
				if (p > 0 && p < 11)
				{
					mne_p_val = mne_val;
					for (int i = 1; i < p; i++) mne_p_val *= mne_val;
				}
				else if (p > 10)
				{
					mne_p_val = pow(mne_val, p);
				}
				multinomial_exp_init(mne_p, -1, 0, &mne_p_val, DTYPE_CONST_COMPLEX);
			}
		}
	}
}

// RESUME: check/fix this!
void compose_multinomial_exp(multinomial_exp **mne1, int dim1, multinomial_exp **mne2, int dim2, multinomial_exp **mne3)
{
	for (int i = 0; i < dim1; i++)
	{
		if (nonzero_test(mne1[i])) {}
		else
		{
			multinomial_exp_init(mne3[i], -1, 0, NULL, DTYPE_CONST_ZERO);
			continue;
		}
		if ((*(mne1[i])).dtype == DTYPE_MULTINOMIAL)
		{
			multinomial_exp_init(mne3[i], -1, 0, NULL, DTYPE_CONST_ZERO);
			multinomial_exp pow_;
			int val = 1;
			multinomial_exp_init(&pow_, -1, 0, &val, DTYPE_CONST_INT);
			multinomial_exp aux_;
			array_voidstar *mne1i_data = (array_voidstar *) (*(mne1[i])).data;
			for (int ii = 0; ii <= (*(mne1[i])).deg; ii++)
			{
				multinomial_exp term_;
				multinomial_exp eval_;
				multinomial_exp *mne1_ii = (multinomial_exp *) (*mne1i_data).e[ii];
				multinomial_exp *eval_addr = &eval_;
				compose_multinomial_exp(&mne1_ii, 1, mne2, dim2, &eval_addr);
				multiply_multinomial_exp(&pow_, &eval_, &term_);
				add2multinomial_exp(mne3[i], &term_);
				multiply_multinomial_exp(mne2[i], &pow_, &aux_);
				multinomial_exp temp_ = pow_;
				pow_ = aux_;
				aux_ = temp_;
			}
		}
		else
		{
			multinomial_exp_init(mne3[i], -1, 0, (*(mne1[i])).data, (*(mne1[i])).dtype);
		}
	}
}

// RESUME: try to simplify this (or ensure that the DTYPE_MULTINOMIAL case can be safely removed)
// 	UPDATE (2/27/2024): consider changing addition to keep track of cancellations (which would
// 	probably be the source of most if not all avoidable/unnecessary complexity.)
char nonzero_test(multinomial_exp *mne)
{
	if (mne == NULL) return 0;
	if (!((*mne).dtype == DTYPE_CONST_ZERO || 
			(*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data == NULL))
	{
		if ((*mne).dtype == DTYPE_MULTINOMIAL)
		{
			array_voidstar *mne_data = (array_voidstar *) (*mne).data;
			if ((*mne_data).len > 0)
			{
				for (int i = 0; i < (*mne_data).len; i++)
				{
					if ((*mne_data).e[i] != NULL) return 1;
				}
			}
			return 0;
		}
		if (((*mne).dtype & DTYPE_CONST) > 0)
		{
			switch ((*mne).dtype)
			{
				case DTYPE_CONST_INT:
					int *ival = (int *) (*mne).data;
					return (*ival) != 0;
				case DTYPE_CONST_FLOAT:
					double *fval = (double *) (*mne).data;
					return (*fval) != 0;
				case DTYPE_CONST_COMPLEX:
					complex *cval = (complex *) (*mne).data;
					return (*cval) != 0;
			}
		}
		return 1;
	}
	else return 0;
}

void polynomial_int2multinomial_exp(int *p, int p_len, multinomial_exp *mne, int dim)
{
	if (p_len > 0)
	{
		if (p_len > 1) {}
		else
		{
			if (p[0] != 0) multinomial_exp_init_int(mne, p[0]);
			else multinomial_exp_init_zero(mne);
			return;
		}
		multinomial_exp_init(mne, dim, p_len - 1, NULL, DTYPE_MULTINOMIAL);
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		array_voidstar_init(mne_data, p_len);
		for (int i = 0; i < p_len; i++)
		{
			if (p[i] != 0)
			{
				multinomial_exp *desc = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init_int(desc, p[i]);
				add2array_voidstar(mne_data, (void *) desc);
			}
			else add2array_voidstar(mne_data, NULL);
		}
	}
	else multinomial_exp_init_zero(mne);
}

void polynomial_double2multinomial_exp(double *p, int p_len, multinomial_exp *mne, int dim)
{
	if (p_len > 0)
	{
		if (p_len > 1) {}
		else
		{
			if (p[0] != 0) multinomial_exp_init_double(mne, p[0]);
			else multinomial_exp_init_zero(mne);
			return;
		}
		multinomial_exp_init(mne, dim, p_len - 1, NULL, DTYPE_MULTINOMIAL);
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		array_voidstar_init(mne_data, p_len);
		for (int i = 0; i < p_len; i++)
		{
			if (p[i] != 0.0)
			{
				multinomial_exp *desc = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init(desc, -1, 0, &(p[i]), DTYPE_CONST_FLOAT);
				add2array_voidstar(mne_data, (void *) desc);
			}
			else add2array_voidstar(mne_data, NULL);
		}
	}
	else
	{
		multinomial_exp_init(mne, -1, 0, NULL, DTYPE_CONST_ZERO);
	}
}

void polynomial_complex2multinomial_exp(complex *p, int p_len, multinomial_exp *mne, int dim)
{
	if (p_len > 0)
	{
		if (p_len > 1) {}
		else
		{
			if (p[0] != 0) multinomial_exp_init_complex(mne, p[0]);
			else multinomial_exp_init_zero(mne);
			return;
		}
		multinomial_exp_init(mne, dim, p_len - 1, NULL, DTYPE_MULTINOMIAL);
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		array_voidstar_init(mne_data, p_len);
		for (int i = 0; i < p_len; i++)
		{
			if (p[i] != 0.0)
			{
				multinomial_exp *desc = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init(desc, -1, 0, &(p[i]), DTYPE_CONST_COMPLEX);
				add2array_voidstar(mne_data, (void *) desc);
			}
			else add2array_voidstar(mne_data, NULL);
		}
	}
	else
	{
		multinomial_exp_init(mne, -1, 0, NULL, DTYPE_CONST_ZERO);
	}
}

double multinomial_exp_eval_double(multinomial_exp *mne, double *vals, int len)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		double rval = multinomial_exp_eval_double((multinomial_exp *) (*mne_data).e[0], vals, len);
		if ((*mne).dim < len)
		{
			double pow_ = vals[(*mne).dim];
			for (int i = 1; i <= (*mne).deg; i++)
			{
				rval += multinomial_exp_eval_double((multinomial_exp *) (*mne_data).e[i], vals, len) * pow_;
				pow_ *= vals[(*mne).dim];
			}
		}
		return rval;
	}
	else
	{
		if (nonzero_test(mne))
		{
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				return (double) *((int *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				return *((double *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
			{
				printf("Error: double evaluation of complex multinomials not yet supported\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			return 0.0;
		}
	}

}

int multinomial_exp_eval_int(multinomial_exp *mne, int *vals, int len)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL && (*mne).dim > 0)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		int rval = multinomial_exp_eval_int((multinomial_exp *) (*mne_data).e[0], vals, len);
		if ((*mne).dim < len)
		{
			int pow_ = vals[(*mne).dim];
			for (int i = 1; i <= (*mne).deg; i++)
			{
				rval += multinomial_exp_eval_int((multinomial_exp *) (*mne_data).e[i], vals, len);
				pow_ *= vals[(*mne).dim];
			}
		}
		return rval;
	}
	else
	{
		if (nonzero_test(mne))
		{
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				return *((int *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				return (int) *((double *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
			{
				printf("Error: real evaluation of complex multinomials not yet supported\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			return 0;
		}
	}
}

complex multinomial_exp_eval_as_complex(multinomial_exp *mne, void *vals, int len, char dtype)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL && (*mne).dim > 0)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		complex rval = multinomial_exp_eval_as_complex((multinomial_exp *) (*mne_data).e[0], vals, len, dtype);
		if ((*mne).dim < len)
		{
			complex v0;
			switch (dtype)
			{
				case DTYPE_FLOAT:
					v0 = (complex) (((double *) vals)[(*mne).dim]);
					break;
				case DTYPE_INT:
					v0 = (complex) (((int *) vals)[(*mne).dim]);
					break;
				case DTYPE_COMPLEX:
					v0 = ((complex *) vals)[(*mne).dim];
			}
			complex pow_ = v0;
			for (int i = 1; i <= (*mne).deg; i++)
			{
				rval += multinomial_exp_eval_as_complex((multinomial_exp *) (*mne_data).e[i], vals, len, dtype) * pow_;
				pow_ *= v0;
			}
		}
		return rval;
	}
	else
	{
		if (nonzero_test(mne))
		{
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				return (complex) *((int *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				return (complex) *((double *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
			{
				return *((complex *) (*mne).data);
			}
		}
		else
		{
			return 0.0;
		}
	}
}

complex multinomial_exp_eval_complex(multinomial_exp *mne, complex *vals, int len)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		complex rval = multinomial_exp_eval_complex((multinomial_exp *) (*mne_data).e[0], vals, len);
		if ((*mne).dim < len)
		{
			complex pow_ = vals[(*mne).dim];
			for (int i = 1; i <= (*mne).deg; i++)
			{
				rval += multinomial_exp_eval_complex((multinomial_exp *) (*mne_data).e[i], vals, len) * pow_;
				pow_ *= vals[(*mne).dim];
			}
		}
		return rval;
	}
	else
	{
		if (nonzero_test(mne))
		{
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				return (complex) *((int *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				return (complex) *((double *) (*mne).data);
			}
			else if ((*mne).dtype == DTYPE_CONST_COMPLEX)
			{
				return *((complex *) (*mne).data);
			}
		}
		else
		{
			return 0.0;
		}
	}
}

// RESUME: Consider changing the 'type' encodings to simplify checks of this sort:
// 	for example, instead of simply assigning a different integer to each type,
// 	types could occupy different bits in the 'data' field (so that CONST types
// 	can be easily distinguished from ARRAY types with relatively simple bit 
// 	operations.)
char has_constant_type(multinomial_exp *mne)
{
	return ((*mne).dtype & DTYPE_CONST) == DTYPE_CONST;
}

complex parse_as_complex(multinomial_exp *mne)
{
	switch ((*mne).dtype)
	{
		case DTYPE_CONST_INT:
			return (complex) *((int *) (*mne).data);
			break;
		case DTYPE_CONST_FLOAT:
			return (complex) *((double *) (*mne).data);
			break;
		case DTYPE_CONST_COMPLEX:
			return *((complex *) (*mne).data);
			break;
	}
}

char multinomial_exp_equals(multinomial_exp *mne1, multinomial_exp *mne2)
{
	if (mne1 != NULL && mne2 != NULL) {}
	else
	{
		return mne1 == mne2;
	}
	char nonzero_test1 = nonzero_test(mne1);
	char nonzero_test2 = nonzero_test(mne2);
	if (nonzero_test1 && nonzero_test2) {}
	else return nonzero_test1 == nonzero_test2;
	if ((*mne1).dtype == (*mne2).dtype)
	{
		switch ((*mne1).dtype)
		{
			case (DTYPE_MULTINOMIAL):
			{
				if ((*mne1).dim == (*mne2).dim) {}
				else return 0;
				if ((*mne1).data != NULL && (*mne2).data != NULL)
				{
					array_voidstar *mne1_data = (array_voidstar *) (*mne1).data;
					array_voidstar *mne2_data = (array_voidstar *) (*mne2).data;
					if ((*mne1_data).len == (*mne2_data).len)
					{
						for (int i = 0; i < (*mne1_data).len; i++)
						{
							if (multinomial_exp_equals((multinomial_exp *) (*mne1_data).e[i], (multinomial_exp *) (*mne2_data).e[i])) {}
							else return 0;
						}
						return 1;
					}
					else return 0;
				}
				else
				{
					return (*mne1).data == (*mne2).data;
				}
				break;
			}
			case (DTYPE_CONST_INT):
			{
				int *val1 = (int *) (*mne1).data;
				int *val2 = (int *) (*mne2).data;
				return (*val1) == (*val2);
				break;
			}
			case (DTYPE_CONST_FLOAT):
			{
				double *val1 = (double *) (*mne1).data;
				double *val2 = (double *) (*mne2).data;
				return (*val1) == (*val2);
			}
			case (DTYPE_CONST_COMPLEX):
			{
				complex *val1 = (complex *) (*mne1).data;
				complex *val2 = (complex *) (*mne2).data;
				return (*val1) == (*val2);
			}
		}
	}
	else
	{
		if (has_constant_type(mne1) && has_constant_type(mne2))
		{
			complex val1 = parse_as_complex(mne1);
			complex val2 = parse_as_complex(mne2);
			return val1 == val2;
		}
		else return 0;
	}
}

int count_indeterminates_rec(multinomial_exp *mne, array_int *sym_list, array_int *sym_, int *n_unspec)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		for (int i = 0; i < (*mne_data).len; i++)
		{
			multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[i];
			count_indeterminates_rec(mne, sym_list, sym_, n_unspec);
		}
	}
	if ((*mne).dtype == DTYPE_INDETERMINATE_CONST)
	{
	       if ((*mne).data != NULL)
	       {
			indeterminate *mne_data = (indeterminate *) (*mne).data;
			if ((*mne_data).sym_id < (*sym_).len) {}
			else
			{
				add_mem_array_int_until(sym_, (*mne_data).sym_id);
				int init_len = (*sym_).len;
				(*sym_).len = (*sym_).mem;
				for (int i = init_len; i < (*sym_).len; i++) (*sym_).e[i] = -1;
			}
			if ((*sym_).e[(*mne_data).sym_id] > -1) {}
			else
			{
				(*sym_).e[(*mne_data).sym_id] = (*sym_list).len;
				add2array_int(sym_list, (*mne_data).sym_id);
			}
	       }
	       else (*n_unspec) += 1;
	}
}

int count_indeterminates(multinomial_exp *mne)
{
	array_int sym_list;
	array_int_init(&sym_list, 1);
	array_int sym_map;
	array_int_init(&sym_map, 256);
	for (int i  = 0; i < 256; i++)
	{
		add2array_int(&sym_map, -1);
	}
	int n_unspec = 0;
	count_indeterminates_rec(mne, &sym_list, &sym_map, &n_unspec);
	int total = n_unspec + sym_list.len;
	free_array_int(&sym_list);
	free_array_int(&sym_map);
	return total;
}

void truncate_multinomial_exp_rec(multinomial_exp *mne, multinomial_exp *tmne, int dim, array_voidstar *tails, array_int *sym_list)
{
	if ((*mne).dim < dim)
	{
		int tmne_deg = (*mne).deg;
		multinomial_exp_init(tmne, (*mne).dim, (*mne).deg, NULL, DTYPE_MULTINOMIAL);
		if ((*mne).dtype == DTYPE_MULTINOMIAL) {}
		else
		{
			printf("Something weird happened!\n");
			exit(EXIT_FAILURE);
		}
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		if (mne_data != NULL)
		{
			(*tmne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
			array_voidstar *tmne_data = (*tmne).data;
			array_voidstar_init(tmne_data, tmne_deg + 1);
			for (int i = 0; i < (*mne_data).len; i++)
			{
				if ((*mne_data).e[i] != NULL)
				{
					multinomial_exp *mne_i = (multinomial_exp *) (*mne_data).e[i];
					multinomial_exp *tmne_i = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
					truncate_multinomial_exp_rec(mne_i, tmne_i, dim, tails, sym_list);
					add2array_voidstar(tmne_data, (void *) tmne_i);
				}
			}
		}
	}
	else if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).deg > 0)
	{
		// Compare (*mne) with other polynomials in (*tails) (RESUME: make sure that comparison works)
		int sym_id = -1;
		for (int i = 0; i < (*tails).len; i++)
		{
			multinomial_exp *tail_i = (multinomial_exp *) (*tails).e[i];
			if (multinomial_exp_equals(tail_i, mne))
			{
				sym_id = i;
				break;
			}
		}
		indeterminate c;
		if (sym_id > -1)
		{
			c.sym_id = sym_id;
			multinomial_exp_init(tmne, -1, 0, &c, DTYPE_INDETERMINATE_CONST);
		}
		else
		{
			c.sym_id = (*tails).len;
			add2array_int(sym_list, (*tails).len);
			multinomial_exp *tail_exp = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
			transcribe_multinomial_exp(mne, tail_exp);
			add2array_voidstar(tails, (void *) tail_exp);
			multinomial_exp_init(tmne, -1, 0, &c, DTYPE_INDETERMINATE_CONST);
		}
	}
	else
	{
		if (nonzero_test(mne))
		{
			transcribe_multinomial_exp(mne, tmne);
		}
		else
		{
			multinomial_exp_init(tmne, -1, 0, NULL, DTYPE_CONST_ZERO);
		}
	}
}

// Generate a 'truncated' multinomial of lower (effective) dimensionality, replacing 
// 	trailing polynomials in low dimensions with one 'indeterminate constant' per
// 	distinct independent term (i.e. if a multinomial is truncated at dimension m,
// 	then all instances of p_i(x_m, x_{m+1}, ...) can be replaced by the same indeterminate constant,
// 	e.g. P_i; the idea is to compute the 'general form' or equivalence class of a 
// 	polynomial that would arise from treating the last few dimensions as [possibly 
// 	complex] indeterminate constants; this may be useful for evaluating/comparing 
// 	different multinomial evaluators.)
void truncate_multinomial_exp(multinomial_exp *mne, multinomial_exp *tmne, int dim)
{
	array_voidstar tails;
	array_voidstar_init(&tails, 1);
	array_int sym_list;
	array_int_init(&sym_list, 1);
	truncate_multinomial_exp_rec(mne, tmne, dim, &tails, &sym_list);
	free_array_voidstar(&tails, free_multinomial_exp);
	free_array_int(&sym_list);
}

void indetermify_multinomial_exp_rec(multinomial_exp *mne, multinomial_exp *imne, int *sym_id)
{
	if (mne != NULL)
	{
		if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
		{
			array_voidstar *mne_data = (array_voidstar *) (*mne).data;
			multinomial_exp_init(imne, (*mne).dim, (*mne).deg, NULL, DTYPE_MULTINOMIAL);
			(*imne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
			array_voidstar *imne_data = (array_voidstar *) (*imne).data;
			array_voidstar_init(imne_data, (*mne_data).len);
			for (int i = 0; i < (*mne_data).len; i++)
			{
				multinomial_exp *mne_i = (multinomial_exp *) (*mne_data).e[i];
				if (mne_i != NULL)
				{
					multinomial_exp *imne_i = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
					indetermify_multinomial_exp_rec(mne_i, imne_i, sym_id);
					add2array_voidstar(imne_data, imne_i);
				}
			}
		}
		else
		{
			if (nonzero_test(mne))
			{
				indeterminate c;
				c.sym_id = (*sym_id);
				multinomial_exp_init(imne, -1, 0, &c, DTYPE_INDETERMINATE_CONST);
				(*sym_id += 1);
			}
		}
	}
}

// Convert a multinomial into an indeterminate version
void indetermify_multinomial_exp(multinomial_exp *mne, multinomial_exp *imne)
{
	int sym_id = 0;
	indetermify_multinomial_exp_rec(mne, imne, &sym_id);
}

int multinomial_exp_complexity(multinomial_exp *mne)
{
	if (mne != NULL) {}
	else return 0;
	int c = 0;
	if ((*mne).dtype == DTYPE_MULTINOMIAL && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		for (int i = 0; i < (*mne_data).len; i++)
		{
			c += 1 + multinomial_exp_complexity((multinomial_exp *) (*mne_data).e[i]);
		}
	}
	return c;
}

void define_root_multinomial_evaluator(multinomial_evaluator *mnev, multinomial_exp *expr, int *perm, int perm_len)
{
	(*mnev).perm = perm;
	(*mnev).perm_len = perm_len;
	(*mnev).expr = expr;
	(*mnev).expr_1 = NULL;
	(*mnev).expr_2 = NULL;
	(*mnev).type = M_EV_ROOT_TYPE;
}

void free_multinomial_evaluator(multinomial_evaluator *mnev)
{
	if ((*mnev).type == M_EV_COMP_TYPE && (*mnev).expr != NULL)
	{
		free_multinomial_exp((*mnev).expr);
		free_multinomial_evaluator((multinomial_evaluator *) (*mnev).expr_1);
		free_multinomial_evaluator((multinomial_evaluator *) (*mnev).expr_2);
	}
}

void multinomial_evaluator_combine(multinomial_evaluator *e1, multinomial_evaluator *e2, char op, multinomial_evaluator *comp)
{
	(*comp).expr = NULL;
	(*comp).expr_1 = (void *) e1;
	(*comp).expr_2 = (void *) e2;
	(*comp).op = op;
	(*comp).type = M_EV_COMP_TYPE;
	(*comp).perm = NULL;
}

void multinomial_evaluator_compute_expr(multinomial_evaluator *mev)
{
	if ((*mev).type == M_EV_COMP_TYPE)
	{
		if ((*mev).expr == NULL) {}
		else 
		{
			printf("Error in multinomial_evaluator_compute_expr! (*mev).expr is not NULL (free (*mev).expr first)\n");
			exit(EXIT_FAILURE);
		}
		(*mev).expr = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mev).expr_2;
		multinomial_evaluator_compute_expr(expr_1);
		multinomial_evaluator_compute_expr(expr_2);
		switch ((*mev).op)
		{
			case MULTINOM_MULT:
				multiply_multinomial_exp((*expr_1).expr, (*expr_2).expr, (*mev).expr);
				break;
			case MULTINOM_ADD:
				add_multinomial_exp((*expr_1).expr, (*expr_2).expr, (*mev).expr);
				break;
			case MULTINOM_SUB:
				sub_multinomial_exp((*expr_1).expr, (*expr_2).expr, (*mev).expr);
				break;
		}
	}
}

void multinomial_evaluator_combine_exp(multinomial_evaluator *e1, multinomial_evaluator *e2, char op, multinomial_evaluator *comp)
{
	(*comp).expr = NULL;
	(*comp).expr_1 = (void *) e1;
	(*comp).expr_2 = (void *) e2;
	(*comp).op = op;
	(*comp).type = M_EV_COMP_TYPE;
	multinomial_evaluator_compute_expr(comp);
}

void multinomial_evaluator_eval_as_complex(multinomial_evaluator *mnev, void *vals, complex *cval, char dtype)
{
	if ((*mnev).type == M_EV_COMP_TYPE)
	{
		complex val1, val2;
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		multinomial_evaluator_eval_as_complex(expr_1, vals, &val1, dtype);
		multinomial_evaluator_eval_as_complex(expr_2, vals, &val2, dtype);
		if ((*mnev).op == MULTINOM_ADD)
		{
			(*cval) = val1 + val2;
		}
		else if ((*mnev).op == MULTINOM_SUB)
		{
			(*cval) = val1 - val2;
		}
		else if ((*mnev).op == MULTINOM_MULT)
		{
			(*cval) = val1 * val2;
		}
	}
	else if ((*mnev).type == M_EV_ROOT_TYPE)
	{
		complex evals[(*mnev).perm_len];
		double *cvals = (double *) vals;
		for (int i = 0; i < (*mnev).perm_len; i++)
		{
			evals[(*mnev).perm[i]] = cvals[i];
		}
		(*cval) = multinomial_exp_eval_as_complex((*mnev).expr, &evals[0], (*mnev).perm_len, dtype);
	}
	else
	{
		printf("Error: unable to evaluate multinomial_evaluator of type %d\n", (*mnev).type);
		exit(EXIT_FAILURE);
	}
}

void multinomial_evaluator_eval_complex(multinomial_evaluator *mnev, complex *vals, complex *cval)
{
	if ((*mnev).type == M_EV_COMP_TYPE)
	{
		complex val1, val2;
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		multinomial_evaluator_eval_complex(expr_1, vals, &val1);
		multinomial_evaluator_eval_complex(expr_2, vals, &val2);
		complex pre;
		if ((*mnev).op == MULTINOM_ADD)
		{
			(*cval) = val1 + val2;
		}
		else if ((*mnev).op == MULTINOM_SUB)
		{
			(*cval) = val1 - val2;
		}
		else if ((*mnev).op == MULTINOM_MULT)
		{
			(*cval) = val1 * val2;
		}
	}
	else if ((*mnev).type == M_EV_ROOT_TYPE)
	{
		complex evals[(*mnev).perm_len];
		for (int i = 0; i < (*mnev).perm_len; i++)
		{
			evals[(*mnev).perm[i]] = vals[i];
		}
		(*cval) = multinomial_exp_eval_complex((*mnev).expr, &(evals[0]), (*mnev).perm_len);
	}
	else
	{
		printf("Error: unable to evaluate multinomial_evaluator of type %d\n", (*mnev).type);
		exit(EXIT_FAILURE);
	}
}

void multinomial_evaluator_eval_double(multinomial_evaluator *mnev, double *vals, double *cval)
{
	if ((*mnev).type == M_EV_COMP_TYPE)
	{
		complex val1, val2;
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		multinomial_evaluator_eval_as_complex(expr_1, vals, &val1, DTYPE_FLOAT);
		multinomial_evaluator_eval_as_complex(expr_2, vals, &val2, DTYPE_FLOAT);
		complex pre;
		if ((*mnev).op == MULTINOM_ADD)
		{
			pre = val1 + val2;
		}
		else if ((*mnev).op == MULTINOM_SUB)
		{
			pre = val1 - val2;
		}
		else if ((*mnev).op == MULTINOM_MULT)
		{
			pre = val1 * val2;
		}
		(*cval) = creal(pre);
	}
	else if ((*mnev).type == M_EV_ROOT_TYPE)
	{
		double evals[(*mnev).perm_len];
		for (int i = 0; i < (*mnev).perm_len; i++)
		{
			evals[(*mnev).perm[i]] = vals[i];
		}
		(*cval) = multinomial_exp_eval_double((*mnev).expr, &evals[0], (*mnev).perm_len);
	}
	else
	{
		printf("Error: unable to evaluate multinomial_evaluator of type %d\n", (*mnev).type);
		exit(EXIT_FAILURE);
	}
}

void multinomial_evaluator_eval_int(multinomial_evaluator *mnev, int *vals, int *cval)
{
	if ((*mnev).type == M_EV_COMP_TYPE)
	{
		complex val1, val2;
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		multinomial_evaluator_eval_as_complex(expr_1, vals, &val1, DTYPE_INT);
		multinomial_evaluator_eval_as_complex(expr_2, vals, &val2, DTYPE_INT);
		if ((*mnev).op == MULTINOM_ADD)
		{
			val1 += val2;
		}
		else if ((*mnev).op == MULTINOM_SUB)
		{
			val1 -= val2;
		}
		else if ((*mnev).op == MULTINOM_MULT)
		{
			val1 *= val2;
		}
		double rpart = creal(val1);
		(*cval) = rpart > 0 ? (int) (rpart + 0.01) : (int) (rpart - 0.01);
	}
	else if ((*mnev).type == M_EV_ROOT_TYPE)
	{
		int evals[(*mnev).perm_len];
		for (int i = 0; i < (*mnev).perm_len; i++)
		{
			evals[(*mnev).perm[i]] = vals[i];
		}
		(*cval) = multinomial_exp_eval_int((*mnev).expr, &evals[0], (*mnev).perm_len);
	}
	else
	{
		printf("Error: unable to evaluate multinomial_evaluator of type %d\n", (*mnev).type);
		exit(EXIT_FAILURE);
	}
}

int multinomial_evaluator_complexity(multinomial_evaluator *mnev)
{
	if ((*mnev).type == M_EV_ROOT_TYPE)
	{
		return multinomial_exp_complexity((*mnev).expr);
	}
	else if ((*mnev).type == M_EV_COMP_TYPE)
	{
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		return multinomial_evaluator_complexity(expr_1) + multinomial_evaluator_complexity(expr_2) + 1;
	}
	else
	{
		printf("Error: unable to compute complexity for multinomial_evaluator of type %d (not %d or %d) \n", (*mnev).type, M_EV_ROOT_TYPE, M_EV_COMP_TYPE);
		exit(EXIT_FAILURE);
	}
}


