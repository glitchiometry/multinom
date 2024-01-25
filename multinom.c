#include "multinom.h"

void multinomial_exp_init(multinomial_exp *mne, int dim, int deg, void *data, int dtype)
{
	(*mne).dtype = dtype;
	(*mne).deg = deg;
	(*mne).dim = dim;
	if ((dtype == DTYPE_ARRAY_VOIDSTAR) && (data != NULL))
	{
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *in_data = (array_voidstar *) data;
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		if ((*in_data).len != deg + 1)
		{
			printf("Warning (multinomial_exp_init): Mismatch between degree + 1 = %d and data length = %d\n", deg + 1, (*in_data).len); 
		}
		array_voidstar_init(mne_data, (*in_data).len);
		(*mne_data).len = (*in_data).len;
		for (int i = 0; i < (*in_data).len; i++)
		{
			if ((*in_data).e[i] != NULL) 
			{
				(*mne_data).e[i] = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				transcribe_multinomial_exp((multinomial_exp *) (*in_data).e[i], (multinomial_exp *) (*mne_data).e[i]);
			}
			else
			{
				(*mne_data).e[i] = NULL;
			}
		}
		for (int i = (*in_data).len; i <= deg; i++)
		{
			add2array_voidstar(mne_data, NULL);
		}
	}
	else 
	{
		if (dtype == DTYPE_ARRAY_VOIDSTAR)
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
	}
}

void free_multinomial_exp(void *mne)
{
	if (mne != NULL)
	{
		multinomial_exp *aux = (multinomial_exp *) mne;
		if (nonzero_test(aux))
		{
			if ((*aux).dtype == DTYPE_ARRAY_VOIDSTAR && (*aux).data != NULL)
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
	return ((*mne).dtype == DTYPE_CONST_INT || (*mne).dtype == DTYPE_CONST_FLOAT) && (*mne).data != NULL;
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
		multinomial_exp_init(mne, 0, 0, NULL, DTYPE_INDETERMINATE_CONST);
		return;
	}
	else if (valid_numeric_test(mne) && valid_numeric_test(sub_mne))
	{
		if ((*mne).dtype == DTYPE_CONST_INT)
		{
			if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *mne_data = (int *) (*mne).data;
				int *sub_mne_data = (int *) (*sub_mne).data;
				(*mne_data) += (*sub_mne_data);
			}
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				int *mne_data = (int *) (*mne).data;
				double *sub_mne_data = (double *) (*sub_mne).data;
				double *new_data = (double *) calloc(1, sizeof(double));
				(*new_data) = (*sub_mne_data) + (*mne_data);
				free(mne_data);
				(*mne).data = (void *) new_data;
				(*mne).dtype == DTYPE_CONST_FLOAT;
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
		}
	}
	else if (!nonzero_test(mne) && nonzero_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, 0, 0, (*sub_mne).data, (*sub_mne).dtype);
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
	if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL)
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
	// RESUME: address indeterminate constant case to allow symbolic arithmetic (although this may be easier to accomplish simply by expanding the space of variables)
}

void add_sub_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, char mode)
{
	if (nonzero_test(mne2)) {}
	else return;
	if ((*mne1).dim > (*mne2).dim)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne1).data;
		if ((*mne_data).len > 0) {}
		else add2array_voidstar(mne_data, NULL);
		if ((*mne_data).e[0] != NULL)
		{
			multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[0];
			//add2multinomial_exp(aux, sub_mne);
			add_sub_multinomial_exp(aux, mne2, mode); // RESUME define this!
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
	else if ((*mne1).dim < (*mne2).dim)
	{
		multinomial_exp prec;
	       	transcribe_multinomial_exp(mne2, &prec);
		if (mode == MULTINOM_SUB) negative_multinomial_exp(&prec);
		add_sub_multinomial_exp(&prec, mne1, mode);
		free_multinomial_exp(mne1);
		(*mne1) = prec;
	}
	else
	{
		if ((*mne1).dim > 0)
		{
			//printf("Combining multinomials of equal dimension of degree %d, %d\n", (*mne).deg, (*sub_mne).deg);
			int min_deg;
			array_voidstar *mne_data = (array_voidstar *) (*mne1).data;
			array_voidstar *sub_mne_data = (array_voidstar *) (*mne2).data;
			if ((*mne1).deg < (*mne2).deg)
			{
				min_deg = (*mne1).deg;
				add_mem_array_voidstar_until(mne_data, (*mne2).deg + 1);
				(*mne_data).len = (*mne_data).mem;
				for (int i = (*mne1).deg + 1; i <= (*mne2).deg; i++)
				{
					if (i < (*sub_mne_data).len && (*sub_mne_data).e[i] != NULL) 
					{
						(*mne_data).e[i] = (void *) ((multinomial_exp *) calloc(1, sizeof(multinomial_exp)));
						transcribe_multinomial_exp((multinomial_exp *) (*sub_mne_data).e[i], (multinomial_exp *) (*mne_data).e[i]);
						if (mode == MULTINOM_ADD) {}
						else if (mode == MULTINOM_SUB)
						{
							negative_multinomial_exp((multinomial_exp *) (*mne_data).e[i]);
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
				}
				else
				{
					if (mne_desc == NULL && sub_mne_desc != NULL)
					{
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
		else
		{
			add_sub_multinomial_data(mne1, mne2, mode);
			// add2multinomial_data(mne, sub_mne); 
		}
	}
}

void add_sub_multinomial_data(multinomial_exp *mne, multinomial_exp *sub_mne, char mode)
{
	// TBD
	if ((*mne).deg == 0 && (*sub_mne).deg == 0) {}
	else printf("Something weird happened! (add_sub_multinomial_data)\n");
	if (indeterminate_const_test(mne)) return;
	else if (indeterminate_const_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, 0, 0, NULL, DTYPE_INDETERMINATE_CONST);
		return;
	}
	else if (valid_numeric_test(mne) && valid_numeric_test(sub_mne))
	{
		if ((*mne).dtype == DTYPE_CONST_INT)
		{
			if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *mne_data = (int *) (*mne).data;
				int *sub_mne_data = (int *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
			}
			if ((*sub_mne).dtype == DTYPE_CONST_FLOAT)
			{
				int *mne_data = (int *) (*mne).data;
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
				(*mne).data = (void *) new_data;
				(*mne).dtype == DTYPE_CONST_FLOAT;
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
			}
			else if ((*sub_mne).dtype == DTYPE_CONST_INT)
			{
				int *sub_mne_data = (int *) (*sub_mne).data;
				if (mode == MULTINOM_ADD) (*mne_data) += (*sub_mne_data);
				else if (mode == MULTINOM_SUB) (*mne_data) -= (*sub_mne_data);
			}
		}
	}
	else if (!nonzero_test(mne) && nonzero_test(sub_mne))
	{
		free_multinomial_exp(mne);
		multinomial_exp_init(mne, 0, 0, (*sub_mne).data, (*sub_mne).dtype);
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
	if ((*mne).deg > 0 && (*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		int k = (*mne).dim - 1;
		int deg_ulim = (*mne).deg;
		while ((*mne_data).e[deg_ulim] == NULL && deg_ulim > 0)
		{
			deg_ulim -= 1;
		}
		add2array_char(mne_str, '(');
		for (int i = 0; i < deg_ulim; i++)
		{
			multinomial_exp *aux = (multinomial_exp *) (*mne_data).e[i];
			if (aux != NULL)
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
		if (aux != NULL)
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
		else 
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

void multiply_multinomial_exp(multinomial_exp *mne1, multinomial_exp *mne2, multinomial_exp *mne3)
{
	if (nonzero_test(mne1) && nonzero_test(mne2)) {}
	else
	{
		multinomial_exp_init(mne3, 0, 0, NULL, DTYPE_CONST_ZERO);
		return;
	}
	if ((*mne1).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne1).dim >= (*mne2).dim && (*mne1).dim > 0)
	{
		int mne3_deg;
		char status = ((*mne1).dim == (*mne2).dim);
		multinomial_exp_init(mne3, 0, 0, NULL, DTYPE_CONST_ZERO);
		array_voidstar *mne1_data = (array_voidstar *) (*mne1).data;
		multinomial_exp mne3_incr;
		multinomial_exp_init(&mne3_incr, 0, 0, NULL, DTYPE_ARRAY_VOIDSTAR);
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
				multinomial_exp aux;
				multinomial_exp *mne1_desc = (multinomial_exp *) (*mne1_data).e[i];
				multiply_multinomial_exp(mne1_desc, mne2, &aux);
				//display_multinomial_exp(&aux);
				if (!status || aux.deg == 0)
				{
					if (nonzero_test(&aux))
					{
						//mne3_incr.deg = i;
						//(*mne3_incr_data).len = i + 1;
						(*mne3_incr_data).e[i] = (void *) &aux;
						//mne3_incr.deg = i;
						add2multinomial_exp(mne3, &mne3_incr);
						//display_multinomial_exp(mne3);
						(*mne3_incr_data).e[i] = NULL;
					}
				}
				else
				{
					// shift aux_data over by i steps
					if (nonzero_test(&aux))
					{
						array_voidstar *aux_data = (array_voidstar *) aux.data;
						for (int ii = 0; ii <= aux.deg; ii++)
						{
							multinomial_exp *aux_di = (multinomial_exp *) (*aux_data).e[ii];
							(*mne3_incr_data).e[i + ii] = (*aux_data).e[ii];
						}
						//(*mne3_incr_data).len = (*aux_data).len + i;
						//mne3_incr.deg = aux.deg + i;
						add2multinomial_exp(mne3, &mne3_incr);
						int ii_max = aux.deg + i;
						for (int ii = i; ii <= ii_max; ii++)
						{
							(*mne3_incr_data).e[ii] = NULL;
						}
					}
				}
				free_multinomial_exp(&aux);
			}
		}
		free_multinomial_exp(&mne3_incr);
	}
	else if ((*mne2).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne2).dim > (*mne1).dim)
	{
		multiply_multinomial_exp(mne2, mne1, mne3);
	}
	else
	{
		int mne3_dtype;
		void *mne3_data;
		if (indeterminate_const_test(mne1) || indeterminate_const_test(mne2))
		{
			mne3_dtype = DTYPE_INDETERMINATE_CONST;
			mne3_data = NULL;
		}
		else if ((*mne1).dtype == DTYPE_CONST_INT)
		{
			int mne1_val = *((int *) (*mne1).data);
			if ((*mne2).dtype == DTYPE_CONST_INT)
			{
				mne3_dtype = DTYPE_CONST_INT;
				int mne2_val = *((int *) (*mne2).data);
				int *aux = (int *) calloc(1, sizeof(int));
				(*aux) = mne1_val * mne2_val;
				mne3_data = (void *) aux;
			}
			else if ((*mne2).dtype == DTYPE_CONST_FLOAT)
			{
				mne3_dtype = DTYPE_CONST_FLOAT;
				double mne2_val = *((double *) (*mne2).data);
				double *aux = (double *) calloc(1, sizeof(int));
				(*aux) = mne1_val * mne2_val;
				mne3_data = (void *) aux;
			}
		}
		else
		{
			mne3_dtype = DTYPE_CONST_FLOAT;
			double mne1_val = *((double *) (*mne1).data);
			double *aux = (double *) calloc(1, sizeof(int));
			if ((*mne2).dtype == DTYPE_CONST_INT)
			{
				int mne2_val = *((int *) (*mne2).data);
				(*aux) = mne1_val * mne2_val;
			}
			else if ((*mne2).dtype == DTYPE_CONST_FLOAT)
			{
				double mne2_val = *((double *) (*mne2).data);
				(*aux) = mne2_val * mne1_val;
			}
			else printf("Something weird happened!\n");
			mne3_data = (void *) aux;
		}
		multinomial_exp_init(mne3, 0, 0, mne3_data, mne3_dtype);
		free(mne3_data);
	}
}

void pow_mutinomial_exp(multinomial_exp *mne, multinomial_exp *mne_p, int p)
{
	if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR)
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
			multinomial_exp_init(mne_p, 0, 0, NULL, DTYPE_INDETERMINATE_CONST);
		}
		else
		{
			if ((*mne).dtype == DTYPE_CONST_INT)
			{
				int mne_val = *((int *) (*mne).data);
				int mne_p_val = mne_val;
				// Consider replacing this with an 'optimized' evaluator for large 'p'
				for (int i = 0; i < p - 1; i++)
				{
					mne_p_val *= mne_val;
				}
				multinomial_exp_init(mne_p, 0, 0, &mne_p_val, DTYPE_CONST_INT);
			}
			if ((*mne).dtype == DTYPE_CONST_FLOAT)
			{
				double mne_val = *((double *) (*mne).data);
				double mne_p_val = mne_val;
				for (int i = 0; i < p - 1; i++) mne_p_val *= mne_val;
				multinomial_exp_init(mne_p, 0, 0, &mne_p_val, DTYPE_CONST_FLOAT);
			}
		}
	}
}

void compose_multinomial_exp(multinomial_exp **mne1, int dim1, multinomial_exp **mne2, int dim2, multinomial_exp **mne3)
{
	for (int i = 0; i < dim1; i++)
	{
		if (nonzero_test(mne1[i])) {}
		else
		{
			multinomial_exp_init(mne3[i], 0, 0, NULL, DTYPE_CONST_ZERO);
			continue;
		}
		if ((*(mne1[i])).dtype == DTYPE_ARRAY_VOIDSTAR)
		{
			multinomial_exp_init(mne3[i], 0, 0, NULL, DTYPE_CONST_ZERO);
			multinomial_exp pow_;
			int val = 1;
			multinomial_exp_init(&pow_, 0, 0, &val, DTYPE_CONST_INT);
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
			multinomial_exp_init(mne3[i], 0, 0, (*(mne1[i])).data, (*(mne1[i])).dtype);
		}
	}
}

char nonzero_test(multinomial_exp *mne)
{
	if (mne == NULL) return 0;
	if (!((*mne).dtype == DTYPE_CONST_ZERO || 
			(*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data == NULL))
	{
		if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR)
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
		return 1;
	}
	else return 0;
}

void polynomial_int2multinomial_exp(array_int *p, multinomial_exp *mne, int dim)
{
	if ((*p).len > 0)
	{
		multinomial_exp_init(mne, dim, (*p).len - 1, NULL, DTYPE_ARRAY_VOIDSTAR);
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		array_voidstar_init(mne_data, (*p).len);
		for (int i = 0; i < (*p).len; i++)
		{
			if ((*p).e[i] != 0)
			{
				multinomial_exp *desc = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init(desc, 0, 0, &((*p).e[i]), DTYPE_CONST_INT);
				add2array_voidstar(mne_data, (void *) desc);
			}
			else add2array_voidstar(mne_data, NULL);
		}
	}
	else
	{
		multinomial_exp_init(mne, 0, 0, NULL, DTYPE_CONST_ZERO);
	}
}

void polynomial_double2multinomial_exp(array_double *p, multinomial_exp *mne, int dim)
{
	if ((*p).len > 0)
	{
		multinomial_exp_init(mne, dim, (*p).len - 1, NULL, DTYPE_ARRAY_VOIDSTAR);
		(*mne).data = (array_voidstar *) calloc(1, sizeof(array_voidstar));
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		array_voidstar_init(mne_data, (*p).len);
		for (int i = 0; i < (*p).len; i++)
		{
			if ((*p).e[i] != 0.0)
			{
				multinomial_exp *desc = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
				multinomial_exp_init(desc, 0, 0, &((*p).e[i]), DTYPE_CONST_FLOAT);
				add2array_voidstar(mne_data, (void *) desc);
			}
			else add2array_voidstar(mne_data, NULL);
		}
	}
	else
	{
		multinomial_exp_init(mne, 0, 0, NULL, DTYPE_CONST_ZERO);
	}
}

double multinomial_exp_eval_double(multinomial_exp *mne, double *vals)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL && (*mne).dim > 0)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		double rval = multinomial_exp_eval_double((multinomial_exp *) (*mne_data).e[0], &(vals[1]));
		double pow_ = vals[0];
		for (int i = 1; i <= (*mne).deg; i++)
		{
			rval += multinomial_exp_eval_double((multinomial_exp *) (*mne_data).e[i], &(vals[1])) * pow_;
			pow_ *= vals[0];
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
		}
		else
		{
			return 0.0;
		}
	}

}

int multinomial_exp_eval_int(multinomial_exp *mne, int *vals)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL && (*mne).dim > 0)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		int rval = multinomial_exp_eval_int((multinomial_exp *) (*mne_data).e[0], &(vals[1]));
		int pow_ = vals[0];
		for (int i = 1; i <= (*mne).deg; i++)
		{
			rval += multinomial_exp_eval_int((multinomial_exp *) (*mne_data).e[i], &(vals[1]));
			pow_ *= vals[0];
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
		}
		else
		{
			return 0;
		}
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
			case (DTYPE_ARRAY_VOIDSTAR):
			{
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
		}
	}
}

int count_indeterminates_rec(multinomial_exp *mne, array_int *sym_list, array_int *sym_, int *n_unspec)
{
	if (mne == NULL) return 0;
	if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL)
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
	if (dim < (*mne).dim)
	{
		int tmne_deg = (*mne).deg;
		multinomial_exp_init(tmne, (*mne).dim, (*mne).deg, NULL, DTYPE_ARRAY_VOIDSTAR);
		if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR) {}
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
	else if ((*mne).dim > 0)
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
			multinomial_exp_init(tmne, 0, 0, &c, DTYPE_INDETERMINATE_CONST);
		}
		else
		{
			c.sym_id = (*tails).len;
			add2array_int(sym_list, (*tails).len);
			multinomial_exp *tail_exp = (multinomial_exp *) calloc(1, sizeof(multinomial_exp));
			transcribe_multinomial_exp(mne, tail_exp);
			add2array_voidstar(tails, (void *) tail_exp);
			multinomial_exp_init(tmne, 0, 0, &c, DTYPE_INDETERMINATE_CONST);
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
			multinomial_exp_init(tmne, 0, 0, NULL, DTYPE_CONST_ZERO);
		}
	}
}

// Generate a 'truncated' multinomial of lower (effective) dimensionality, replacing 
// 	trailing polynomials in low dimensions with one 'indeterminate constant' per
// 	distinct 'equivalence class' (i.e. if a multinomial is truncated at dimension 1,
// 	then all instances of p_i(x_0) can be replaced by the same indeterminate constant,
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
		if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL)
		{
			array_voidstar *mne_data = (array_voidstar *) (*mne).data;
			multinomial_exp_init(imne, (*mne).dim, (*mne).deg, NULL, DTYPE_ARRAY_VOIDSTAR);
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
				multinomial_exp_init(imne, 0, 0, &c, DTYPE_INDETERMINATE_CONST);
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
	if ((*mne).dtype == DTYPE_ARRAY_VOIDSTAR && (*mne).data != NULL)
	{
		array_voidstar *mne_data = (array_voidstar *) (*mne).data;
		for (int i = 0; i < (*mne_data).len; i++)
		{
			c += 1 + multinomial_exp_complexity((multinomial_exp *) (*mne_data).e[i]);
		}
	}
	return c;
}

// RESUME

void multinomial_evaluator_init(multinomial_evaluator *mnev, multinomial_exp *expr)
{
	(*mnev).expr = expr;
	(*mnev).expr_1 = NULL;
	(*mnev).expr_2 = NULL;
	(*mnev).type = M_EV_ROOT_TYPE;
}

void free_multinomial_evaluator(multinomial_evaluator *mnev)
{
	// TBD
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

void multinomial_evaluator_eval_double(multinomial_evaluator *mnev, double *vals, double *cval)
{
	if ((*mnev).type == M_EV_COMP_TYPE)
	{
		double val1, val2;
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		multinomial_evaluator_eval_double(expr_1, vals, &val1);
		multinomial_evaluator_eval_double(expr_2, vals, &val2);
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
		(*cval) = multinomial_exp_eval_double((*mnev).expr, vals);
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
		int val1, val2;
		multinomial_evaluator *expr_1 = (multinomial_evaluator *) (*mnev).expr_1;
		multinomial_evaluator *expr_2 = (multinomial_evaluator *) (*mnev).expr_2;
		multinomial_evaluator_eval_int(expr_1, vals, &val1);
		multinomial_evaluator_eval_int(expr_2, vals, &val2);
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
		(*cval) = multinomial_exp_eval_int((*mnev).expr, vals);
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
