// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "lda-model.h"

/*
 * compute MLE lda model from sufficient statistics
 *
 */

void lda_mle(lda_model* model, lda_suffstats* ss, int estimate_alpha)
{
    int k; int w;

    for (k = 0; k < model->num_topics; k++)
    {
        for (w = 0; w < model->num_locations; w++)
        {
            if (ss->class_word[k][w] > 0)
            {
                model->log_prob_w[k][w] =
                    log(ss->class_word[k][w]) -
                    log(ss->class_totalw[k]);
            }
            else
                model->log_prob_w[k][w] = -100;
        }
    }
    for (k = 0; k < model->num_locations; k++)
    {
        for (w = 0; w < model->num_angles; w++)
        {
            if (ss->class_angle[k][w] > 0)
            {
                model->log_prob_a[k][w] =
                    log(ss->class_angle[k][w]) -
                    log(ss->class_totala[k]);
            }
            else
                model->log_prob_a[k][w] = -100;
        }
    }
    if (estimate_alpha == 1)
    {
        model->alpha = opt_alpha(ss->alpha_suffstats,
                                 ss->num_docs,
                                 model->num_topics);

        printf("new alpha = %5.5f\n", model->alpha);
    }
    // does eta need?
}

/*
 * allocate sufficient statistics
 *
 */

lda_suffstats* new_lda_suffstats(lda_model* model)
{
    int num_topics = model->num_topics;
    int num_locations = model->num_locations;
    int num_angles = model->num_angles;
    int i,j,k;

    lda_suffstats* ss = malloc(sizeof(lda_suffstats));
    ss->class_totalw = malloc(sizeof(double)*num_topics);
    ss->class_word = malloc(sizeof(double*)*num_topics);
    ss->class_totala = malloc(sizeof(double)*num_locations);
    ss->class_angle = malloc(sizeof(double*)*num_locations);
    for (i = 0; i < num_topics; i++)
    {
        ss->class_totalw[i] = 0;
        ss->class_word[i] = malloc(sizeof(double)*num_locations);
        for (j = 0; j < num_locations; j++)
        {
            ss->class_word[i][j] = 0;
        }
    }

    for(i = 0; i< num_locations; i++)
    {
        ss->class_totala[i]=0;
        ss->class_angle[i]=malloc(sizeof(double)*num_angles);
        for(k=0;k<num_angles;k++)
        {
	    ss->class_angle[i][k]=0;
        }
    }

    return(ss);
}


/*
 * various intializations for the sufficient statistics
 *
 */

void zero_initialize_ss(lda_suffstats* ss, lda_model* model)
{
    int k, w;
    for (k = 0; k < model->num_topics; k++)
    {
        ss->class_totalw[k] = 0;
        for (w = 0; w < model->num_locations; w++)
        {
            ss->class_word[k][w] = 0;
        }
    }
    for (k = 0; k < model->num_locations; k++)
    {
        ss->class_totala[k] = 0;
        for (w = 0; w < model->num_angles; w++)
        {
            ss->class_angle[k][w] = 0;
        }
    }
    ss->num_docs = 0;
    ss->alpha_suffstats = 0;
}


void random_initialize_ss(lda_suffstats* ss, lda_model* model)
{
    int num_topics = model->num_topics;
    int num_locations = model->num_locations;
    int num_angles= model->num_angles;
    int k, n;
    for (k = 0; k < num_topics; k++)
    {
        for (n = 0; n < num_locations; n++)
        {
            ss->class_word[k][n] += 1.0/num_locations + myrand();
            ss->class_totalw[k] += ss->class_word[k][n];
        }
    }
    for (k = 0; k < num_locations; k++)
    {
        for (n = 0; n < num_angles; n++)
        {
            ss->class_angle[k][n] += 1.0/num_angles + myrand();
            ss->class_totala[k] += ss->class_angle[k][n];
        }
    }
}


void corpus_initialize_ss(lda_suffstats* ss, lda_model* model, corpus* c)
{
    int num_topics = model->num_topics;
    int num_angles = model->num_angles;
    int num_locations = model->num_locations;
    int i, k, d, n,m;
    document* doc;

    for (k = 0; k < num_topics; k++)
    {
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(myrand() * c->num_docs);
            printf("initialized with document %d\n", d);
            doc = &(c->docs[d]);
            for (n = 0; n < doc->length; n++)
            {
                ss->class_word[k][doc->words[n]/num_angles] += doc->counts[n];
                ss->class_angle[k][doc->words[n]%num_angles]+=doc->counts[n];
            }
        }
        for (n = 0; n < model->num_locations; n++)
        {
            ss->class_word[k][n] += 1.0;
            ss->class_totalw[k] = ss->class_totalw[k] + ss->class_word[k][n];
        }
    }
    for (k = 0; k < num_locations; k++)
    {
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(myrand() * c->num_docs);
            printf("initialized with document %d\n", d);
            doc = &(c->docs[d]);
            for (n = 0; n < doc->length; n++)
            {
                ss->class_angle[k][doc->words[n]%num_angles]+=doc->counts[n];
            }
        }
        for(m = 0;m < model->num_angles; m++)
        {
            ss->class_angle[k][m]+=1.0;
            ss->class_totala[k]=ss->class_totala[k]+ss->class_angle[k][m];
        }
    }
}

/*
 * allocate new lda model
 *
 */

lda_model* new_lda_model(int num_locations, int num_topics,int num_angles)
{
    int i,j;
    lda_model* model;

    model = malloc(sizeof(lda_model));
    model->num_topics = num_topics;
    model->num_locations = num_locations;
    model->num_angles = num_angles;
    model->alpha = 1.0;
    model->log_prob_w = malloc(sizeof(double*)*num_topics);
    model->log_prob_a = malloc(sizeof(double*)*num_locations);
    for (i = 0; i < num_topics; i++)
    {
	model->log_prob_w[i] = malloc(sizeof(double)*num_locations);
	for (j = 0; j < num_locations; j++)
	    model->log_prob_w[i][j] = 0;
    }
    for (i = 0; i< num_locations; i++)
    {
    model->log_prob_a[i] = malloc(sizeof(double)*num_angles);
	for(j=0; j<num_angles;j++)
        model->log_prob_a[i][j]=0;
    }
    return(model);
}


/*
 * deallocate new lda model
 *
 */

void free_lda_model(lda_model* model)
{
    int i;

    for (i = 0; i < model->num_topics; i++)
    {
	free(model->log_prob_w[i]);
    }
    free(model->log_prob_w);
}


/*
 * save an lda model
 *
 */

void save_lda_model(lda_model* model, char* model_root)
{
    char filename[100];
    FILE* fileptr;
    int i, j;

    sprintf(filename, "%s.beta", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics; i++)
    {
	for (j = 0; j < model->num_locations; j++)
	{
	    fprintf(fileptr, " %5.10f", model->log_prob_w[i][j]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.eta", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_locations; i++)
    {
	for (j = 0; j < model->num_angles; j++)
	{
	    fprintf(fileptr, " %5.10f", model->log_prob_a[i][j]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->num_topics);
    fprintf(fileptr, "num_terms %d\n", model->num_terms);
    fprintf(fileptr, "num_angles %d\n", model->num_angles);
    fprintf(fileptr, "num_locations %d\n", model->num_locations);
    fprintf(fileptr, "alpha %5.10f\n", model->alpha);
    fclose(fileptr);
}


lda_model* load_lda_model(char* model_root)
{
    char filename[100];
    FILE* fileptr;
    int i, j, num_terms, num_topics,num_angles,num_locations;
    float x, alpha;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
    fscanf(fileptr, "num_angles %d\n", &num_angles);
    fscanf(fileptr, "num_locations %d\n", &num_locations);
    fscanf(fileptr, "alpha %f\n", &alpha);
    fclose(fileptr);

    lda_model* model = new_lda_model(num_locations, num_topics,num_angles);
    model->alpha = alpha;

    sprintf(filename, "%s.beta", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_locations; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_w[i][j] = x;
        }
    }

    sprintf(filename, "%s.eta", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_locations; i++)
    {
        for (j = 0; j < num_angles; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_a[i][j] = x;
            //printf("%f",x);
        }
    }
    fclose(fileptr);
    return(model);
}
