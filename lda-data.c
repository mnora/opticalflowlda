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

#include "lda-data.h"

corpus* read_data(char* data_filename, int num_angles)
{
    FILE *fileptr;
    int length, count, word, n, nd, nw,index,direction,location_length;
    corpus* c;
    index = -1;//location 计数
    printf("reading data from %s\n", data_filename);
    c = malloc(sizeof(corpus));
    c->docs = 0;
    c->num_terms = 0;
    c->num_docs = 0;
    fileptr = fopen(data_filename, "r");
    nd = 0; nw = 0;
    while ((fscanf(fileptr, "%10d", &length) != EOF))
    {
	c->docs = (document*) realloc(c->docs, sizeof(document)*(nd+1));
	c->docs[nd].length = length;
	c->docs[nd].total = 0;
	c->docs[nd].nlocation = 0;
	c->docs[nd].words = malloc(sizeof(int)*length);
	c->docs[nd].counts = malloc(sizeof(int)*length);
	c->docs[nd].locations = malloc(sizeof(int)*length);//in case word's length == location's length
	c->docs[nd].angles = malloc(sizeof(int)*length);
	c->docs[nd].count_locations = malloc(sizeof(int)*length);
	for (n = 0; n < length; n++)
	{
	    fscanf(fileptr, "%10d:%10d", &word, &count);
	    word = word - OFFSET;
	    c->docs[nd].words[n] = word;
	    c->docs[nd].counts[n] = count;
	    c->docs[nd].total += count;
	    if (word/num_angles!=index)
	    {
	        index = word/num_angles;//some location
            direction=word%num_angles;
	        location_length=(++c->docs[nd].nlocation);
	        c->docs[nd].locations[location_length-1]=index;//start from 0
	        c->docs[nd].angles[location_length-1]=direction; //only choose the first direction ex: 8:0 9:3 10:4 11:2 choose 9%4==1 as the location's angle
            c->docs[nd].count_locations[location_length-1]=count;
            c->docs[nd].total_locations+=count;
        }
	    if (word >= nw) { nw = word + 1; }
	}
	nd++;
    }
    fclose(fileptr);
    c->num_docs = nd;
    c->num_terms = nw;
    c->num_locations = nw/num_angles;
    printf("number of docs    : %d\n", nd);
    printf("number of terms   : %d\n", nw);
    printf("number of locations  : %d\n", c->num_locations);
    return(c);
}

int max_corpus_length(corpus* c)
{
    int n, max = 0;
    for (n = 0; n < c->num_docs; n++)
	if (c->docs[n].nlocation > max) max = c->docs[n].nlocation;
    return(max);
}
