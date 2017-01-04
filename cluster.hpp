/* Copyright (C) 2016  Whansun Kim
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include "common.h"

/* Trigger information structure */
typedef struct  
{
	FLOAT id;
	FLOAT start_index;
	FLOAT end_index;
	FLOAT peak_index;
	FLOAT amplitude;
	FLOAT frequency;
	FLOAT fmin;
	FLOAT fmax;
	FLOAT snr;
}trginfo;

typedef struct 
{
	FLOAT start_index;   	//cluster start time
	FLOAT end_index;    	//clustere end time
	FLOAT fmin;             //clustered min frequency
	FLOAT fmax;             //clustered max frequency
	FLOAT c_index;          // sum( peak_index * snr) / total_snr;
	FLOAT c_amp;            // sum( peak_amplitude * snr) / total_snr;
	FLOAT c_freq;           // sum( peak_frequency * snr) / total_snr;
	FLOAT p_index;          // peak index
	FLOAT p_amp;            // peak amplitude
	FLOAT p_freq;           // peak frequency
	FLOAT p_imfindex;       // peak imf id
	FLOAT p_snr;            // peak snr
	FLOAT nptr;             // number of triggers
	FLOAT snr;	            // total snr;
}cltinfo;

/* Linked list of Same Cluster */
class trgLink
{
public:
	trginfo *trg;
	trgLink *next;

	trgLink();
	~trgLink();
	int show();
	void merge(trgLink *data);
	void add(trginfo *t);
	int matching(trginfo *t, FLOAT alpha, FLOAT beta);
	int numlink();
};

/* Linked list of Clusters */
class cluster
{
public:
	/* Representive information of this Cluster */
	FLOAT start_index;
	FLOAT end_index;
	FLOAT fmin;
	FLOAT fmax;

	/* Flag of last link */
	int last;

	trgLink *data;
	cluster *next;

	cluster();
	~cluster();
	void add(trginfo *trg);
	void merge(cluster *clt);
	cluster* clustering(trginfo *trg, FLOAT alpha, FLOAT beta);
	cluster* clustering(cluster *clt);
	int show();
	int numCluster();
	void info(cltinfo*, FLOAT, FLOAT*, int, int, FLOAT, FLOAT);
	trgLink* find(FLOAT, FLOAT, FLOAT, FLOAT);
};

class triggerCluster
{
private:
	cltinfo *clt;
	cluster *root;
	FLOAT alpha;
	FLOAT beta;
	FLOAT median;
	FLOAT* imf;
	int imf_num;
	int data_size;
	FLOAT start_time;
	FLOAT fsr;

public:
	int numofCluster;
	triggerCluster();
	~triggerCluster();

	void feed(trginfo*);
	void feed(cluster*);
	void show();
	void set_param(FLOAT, FLOAT, FLOAT, FLOAT*, int, int, FLOAT, FLOAT);
	cltinfo* getClusteredTrigger(FLOAT);
	int getWaveform(int index, FLOAT**);
};

#endif
