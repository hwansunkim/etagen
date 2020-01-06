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

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <unistd.h>
#include "cluster.hpp"

trgLink::trgLink()
{
	trg = NULL;
	next = NULL;
}

trgLink::~trgLink()
{
	if(next != NULL)
	{
		delete next;
		next = NULL;
	}
}
void trgLink::merge(trgLink *data)
{
	if(next != NULL)
	{
		next->merge(data);
	}
	else
	{
		next = data;
	}
}
void trgLink::add(trginfo *t)
{
	if(next != NULL)
	{
		next->add(t);
	}
	else
	{
		trgLink *tmp = new trgLink();
		tmp->trg = t;
		next=tmp;
	}
}

int trgLink::matching(trginfo *t, FLOAT alpha, FLOAT beta)
{
	int ret = -1;
	/* Return 1 when matching succesively*/
	if((t->start_index - alpha)<= (trg->end_index + alpha) && (t->end_index + alpha) >= (trg->start_index - alpha) && t->fmin / (1.0 + beta) <= trg->fmax * (1.0 + beta) && t->fmax * (1.0 + beta) >= trg->fmin / (1.0 + beta))
	{
		return 1;
	}
	/* if not matched next trigger information matching start */
	if(next != NULL)
	{
		return next->matching(t, alpha, beta);
	}
	return ret;
}

int trgLink::numlink()
{
	int ret = 1;
	if(next != NULL)
	{
		return ret + next->numlink();
	}
	return ret;
}

int trgLink::show()
{
	int ret = 1;
	#if 0 
	printf("id : %.2f, ", trg->id);
	printf("st : %.2f, ", trg->start_index);
	printf("et : %.2f, ", trg->end_index);
	printf("pt : %.2f, ", trg->peak_index);
	printf("am : %.2f, ", trg->amplitude);
	printf("fr : %.2f, ", trg->frequency);
	printf("fmin : %.2f, ", trg->fmin);
	printf("fmax : %.2f, ", trg->fmax);
	printf("snr : %.2f\n", trg->snr);
	#else
	printf("id : %d, ", (int)trg->id);
	printf("st : %lf, ", trg->start_index);
	printf("et : %lf, ", trg->end_index);
	printf("fmin : %.2f, ", trg->fmin);
	printf("fmax : %.2f\n ", trg->fmax);
	#endif
 
	if(next != NULL)
	{
		ret += next->show();
	}
	return ret;
}

cluster::cluster()
{
	last = 1;	//tail node : 1, first or middle : 0;
	next = NULL;
	data = NULL;
}

cluster::~cluster()
{
	if(next != NULL)
	{
		delete next;
		next = NULL;
	}	
}

void cluster::merge(cluster *clt)
{
	/* representive  information update */
	start_index = std::fmin(start_index, clt->start_index);
	end_index = std::fmax(end_index, clt->end_index);
	fmin = std::fmin(fmin, clt->fmin);
	fmax = std::fmax(fmax, clt->fmax);

	data->merge(clt->data);	
}

int cluster::show()
{
	printf("cluster : %d\n", last);
	int ret = 0;
	int n = 0;
	if(data != NULL)
	{
		n = data->show();
	}
	if(next != NULL)
	{
		ret += next->show();
		printf("total : %d\n", n);
	}
	return ret + n;
}

void cluster::add(trginfo *t)
{
	if(next != NULL)
	{
		next->add(t);
	}
	else
	{
		last = 0; 
		cluster *tmp = new cluster();
		tmp->data = new trgLink();
		tmp->data->trg = t;
		tmp->start_index = t->start_index;
		tmp->end_index = t->end_index;
		tmp->fmin = t->fmin;
		tmp->fmax = t->fmax;

		next=tmp;
	}
}

cluster* cluster::clustering(trginfo *trg, FLOAT alpha, FLOAT beta)
{
	cluster *clt = NULL;
	int ret = -1;

	/* At first, It checks whether given trigger information matches the own information in the cluster. */
	if((start_index - alpha) <= (trg->end_index + alpha) && ( end_index + alpha )>= (trg->start_index - alpha) && fmin / (1.0 + beta) <= trg->fmax * (1.0 + beta) && fmax * (1.0 + beta) >= trg->fmin / (1.0 + beta))
	{
		ret = data->matching(trg, alpha, beta);
		if(ret) 
		{
			clt = this;
		}
	}
	
	/* If there are other clusters, It has to check whether there are matched clusters. */
	if(next != NULL)
	{
		cluster *tmp = next->clustering(trg, alpha, beta);
		if(tmp != NULL)
		{
			if(ret == 1)
			{
				//merge

				merge(tmp);
				if(tmp->last == 1)
				{
					last = 1;
					next = NULL;
				}
				if(tmp->next != NULL)
				{
					next = tmp->next;
					tmp->next = NULL;
					delete tmp;
				}

				clt = this;
			}
			else	//bypass
			{
				if(tmp->last == 1)
				{
					last = 1;
					tmp->last = 0;
					next = NULL;
				}
				if(tmp->next != NULL)
				{
					next = tmp->next;
					tmp->next = NULL;
				}
				clt = tmp;
			}		
		}
	}	
	return clt;
}
cluster* cluster::clustering(cluster *trg)
{
	return NULL;
}

void cluster::info(cltinfo *clt, FLOAT m, FLOAT* imf, int imf_num, int data_size, FLOAT start_time, FLOAT fsr)
{
	clt->start_index = start_index;
	clt->end_index = end_index;
	clt->fmin = fmin;
	clt->fmax = fmax;
	clt->nptr = data->numlink();
	clt->snr = 0.0;
	clt->snr_rss = 0.0;
	clt->p_snr = 0.0;
	clt->c_index = 0.0;
	clt->c_amp = 0.0;
	clt->c_freq = 0.0;

	trgLink *tmp = data;

	int n = (int)(end_index*fsr) - (int)(start_index * fsr) + 1;

	FLOAT *wave = new FLOAT[n];
	std::memset(wave, 0.0, sizeof(FLOAT)*n);
	int start = (start_index - start_time) * fsr;
	
	while(tmp != NULL)
	{
		if(tmp->trg->snr > clt->p_snr)
		{
			clt->p_index = tmp->trg->peak_index;
			clt->p_amp = tmp->trg->amplitude;
			clt->p_freq = tmp->trg->frequency;
			clt->p_imfindex = tmp->trg->id;
			clt->p_snr = tmp->trg->snr;
		}
		clt->c_index += tmp->trg->peak_index * tmp->trg->snr * tmp->trg->snr;
		clt->c_amp += tmp->trg->amplitude * tmp->trg->snr * tmp->trg->snr;
		clt->c_freq += tmp->trg->frequency * tmp->trg->snr * tmp->trg->snr;
		clt->snr_rss += tmp->trg->snr * tmp->trg->snr;

		trginfo *t = tmp->trg;
		int s = (int)(t->start_index*fsr) - (int)(start_time * fsr);
		int e = (int)(t->end_index*fsr) - (int)(start_time * fsr);
		int id = (int)t->id;

		for(int i= 0; i < e-s; i++)
		{

			wave[s-start + i] += *(imf + data_size * id + s +i); 
		}
		
		tmp = tmp->next;
	}
	//start snr_rss
	clt->c_index /= clt->snr_rss;
	clt->c_amp /= clt->snr_rss;
	clt->c_freq /= clt->snr_rss;
	clt->snr_rss = sqrt(clt->snr_rss);
	//end snr_rss

	for (int j=0; j < n; j++)
	{
		clt->snr += wave[j] * wave[j];
	}
	clt->snr = sqrt(clt->snr) / m / MADFACTOR;
	delete [] wave;
}

int cluster::numCluster()
{
	int ret = 1;
	if(next != NULL)
	{
		return ret + next->numCluster();
	}
	return ret;
}

trgLink* cluster::find(FLOAT s, FLOAT e, FLOAT fl, FLOAT fh)
{
	if (start_index == s && end_index == e && fmin == fl && fmax == fh)
	{
		return data;
	}
	else if(next != NULL)
	{
		return next->find(s,e,fl,fh);
	}

	return NULL;
}
triggerCluster::triggerCluster()
{
	root = NULL;
	alpha = 0.0;
	beta = 0.0;
}
triggerCluster::~triggerCluster()
{
	if(root != NULL)
	{
		delete root;
	}
}

void triggerCluster::feed(trginfo *trg)
{
	if(root == NULL)
	{
		root = new cluster();

		root->data = new trgLink();

		root->data->trg = trg;
		root->start_index = trg->start_index;
		root->end_index = trg->end_index;
		root->fmin = trg->fmin;
		root->fmax = trg->fmax;
	}
	else
	{
		cluster *tmp = root->clustering(trg, alpha, beta);
		if(tmp == NULL)
		{
			root->add(trg);
		}
		else
		{
			tmp->data->add(trg);
			tmp->start_index = std::fmin(tmp->start_index, trg->start_index);
			tmp->end_index = std::fmax(tmp->end_index, trg->end_index);
			tmp->fmin = std::fmin(tmp->fmin, trg->fmin);
			tmp->fmax = std::fmax(tmp->fmax, trg->fmax);
			if(tmp != root && tmp->next == NULL)
			{
				tmp->next = root;
				root = tmp;
			}
		}
	}

}

void triggerCluster::feed(cluster *clt)
{
	if(root == NULL)
	{
		root = clt;
	}
	else
	{
		cluster *tmp = root->clustering(clt);
		if(tmp == NULL)
		{
			clt->next = root;
			root = clt;
		}
		else
		{
			//merge clt and tmp
		}
	}

}

void triggerCluster::show()
{
	if(root != NULL)
	{
		int n = root->show();
		printf("Total cluster : %d\n", n);
	}
}

void triggerCluster::set_param(FLOAT a, FLOAT b, FLOAT m, FLOAT* imf_, int imf_num_, int data_size_, FLOAT start_time_, FLOAT fsr_)
{
	alpha = a;
	beta = b;
	median = m;
	imf = imf_;
	imf_num = imf_num_;
	data_size = data_size_;
	start_time = start_time_;
	fsr = fsr_;
}

cltinfo* triggerCluster::getClusteredTrigger(FLOAT th_snr)
{
	numofCluster = root->numCluster();	//maximum number of Clustered Trigger
	
	clt = new cltinfo[numofCluster];
	cluster *tmp = root;
	int i = 0;
	while(tmp != NULL)
	{
		tmp->info(&clt[i], median, imf, imf_num, data_size, start_time, fsr);

		if(clt[i].snr >= th_snr)
		{
			i++;
		}
		tmp = tmp->next;
	}
	numofCluster = i;
	printf("total Clusters : %d > %lf\n", numofCluster, th_snr);
	return clt;
}

int triggerCluster::getWaveform(int index, FLOAT **ret)
{
	long indx = (long)index;
	return getWaveform(indx, ret);
}

int triggerCluster::getWaveform(long index, FLOAT **ret)
{
	trgLink *trg = root->find(clt[index].start_index, clt[index].end_index, clt[index].fmin, clt[index].fmax);
	if (trg ==NULL)
		return 0;
	int n = (clt[index].end_index - clt[index].start_index) * fsr;
	FLOAT *wave = NULL;
	
	try{
		wave = new FLOAT[n];
	} catch (std::bad_alloc&) {
		throw; 
	}
	std::memset(wave, 0.0, sizeof(FLOAT)*n);
	int start = (clt[index].start_index - start_time) * fsr;
	while(trg != NULL)
	{
		trginfo *t = trg->trg;
		int s = (t->start_index - start_time) * fsr;
		int e = (t->end_index - start_time) * fsr;
		int id = (int)t->id;
		for(int i= 0; i < e-s; i++)
		{

			wave[s-start + i] += *(imf + data_size * id + s +i); 
		}

		trg = trg->next;
	}
	*ret = wave;
	return n;
}
