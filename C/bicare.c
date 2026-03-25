#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <float.h>

void echange (double *t, int *tab, int i, int j)   
{
    int tampon1 = tab[i];
    double tampon2 = t[i];
    tab[i] = tab[j];
    tab[j] = tampon1;
    t[i] = t[j];
    t[j] = tampon2;
}

void tri(double *t, int *tab, int G, int D)    
{
    int g,d;
    double val;
    
    if(D <= G)
 	return;
    val = t[D];
    g = G - 1;
    d = D;
    do
	{
	    while ((t[++g]>val));
	    while((t[--d]<val) && (d>G));
	    
	    if (g < d) echange(t, tab, g, d);
	} while (g < d);
    echange(t, tab, g, D);
    tri(t, tab, G, g-1);
    tri(t, tab, g+1, D);
}

int count_row_col(int z, int nrowColData, int *bicRowCol)
{
    int N = 0, i;
    for (i = 0 ; i < nrowColData ; i++) N += bicRowCol[nrowColData * z + i];
    return(N);
}

void sum(int z, int nrowData, int ncolData, double *Data, int *bicRow, int *bicCol, double *sumBic, double *sumRow, double *sumCol)
{
    int i, j, first = 0;
    
    sumBic[z] = 0;                                                  
    
    for (i = 0 ; i < nrowData ; i++)
	{
	    if (bicRow[nrowData * z + i])                             
		{
		    sumRow[nrowData * z + i] = 0;   
		    for (j = 0 ; j < ncolData ; j++)
			{
			    if (bicCol[ncolData * z + j])
				{
				    if (first == 0) sumCol[ncolData * z + j] = 0;  
				    sumRow[z * nrowData + i] += Data[ncolData * i + j];
				    sumCol[z * ncolData + j] += Data[ncolData * i + j];
				}
			}
		    sumBic[z] += sumRow[z * nrowData + i];
		    first = 1;
		}
	}
}

double residu(int z, int nrowData, int ncolData, double *Data, int *bicRow, int *bicCol, double *sumBic, double *sumRow, double *sumCol)
{
    int i, j;
    int nrowBic = count_row_col(z, nrowData, bicRow);
    int ncolBic = count_row_col(z, ncolData, bicCol);
    double res = 0, reselt;
    
    double const invrow = 1/(double)nrowBic;
    double const invcol = 1/(double)ncolBic;
    double const invvol = 1/(double)(nrowBic * ncolBic);
    double const meanbic = sumBic[z] * invvol;
    
    for (i = 0 ; i < nrowData ; i++)
	{
	    if (bicRow[nrowData * z + i])
		{ 
		    for (j = 0 ; j < ncolData ; j++)
			{
			    if (bicCol[ncolData * z + j])
				{
				    reselt = Data[ncolData * i + j]  - sumRow[z * nrowData + i] * invcol - sumCol[z * ncolData + j] * invrow + meanbic;   
				    res += reselt * reselt;
				}
			}
		}
	}
    res = res * invvol;
    return(res);
}

void bestgain(int k, double r, int nrowData, int ncolData, double *Data, int *bicRow, int *bicCol, int *bicRow2, int *bicCol2, double *sumBic, double *sumRow, double *sumCol, double *sumBic2, double *sumRow2, double *sumCol2, double *vecBestGain, int *vecBestBic, double *inv2R, double *vecResvolBic, int N, int M, int *vecBlocGene, int *vecBlocSample)
{
    int z, i, j, u;
    int nrowBic, ncolBic, volBic, volBic2;
    double invvol, invr2onres;
    double valGain, res1, res2;
    double gainMax = -DBL_MAX, gainMin = DBL_MAX;
    
    for (z = 0 ; z < k ; z++)              
	{
	    nrowBic = vecResvolBic[z * 4 + 2];
	    ncolBic = vecResvolBic[z * 4 + 3];
	    volBic = nrowBic * ncolBic;
	    invvol = 1 / (double)volBic;
	    
	    res1 = vecResvolBic[z * 4];
	    invr2onres = 1 / (r * r / res1);
	    
	    for (i = 0 ; i < nrowData ; i++) 
		{
		    if (!vecBlocGene[nrowData * z + i]){
			bicRow2[nrowData * z + i] = 1 - bicRow2[nrowData * z + i];
			u = 2 * bicRow2[nrowData * z + i] - 1; 
			if ((nrowBic + u) >= N)
			    { 
				sumBic2[z] = sumBic[z]; 
				sumRow2[nrowData * z + i] = 0;
				
				for (j = 0 ; j < ncolData ; j++)
				    {
					if (bicCol[ncolData * z + j] == 1)
					    {
						sumRow2[nrowData * z + i] += Data[ncolData * i + j];
						sumCol2[ncolData * z + j] = sumCol[ncolData * z +j] + u * Data[ncolData * i +j];
					    } 
				    }
				sumBic2[z] += u * sumRow2[nrowData * z + i];
				volBic2 = (nrowBic + u) * ncolBic;
				res2 = residu(z, nrowData, ncolData, Data, bicRow2, bicCol, sumBic2, sumRow2, sumCol2);
				valGain = (res1 - res2) * invr2onres + (volBic2 - volBic) * invvol;
			    }
			else valGain = -DBL_MAX;
		    }
		    
		    else valGain = -DBL_MAX;

		    if (z == 0)
			{
			    vecBestGain[i] = valGain;
			    vecBestBic[i] = 0;
			}
		    else if (valGain >= vecBestGain[i])
			{
			    vecBestGain[i] = valGain;
			    vecBestBic[i] = z;
			}
		    
		    if (vecBestGain[i] > gainMax) gainMax = vecBestGain[i];
		    if (vecBestGain[i] < gainMin) gainMin = vecBestGain[i];
		    bicRow2[nrowData * z + i] = bicRow[nrowData * z + i];
		}
	    
	    for (j = 0 ; j < ncolData ; j++)
		{    
		    if (bicCol[ncolData * z + j] == 1) sumCol2[ncolData * z + j] = sumCol[ncolData * z + j]; 
		}
	    sum(z, nrowData, ncolData, Data, bicRow, bicCol, sumBic2, sumRow2, sumCol2);
	    
	    for (j = 0; j < ncolData; j++)
		{
		    if (!vecBlocSample[ncolData *z + j]){
			bicCol2[ncolData * z + j] = 1 - bicCol2[ncolData * z + j];
			u = 2 * bicCol2[ncolData * z + j] - 1;
			
			if ((ncolBic + u) >= M)
			    {
				sumBic2[z] = sumBic[z];
				sumCol2[ncolData * z + j] = 0;
				for (i = 0; i < nrowData; i++)
				    {
					if (bicRow[nrowData * z + i] == 1)
					    {
						sumCol2[ncolData * z + j] += Data[ncolData * i + j];
						sumRow2[nrowData * z + i] = sumRow[nrowData * z + i] + u * Data[ncolData * i + j];
					    }
				    }
				sumBic2[z] += u * sumCol2[ncolData * z + j];
				volBic2 = nrowBic * (ncolBic + u);
				res2 = residu(z, nrowData, ncolData, Data, bicRow, bicCol2, sumBic2, sumRow2, sumCol2);
				valGain = (res1 - res2) * invr2onres + (volBic2 - volBic) * invvol;
			    }
			else valGain = -DBL_MAX;
		    }
		    else valGain= -DBL_MAX;
		    
		    if (z == 0)
			{
			    vecBestGain[j + nrowData] = valGain;
			    vecBestBic[j + nrowData] = 0;
			}
		    else if (valGain >= vecBestGain[j + nrowData])
			{
			    vecBestGain[j + nrowData] = valGain;
			    vecBestBic[j + nrowData] = z;
			}
		    if (vecBestGain[j + nrowData] > gainMax) gainMax = vecBestGain[i];
		    if (vecBestGain[j + nrowData] < gainMin) gainMin = vecBestGain[i];
		    bicCol2[ncolData * z +j] = bicCol[ncolData * z +j]; 
		    
		}

	    for (i = 0 ; i < nrowData ; i++)
		{
		    if (bicRow[nrowData * z + i] == 1) sumRow2[nrowData * z + i] = sumRow[nrowData * z + i];
		}
	    
	    vecResvolBic[z * 4] = res1;
	    vecResvolBic[z * 4 + 1] = volBic;
	    vecResvolBic[z * 4 + 2] = nrowBic;
	    vecResvolBic[z * 4 + 3] = ncolBic;
	}
    *inv2R = 1/(gainMax - gainMin);
}

void order(double *inv2R,int nrowData,int ncolData,double *vecBestGain, int *vecOrder)
{
    int i, rand1, rand2, rand1b, rand2b; 
    double proba, pij;
    
    for (i = 0 ; i < 2 * (nrowData + ncolData) ; i++) 
	{
	    rand1 = rand() % (nrowData + ncolData);
        rand2 = rand() % (nrowData + ncolData);
        rand1b = vecOrder[rand1];
        rand2b = vecOrder[rand2];
	    
	    pij = 0.5 + (vecBestGain[rand2b] - vecBestGain[rand1b]) * *inv2R;
	    
	    proba = (double)rand() / (RAND_MAX + 1.0);
	    
	    if (pij >= proba)  
		{
		    vecOrder[rand1] = rand2b;
		    vecOrder[rand2] = rand1b;
		}
	} 
}

void action(int k, int nrowData, int ncolData, double *Data, int *vecOrder, int *vecBestBic, int *bicRow, int *bicCol, int *bicRow2, int *bicCol2, double r, int *valbreak, double *vecResvolBic, double *sumBic, double *sumRow, double *sumCol, double *sumBic2, double *sumRow2, double *sumCol2, int N, int M, int zz, int *aa, int *vecBlocGene, int *vecBlocSample)
{
    int x, b, z, i, j, u, nrowBic, ncolBic, nrowBic2, ncolBic2, volBic, volBic2, amelio, contrainte;
    double res1, res2;
    
    for (x = 0 ; x < (nrowData + ncolData) ; x++)  
	{
	    b = vecOrder[x]; 
	    z = vecBestBic[b]; 
	    amelio = contrainte = 0;
	    
	    res1 = vecResvolBic[z * 4];
	    nrowBic = vecResvolBic[z * 4 + 2];
	    ncolBic = vecResvolBic[z * 4 + 3];
	    volBic = vecResvolBic[z * 4 + 1];
	    
	    sumBic2[z] = sumBic[z];
	    
	    if (b < nrowData) 
		{
		    i = b;
		    bicRow2[nrowData * z + i] = 1 - bicRow[nrowData * z + i];
		    u = 2 * bicRow2[nrowData * z + i] - 1;
		    nrowBic2 = nrowBic + u;
		    ncolBic2 = ncolBic;
		    
		    if(nrowBic2 >= N && !vecBlocGene[nrowData * z + i]) 
			{  
			    contrainte = 1;
			    sumRow2[nrowData * z + i] = 0;
			    for (j = 0 ; j < ncolData ; j++)
				{
				    if (bicCol[ncolData * z + j] == 1)
					{
					    sumRow2[nrowData * z + i] += Data[ncolData * i + j];
					    sumCol2[ncolData * z + j] += u * Data[ncolData * i + j];
					}
				}
			    sumBic2[z] += u * sumRow2[nrowData * z + i];
			}
		}
	    else  
		{
		    j = b - nrowData;
		    bicCol2[ncolData * z + j] = 1 - bicCol[ncolData * z + j];
		    u = 2 * bicCol2[ncolData * z + j] - 1;
		    ncolBic2 = ncolBic + u;
		    nrowBic2 = nrowBic;
		    
		    if (ncolBic2 >= M && !vecBlocSample[ncolData * z + j])
			{
			    contrainte = 1;
			    sumCol2[ncolData * z + j] = 0;
			    for (i = 0; i < nrowData; i++)
				{
				    if (bicRow[nrowData * z + i] == 1)
					{
					    sumCol2[ncolData * z + j] += Data[ncolData * i + j];
					    sumRow2[nrowData * z + i] += u * Data[ncolData * i + j];
					}
				}
			    sumBic2[z] += u * sumCol2[ncolData * z + j];
			}
		}
	    
	    if (contrainte)
		{
		    res2 = residu(z, nrowData, ncolData, Data, bicRow2, bicCol2, sumBic2, sumRow2, sumCol2); 
		    volBic2 = nrowBic2 * ncolBic2;
		    
		    if ((res2 < res1) || ((res2 < r) && (volBic2 > volBic))) 
			{
                            *aa += 1;
			    
			    amelio = 1;           
			    vecResvolBic[z * 4] = res2;
			    vecResvolBic[z * 4 + 1] = volBic2;
			    vecResvolBic[z * 4 + 2] = nrowBic2;
			    vecResvolBic[z * 4 + 3] = ncolBic2;
			    *valbreak = 1;
			    sumBic[z] = sumBic2[z];
			    
			    if (b < nrowData)  
				{
				    bicRow[nrowData * z + i] = bicRow2[nrowData * z + i];
				    sumRow[nrowData * z + i] = sumRow2[nrowData * z + i];
				    for (j = 0 ; j < ncolData ; j++)
					{
					    if (bicCol[ncolData * z + j] == 1)
						{
						    sumCol[ncolData * z + j] = sumCol2[ncolData * z + j];
						}
					}
				}
			    else 
				{
				    bicCol[ncolData * z + j] = bicCol2[ncolData * z +j];
				    sumCol[ncolData * z + j] = sumCol2[ncolData * z + j];
				    for(i = 0 ; i < nrowData ; i++)
					{
					    if (bicRow[nrowData * z + i] == 1)
						{
						    sumRow[nrowData * z + i] = sumRow2[nrowData * z + i];
						}
					}
				}
			}
		}	     
	    if(!amelio)
		{
		    if (b < nrowData)
			{
			    bicRow2[nrowData * z + i] = bicRow[nrowData * z + i];
			    for (j=0; j<ncolData; j++)
				{
				    if (bicCol[ncolData * z + j] == 1)
					{
					    sumCol2[ncolData * z + j] = sumCol[ncolData * z + j];
					}
				}
			}
		    else
			{
			    bicCol2[ncolData * z + j] = bicCol[ncolData * z +j];
			    for(i = 0 ; i < nrowData ; i++)
				{
				    if (bicRow[nrowData * z + i] == 1)
					{
					    sumRow2[nrowData * z + i] = sumRow[nrowData * z + i];
					}
				}
			}
		}
	}
}

void floc(double *Data, int *nrowData, int *ncolData, int *bicRow, int *bicCol, double *vecResvolBic, double *r, int *k, int *N, int *M, int *t, int *vecBlocGene, int *vecBlocSample)
{
	
	srand(time(NULL));
	
    int zz, i, valbreak,j;    
    int total = *nrowData + *ncolData;
    double inv2R;
    int *bicRow2 = malloc(*k * *nrowData * sizeof(int));
    int *bicCol2 = malloc(*k * *ncolData * sizeof(int));
    int *vecOrder = malloc(total * sizeof(int));
    double *vecBestGain = malloc(total * sizeof(double));
    int *vecBestBic = malloc(total * sizeof(int));
    double *sumRow = malloc(*k * *nrowData * sizeof(double));
    double *sumCol = malloc(*k * *ncolData * sizeof(double));
    double *sumBic = malloc(*k * sizeof(double));
    double *sumRow2 = malloc(*k * *nrowData * sizeof(double));
    double *sumCol2 = malloc(*k * *ncolData * sizeof(double));
    double *sumBic2 = malloc(*k * sizeof(double));
    
    double invk = 1/(double)*k;
    
    memcpy(bicRow2, bicRow, *k * *nrowData * sizeof(int));
    memcpy(bicCol2, bicCol, *k * *ncolData * sizeof(int));

    for (zz = 0 ; zz < *k ; zz++)
	{
	    vecResvolBic[zz * 4 + 2] = count_row_col(zz, *nrowData, bicRow);
	    vecResvolBic[zz * 4 + 3] = count_row_col(zz, *ncolData, bicCol);
	    sum(zz, *nrowData, *ncolData, Data, bicRow, bicCol, sumBic, sumRow, sumCol);  
	    sum(zz, *nrowData, *ncolData, Data, bicRow, bicCol, sumBic2, sumRow2, sumCol2);
	    vecResvolBic[zz * 4] = residu(zz, *nrowData, *ncolData, Data, bicRow, bicCol, sumBic, sumRow, sumCol);
	}
    j = 0;
    for (zz = 0 ; zz < *t ; zz++)
	{
	    valbreak = 0;
	    
            bestgain(*k, *r, *nrowData, *ncolData, Data, bicRow, bicCol, bicRow2, bicCol2, sumBic, sumRow, sumCol, sumBic2, sumRow2, sumCol2, vecBestGain, vecBestBic, &inv2R, vecResvolBic, *N, *M, vecBlocGene, vecBlocSample);
	    
	    for (i = 0 ; i < total ; i++) vecOrder[i] = i; 
	    
	    tri(vecBestGain, vecOrder, 0, total-1);
	    order(&inv2R, *nrowData, *ncolData, vecBestGain, vecOrder);
	    
	    action(*k, *nrowData, *ncolData, Data, vecOrder, vecBestBic, bicRow, bicCol, bicRow2, bicCol2, *r, &valbreak, vecResvolBic, sumBic, sumRow, sumCol, sumBic2, sumRow2, sumCol2, *N, *M, zz, &j, vecBlocGene, vecBlocSample);
	    
	    if (valbreak == 0)
		{
		    printf("\n STOP\n ");
		    break;
		}   
	}
    
    free(bicRow2);
    free(bicCol2);
    free(vecOrder);
    free(vecBestGain);
    free(vecBestBic);
    free(sumRow);
    free(sumCol);
    free(sumBic);
    free(sumRow2);
    free(sumCol2);
    free(sumBic2);
    
}

void printres(int *nrowData, int *ncolData, double *Data, int *bicRow, int *bicCol, double *res)
{
    double *sumRow = malloc(*nrowData * sizeof(double));
    double *sumCol = malloc(*ncolData * sizeof(double));
    double *sumBic = malloc(sizeof(double));
    
    sum(0, *nrowData, *ncolData, Data, bicRow, bicCol, sumBic, sumRow, sumCol);
    *res = residu(0, *nrowData, *ncolData, Data, bicRow, bicCol, sumBic, sumRow, sumCol);
    
    free(sumRow);
    free(sumCol);
    free(sumBic);
}