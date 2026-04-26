#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <float.h>

#define IDX(i,j,k,l) ((((l)*ncondData + (k))*ncolData + (j))*nrowData + (i))

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

void sum4d(int z,
           int nrowData, int ncolData, int ncondData, int ntimeData,
           double *Data,
           int *bicRow, int *bicCol, int *bicCond, int *bicTime,
           double *sumBic,
           double *sumRow, double *sumCol,
           double *sumCond, double *sumTime)
{
    int i, j, k, l;

    sumBic[z] = 0;

    /* reset sums */
    for (i = 0; i < nrowData; i++)
        sumRow[z * nrowData + i] = 0;

    for (j = 0; j < ncolData; j++)
        sumCol[z * ncolData + j] = 0;

    for (k = 0; k < ncondData; k++)
        sumCond[z * ncondData + k] = 0;

    for (l = 0; l < ntimeData; l++)
        sumTime[z * ntimeData + l] = 0;

    /* accumulate */
    for (i = 0; i < nrowData; i++)
    {
        if (!bicRow[z * nrowData + i]) continue;

        for (j = 0; j < ncolData; j++)
        {
            if (!bicCol[z * ncolData + j]) continue;

            for (k = 0; k < ncondData; k++)
            {
                if (!bicCond[z * ncondData + k]) continue;

                for (l = 0; l < ntimeData; l++)
                {
                    if (!bicTime[z * ntimeData + l]) continue;

                    double val = Data[IDX(i,j,k,l)];

                    sumRow[z * nrowData + i] += val;
                    sumCol[z * ncolData + j] += val;
                    sumCond[z * ncondData + k] += val;
                    sumTime[z * ntimeData + l] += val;

                    sumBic[z] += val;
                }
            }
        }
    }
}

double residu4d(int z,
                int nrowData, int ncolData, int ncondData, int ntimeData,
                double *Data,
                int *bicRow, int *bicCol, int *bicCond, int *bicTime,
                double *sumBic,
                double *sumRow, double *sumCol,
                double *sumCond, double *sumTime)
{
    int i, j, k, l;

    int nrowBic  = count_row_col(z, nrowData,  bicRow);
    int ncolBic  = count_row_col(z, ncolData,  bicCol);
    int ncondBic = count_row_col(z, ncondData, bicCond);
    int ntimeBic = count_row_col(z, ntimeData, bicTime);

    double res = 0, reselt;

    double invrow  = 1.0 / nrowBic;
    double invcol  = 1.0 / ncolBic;
    double invcond = 1.0 / ncondBic;
    double invtime = 1.0 / ntimeBic;

    double invvol = 1.0 / (nrowBic * ncolBic * ncondBic * ntimeBic);

    double meanbic = sumBic[z] * invvol;

    for (i = 0; i < nrowData; i++)
    {
        if (!bicRow[z * nrowData + i]) continue;

        for (j = 0; j < ncolData; j++)
        {
            if (!bicCol[z * ncolData + j]) continue;

            for (k = 0; k < ncondData; k++)
            {
                if (!bicCond[z * ncondData + k]) continue;

                for (l = 0; l < ntimeData; l++)
                {
                    if (!bicTime[z * ntimeData + l]) continue;

                    double x = Data[IDX(i,j,k,l)];

                    double rowMean  = sumRow[z * nrowData + i] * (invcol * invcond * invtime);
                    double colMean  = sumCol[z * ncolData + j] * (invrow * invcond * invtime);
                    double condMean = sumCond[z * ncondData + k] * (invrow * invcol * invtime);
                    double timeMean = sumTime[z * ntimeData + l] * (invrow * invcol * invcond);

                    reselt = x
                        - rowMean
                        - colMean
                        - condMean
                        - timeMean
                        + 3.0 * meanbic;

                    res += reselt * reselt;
                }
            }
        }
    }

    res *= invvol;
    return res;
}

void bestgain4d(
    int k, double r,
    int nrowData, int ncolData, int ncondData, int ntimeData,
    double *Data,
    int *bicRow, int *bicCol, int *bicCond, int *bicTime,
    int *bicRow2, int *bicCol2, int *bicCond2, int *bicTime2,
    double *sumBic,
    double *sumRow, double *sumCol, double *sumCond, double *sumTime,
    double *sumBic2,
    double *sumRow2, double *sumCol2, double *sumCond2, double *sumTime2,
    double *vecBestGain, int *vecBestBic,
    double *inv2R,
    double *vecResvolBic,
    int N, int M, int P, int Q,
    int *vecBlocGene, int *vecBlocSample,
    int *vecBlocCond, int *vecBlocTime)
{
    int z, i, j, k2, l, u;
    double gainMax = -DBL_MAX, gainMin = DBL_MAX;

    int total = nrowData + ncolData + ncondData + ntimeData;

    for (z = 0; z < k; z++)
    {
        double res1 = vecResvolBic[z * 6];
        double invr2onres = 1 / (r * r / res1);

        int nrowBic  = vecResvolBic[z * 6 + 2];
        int ncolBic  = vecResvolBic[z * 6 + 3];
        int ncondBic = vecResvolBic[z * 6 + 4];
        int ntimeBic = vecResvolBic[z * 6 + 5];

        int volBic = nrowBic * ncolBic * ncondBic * ntimeBic;
        double invvol = 1.0 / volBic;

        /* iterate over ALL dimensions */
        for (i = 0; i < total; i++)
        {
            double valGain = -DBL_MAX;

            /* determine dimension */
            if (i < nrowData)  /* gene */
            {
                if (!vecBlocGene[nrowData * z + i])
                {
                    bicRow2[nrowData * z + i] = 1 - bicRow[nrowData * z + i];
                    u = 2 * bicRow2[nrowData * z + i] - 1;

                    if ((nrowBic + u) >= N)
                    {
                        double res2 = residu4d(z,
                            nrowData, ncolData, ncondData, ntimeData,
                            Data,
                            bicRow2, bicCol, bicCond, bicTime,
                            sumBic, sumRow, sumCol, sumCond, sumTime);

                        int volBic2 = (nrowBic + u) * ncolBic * ncondBic * ntimeBic;

                        valGain = (res1 - res2) * invr2onres +
                                  (volBic2 - volBic) * invvol;
                    }

                    bicRow2[nrowData * z + i] = bicRow[nrowData * z + i];
                }
            }
            else if (i < nrowData + ncolData)  /* sample */
            {
                j = i - nrowData;

                if (!vecBlocSample[ncolData * z + j])
                {
                    bicCol2[ncolData * z + j] = 1 - bicCol[ncolData * z + j];
                    u = 2 * bicCol2[ncolData * z + j] - 1;

                    if ((ncolBic + u) >= M)
                    {
                        double res2 = residu4d(z,
                            nrowData, ncolData, ncondData, ntimeData,
                            Data,
                            bicRow, bicCol2, bicCond, bicTime,
                            sumBic, sumRow, sumCol, sumCond, sumTime);

                        int volBic2 = nrowBic * (ncolBic + u) * ncondBic * ntimeBic;

                        valGain = (res1 - res2) * invr2onres +
                                  (volBic2 - volBic) * invvol;
                    }

                    bicCol2[ncolData * z + j] = bicCol[ncolData * z + j];
                }
            }
            else if (i < nrowData + ncolData + ncondData)  /* condition */
            {
                k2 = i - nrowData - ncolData;

                if (!vecBlocCond[ncondData * z + k2])
                {
                    bicCond2[ncondData * z + k2] = 1 - bicCond[ncondData * z + k2];
                    u = 2 * bicCond2[ncondData * z + k2] - 1;

                    if ((ncondBic + u) >= P)
                    {
                        double res2 = residu4d(z,
                            nrowData, ncolData, ncondData, ntimeData,
                            Data,
                            bicRow, bicCol, bicCond2, bicTime,
                            sumBic, sumRow, sumCol, sumCond, sumTime);

                        int volBic2 = nrowBic * ncolBic * (ncondBic + u) * ntimeBic;

                        valGain = (res1 - res2) * invr2onres +
                                  (volBic2 - volBic) * invvol;
                    }

                    bicCond2[ncondData * z + k2] = bicCond[ncondData * z + k2];
                }
            }
            else  /* time */
            {
                l = i - nrowData - ncolData - ncondData;

                if (!vecBlocTime[ntimeData * z + l])
                {
                    bicTime2[ntimeData * z + l] = 1 - bicTime[ntimeData * z + l];
                    u = 2 * bicTime2[ntimeData * z + l] - 1;

                    if ((ntimeBic + u) >= Q)
                    {
                        double res2 = residu4d(z,
                            nrowData, ncolData, ncondData, ntimeData,
                            Data,
                            bicRow, bicCol, bicCond, bicTime2,
                            sumBic, sumRow, sumCol, sumCond, sumTime);

                        int volBic2 = nrowBic * ncolBic * ncondBic * (ntimeBic + u);

                        valGain = (res1 - res2) * invr2onres +
                                  (volBic2 - volBic) * invvol;
                    }

                    bicTime2[ntimeData * z + l] = bicTime[ntimeData * z + l];
                }
            }

            if (z == 0 || valGain >= vecBestGain[i])
            {
                vecBestGain[i] = valGain;
                vecBestBic[i] = z;
            }

            if (vecBestGain[i] > gainMax) gainMax = vecBestGain[i];
            if (vecBestGain[i] < gainMin) gainMin = vecBestGain[i];
        }
    }

    *inv2R = 1.0 / (gainMax - gainMin);
}

void order4d(double *inv2R, int total,
             double *vecBestGain, int *vecOrder)
{
    int i, rand1, rand2, rand1b, rand2b;
    double proba, pij;

    for (i = 0; i < 2 * total; i++)
    {
        rand1 = rand() % total;
        rand2 = rand() % total;

        rand1b = vecOrder[rand1];
        rand2b = vecOrder[rand2];

        pij = 0.5 + (vecBestGain[rand2b] - vecBestGain[rand1b]) * (*inv2R);

        proba = (double)rand() / (RAND_MAX + 1.0);

        if (pij >= proba)
        {
            vecOrder[rand1] = rand2b;
            vecOrder[rand2] = rand1b;
        }
    }
}

void action4d(
    int k,
    int nrowData, int ncolData, int ncondData, int ntimeData,
    double *Data,
    int *vecOrder, int *vecBestBic,
    int *bicRow, int *bicCol, int *bicCond, int *bicTime,
    int *bicRow2, int *bicCol2, int *bicCond2, int *bicTime2,
    double r, int *valbreak,
    double *vecResvolBic,
    double *sumBic,
    double *sumRow, double *sumCol, double *sumCond, double *sumTime,
    double *sumBic2,
    double *sumRow2, double *sumCol2, double *sumCond2, double *sumTime2,
    int N, int M, int P, int Q,
    int zz, int *aa,
    int *vecBlocGene, int *vecBlocSample,
    int *vecBlocCond, int *vecBlocTime)
{
    int x, b, z;

    int total = nrowData + ncolData + ncondData + ntimeData;

    for (x = 0; x < total; x++)
    {
        b = vecOrder[x];
        z = vecBestBic[b];

        double res1 = vecResvolBic[z * 6];

        int nrowBic  = vecResvolBic[z * 6 + 2];
        int ncolBic  = vecResvolBic[z * 6 + 3];
        int ncondBic = vecResvolBic[z * 6 + 4];
        int ntimeBic = vecResvolBic[z * 6 + 5];

        int volBic = nrowBic * ncolBic * ncondBic * ntimeBic;

        /* determine dimension */
        int changed = 0;

        if (b < nrowData)
        {
            int i = b;
            if (!vecBlocGene[nrowData * z + i])
            {
                bicRow2[nrowData * z + i] = 1 - bicRow[nrowData * z + i];
                changed = 1;
            }
        }
        else if (b < nrowData + ncolData)
        {
            int j = b - nrowData;
            if (!vecBlocSample[ncolData * z + j])
            {
                bicCol2[ncolData * z + j] = 1 - bicCol[ncolData * z + j];
                changed = 1;
            }
        }
        else if (b < nrowData + ncolData + ncondData)
        {
            int k2 = b - nrowData - ncolData;
            if (!vecBlocCond[ncondData * z + k2])
            {
                bicCond2[ncondData * z + k2] = 1 - bicCond[ncondData * z + k2];
                changed = 1;
            }
        }
        else
        {
            int l = b - nrowData - ncolData - ncondData;
            if (!vecBlocTime[ntimeData * z + l])
            {
                bicTime2[ntimeData * z + l] = 1 - bicTime[ntimeData * z + l];
                changed = 1;
            }
        }

        if (!changed) continue;

        double res2 = residu4d(z,
            nrowData, ncolData, ncondData, ntimeData,
            Data,
            bicRow2, bicCol2, bicCond2, bicTime2,
            sumBic, sumRow, sumCol, sumCond, sumTime);

        if (res2 < res1)
        {
            (*aa)++;
            *valbreak = 1;

            vecResvolBic[z * 6] = res2;

            memcpy(bicRow,  bicRow2,  k * nrowData  * sizeof(int));
            memcpy(bicCol,  bicCol2,  k * ncolData  * sizeof(int));
            memcpy(bicCond, bicCond2, k * ncondData * sizeof(int));
            memcpy(bicTime, bicTime2, k * ntimeData * sizeof(int));
        }
    }
}

void floc(double *Data, int *nrowData, int *ncolData, int *ncondData, int *ntimeData, int *bicRow, int *bicCol, 
		  int *bicCond, int *bicTime, double *vecResvolBic, double *r, int *k, int *N, int *M, int *P, int *Q, 
		  int *t, int *vecBlocGene, int *vecBlocSample, int *vecBlocCond, int *vecBlocTime)
{
	
	srand(time(NULL));
	
    int zz, i, valbreak,j;
	
    int total = *nrowData + *ncolData + *ncondData + *ntimeData;
	
    double inv2R;
	
    int *bicRow2  = malloc(*k * *nrowData  * sizeof(int));
    int *bicCol2  = malloc(*k * *ncolData  * sizeof(int));
	int *bicCond2 = malloc(*k * *ncondData * sizeof(int));
	int *bicTime2 = malloc(*k * *ntimeData * sizeof(int));
	
    int *vecOrder = malloc(total * sizeof(int));
    double *vecBestGain = malloc(total * sizeof(double));
    int *vecBestBic = malloc(total * sizeof(int));
	
    double *sumRow  = malloc(*k * *nrowData  * sizeof(double));
    double *sumCol  = malloc(*k * *ncolData  * sizeof(double));
	double *sumCond = malloc(*k * *ncondData * sizeof(double));
	double *sumTime = malloc(*k * *ntimeData * sizeof(double));
	
    double *sumBic = malloc(*k * sizeof(double));
	
    double *sumRow2  = malloc(*k * *nrowData  * sizeof(double));
    double *sumCol2  = malloc(*k * *ncolData  * sizeof(double));
	double *sumCond2 = malloc(*k * *ncondData * sizeof(double));
	double *sumTime2 = malloc(*k * *ntimeData * sizeof(double));
	
    double *sumBic2 = malloc(*k * sizeof(double));
    
    double invk = 1/(double)*k;
    
    memcpy(bicRow2, bicRow,   *k * *nrowData  * sizeof(int));
    memcpy(bicCol2, bicCol,   *k * *ncolData  * sizeof(int));
	memcpy(bicCond2, bicCond, *k * *ncondData * sizeof(int));
	memcpy(bicTime2, bicTime, *k * *ntimeData * sizeof(int));

    for (zz = 0 ; zz < *k ; zz++)
	{
		int nrowBic  = count_row_col(zz, *nrowData,  bicRow);
        int ncolBic  = count_row_col(zz, *ncolData,  bicCol);
        int ncondBic = count_row_col(zz, *ncondData, bicCond);
        int ntimeBic = count_row_col(zz, *ntimeData, bicTime);

        vecResvolBic[zz * 6 + 2] = nrowBic;
        vecResvolBic[zz * 6 + 3] = ncolBic;
        vecResvolBic[zz * 6 + 4] = ncondBic;
        vecResvolBic[zz * 6 + 5] = ntimeBic;
		
	    sum4d(zz, *nrowData, *ncolData, *ncondData, *ntimeData,
              Data, bicRow, bicCol, bicCond, bicTime, sumBic, 
			  sumRow, sumCol, sumCond, sumTime);

        sum4d(zz, *nrowData, *ncolData, *ncondData, *ntimeData,
              Data, bicRow, bicCol, bicCond, bicTime, sumBic2, 
			  sumRow2, sumCol2, sumCond2, sumTime2);

        vecResvolBic[zz * 6] = residu4d(zz, *nrowData, *ncolData, *ncondData, *ntimeData,
										Data, bicRow, bicCol, bicCond, bicTime, sumBic, 
										sumRow, sumCol, sumCond, sumTime);

        vecResvolBic[zz * 6 + 1] = nrowBic * ncolBic * ncondBic * ntimeBic;
	}
	
    j = 0;
	
    for (zz = 0 ; zz < *t ; zz++)
	{
	    valbreak = 0;
	    
		bestgain4d(
            *k, *r,
            *nrowData, *ncolData, *ncondData, *ntimeData,
            Data,
            bicRow, bicCol, bicCond, bicTime,
            bicRow2, bicCol2, bicCond2, bicTime2,
            sumBic, sumRow, sumCol, sumCond, sumTime,
            sumBic2, sumRow2, sumCol2, sumCond2, sumTime2,
            vecBestGain, vecBestBic,
            &inv2R,
            vecResvolBic,
            *N, *M, *P, *Q,
            vecBlocGene, vecBlocSample, vecBlocCond, vecBlocTime
        );
        
	    for (i = 0; i < total; i++)
            vecOrder[i] = i;

        tri(vecBestGain, vecOrder, 0, total - 1);

        order4d(&inv2R, total, vecBestGain, vecOrder);

        action4d(
            *k,
            *nrowData, *ncolData, *ncondData, *ntimeData,
            Data,
            vecOrder, vecBestBic,
            bicRow, bicCol, bicCond, bicTime,
            bicRow2, bicCol2, bicCond2, bicTime2,
            *r, &valbreak,
            vecResvolBic,
            sumBic, sumRow, sumCol, sumCond, sumTime,
            sumBic2, sumRow2, sumCol2, sumCond2, sumTime2,
            *N, *M, *P, *Q,
            zz, &j,
            vecBlocGene, vecBlocSample, vecBlocCond, vecBlocTime
        );

        if (valbreak == 0)
        {
            printf("\n STOP\n ");
            break;
        }
	}
    
    free(bicRow2);
    free(bicCol2);
    free(bicCond2);
    free(bicTime2);

    free(vecOrder);
    free(vecBestGain);
    free(vecBestBic);

    free(sumRow);
    free(sumCol);
    free(sumCond);
    free(sumTime);

    free(sumBic);

    free(sumRow2);
    free(sumCol2);
    free(sumCond2);
    free(sumTime2);

    free(sumBic2);
    
}