## Szükséges könyvtárak betöltése (numpy, pandas -> adatszerkezetek; matplotlib -> ábrázolás; 
## ctypes -> C kód futtatása; scipy, statsmodels -> statisztikai tesztek)

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from datetime import datetime
import ctypes

## C kód beolvasása (először dll formátumra kell alakítani)
clib = ctypes.CDLL('C/bicare.dll')

## C-kompatibilis típusok beállítása a paraméterekre
clib.floc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
]

## C függvény visszatérési típusának beállítása
clib.floc.restype = None

#%%

## Eredmények értelmezésére és ábrázolására használt osztályok
class ExpressionSet:
    def __init__(self, exprs, pheno_data=None, feature_data=None, annotation=None):
        self.exprs = pd.DataFrame(exprs)
        self.pheno_data = pd.DataFrame(pheno_data) if pheno_data is not None else None
        self.feature_data = pd.DataFrame(feature_data) if feature_data is not None else None
        self.annotation = annotation

class biclustering:
    def __init__(self, Call, eset, param, bicRow, bicCol, mat_resvol_bic):
        self.Call = Call
        self.eset = eset
        self.param = param
        self.bicRow = bicRow
        self.bicCol = bicCol
        self.mat_resvol_bic = mat_resvol_bic
    
    def __repr__(self):
        out = []
        out.append("Biclustering Result\n")
        out.append("Call:\n")
        out.append(f"{self.Call}\n")

        out.append("\nParameters:\n")
        for name, value in self.param:
            out.append(f"  {name}: {value}")

        out.append("\n\nBiclusters summary:\n")
        header = ["Residue", "Volume", "Genes", "Conditions", "RowVar"]
        out.append("  " + "  ".join(header))

        for row in self.mat_resvol_bic:
            out.append("  " + "  ".join(f"{x:.4f}" for x in row))

        return "\n".join(out)

    def summary(self):
        return self.mat_resvol_bic


## Adatmátrix reziduális értékét kiszámító függvény
def residue(data):
    row_means = data.mean(axis=1, keepdims=True)
    col_means = data.mean(axis=0, keepdims=True)
    overall_mean = data.mean()
    
    resid = data - row_means - col_means + overall_mean
    return np.mean(resid ** 2)


## Adatmátrix reziduális értékét kiszámító függvény
def compute_residue_submatrix(data):
    row_means = data.mean(axis=1, keepdims=True)
    col_means = data.mean(axis=0, keepdims=True)
    overall_mean = data.mean()
    return np.mean((data - row_means - col_means + overall_mean) ** 2)


## Központi biklaszterező függvény, amely meghívja a C kódot a beállított paraméterek szerint, 
## majd tárolja az eredményeket
def FLOC(eset, k=15, pGene=0.3, pSample=0.6, 
         r=None, N=10, M=8, t=200, random_state=None,
         blocGene=None, blocSample=None):
    
    
    ## Seed generálása reprodukálhatóság céljából
    rng = np.random.default_rng(random_state)
    
    if pSample is None:
        pSample = pGene
    
    ## Génexpressziós értékek kinyerése
    data = np.ascontiguousarray(eset.exprs, dtype=np.float64)

    n_genes, n_samples = data.shape

    if r is None:
        r = residue(data) / 10

    ## Bemeneti vektor C-kompatibilissé alakítása lapítással
    vecData = data.T.flatten(order="F").astype(np.float64)

    vecBicRow = np.zeros((n_genes, k), dtype=np.int32)
    vecBicCol = np.zeros((n_samples, k), dtype=np.int32)

    ## Gének és minták inicializálásának ellenőrzése
    if blocGene is not None:
        k = max(k, blocGene.shape[1])
        vecBicRow[:, :blocGene.shape[1]] = blocGene

    if blocSample is not None:
        k = max(k, blocSample.shape[1])
        vecBicCol[:, :blocSample.shape[1]] = blocSample

    rand1 = rng.random(k * n_genes)
    rand2 = rng.random(k * n_samples)
    
    ## Bemeneti vektorok C-kompatibilissé alakítása lapítással
    vecBicRow = vecBicRow.flatten(order="F")
    vecBicCol = vecBicCol.flatten(order="F")
    
    vecBlocGene = vecBicRow.copy()
    vecBlocSample = vecBicCol.copy()
    
    vecBicRow[rand1 < pGene] = 1
    vecBicCol[rand2 < pSample] = 1
    
    vec_resvol_bic = np.zeros(k * 4, dtype=np.float64)
    
    ## Paraméterek C-kompatibilissé alakítása
    n_genes_c = ctypes.c_int(n_genes)
    n_samples_c = ctypes.c_int(n_samples)
    r_c = ctypes.c_double(r)
    k_c = ctypes.c_int(k)
    N_c = ctypes.c_int(N)
    M_c = ctypes.c_int(M)
    t_c = ctypes.c_int(t)
    
    ## C függvény meghívása
    clib.floc(
        vecData,
        ctypes.byref(n_genes_c),
        ctypes.byref(n_samples_c),
        vecBicRow,
        vecBicCol,
        vec_resvol_bic,
        ctypes.byref(r_c),
        ctypes.byref(k_c),
        ctypes.byref(N_c),
        ctypes.byref(M_c),
        ctypes.byref(t_c),
        vecBlocGene,
        vecBlocSample
    )
    
    ## C kód által módosított vektorok visszaalakítása eredeti struktúrára
    bicRow = vecBicRow.reshape((n_genes, k), order="F").T
    bicCol = vecBicCol.reshape((n_samples, k), order="F").T

    mat_resvol_bic = np.zeros((k, 5))
    mat_resvol_bic[:, :4] = vec_resvol_bic.reshape((k, 4), order="C")

    for iteration in range(t):
        for i in range(k):
            rows = np.where(bicRow[i] == 1)[0]
            cols = np.where(bicCol[i] == 1)[0]
            
            if len(rows) > 0 and len(cols) > 0:
                sub = data[np.ix_(rows, cols)]
                rowvars = np.var(sub, axis=1, ddof=1)
                mat_resvol_bic[i, 4] = np.mean(rowvars)
    
    param = [
        ("Number of biclusters", k),
        ("Residue threshold", r),
        ("Genes initial probability", pGene),
        ("Samples initial probability", pSample),
        ("Number of iterations", t),
        ("Date", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
    ]

    Call = "floc_wrapper(...)"

    ## Osztály példányosítás
    return biclustering(
        Call=Call,
        eset=eset,
        param=param,
        bicRow=bicRow,
        bicCol=bicCol,
        mat_resvol_bic=mat_resvol_bic
    )


## Az eredményekből tetszőleges biklasztert kinyerő függvény
def bicluster(res, data, k):
    bicRow = res.bicRow
    bicCol = res.bicCol

    rows = np.where(bicRow[k] == 1)[0]
    cols = np.where(bicCol[k] == 1)[0]

    if len(rows) == 0 or len(cols) == 0:
        return np.array([])

    return res.eset.exprs.values[np.ix_(rows, cols)]


## Minta<->Kovariáns teszt
def testAnnot(res, covariates="all"):
    
    pheno = res.eset.pheno_data
    
    if covariates == "all":
        covariates = pheno.columns
    
    pvalues = pd.DataFrame(index=range(res.param[0][1]), columns=covariates)
    
    for cov in covariates:
        full_counts = pheno[cov].value_counts()
        
        for j in range(res.param[0][1]):
            selected = res.bicCol[j]
            sample_mask = selected == 1
            samples = res.eset.exprs.columns[sample_mask]
            
            if len(samples) == 0:
                continue
            
            subset_counts = pheno.loc[samples, cov].value_counts()
            obs = subset_counts.reindex(full_counts.index, fill_value=0)
            
            chi2, p, _, _ = chi2_contingency(
                [obs.values, full_counts.values],
                correction=False
            )
            pvalues.loc[j, cov] = p
    
    adj = pvalues.copy()
    for cov in covariates:
        mask = pvalues[cov].notna()
        adj.loc[mask, cov] = multipletests(
            pvalues.loc[mask, cov],
            method="fdr_by"
        )[1]
    
    res.covar = {
        "pvalues": pvalues,
        "adjpvalues": adj
    }
    
    return res


## Génkészlet felülreprezentációs teszt
def testSet(resBic, gene_sets):

    expr = resBic.eset.exprs
    all_genes = set(expr.index)

    nb_bic = resBic.bicRow.shape[0]
    nb_sets = len(gene_sets)

    pvals = np.full((nb_bic, nb_sets), np.nan)
    inter = np.zeros((nb_bic, nb_sets))

    gene_sets_names = list(gene_sets.keys())

    for i, (set_name, gene_set) in enumerate(gene_sets.items()):

        gene_set = set(gene_set) & all_genes
        set_size = len(gene_set)

        for j in range(nb_bic):

            bic_genes = set(expr.index[np.where(resBic.bicRow[j] == 1)[0]])
            bic_size = len(bic_genes)

            overlap = len(bic_genes & gene_set)

            inter[j, i] = overlap

            if overlap > 0:

                pvals[j, i] = hypergeom.sf(
                    overlap - 1,
                    len(all_genes),
                    set_size,
                    bic_size
                )

    return pd.DataFrame(
        pvals,
        columns=gene_sets_names
    )


## Biklaszterek ábrázolására használt függvény
def plot_bicluster(res, data, k):

    if res.size == 0:
        print("Empty bicluster")
        return
    res_t = res.T

    plt.figure()
    plt.plot(res_t)
    plt.xlabel("Conditions")
    plt.ylabel("Expression Level")
    plt.title(f"Bicluster {k}")
    plt.show()

#%%

## Bemeneti adatok beolvasása és létrehozása
exprs = pd.read_csv("exprs.csv", index_col=0)
pheno = pd.read_csv("pheno.csv", index_col=0)

sample_bicData = ExpressionSet(
    exprs=exprs,
    pheno_data=pheno,
    annotation="hgu95av2"
)

#%%

## Biklaszterek meghatározása, kinyerése és ábrázolása
res = FLOC(sample_bicData,
           k=15,
           pGene=0.3,
           pSample=0.6,
           r=0.01,
           t=200,
           random_state=1)

print(res)

bic6 = bicluster(res, sample_bicData.exprs.values, 4)

plot_bicluster(bic6, sample_bicData.exprs.values, 4)

#%%

## Gén- és minta inicializálás
init_genes = np.zeros((352, 5), dtype=int)
init_samples = np.zeros((26, 5), dtype=int)
init_genes[0:10, 0] = 1
init_genes[19:31, 1] = 1
init_genes[49:61, 2] = 1
init_samples[0:5, 2] = 1
init_samples[0:5, 3] = 1
init_samples[9:15, 4] = 1

res2 = FLOC(sample_bicData, 
            k=15, 
            pGene=0.3, 
            pSample=0.6, 
            r=0.01, 
            t=200, 
            random_state=1, 
            blocGene=init_genes,
            blocSample=init_samples)

print(res2)

bic5 = bicluster(res2, sample_bicData.exprs.values, 4)

plot_bicluster(bic5, sample_bicData.exprs.values, 4)

#%%

## Génkészlet felülreprezentációs teszt elvégzése
expr = sample_bicData.exprs

gene_sets = {
    "GO_set1": list(expr.index[:20]),
    "GO_set2": list(expr.index[30:50]),
}

res_test = testSet(res, gene_sets)
print(res_test)

#%%

## Minta<->Kovariáns teszt elvégzése
resTA = testAnnot(res, covariates=["sex", "type"])

print(resTA.covar)
