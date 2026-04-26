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
clib = ctypes.CDLL('C/bicare4d.dll')

## C-kompatibilis típusok beállítása a paraméterekre
clib.floc.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
]

## C függvény visszatérési típusának beállítása
clib.floc.restype = None

#%%

## Eredmények értelmezésére és ábrázolására használt osztály
class clustering:
    def __init__(self, Call, eset, param, bicRow, bicCol, bicCond, bicTime, mat_resvol_bic):
        self.Call = Call
        self.eset = eset
        self.param = param
        self.bicRow = bicRow
        self.bicCol = bicCol
        self.bicCond = bicCond
        self.bicTime = bicTime
        self.mat_resvol_bic = mat_resvol_bic
    
    def __repr__(self):
        out = []
        out.append("Clustering Result\n")
        out.append("Call:\n")
        out.append(f"{self.Call}\n")

        out.append("\nParameters:\n")
        for name, value in self.param:
            out.append(f"  {name}: {value}")

        out.append("\n\nBiclusters summary:\n")
        header = ["Residue", "Volume", "Genes", "Samples", "Conditions", "Times", "RowVar"]
        out.append("  " + "  ".join(header))

        for row in self.mat_resvol_bic:
            out.append("  " + "  ".join(f"{x:.4f}" for x in row))

        return "\n".join(out)

    def summary(self):
        return self.mat_resvol_bic


## Adatmátrix reziduális értékét kiszámító függvény
def residue_4d(X):

    mi = X.mean(axis=(1,2,3), keepdims=True)
    mj = X.mean(axis=(0,2,3), keepdims=True)
    mk = X.mean(axis=(0,1,3), keepdims=True)
    ml = X.mean(axis=(0,1,2), keepdims=True)

    mij = X.mean(axis=(2,3), keepdims=True)
    mik = X.mean(axis=(1,3), keepdims=True)
    mil = X.mean(axis=(1,2), keepdims=True)
    mjk = X.mean(axis=(0,3), keepdims=True)
    mjl = X.mean(axis=(0,2), keepdims=True)
    mkl = X.mean(axis=(0,1), keepdims=True)

    mijk = X.mean(axis=3, keepdims=True)
    mijl = X.mean(axis=2, keepdims=True)
    mikl = X.mean(axis=1, keepdims=True)
    mjkl = X.mean(axis=0, keepdims=True)

    overall = X.mean()

    resid = (X - mi - mj - mk - ml + mij + mik + mil + mjk + mjl + mkl - mijk - mijl - mikl - mjkl + overall)

    return np.mean(resid**2)


## Központi klaszterező függvény, amely meghívja a C kódot a beállított paraméterek szerint, 
## majd tárolja az eredményeket
def FLOC(eset, k=15, pGene=0.3, pSample=0.6, pCond=0.5, pTime=0.4, 
         r=None, N=10, M=8, P=2, Q=2, t=200, random_state=None,
         blocGene=None, blocSample=None, blocCond=None, blocTime=None):
    
    ## Seed generálása reprodukálhatóság céljából
    rng = np.random.default_rng(random_state)
    
    if pSample is None:
        pSample = pGene
    
    ## Génexpressziós értékek kinyerése
    data = np.ascontiguousarray(eset, dtype=np.float64)

    n_genes, n_samples, n_cond, n_time = data.shape

    if r is None:
        r = residue_4d(data) / 10

    ## Bemeneti vektor C-kompatibilissé alakítása lapítással
    vecData = data.T.flatten(order="F").astype(np.float64)

    vecBicRow = np.zeros((n_genes, k), dtype=np.int32)
    vecBicCol = np.zeros((n_samples, k), dtype=np.int32)
    vecBicCond = np.zeros((n_cond, k), dtype=np.int32)
    vecBicTime = np.zeros((n_time, k), dtype=np.int32)

    ## Gének és minták inicializálásának ellenőrzése
    if blocGene is not None:
        k = max(k, blocGene.shape[1])
        vecBicRow[:, :blocGene.shape[1]] = blocGene

    if blocSample is not None:
        k = max(k, blocSample.shape[1])
        vecBicCol[:, :blocSample.shape[1]] = blocSample

    if blocCond is not None:
        k = max(k, blocCond.shape[1])
        vecBicCond[:, :blocCond.shape[1]] = blocCond

    if blocTime is not None:
        k = max(k, blocTime.shape[1])
        vecBicTime[:, :blocTime.shape[1]] = blocTime

    rand1 = rng.random(k * n_genes)
    rand2 = rng.random(k * n_samples)
    rand3 = rng.random(k * n_cond)
    rand4 = rng.random(k * n_time)
    
    ## Bemeneti vektorok C-kompatibilissé alakítása lapítással
    vecBicRow = vecBicRow.flatten(order="F")
    vecBicCol = vecBicCol.flatten(order="F")
    vecBicCond = vecBicCond.flatten(order="F")
    vecBicTime = vecBicTime.flatten(order="F")
    
    vecBlocGene = vecBicRow.copy()
    vecBlocSample = vecBicCol.copy()
    vecBlocCond = vecBicCond.copy()
    vecBlocTime = vecBicTime.copy()
    
    vecBicRow[rand1 < pGene] = 1
    vecBicCol[rand2 < pSample] = 1
    vecBicCond[rand3 < pCond] = 1
    vecBicTime[rand4 < pTime] = 1
    
    vec_resvol_bic = np.zeros(k * 6, dtype=np.float64)
    
    ## Paraméterek C-kompatibilissé alakítása
    n_genes_c = ctypes.c_int(n_genes)
    n_samples_c = ctypes.c_int(n_samples)
    n_cond_c = ctypes.c_int(n_cond)
    n_time_c = ctypes.c_int(n_time)
    r_c = ctypes.c_double(r)
    k_c = ctypes.c_int(k)
    N_c = ctypes.c_int(N)
    M_c = ctypes.c_int(M)
    P_c = ctypes.c_int(P)
    Q_c = ctypes.c_int(Q)
    t_c = ctypes.c_int(t)
    
    ## C függvény meghívása
    clib.floc(
        vecData,
        ctypes.byref(n_genes_c),
        ctypes.byref(n_samples_c),
        ctypes.byref(n_cond_c),
        ctypes.byref(n_time_c),
        vecBicRow,
        vecBicCol,
        vecBicCond,
        vecBicTime,
        vec_resvol_bic,
        ctypes.byref(r_c),
        ctypes.byref(k_c),
        ctypes.byref(N_c),
        ctypes.byref(M_c),
        ctypes.byref(P_c),
        ctypes.byref(Q_c),
        ctypes.byref(t_c),
        vecBlocGene,
        vecBlocSample,
        vecBlocCond,
        vecBlocTime
    )
    
    ## C kód által módosított vektorok visszaalakítása eredeti struktúrára
    bicRow = vecBicRow.reshape((n_genes, k), order="F").T
    bicCol = vecBicCol.reshape((n_samples, k), order="F").T
    bicCond = vecBicCond.reshape((n_cond, k), order="F").T
    bicTime = vecBicTime.reshape((n_time, k), order="F").T
    
    mat_resvol_bic = np.zeros((k, 7))
    mat_resvol_bic[:, :6] = vec_resvol_bic.reshape((k, 6), order="C")

    for iteration in range(t):
        for i in range(k):
            rows = np.where(bicRow[i] == 1)[0]
            cols = np.where(bicCol[i] == 1)[0]
            conds = np.where(bicCond[i] == 1)[0]
            times = np.where(bicTime[i] == 1)[0]
            
            if len(rows) > 0 and len(cols) > 0 and len(conds) > 0 and len(times) > 0:
                sub = data[np.ix_(rows, cols, conds, times)]
                rowvars = np.var(sub, axis=1, ddof=1)
                mat_resvol_bic[i, 6] = np.mean(rowvars)
    
    param = [
        ("Number of biclusters", k),
        ("Residue threshold", r),
        ("Genes initial probability", pGene),
        ("Samples initial probability", pSample),
        ("Conditions initial probability", pCond),
        ("Times initial probability", pTime),
        ("Number of iterations", t),
        ("Date", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
    ]

    Call = "floc_wrapper(...)"

    ## Osztály példányosítás
    return clustering(
        Call=Call,
        eset=eset,
        param=param,
        bicRow=bicRow,
        bicCol=bicCol,
        bicCond=bicCond,
        bicTime=bicTime,
        mat_resvol_bic=mat_resvol_bic
    )


## Az eredményekből tetszőleges klasztert kinyerő függvény
def bicluster_4d(res, k):
    rows = np.where(res.bicRow[k] == 1)[0]
    cols = np.where(res.bicCol[k] == 1)[0]
    conds = np.where(res.bicCond[k] == 1)[0]
    times = np.where(res.bicTime[k] == 1)[0]

    return rows, cols, conds, times

#%%

## Kezdeti négydimenziós adattér létrehozása
def generate_base_4d(n_genes, n_samples, n_cond, n_time, seed=0):
    rng = np.random.default_rng(seed)

    data = rng.normal(
        loc=7.16,
        scale=1.5,
        size=(n_genes, n_samples, n_cond, n_time)
    )

    data = np.clip(data, 1.95, 12.6)

    return data


## Kezdeti klaszterek hozzáadása a szintetikus adatokhoz
def add_quadcluster(data,
                   gene_idx,
                   sample_idx,
                   cond_idx,
                   time_idx,
                   pattern="constant",
                   strength=2.0):
    
    sub = data[np.ix_(gene_idx, sample_idx, cond_idx, time_idx)]
    
    if pattern == "constant":
        sub += strength
    
    elif pattern == "additive":
        g_effect = np.linspace(-1, 1, len(gene_idx))[:, None, None, None]
        s_effect = np.linspace(-1, 1, len(sample_idx))[None, :, None, None]
        sub += strength * (g_effect + s_effect)
    
    elif pattern == "multiplicative":
        sub *= (1 + 0.2 * strength)
    
    data[np.ix_(gene_idx, sample_idx, cond_idx, time_idx)] = sub
    
    return data


## Fő szintetikus adatot generáló függvény
def generate_synthetic_4d(
    n_genes=20,
    n_samples=20,
    n_cond=20,
    n_time=20,
    n_clusters=5,
    seed=1
):
    rng = np.random.default_rng(seed)

    data = generate_base_4d(n_genes, n_samples, n_cond, n_time, seed)

    clusters = []

    for _ in range(n_clusters):

        genes = rng.choice(n_genes, size=rng.integers(5, 15), replace=False)
        samples = rng.choice(n_samples, size=rng.integers(5, 15), replace=False)
        conds = rng.choice(n_cond, size=rng.integers(1, n_cond), replace=False)
        times = rng.choice(n_time, size=rng.integers(1, n_time), replace=False)

        pattern = rng.choice(["constant", "additive", "multiplicative"])

        data = add_quadcluster(
            data,
            genes,
            samples,
            conds,
            times,
            pattern=pattern,
            strength=rng.uniform(1.5, 3.0)
        )

        clusters.append({
            "genes": genes,
            "samples": samples,
            "conditions": conds,
            "times": times,
            "pattern": pattern
        })

    return data, clusters

## Klaszterek ábrázolására használt függvény
def plot_bicluster_4d(res, data4d, k, cond_idx=0, time_idx=0):

    rows, cols, conds, times = bicluster_4d(res, k)

    if len(rows) == 0 or len(cols) == 0:
        print("Empty bicluster")
        return

    if cond_idx >= len(conds) or time_idx >= len(times):
        print("Invalid condition/time index")
        return

    c = conds[cond_idx]
    t = times[time_idx]

    sub = data4d[np.ix_(rows, cols, [c], [t])].squeeze()

    res_t = sub.T

    plt.figure()
    plt.plot(res_t)
    plt.xlabel("Samples")
    plt.ylabel("Expression Level")
    plt.title(f"Bicluster {k} (Cond={c}, Time={t})")
    plt.show()
    

def plot_bicluster(res, k):

    rows = np.where(res.bicRow[k] == 1)[0]
    cols = np.where(res.bicCol[k] == 1)[0]

    if len(rows) == 0 or len(cols) == 0:
        print("Empty bicluster")
        return

    data = res.eset

    sub = data.iloc[rows, cols]
    res_t = sub.T

    plt.figure()
    for gene in res_t.columns:
        plt.plot(res_t.index, res_t[gene])

    plt.xlabel("Samples")
    plt.ylabel("Expression Level")
    plt.title(f"Bicluster {k}")
    plt.xticks(rotation=45)
    plt.show()

#%%

## Szintetikus adathalmaz generálása
data4d, true_clusters = generate_synthetic_4d()

## Biklaszterek meghatározása, kinyerése és ábrázolása
res = FLOC(data4d,
           k=10,
           pGene=0.3,
           pSample=0.4,
           r=0.01,
           t=500,
           random_state=1)

print(res)

plot_bicluster(res, 2)

plot_bicluster_4d(res, data4d, 2, 0, 2)

plot_bicluster_4d(res, data4d, 2, 1, 3)

#%%

## Gén- és minta inicializálás
init_genes = np.zeros((20, 5), dtype=int)
init_samples = np.zeros((20, 5), dtype=int)
init_genes[0:10, 0] = 1
init_genes[19:31, 1] = 1
init_genes[49:61, 2] = 1
init_samples[0:5, 2] = 1
init_samples[0:5, 3] = 1
init_samples[9:15, 4] = 1

res2 = FLOC(data4d, 
            k=15, 
            pGene=0.3, 
            pSample=0.6, 
            r=0.01, 
            t=200, 
            random_state=1, 
            blocGene=init_genes,
            blocSample=init_samples)

print(res2)

plot_bicluster_4d(res, data4d, 2, 0, 2)

plot_bicluster_4d(res, data4d, 2, 1, 3)
