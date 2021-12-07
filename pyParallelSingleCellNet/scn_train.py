import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.ensemble import RandomForestClassifier
import warnings
from .utils import *
from .tsp_rf import *

def randomize(expDat, num=50):
    temp=expDat.to_numpy()
    temp=np.array([np.random.choice(x, len(x), replace=False) for x in temp])
    temp=temp.T
    temp=np.array([np.random.choice(x, len(x), replace=False) for x in temp]).T
    return pd.DataFrame(data=temp, columns=expDat.columns).iloc[0:num,:]

def sc_trans_rnaseq(aDat,total = 10000 ):
    sc.pp.normalize_per_cell(aDat, counts_per_cell_after=total)
    sc.pp.log1p(aDat)
    sc.pp.scale(aDat, max_value=10)
    # return np.log(1+expCountDnW)
    return aDat

def sc_makeClassifier(expTrain, genes, groups, nRand=70, ntrees=2000, stratify=False, n_procs=1):
    randDat = randomize(expTrain, num=nRand)
    expT = pd.concat([expTrain, randDat])
    allgenes = expT.columns.values
    # missingGenes = np.setdiff1d(np.unique(genes), allgenes)
    ggenes = tuple(set(genes) & set(allgenes))
    clf_params = dict(
        n_estimators=ntrees,
        n_jobs=n_procs,
        random_state=100
    )
    if not stratify:
        clf = RandomForestClassifier(**clf_params)
    else:
        clf = RandomForestClassifier(class_weight="balanced", **clf_params)
    ggroups = np.append(np.array(groups), np.repeat("rand", nRand)).flatten()
    clf.fit(expT.loc[:, ggenes].to_numpy(), ggroups)
    return clf

def scn_train(
    adata,
    dLevel,
    nTopGenes=100,
    nTopGenePairs=100,
    nRand=100,
    nTrees=1000,
    stratify=False,
    counts_per_cell_after=1e4,
    scaleMax=10,
    limitToHVG=False,
    n_procs=1
    ):
    warnings.filterwarnings('ignore')
    aTrain = adata.copy()

    sc.pp.normalize_per_cell(aTrain, counts_per_cell_after=counts_per_cell_after)
    sc.pp.log1p(aTrain)
    print("Matrix normalized")
    if limitToHVG:
        sc.pp.highly_variable_genes(aTrain, min_mean=0.0125, max_mean=4, min_disp=0.5, subset=True)
        print("Subset HVG")
    sc.pp.scale(aTrain, max_value=scaleMax)
    print("Matrix scaled\n")

    aTrain = aTrain[adata.obs.index.values, :].copy()
    expTrain = pd.DataFrame(data=aTrain.X, index=aTrain.obs.index.values, columns=aTrain.var.index.values)
    del(aTrain)
    cgenesA, grps, cgenes_list = findClassyGenes(expTrain, adata.obs[dLevel], topX=nTopGenes)
    print(f"There are {len(cgenesA)} classification genes\n")

    xpairs = ptGetTop(expTrain.loc[:, cgenesA], grps, cgenes_list, topX=nTopGenePairs, sliceSize=5000, n_procs=n_procs)
    del(expTrain)
    print(f"There are {len(xpairs)} top gene pairs\n")
    
    aRaw = adata[adata.obs.index.values, :].copy()
    expRaw = pd.DataFrame(data=aRaw.X.toarray(), index=aRaw.obs.index.values, columns=aRaw.var.index.values)
    del(aRaw)
    pdTrain = query_transform(expRaw.loc[:, cgenesA], xpairs)
    del(expRaw)
    print("Finished pair transforming the data\n")
    tspRF = sc_makeClassifier(
        pdTrain.loc[:, xpairs],
        genes=xpairs,
        groups=grps,
        nRand=nRand,
        ntrees=nTrees,
        stratify=stratify,
        n_procs=n_procs
        )
    print("Classifier trained")
    return [cgenesA, xpairs, tspRF]

def scn_classify(adata, cgenes, xpairs, rf_tsp, nrand=0):
    classRes = scn_predict(cgenes, xpairs, rf_tsp, adata, nrand=nrand)
    categories = classRes.columns.values
    adNew = ad.AnnData(classRes, obs=adata.obs, var=pd.DataFrame(index=categories))
    adNew.obs['SCN_class'] = classRes.idxmax(axis=1)
    return adNew

def add_classRes(adata: AnnData, adClassRes, copy=False) -> AnnData:
    cNames = adClassRes.var_names
    for cname in cNames:
        adata.obs[cname] = adClassRes[:,cname].X
    # adata.obs['category'] = adClassRes.obs['category']
    adata.obs['SCN_class'] = adClassRes.obs['SCN_class']
    return adata if copy else None

def check_adX(adata: AnnData) -> AnnData:
    from scipy import sparse
    if(isinstance(adata.X, np.ndarray)):
        adata.X = sparse.csr_matrix(adata.X)

def scn_predict(cgenes, xpairs, rf_tsp, aDat, nrand=2):
    if isinstance(aDat.X, np.ndarray):
        # in the case of aDat.X is a numpy array
        aDat.X = ad._core.views.ArrayView(aDat.X)
    expDat = pd.DataFrame(data=aDat.X.toarray(), index=aDat.obs.index.values, columns=aDat.var.index.values)
    expValTrans = query_transform(expDat.reindex(labels=cgenes, axis='columns', fill_value=0), xpairs)
    classRes_val = rf_classPredict(rf_tsp, expValTrans, numRand=nrand)
    return classRes_val

def rf_classPredict(rfObj,expQuery,numRand=50):
    if numRand > 0 :
        randDat=randomize(expQuery, num=numRand)
        expQuery=pd.concat([expQuery, randDat])
    xpreds= pd.DataFrame(rfObj.predict_proba(expQuery), columns= rfObj.classes_, index=expQuery.index)
    return xpreds
