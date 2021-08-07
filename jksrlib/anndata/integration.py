import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import scanorama
import harmonypy
import mnnpy

def var_subset_selection(adatas, var_subset=None):
    comm_vars = adatas[0].var_names
    for adata in adatas[1:]:
        comm_vars = comm_vars[comm_vars.isin(adata.var_names)]
    if var_subset is not None:
        comm_vars = comm_vars[comm_vars.isin(var_subset)]
    
    return [ adata[:,comm_vars].copy() for adata in adatas]


def integrate_scanorama(*adatas, var_subset=None):
    tmpadatas = var_subset_selection(adatas, var_subset)
    integrated = scanorama.integrate_scanpy(tmpadatas)

    rtns = []
    for i in range(len(integrated)):
        iadata = anndata.AnnData(integrated[i])
        iadata.obs = tmpadatas[i].obs
        rtns.append(iadata)
    return rtns


def integrate_harmony(*adatas, var_subset=None):
    tmpadatas = var_subset_selection(adatas, var_subset)

    meta = []
    pca = []
    for i in range(len(tmpadatas)):
        tmpadatas[i].obs['tmptmp'] = f'tmptmp_{i}'
        meta.append(tmpadatas[i].obs)

        sc.tl.pca(tmpadatas[i])
        pca.append(tmpadatas[i].obsm['X_pca'])

    meta = pd.concat(meta)
    df = pd.DataFrame(np.vstack(pca))
    df.index = meta.index

    ho = harmonypy.run_harmony(df, meta, ['tmptmp'])
    iadatas = anndata.AnnData(ho.Z_corr.T)

    rtns = []
    start = 0
    for i in range(len(adatas)):
        end = start+len(adatas[i])
        rtns.append(iadatas[start:end])
        start = end

    for i in range(len(adatas)):
        rtns[i].obs=adatas[i].obs

    return rtns


def integrate_mnn(*adatas, var_subset=None):
    tmpadatas = var_subset_selection(adatas, var_subset)
    
    corrected = mnnpy.mnn_correct(tmpadatas)
    corrected = corrected[0]
    
    iadatas = anndata.AnnData(np.vstack([x.X for x in corrected]))

    rtns = []
    #start = 0
    #for i in range(len(adatas)):
    #    end = start+len(adatas[i])
    #    rtns.append(iadatas[start:end])
    #    start = end
    #    
    #for i in range(len(adatas)):
    #    rtns[i].obs=adatas[i].obs

    rtns = [ iadata.copy() for iadata in iadatas ]
    return rtns

