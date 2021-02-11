def grouped_obs_mean(adata, group_key, mean_func=None, layer=None, gene_symbols=None):
    # https://github.com/theislab/scanpy/issues/181#issuecomment-534867254
    import pandas as pd
    import numpy as np
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )


    if mean_func is None:
        mean_func = np.mean

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        if len(X.shape)<2:
            X = np.expand_dims(X,0)
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out
