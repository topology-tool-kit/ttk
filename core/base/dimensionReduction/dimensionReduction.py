def doIt(X, method, ncomponents, nneighbors, njobs, rstate, params):

    import importlib

    # check if numpy is installed
    loader = importlib.find_loader('numpy')
    found = loader is not None
    if found:
        print("[DimensionReduction] Python: numpy module found.")
    else:
        print("[DimensionReduction] Python error: numpy module not found.")
        return 0

    # check if scipy is installed
    loader = importlib.find_loader('scipy')
    found = loader is not None
    if found:
        print("[DimensionReduction] Python: scipy module found.")
    else:
        print("[DimensionReduction] Python error: scipy module not found.")
        return 0

    # check if scikit-learn is installed
    loader = importlib.find_loader('sklearn')
    found = loader is not None
    if found:
        print("[DimensionReduction] Python: sklearn module found.")
    else:
        print("[DimensionReduction] Python error: sklearn module not found.")
        return 0

    from sklearn import manifold
    from sklearn import decomposition
    import numpy as np
    from sys import platform
    
    if platform == "darwin":
        import sklearn
        sklearn.utils.parallel_backend('threading')

    if rstate > 0:
        np.random.seed(0)

    try:
        if method == 0:
            seParams = params[0]
            if seParams[2] == "None":
                se = manifold.SpectralEmbedding(n_components=ncomponents,
                                                affinity=seParams[0],
                                                gamma=seParams[1],
                                                eigen_solver=None,
                                                n_neighbors=nneighbors,
                                                n_jobs=njobs)
            else:
                se = manifold.SpectralEmbedding(n_components=ncomponents,
                                                affinity=seParams[0],
                                                gamma=seParams[1],
                                                eigen_solver=seParams[2],
                                                n_neighbors=nneighbors,
                                                n_jobs=njobs)
            Y = se.fit_transform(X)
        elif method == 1:
            lleParams = params[1]
            lle = manifold.LocallyLinearEmbedding(n_neighbors=nneighbors,
                                                  n_components=ncomponents,
                                                  reg=lleParams[0],
                                                  eigen_solver=lleParams[1],
                                                  tol=lleParams[2],
                                                  max_iter=lleParams[3],
                                                  method=lleParams[4],
                                                  hessian_tol=lleParams[5],
                                                  modified_tol=lleParams[6],
                                                  neighbors_algorithm=lleParams[7],
                                                  n_jobs=njobs)
            Y = lle.fit_transform(X)
        elif method == 2:
            mdsParams = params[2]
            mds = manifold.MDS(n_components=ncomponents,
                               metric=mdsParams[0],
                               n_init=mdsParams[1],
                               max_iter=mdsParams[2],
                               verbose=mdsParams[3],
                               eps=mdsParams[4],
                               dissimilarity=mdsParams[5],
                               n_jobs=njobs)
            Y = mds.fit_transform(X)
        elif method == 3:
            tsneParams = params[3]
            tsne = manifold.TSNE(n_components=ncomponents,
                                 perplexity=tsneParams[0],
                                 early_exaggeration=tsneParams[1],
                                 learning_rate=tsneParams[2],
                                 n_iter=tsneParams[3],
                                 n_iter_without_progress=tsneParams[4],
                                 min_grad_norm=tsneParams[5],
                                 metric=tsneParams[6],
                                 init=tsneParams[7],
                                 verbose=tsneParams[8],
                                 method=tsneParams[9],
                                 angle=tsneParams[10])
            Y = tsne.fit_transform(X)
        elif method == 4:
            isoParams = params[4]
            iso = manifold.Isomap(n_neighbors=nneighbors,
                                  n_components=ncomponents,
                                  eigen_solver=isoParams[0],
                                  tol=isoParams[1],
                                  max_iter=isoParams[2],
                                  path_method=isoParams[3],
                                  neighbors_algorithm=isoParams[4],
                                  n_jobs=njobs)
            Y = iso.fit_transform(X)
        elif method == 5:
            pcaParams = params[5]
            if pcaParams[4] == "auto":
                pca = decomposition.PCA(n_components=ncomponents,
                                        copy=pcaParams[0],
                                        whiten=pcaParams[1],
                                        svd_solver=pcaParams[2],
                                        tol=pcaParams[3],
                                        iterated_power="auto")
            else:
                pca = decomposition.PCA(n_components=ncomponents,
                                        copy=pcaParams[0],
                                        whiten=pcaParams[1],
                                        svd_solver=pcaParams[2],
                                        tol=pcaParams[3],
                                        iterated_power=int(pcaParams[4]))
            Y = pca.fit_transform(X)

        L = [Y.shape[0], Y.shape[1], np.ravel(Y, 'F')]
        return L
    except Exception as inst:
        print('[DimensionReduction] Error: unexpected behaviour detected in the python script.')
        print(type(inst))    # the exception instance
        print(inst.args)     # arguments stored in .args
        print(inst)
        return []
