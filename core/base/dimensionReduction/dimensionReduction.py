def doIt(X, method, ncomponents, nneighbors, njobs, rstate):
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

    from sklearn import manifold
    from sklearn import decomposition
    import numpy as np

    if rstate > 0:
        np.random.seed(rstate)

    if method == 0:
        se = manifold.SpectralEmbedding(n_components=ncomponents, n_neighbors=nneighbors, n_jobs=njobs)
        Y = se.fit_transform(X)
    elif method == 1:
        lle = manifold.LocallyLinearEmbedding(n_components=ncomponents, n_neighbors=nneighbors, eigen_solver='auto', method='standard', n_jobs=njobs)
        Y = lle.fit_transform(X)
    elif method == 2:
        mds = manifold.MDS(n_components=ncomponents, max_iter=100, n_init=1, n_jobs=njobs)
        Y = mds.fit_transform(X)
    elif method == 3:
        tsne = manifold.TSNE(n_components=ncomponents, init='pca')
        Y = tsne.fit_transform(X)
    elif method == 4:
        iso = manifold.Isomap(n_components=ncomponents, n_neighbors=nneighbors, n_jobs=njobs)
        Y = iso.fit_transform(X)
    elif method == 5:
        pca = decomposition.PCA(n_components=ncomponents)
        Y = pca.fit_transform(X)

    L = [Y.shape[0], Y.shape[1], np.ravel(Y, 'F')]

    return L
