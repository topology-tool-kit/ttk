def doIt(X, method, ncomponents, nneighbors, njobs):
    import importlib

    # check if numpy is installed
    loader = importlib.find_loader('numpy')
    found = loader is not None
    if found:
        print("[ManifoldLearning] Python: numpy module found.")
    else:
        print("[ManifoldLearning] Python error: numpy module not found.")
        return 0

    # check if scipy is installed
    loader = importlib.find_loader('scipy')
    found = loader is not None
    if found:
        print("[ManifoldLearning] Python: scipy module found.")
    else:
        print("[ManifoldLearning] Python error: scipy module not found.")
        return 0

    # check if scikit-learn is installed
    loader = importlib.find_loader('sklearn')
    found = loader is not None
    if found:
        print("[ManifoldLearning] Python: sklearn module found.")
    else:
        print("[ManifoldLearning] Python error: sklearn module not found.")

    from sklearn import manifold
    import numpy

    if method == 0:
        se = manifold.SpectralEmbedding(n_components=ncomponents, n_neighbors=nneighbors)
        Y = se.fit_transform(X)
    elif method == 1:
        lle = manifold.LocallyLinearEmbedding(nneighbors, ncomponents, eigen_solver='auto', method='standard')
        Y = lle.fit_transform(X)
    elif method == 2:
        mds = manifold.MDS(ncomponents, max_iter=100, n_init=1)
        Y = mds.fit_transform(X)
    elif method == 3:
        tsne = manifold.TSNE(n_components=ncomponents, init='pca', random_state=0)
        Y = tsne.fit_transform(X)
    elif method == 4:
        iso = manifold.Isomap(nneighbors, ncomponents)
        Y = iso.fit_transform(X)

    L = list()
    for i in range(ncomponents):
        L.append(numpy.copy(Y[:, i]))

    return L
