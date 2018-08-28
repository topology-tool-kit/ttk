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

    return 0
