# from sklearn.base import TransformerMixin, BaseEstimator
from sklearn.cluster import MiniBatchKMeans
from collections import defaultdict
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import StandardScaler

# class KMeansFeatureSelection(TransformerMixin, BaseEstimator):
class KMeansFeatureSelection(object):
    def __init__(self, n_features=10, standardized=False):
        self.n_features = n_features
        self.standardized = standardized
        self.indices_ = []

    def fit(self, X, y=None):
        if not self.standardized:
            X = StandardScaler().fit_transform(X)
        
        X = X.T
        kmeans = MiniBatchKMeans(n_clusters=self.n_features).fit(X)
        clusters = kmeans.predict(X)
        cluster_centers = kmeans.cluster_centers_

        dists = defaultdict(list)
        for i, c in enumerate(clusters):
            dist = euclidean_distances([X[i, :]], [cluster_centers[c, :]])[0][0]
            dists[c].append((i, dist))

        self.indices_ = [sorted(f, key=lambda x: x[1])[0][0] for f in dists.values()]
        
        return self

    def transform(self, X, y=None):
        return X[:, self.indices_]
    
    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)
