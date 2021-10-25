import numpy as np
import scipy
import cvxopt
import itertools
from sklearn.cluster import KMeans


def topic_score(word_doc, K_topic, K_cand, n_cluster, word_freq_quantile = 1, kmeans_start = 200):
    n, p = word_doc.shape
    word_freq = np.mean(word_doc, axis = 1)
    word_freq_trunk = np.maximum(word_freq, np.quantile(word_freq, word_freq_quantile))

    word_doc = (word_doc.T / word_freq_trunk).T
    svd = scipy.sparse.linalg.svds(word_doc, k = K_topic)
    Xi = svd[0]

    # Step 1
    Xi[:,0] = abs(Xi[:,0])
    R = (Xi[:,1:].T / Xi[:,0]).T

    # Step 2
    vertices = vertices_estimate(R, K_cand, n_cluster, kmeans_start = kmeans_start)
    V = vertices[0]
    theta = vertices[1]

    # Step 3
    R_extend = np.concatenate((np.ones((n,1)), R), axis = 1)
    V_extend = np.concatenate((np.ones((K_topic,1)), V), axis = 1)
    Pi = R_extend @ np.linalg.inv(V_extend, )
    Pi = np.maximum(Pi, 0)
    rowsum = Pi.sum(axis = 1)
    Pi = Pi.T / rowsum

    # Step 4
    A_hat = np.sqrt(word_freq_trunk) * Xi[:,0] * Pi
    Pi = Pi.T
    A_hat = A_hat.T
    colsum = A_hat.sum(axis = 0)
    A_hat = A_hat / colsum

    return A_hat, R, V, Pi, theta


def vertices_estimate(R, K_cand, n_cluster, kmeans_start):
    n, K = R.shape
    K = K + 1 # K = K_topic

    # Step 2a
    obj = KMeans(n_clusters = n_cluster, max_iter = 300, n_init = kmeans_start).fit(R)
    theta = obj.cluster_centers_
    theta_original = theta.copy()

    # Step 2b'
    inner = theta @ theta.T
    d = np.diag(inner).reshape((n_cluster,1))
    one = np.ones((n_cluster,1))
    distance = d @ one.T + one @ d.T - 2 * inner
    top2 = np.where(distance == np.max(distance))
    top2 = np.array([top2[0][0], top2[1][0]])
    
    theta0 = theta[top2,:]
    theta = theta[np.setdiff1d(np.arange(n_cluster),top2),:]

    if K_cand > 2:
        for k0 in range(3, K_cand + 1):
            inner = theta @ theta.T
            n_k0 = theta.shape[0]
            d = np.diag(inner).reshape((n_k0,1))
            distance = np.ones((k0 - 1,1)) @ d.T - 2 * theta0 @ theta.T
            avg_dist = distance.mean(axis = 0)
            max_ix = np.argmax(avg_dist)
            
            theta0 = np.concatenate((theta0, theta[max_ix:max_ix + 1,:]), axis = 0)
            theta = theta[np.setdiff1d(np.arange(n_k0),max_ix),:]
        theta = theta0

    # Step 2b
    index_matrix = np.array(list(itertools.combinations(np.arange(K_cand),K))) # (K_cand choose K_topic) X K_topic array
    max_outlier = np.zeros(index_matrix.shape[0])
    for i in range(index_matrix.shape[0]):
        for j in range(K_cand):
            max_outlier[i] = max(max_outlier[i], simplex_distance(theta[j,:], theta[index_matrix[i],:]))
    min_ix = np.argmin(max_outlier)

    V = theta[index_matrix[min_ix],:]
    return V, theta_original

def simplex_distance(x, V):
    n_v, p = V.shape
    
    G = np.concatenate((np.eye(n_v - 1), -np.ones((n_v - 1, 1))), axis = 1)
    VV = G @ V
    P = VV @ VV.T
    q = VV @ (x - V[n_v - 1, :]).reshape((p,1))
    h = np.zeros(n_v)
    h[n_v - 1] = -1

    obj = qp_const(P = P, q = q, G = G.T, h = h)
    return np.sum((x - V[n_v - 1,:])**2) + 2 * obj

    
def qp_const(P,q,E = None,f = None,G = None,h = None,obj_only = True):
    '''
    solve a quadratic programming with constraints:

    min 1/2 * x^T P x - q^T x
    s.t. Ex = f
         Gx >= h
    -----------------
    Returns:
        x: solution
        obj: objective function at x
        only return obj by default
    '''
    args = [cvxopt.matrix(P), cvxopt.matrix(q)]
    if G is not None:
        args.extend([-cvxopt.matrix(G), -cvxopt.matrix(h)])
    if E is not None:
        f = np.array(f * 1.).reshape((-1,1))
        f_d = f.shape[0]
        E = np.array(E).reshape((f_d,-1))    
        args.extend([cvxopt.matrix(E), cvxopt.matrix(f)])
        
    cvxopt.solvers.options['show_progress'] = False
    try:
        sol = cvxopt.solvers.qp(*args)
        if 'optimal' not in sol['status'] or not sol['y']:
            # not sure if this is the best way to handle all non-optimal cases
            obj = np.nan
            x = None
        else:
            x = np.array(sol['x']).reshape((P.shape[1],))
        #     obj = np.array(sol['y'])[0,0]
            obj = sol['y'][0]
    except:
        # e.g., P is not positive-semidefinite
        obj = np.nan
        x = None
    if obj_only:
        result = obj
    else:
        result = [x,obj]
    return result