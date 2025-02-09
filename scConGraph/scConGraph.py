import networkx as nx
import matplotlib.pyplot as plt
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import distance
from scipy import sparse
import random
# import community
import time
import scanpy as sc
from sklearn.preprocessing import normalize
from matplotlib.sankey import Sankey
from matplotlib.gridspec import GridSpec
from sklearn import metrics
from collections import Counter
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, rgb2hex
from matplotlib.patches import Rectangle
import subprocess
import os
import scanpy.external as sce

def createScConGraphObj(sc_obj_comb, ctrl_colors=None, treat_colors=None, key = 'conditions', ctrl_name = 'Ctrl', treat_name = 'Gemc', 
                        runLabel = "Sample", resultPath = ""):
                    
    '''
    Create an scConGraph object based on the scRNA-seq data under control and treated conditions. 
    This function also adds the meta-information into the created object.
    
    Parameters
    ----------
    sc_obj_comb : AnnData
        A Scanpy AnnData object that includes the scRNA-seq data for both control and treated conditions.
    key : str
        The name of the column in `sc_obj_comb.obs` that records the sample source, indicating whether the data is from control or treated conditions. Default is "conditions".
    ctrl_name : str
        The identifier for the control condition in the `key` column. Default is 'Ctrl'.
    treat_name : str
        The identifier for the treated condition in the `key` column. Default is 'Treat'.
    ctrl_colors : list, optional
        A list of colors to map the control clusters. If not provided, default colors will be used.
    treat_colors : list, optional
        A list of colors to map the treated clusters. If not provided, default colors will be used.
    runLabel : str, optional
        A label for the analysis run. This label will be used as part of the filenames for saved analysis results. Default is "Sample".
    resultPath : str, optional
        A path to save the analysis results. By default, results will be saved in the current directory.
    
    Returns
    -------
    scConGraph
        An initialized scConGraph object with the input scRNA-seq data and additional metadata.
        
    Notes
    -----
    - The scConGraph object facilitates comparative analysis of scRNA-seq data between control and treated conditions, including preprocessing, clustering, and visualization.
    - It's important to ensure that the `key`, `ctrl_name`, and `treat_name` parameters accurately reflect the structure of the input `sc_obj_comb` to correctly separate
    '''
    data_ctrl, data_treat = splitScObjects(sc_obj_comb, sam_key=key, sam_values=[ctrl_name, treat_name])
    scg = scConGraph(sc_obj_comb, data_ctrl, data_treat, key, ctrl_name, treat_name,  ctrl_colors, treat_colors,
                runLabel, resultPath)
    return (scg)



def splitScObjects(scobj, sam_key, sam_values):
    data_split = [scobj[scobj.obs[sam_key] == v,:] for v in sam_values]
    return(data_split)



def runScanpyProcess(scobj, n_top_genes = 2000, regress_out_variables = ['nCount_RNA', 'percent_mito', 'percent_ribo'], 
                     n_pcs = 40, harmony = False, condition = 'conditions' ):
    '''
    Preprocess and analyze single-cell RNA-seq data using Scanpy. This includes selecting highly variable genes,
    regressing out specific variables, scaling, running PCA, and optionally integrating data with Harmony for batch correction.
    
    Parameters
    ----------
    scobj : AnnData
        The single-cell RNA-seq dataset wrapped in an AnnData object.
    n_top_genes : int, optional
        The number of highly variable genes to retain for further analysis. Default is 2000.
    regress_out_variables : list of str, optional
        A list of variables to regress out during preprocessing. These might include technical factors like 
        total counts per cell (`nCount_RNA`), mitochondrial gene expression percent (`percent_mito`), 
        and ribosomal gene expression percent (`percent_ribo`). Default is ['nCount_RNA', 'percent_mito', 'percent_ribo'].
    n_pcs : int, optional
        The number of principal components to compute during PCA. Default is 40.
    harmony : bool, optional
        Whether to use Harmony for batch correction after PCA. Default is False.
    condition : str, optional
        The key in `scobj.obs` used for batch information if Harmony is applied. Only relevant if `harmony` is True. Default is 'conditions'.
    
    Returns
    -------
    AnnData
        The processed AnnData object.    
    Notes
    -----
    - This function is designed to work as part of a larger single-cell RNA-seq analysis pipeline, 
      where preprocessing steps like cell filtering and normalization are assumed to have been done beforehand.
    - Using Harmony for batch correction is optional and should be based on whether batch effects are a concern in the dataset.
    - The choice of variables to regress out can significantly impact downstream analysis and should be considered carefully.
    '''
    sc.pp.highly_variable_genes(scobj, n_top_genes = n_top_genes)
    scobj.raw = scobj
    scobj = scobj[:, scobj.var.highly_variable]

    if regress_out_variables is not None:
        sc.pp.regress_out(scobj, regress_out_variables)
    sc.pp.scale(scobj)

    sc.tl.pca(scobj, svd_solver='arpack')
    if harmony:
        sce.pp.harmony_integrate(scobj,  condition)
        scobj.obsm['X_pca'] = scobj.obsm['X_pca_harmony'].copy()
        
    sc.pp.neighbors(scobj, n_neighbors=10, n_pcs=n_pcs)
    sc.tl.umap(scobj)
    return (scobj)




def getSimilarity(PC_data1, PC_data2, k, alpha):
    '''
    Calculate the similarity between two sets of data points based on their Principal Component Analysis (PCA) results. 
    The similarity is computed using an alpha-decaying generalized Gaussian kernel based on the euclidean distance, a methods in PHATE (Moon, K.R. et al, Nature Biotechnology, 2019)

    Parameters
    ----------
    PC_data1 : ndarray
        An array of shape (n_samples_1, n_features) containing the PC representations of the first dataset.
    PC_data2 : ndarray
        An array of shape (n_samples_2, n_features) containing the PC representations of the second dataset.
    k : int
        The number of nearest neighbors to consider for each point when computing similarity.
    alpha : float
        The scaling factor for the Gaussian kernel, controlling the width of the kernel and hence the decay of similarity with distance.

    Returns
    -------
    ndarray
        A similarity matrix of shape (n_samples_1, n_samples_2). Each element in the matrix represents the similarity 
        between a pair of samples from the first and second datasets, respectively.

    Notes
    -----
    - The function first computes the square Euclidean distance between all pairs of points across the two datasets.
    - For each point, it identifies its k-nearest neighbors and calculates the Gaussian kernel-based similarity with a decay factor determined by alpha.
    - The final similarity score is normalized such that similarities of a point to its k-nearest neighbors sum up to 1.
    - This method provides a way to assess the pairwise similarity between cells from different conditions or time points based on their transcriptional profiles.
    '''

    dist = distance.cdist(PC_data1, PC_data2, 'sqeuclidean')
    dist_ord = dist.argsort()
    k_dist = [dist[i, j] for a, (i, j) in enumerate(zip(np.array(range(0, dist.shape[0])), dist_ord[:, k]))]
    k_dist = np.array(k_dist).reshape(-1, 1)
    affinity = np.exp(- (dist / k_dist) ** alpha)

    dist2 = distance.cdist(PC_data2, PC_data1, 'sqeuclidean')
    dist_ord2 = dist2.argsort()
    k_dist2 = [dist2[i, j] for a, (i, j) in enumerate(zip(np.array(range(0, dist2.shape[0])), dist_ord2[:, k]))]
    k_dist2 = np.array(k_dist2).reshape(-1, 1)
    affinity2 = np.exp(- (dist2 / k_dist2) ** alpha)

    affinity = affinity / 2 + affinity2.T / 2
    if (PC_data1.shape == PC_data2.shape):
        # print('TRUE')
        affinity[np.eye(affinity.shape[0], dtype=np.bool_)] = 0
        diff_prob = affinity
    else:
        diff_prob = affinity
    return (diff_prob)




def runLINE(path, edges, output, order=1, size=100, negative=5, samples=500, threads = 20):
    '''
    Runs the LINE (Large-scale Information Network Embedding) algorithm to generate embeddings.

    Parameters
    ----------
    - path (str): The file path to the LINE executable.
    - edges (str): The file path for the input graph in the form of edges (source target wight).
    - output (str): The file path where the embeddings will be saved.
    - order (int, optional): The order of the proximity used in LINE. Can be 1 or 2. Default is 1.
    - size (int, optional): The dimension of the embeddings. Default is 100.
    - negative (int, optional): The number of negative samples used in negative sampling. Default is 5.
    - samples (int, optional): The total number of samples used in stochastic gradient descent. Default is 500.
    - threads (int, optional): The number of parallel threads to use. Default is 20.

    This function calls the LINE command line tool with the provided parameters to generate
    embeddings for the input graph. It captures the output and errors of the LINE tool and
    prints the output to the console.
    '''
    
    LINE = path.strip()  
    cmd = [
        LINE,
        '-train', edges,
        '-output', output,
        '-binary', '0',
        '-size', str(size),
        '-order', str(order),
        '-negative', str(negative),
        '-samples', str(samples),
        '-threads', str(threads) # 20 in previous version
    ]
    print('Command line - Run LINE:', ' '.join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    last_print_time = time.time()
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break  
        
        current_time = time.time()
        if current_time - last_print_time >= 10:
            print(line.strip())
            last_print_time = current_time
    process.wait()  
    



def norm_embedding(embed_1st = None, embed_2nd_1 = None, embed_2nd_2 = None, mode = 0, weight = [0.5, 0.5]):
    '''
    Normalize embeddings derived from first-order and/or second-order proximity. This function supports
    multiple modes of normalization to accommodate different types of data embeddings (e.g., time-series
    or perturbation data).

    Parameters
    ----------
    embed_1st : DataFrame 
        DataFrame containing the first-order embeddings. 
    embed_2nd_1 : DataFrame 
        DataFrame containing the second-order embeddings.
    embed_2nd_2 : DataFrame
        DataFrame containing the second set of second-order embeddings, used only in time-series data (mode=1).
    mode : int, optional
        Determines the mode of normalization based on the type of data:
            - 0: Perturbation data, combining first and second-order embeddings.
            - 1: Time-series data, combining first and two sets of second-order embeddings.
            - 2: Uses only first-order embeddings.
            - 3: Uses only second-order embeddings (requires embed_2nd_1).
        Default is 0.
    weight : list of float, optional
        Weights to apply to first and second-order embeddings during combination. Only relevant for modes 0 and 1.
        Default is [0.5, 0.5], indicating equal weighting.

    Returns
    -------
    ndarray
        Normalized embeddings, where each embedding is L2-normalized. For modes 0 and 1, the embeddings are
        concatenated according to specified weights before normalization.

    Notes
    -----
    - This function facilitates the integration and comparison of embeddings from different sources or generated
      through different methods (e.g., LINE for first-order and second-order proximity).
    - The choice of mode and weights allows flexibility in handling a variety of single-cell RNA-seq data analysis
      scenarios, including differential expression analysis under perturbation or the study of dynamic processes.
    '''
    if mode == 0:
        print('Mode 0: perturbed data')
        # normalize for sample, like tsne
        embed_1st_norm = pd.DataFrame(normalize(embed_1st, norm='l2'), index = embed_1st.index)
        embed_2nd_norm = pd.DataFrame(normalize(embed_2nd_1, norm='l2'), index = embed_2nd_1.index)
        embed_cat = pd.concat([weight[0]*embed_1st_norm, weight[1]*embed_2nd_norm], axis=1, join='inner', ignore_index=True)
        embed_norm = normalize(embed_cat, norm='l2')
    elif mode == 1:
        print('Mode 1: time-series data')
        embed_1st_norm = pd.DataFrame(normalize(embed_1st, norm='l2'), index = embed_1st.index)
        embed_2nd_norm1 = pd.DataFrame(normalize(embed_2nd_1, norm='l2'), index = embed_2nd_1.index)
        embed_2nd_norm2 = pd.DataFrame(normalize(embed_2nd_2, norm='l2'), index = embed_2nd_1.index)
        embed_cat = pd.concat([weight[0]*embed_1st_norm, weight[1]*embed_2nd_norm1, weight[2]*embed_2nd_norm2],
                              axis=1, join='inner', ignore_index=True)
        embed_norm = normalize(embed_cat, norm='l2')
    elif mode == 2:
        print('Mode2: only first-order')
        embed_norm = normalize(embed_1st, norm='l2')
    elif mode == 3:
        print('Mode3: only second-order')
        embed_norm = normalize(embed_2nd_1, norm='l2')
    return embed_norm




def runLINEClustering(obj, embed, key = 'SCG', resolution = 0.5):
    '''
    Apply clustering on an AnnData object using embeddings generated by the LINE algorithm. 

    Parameters
    ----------
    obj : AnnData
        The AnnData object.
    embed : ndarray
        The embedding matrix to use for clustering, where rows correspond to cells and columns to the
        dimensions of the embedding space.
    key : str, optional
        The key under which to store the clustering results in the AnnData object's `.obs` attribute.
        Default is 'SCG', standing for scConGraph clustering.
    resolution : float, optional
        The resolution parameter for the Leiden clustering algorithm, controlling the granularity of
        clustering. Higher values lead to more clusters. Default is 0.5.

    Returns
    -------
    AnnData
        The input AnnData object, now augmented with the clustering results stored under `obj.obs[key+'_cluster']`.
        The function also updates the `.obsm` attribute with UMAP coordinates based on the provided embeddings.

    Notes
    -----
    - This function internally uses Scanpy's neighbor graph construction, UMAP dimensionality reduction,
      and Leiden clustering algorithms to process the scCobGraph clustering.
    '''
    obj.obsm[key] = embed
    sc.pp.neighbors(obj, use_rep = key)
    sc.tl.umap(obj)
    sc.tl.leiden(obj, key_added = key+'_cluster', resolution = resolution)
    return obj




# Cluster-Cluster transition probability matrix
def cls2cls_prob(data_ctrl, data_treat, weight_df, thre_cell=0.05, thre_cls=0.05):
    '''
    Calculate cluster-to-cluster transition probabilities between control and treated samples based on
    a provided similarity matrix. 

    Parameters
    ----------
    data_ctrl : AnnData
        The AnnData object containing the control sample data, with clusters already identified.
    data_treat : AnnData
        The AnnData object containing the treated sample data, with clusters already identified.
    weight_df : ndarray or DataFrame
        A matrix of similarities between individual cells in the control and treated datasets.
    thre_cell : float, optional
        The threshold to apply to the cell-to-cluster transition probabilities. Values below this threshold
        will be set to zero. Default is 0.05.
    thre_cls : float, optional
        The threshold to apply to the cluster-to-cluster transition probabilities. Values below this threshold
        will be set to zero. Default is 0.05.

    Returns
    -------
    ndarray
        A matrix of cluster-to-cluster transition probabilities, where the element at index (i, j) represents
        the probability of transition from cluster i in the control sample to cluster j in the treated sample.

    Notes
    -----
    - This function facilitates alignment of clusters across conditions.
    '''
    # probability of cell to cell
    n_clus = np.array([len(data_ctrl.obs["SCG_cluster"].unique()), len(data_treat.obs["SCG_cluster"].unique())])
    n_cells = np.array([data_ctrl.obs.shape[0], data_treat.obs.shape[0]])
    cell2cell = normalize(np.array(weight_df), norm='l1')

    # probability of cell to cluster
    cell_cls_ctrl = data_ctrl.obs['SCG_cluster'] 
    cell_cls_treat = data_treat.obs['SCG_cluster']
    cls_cellnum_treat = [np.sum(cell_cls_treat == str(j)) for j in range(n_clus[1])]

    cell2cls = [np.sum(cell2cell[:, np.where(cell_cls_treat == str(j))[0]], axis=1) / cls_cellnum_treat[j] for j in
                range(n_clus[1])]
    cell2cls = np.transpose(np.array(cell2cls))
    cell2cls = normalize(np.array(cell2cls), norm='l1')
    cell2cls[cell2cls < thre_cell] = 0
    cell2cls = normalize(np.array(cell2cls), norm='l1')

    # probability of cluster to cluster
    cls2cls = [np.sum(cell2cls[np.where(cell_cls_ctrl == str(i))[0], :], axis=0) for i in range(n_clus[0])]
    cls2cls = normalize(np.array(cls2cls), norm='l1')
    cls2cls[cls2cls < thre_cls] = 0
    cls2cls = normalize(cls2cls, norm='l1')
    return cls2cls


def rapid_QC(transition_TPM, origin_TPM, threshold = 0.1):
    '''
    Perform a rapid quality control (QC) check on transition and origin probability matrices to ensure
    that low-probability transitions are filtered out, improving the interpretability and reliability of
    cluster-to-cluster transition analyses.

    Parameters
    ----------
    transition_TPM : ndarray
        The transition probability matrix, where each element (i, j) represents the probability of a cell
        transitioning from cluster i in the control condition to cluster j in the treated condition.
    origin_TPM : ndarray
        The origin probability matrix, where each element (i, j) represents the probability of a cell
        in cluster j in the treated condition originating from cluster i in the control condition.

    Returns
    -------
    tuple of ndarray
        A tuple containing the updated transition and origin probability matrices after QC. Low-probability
        transitions are set to zero, and the matrices are renormalized to ensure that probabilities sum to 1.
    '''
    n_clus = transition_TPM.shape
    for i in range(n_clus[0]):
        for j in range(n_clus[1]):
            if transition_TPM[i, j] + origin_TPM[i, j] < threshold:
                transition_TPM[i, j] = 0
                origin_TPM[i, j] = 0
                transition_TPM = normalize(transition_TPM, norm='l1', axis=1)
                origin_TPM = normalize(origin_TPM, norm='l1', axis=0)
    return transition_TPM, origin_TPM





def plotUMAP(sc_obj, key = 'X_umap', group_by = 'leiden', title = '',  size = 8, colors = None):
    '''
    Visualize single-cell data using UMAP embeddings. This function plots the UMAP embeddings of the cells, 
    optionally coloring them by a specified cluster or group to highlight the cellular heterogeneity.

    Parameters
    ----------
    sc_obj : AnnData
        The AnnData object containing the single-cell data with UMAP coordinates stored in `.obsm[key]`.
    key : str, optional
        The key in the `.obsm` attribute of `sc_obj` that contains the UMAP coordinates. Default is 'X_umap'.
    group_by : str, optional
        The key in the `.obs` attribute of `sc_obj` used for coloring the cells in the UMAP plot. Typically,
        this will be a clustering label such as 'leiden' or any other categorical annotation. Default is 'leiden'.
    title : str, optional
        The title for the plot. Default is an empty string (no title).
    size : int, optional
        The size of the points in the plot. Default is 8.
    colors : list of str, optional
        A list of colors used for plotting the different groups identified by `group_by`. If None, a default
        color palette is used. Default is None.

    Returns
    -------
    The matplotlib Axes object with the UMAP plot.
    '''
    plot_data = pd.DataFrame({'UMAP1':sc_obj.obsm[key][:,0],
                              'UMAP2':sc_obj.obsm[key][:,1],
                              'Cluster':list(sc_obj.obs[ group_by ])})
    plot_data = plot_data.sort_values(by='Cluster', axis=0, ascending=True)
    
    palette = colors[:len(plot_data['Cluster'].unique())]
    ax = sns.scatterplot(data=plot_data, x="UMAP1", y="UMAP2", hue="Cluster", linewidth=0, s=size, palette=palette)
    ax.legend(loc = 2, bbox_to_anchor = (1.01,1), ncol=1)
    ax.set(xlabel = key+'1', ylabel = key+'2')
    plt.title(label = title, fontsize= 14)  
    return(ax)




def plotScore(data, title, ax, score_type =  'silhouette_score', show_legend=False ):
    '''
    Plot scoring metrics (e.g., silhouette score or Davies-Bouldin index) for different combinations of clustering and embeddings.

    Parameters
    ----------
    data : AnnData
        The AnnData object containing the single-cell data with clustering results stored in `.obs`.
    title : str
        The title of the plot, usually indicating the type of data or the specific comparison being made.
    ax : matplotlib.axes._subplots.AxesSubplot
        The matplotlib Axes object where the plot will be drawn. This allows for integration into subplots or
        more complex figures.
    score_type : str, optional
        The type of score to plot. Supported values include 'silhouette_score' and 'davies_bouldin_score'.
        Default is 'silhouette_score'.
    show_legend : bool, optional
        Whether to display a legend on the plot. Default is False.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the calculated scores for three clustering method (scConGraph, single, combined).
    '''

    size = 3
    cluster_comb = data.obs['leiden_comb']
    cluster_single = data.obs['leiden_single']
    cluster_SCG = data.obs['SCG_cluster']

    umap = ['comb_UMAP', 'single_UMAP', 'SCG_UMAP']
    
    a = b = c = np.zeros(len(umap))
    if score_type == 'silhouette_score':
        a = np.array([metrics.silhouette_score(data.obsm[i], cluster_comb) for i in umap])
        b = np.array([metrics.silhouette_score(data.obsm[i], cluster_single) for i in umap])
        c = np.array([metrics.silhouette_score(data.obsm[i], cluster_SCG) for i in umap])
        #print(a, b, c)
    elif score_type == 'davies_bouldin_score':
        a = np.array([metrics.davies_bouldin_score(data.obsm[i], cluster_comb) for i in umap])
        b = np.array([metrics.davies_bouldin_score(data.obsm[i], cluster_single) for i in umap])
        c = np.array([metrics.davies_bouldin_score(data.obsm[i], cluster_SCG) for i in umap])
        #print(a, b, c)
    elif score_type == 'calinski_harabasz_score':
        a = np.array([metrics.calinski_harabasz_score(data.obsm[i], cluster_comb) for i in umap])
        b = np.array([metrics.calinski_harabasz_score(data.obsm[i], cluster_single) for i in umap])
        c = np.array([metrics.calinski_harabasz_score(data.obsm[i], cluster_SCG) for i in umap])
        #print(a, b, c)
    else:
        print("Invalid score_type. Please choose 'silhouette_score', 'calinski_harabasz_score' or 'davies_bouldin_score'.")
        return None  # 返回None或适当的默认值
    score_df = pd.DataFrame([a, b, c], columns =umap, index = ['leiden_comb', 'leiden_single', 'SCG_cluster'])
    
    score_df.T.plot(kind='bar', ax=ax, fontsize=12)
    ax.set_title(title+': '+score_type, fontsize=14) 
    ax.set_ylabel(score_type, fontsize=14) 
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=12)     
    if show_legend:
        ax.legend(title="Clusters", title_fontsize='13', fontsize='12', loc='center left', framealpha=0.4, bbox_to_anchor=(1, 0.5))
    else:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
    
    for p in ax.patches:
        ax.annotate('{:.4f}'.format(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', fontsize=10, color='black', rotation=0, xytext=(0, 5),
                    textcoords='offset points')
    return score_df
    
    
    

def getClusterExpr(scobj, adata, key):
    expr = adata[scobj.obs.index,:].copy()
    expr_data = np.array(expr.X)
    clus = np.unique(list(scobj.obs[key]))
    clus_pc_data = [np.average(expr_data[np.where(scobj.obs[key] == x)[0], :], axis=0) for x in clus]
    clus_pc_data = np.array(clus_pc_data)
    return (clus_pc_data, clus)





def getGraph(cls2cls_sim, cls2cls_sim_Ctrl, cls2cls_sim_Treat, transition_TPM, origin_TPM, graph_comb_QC = True):
    G = nx.Graph()
    weight_list=[]
    set1 = []
    set2 = []
    set_pair = []

    weight_list1 = []
    weight_list2= []
    weight_list_pair = []

    for i in range(origin_TPM.shape[0]):
        for j in range(origin_TPM.shape[1]):
            if graph_comb_QC:
                if (origin_TPM[i, j] > 0  or transition_TPM[i, j] > 0):
                    nodes1 = 'C'+str(i)
                    nodes2 = 'T'+str(j)
                    weight = cls2cls_sim[i, j]
                    G.add_edge(nodes1, nodes2, weight = weight)
                    set_pair.append((nodes1, nodes2))
                    weight_list_pair.append(weight)
            else:
                nodes1 = 'C'+str(i)
                nodes2 = 'T'+str(j)
                weight = cls2cls_sim[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                set_pair.append((nodes1, nodes2))
                weight_list_pair.append(weight)

    for i in range(origin_TPM.shape[0]):
        for j in range(origin_TPM.shape[0]):
            if (i != j and cls2cls_sim_Ctrl[i, j] > 0):
                nodes1 = 'C'+str(i)
                nodes2 = 'C'+str(j)
                weight = cls2cls_sim_Ctrl[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                set1.append((nodes1, nodes2))
                weight_list1.append(weight)

    for i in range(origin_TPM.shape[1]):
        for j in range(origin_TPM.shape[1]):
            if (i != j and cls2cls_sim_Treat[i, j] > 0):
                nodes1 = 'T'+str(i)
                nodes2 = 'T'+str(j)
                weight = cls2cls_sim_Treat[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                set2.append((nodes1, nodes2))
                weight_list2.append(weight)
    return G, set1, set2, set_pair, weight_list1, weight_list2, weight_list_pair





def plotGraph(G, origin_TPM, set1, set2, set_pair,  weight_list1, weight_list2, weight_list_pair):
    fig = plt.figure(figsize=(9, 4), dpi=300)
    pos2 = nx.spring_layout(G, iterations=200, seed=2023)

    options = { "node_size": 800, "alpha": 1}
    nx.draw_networkx_nodes(G, pos2, nodelist=['C'+str(i) for i in range(origin_TPM.shape[0])], node_color="tab:blue", **options)
    nx.draw_networkx_nodes(G, pos2, nodelist=['T'+str(i) for i in range(origin_TPM.shape[1])], node_color="tab:red", **options)

    a = nx.draw_networkx_edges(
        G,
        pos2,
        edgelist=set2,
        width= 2,
        alpha=0.6,
        edge_color = weight_list2, 
        edge_cmap = plt.cm.YlGn, 
        edge_vmin = 0.05, edge_vmax = 1
    )


    b = nx.draw_networkx_edges(
        G,
        pos2,
        edgelist= set1,
        width=2,
        alpha=0.6,
        edge_color = weight_list1, 
        edge_cmap = plt.cm.YlGn, 
        edge_vmin = 0.05, edge_vmax = 1
    )
    
    c= nx.draw_networkx_edges(G,
        pos2,
        edgelist= set_pair,
        width= 2.5, 
        alpha=0.5,
        edge_color = weight_list_pair, 
        edge_cmap = plt.cm.YlGn, 
        edge_vmin = 0.05, edge_vmax = 1
    ) 
    plt.colorbar(a)

    #nx.draw_networkx_labels(G, pos2, font_size=4, font_color="whitesmoke")
    nx.draw_networkx_labels(G, pos2, font_size=16, font_color="white")

    plt.tight_layout()
    plt.axis("off")





def getDiGraph(cls2cls_sim, cls2cls_sim_Ctrl, cls2cls_sim_Treat, transition_TPM, origin_TPM):

    G = nx.DiGraph()

    for i in range(origin_TPM.shape[0]):
        for j in range(origin_TPM.shape[1]):
            if True:
                nodes1 = 'C'+str(i)
                nodes2 = 'T'+str(j)
                weight = cls2cls_sim[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)

    for i in range(origin_TPM.shape[0]):
        for j in range(origin_TPM.shape[0]):
            if (i != j ):
                nodes1 = 'C'+str(i)
                nodes2 = 'C'+str(j)
                weight = cls2cls_sim_Ctrl[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                G.add_edge(nodes2, nodes1, weight = weight)

    for i in range(origin_TPM.shape[1]):
        for j in range(origin_TPM.shape[1]):
            if (i != j ):
                nodes1 = 'T'+str(i)
                nodes2 = 'T'+str(j)
                weight = cls2cls_sim_Treat[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                G.add_edge(nodes2, nodes1, weight = weight)
    return G


def sigmoid(x):
    return 1 / (1 + np.exp(-x))
    
    
def evaluateCertainty(P):
    P = np.where(P == 0, 1, P)
    entropy = -np.sum(P * np.log2(P), axis=1)
    average_entropy = np.mean(entropy)
    
    P_ref = np.full((P.shape[0], P.shape[1]), fill_value=1/P.shape[1])
    entropy_ref = -np.sum(P_ref * np.log2(P_ref), axis=1)
    average_entropy_ref = np.mean(entropy_ref)
    
    Certainty = 1 - average_entropy/average_entropy_ref
    return Certainty




def assignFlowType(DrugResponse,  thre_abundance = 0.9, thre_STC = 0.2544, Vol_Ratio = True):
    if Vol_Ratio:
        if  DrugResponse['STC'] >= 0.2544:
            return 'Acquired resistance'
        elif DrugResponse['AAC'] >= 0.9:
            return 'Absolutely intrinsic resistance'
        elif DrugResponse['RAC'] >= 0.9:
            return 'Relatively intrinsic resistance'
        else:
            return 'Sensitivity'
    else:
        if  DrugResponse['STC'] >= 0.2544:
            return 'Acquired resistance'
        elif DrugResponse['RAC'] >= 0.9:
            return 'Intrinsic resistance'
        else:
            return 'Sensitivity'
    
    
    
    
    
class scConGraph:

    def __init__(self, data_comb, data_ctrl, data_treat, key, ctrl_name, treat_name, ctrl_colors=None,
                 treat_colors=None,  runLabel = 'Sample', resultPath = ''):
        self.data_comb = data_comb
        self.data_ctrl = data_ctrl
        self.data_treat = data_treat
        self.key = key
        self.ctrl_name = ctrl_name
        self.treat_name = treat_name
        self.runLabel = runLabel
        self.resultPath = resultPath
        self.ctrl_colors = ctrl_colors
        self.treat_colors = treat_colors
        self.transition_TPM = None  
        self.origin_TPM = None 
        
            
        if self.ctrl_colors is None:
            self.ctrl_colors = ["#FEC643", "#437BFE", "#43FE69", "#FE6943", "#E78AC3",
                               "#43D9FE", "#FFEC1A", "#E5C494", "#A6D854", "#33AEB1",
                               "#EA6F5D", "#FEE8C3", "#3DBA79", "#D6EBF9", "#7F6699",
                               "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb"]
        if self.treat_colors is None:
            self.treat_colors = ['#EE854A', '#4878D0', '#6ACC64', '#D65F5F', '#956CB4',
                               '#82C6E2', '#D5BB67', '#8C613C', '#DC7EC0', '#797979',
                               "#333399", "#679966", "#c12e34", "#66d5a5", "#5599ff"]

        

        
        
    def preProcess(self, regress_out = None, npcs=40, n_top_genes_single = 2000, n_top_genes_comb=3000, process_comb = True, harmony = False):
        '''
        Preprocess the single-cell RNA-seq data contained within the scConGraph object. This method applies a series
        of preprocessing steps including the selection of highly variable genes, regression out of unwanted sources
        of variation, and data scaling. It also performs PCA to reduce dimensionality.

        Parameters
        ----------
        regress_out : list of str or None, optional
        A list of variables (in adata.obs) to regress out of the analysis, typically technical factors that are not of primary
        biological interest. Common examples include the total counts per cell (nCount_RNA), the percentage of
        mitochondrial gene expression (percent_mito), and the percentage of ribosomal gene expression (percent_ribo).
        If None, it defaults to ['nCount_RNA', 'percent_mito', 'percent_ribo']. Default is None.
        npcs : int, optional
            The number of principal components to retain in the PCA step. This parameter controls the dimensionality
            reduction aspect of the preprocessing. Default is 40.

        '''
        print('Preprocess control sample')
        self.data_ctrl = runScanpyProcess(self.data_ctrl, n_top_genes = n_top_genes_single, regress_out_variables = regress_out )
        print('Preprocess treated sample')
        self.data_treat = runScanpyProcess(self.data_treat, n_top_genes = n_top_genes_single, regress_out_variables = regress_out)
        if process_comb:
            print('Preprocess combined sample')
            self.data_comb= runScanpyProcess(self.data_comb, n_top_genes = n_top_genes_comb, regress_out_variables = regress_out, harmony = harmony, condition = self.key)




    def runClustering(self, resolution_comb = 0.5, resolution_ctrl = 0.5, resolution_treat = 0.5):
        '''
        Perform clustering on single-cell RNA-seq data within the scConGraph object. This method applies clustering
        separately to control and treated datasets, as well as to the combined dataset, using the Leiden algorithm
        in Scanpy. The goal is to identify distinct cell populations within and across different conditions.

        Parameters
        ----------
        resolution_comb : float, optional
            The resolution parameter for the Leiden algorithm applied to the combined dataset. Higher values lead
            to a finer clustering with more clusters. Default is 0.5.
        resolution_ctrl : float, optional
            The resolution parameter for the Leiden algorithm applied to the control dataset. Default is 0.5.
        resolution_treat : float, optional
            The resolution parameter for the Leiden algorithm applied to the treated dataset. Default is 0.5.

        Returns
        -------
        None
        The method updates the scConGraph object in-place, adding clustering labels to the `.obs` attribute
        of the control, treated, and combined AnnData objects contained within it.

        '''
        sc.tl.leiden(self.data_ctrl, key_added='leiden_single', resolution=resolution_ctrl)
        sc.tl.leiden(self.data_treat, key_added='leiden_single', resolution=resolution_treat)
        sc.tl.leiden(self.data_comb, key_added='leiden_comb', resolution=resolution_comb)

        self.data_ctrl.obs['leiden_comb'] = self.data_comb.obs['leiden_comb'][self.data_comb.obs[self.key] == self.ctrl_name].copy()
        self.data_treat.obs['leiden_comb'] = self.data_comb.obs['leiden_comb'][self.data_comb.obs[self.key] == self.treat_name].copy()
        self.data_comb.obs['leiden_single'] = pd.concat([self.data_ctrl.obs['leiden_single'],
                                                          self.data_treat.obs['leiden_single']], axis=0).copy()
        
        
        self.data_comb.obsm['comb_UMAP'] = self.data_comb.obsm['X_umap'].copy()
        self.data_ctrl.obsm['comb_UMAP'] = self.data_comb.obsm['X_umap'][self.data_comb.obs[self.key] == self.ctrl_name, :].copy()
        self.data_treat.obsm['comb_UMAP'] = self.data_comb.obsm['X_umap'][self.data_comb.obs[self.key] == self.treat_name, :].copy()
        
        self.data_ctrl.obsm['single_UMAP'] = self.data_ctrl.obsm['X_umap'].copy()
        self.data_treat.obsm['single_UMAP'] = self.data_treat.obsm['X_umap'].copy()



    def calculateAffinity(self, k=3, alpha= [20, 50], thre_edge_1st = None, thre_edge_2nd = None, resultPath = None):
        '''
        Calculate cell-to-cell affinity based on Principal Component (PC) analysis results, using a k-nearest
        neighbors approach combined with an alpha-decaying function. This method generates affinity matrices
        for control, treated, and combined datasets within the scConGraph object, facilitating the understanding
        of cellular relationships and transitions.

        Parameters
        ----------
        k : int, optional
            The number of nearest neighbors to consider in the affinity calculation. Default is 3.
        alpha : list of int, optional
            A list containing two alpha values for controlling the decay rate of the affinity function, where the
            first value is used for within-condition comparisons and the second for cross-condition comparisons.
            Higher alpha values result in a sharper decay of affinity with distance. Default is [20, 50].
        resultPath : str or None, optional
            The file path where the affinity matrices should be saved. If None, the matrices are not saved to disk.
            Default is None.

        Returns
        -------
        tuple of ndarray
            A tuple containing three affinity matrices: one for the control dataset, one for the treated dataset,
            and one for the combined dataset, corresponding to the within-condition and cross-condition cellular
            affinities calculated based on the method's parameters.
        '''
        runLabel = self.runLabel
        if resultPath is None:
            resultPath = self.resultPath
        
        PC_data_ctrl = self.data_ctrl.obsm['X_pca']
        PC_data_treat = self.data_treat.obsm['X_pca']
        PC_data_comb = pd.DataFrame(self.data_comb.obsm['X_pca'], index=self.data_comb.obs_names)
        n_cells = np.array([self.data_ctrl.shape[0], self.data_treat.obs.shape[0]])
     
        #PC_data_comb_ctrl = PC_data_comb[0: n_cells[0], :]
        #PC_data_comb_treat = PC_data_comb[n_cells[0]: (n_cells[0] + n_cells[1]), :]

        indices_ctrl = self.data_ctrl.obs_names
        indices_treat = self.data_treat.obs_names

        PC_data_comb_ctrl = PC_data_comb.loc[indices_ctrl, :]
        PC_data_comb_ctrl = PC_data_comb_ctrl.to_numpy()
        PC_data_comb_treat = PC_data_comb.loc[indices_treat, :]
        PC_data_comb_treat = PC_data_comb_treat.to_numpy()

        
        sim_Ctrl = getSimilarity(PC_data_ctrl, PC_data_ctrl, k=k, alpha=alpha[0])
        sim_Treat = getSimilarity(PC_data_treat, PC_data_treat, k=k, alpha=alpha[0])
        sim_Ctrl_Treat = getSimilarity(PC_data_comb_ctrl, PC_data_comb_treat, k=k, alpha=alpha[0])

        np.savetxt(resultPath + runLabel + '_Comb_similarity_alpha' + str(alpha[0]) + '_raw.txt', sim_Ctrl_Treat)
        # np.savetxt(resultPath + runLabel + '_Ctrl_similarity_alpha' + str(alpha[0]) + '_raw.txt', sim_Ctrl)
        # np.savetxt(resultPath + runLabel + '_Treat_similarity_alpha' + str(alpha[0]) + '_raw.txt', sim_Treat)

        rows, cols = np.where(sim_Ctrl_Treat != 0)
        # edge_weight_comb = [(f'C{i}', f'T{j}', sim_Ctrl_Treat[i, j]) for i, j in zip(rows, cols)]
        edge_weight_comb = [(f'C{i}', f'T{j}', sim_Ctrl_Treat[i, j]) for i, j in zip(rows, cols)] + [(f'T{j}', f'C{i}', sim_Ctrl_Treat[i, j]) for i, j in zip(rows, cols)] 

        rows, cols = np.where(np.triu(sim_Ctrl, k=1) != 0)
        edge_weight_ctrl = [(f'C{i}', f'C{j}', sim_Ctrl[i, j]) for i, j in zip(rows, cols)] + [(f'C{j}', f'C{i}', sim_Ctrl[i, j]) for i, j in zip(rows, cols)]

        rows, cols = np.where(np.triu(sim_Treat, k=1) != 0)
        edge_weight_treat = [(f'T{i}', f'T{j}', sim_Treat[i, j]) for i, j in zip(rows, cols)] +  [(f'T{j}', f'T{i}', sim_Treat[i, j]) for i, j in zip(rows, cols)] 

        if thre_edge_1st is not None:
            edge_weight_ctrl = [(u, v, w) for u, v, w in edge_weight_ctrl if w >= thre_edge_1st]
            edge_weight_treat = [(u, v, w) for u, v, w in edge_weight_treat if w >= thre_edge_1st]

        if thre_edge_2nd is not None:
            edge_weight_comb = [(u, v, w) for u, v, w in edge_weight_comb if w >= thre_edge_2nd]


        np.savetxt(resultPath + runLabel + '_Comb_edge_weight_alpha' + str(alpha[0]) + '_raw.txt', edge_weight_comb,
                   fmt='%s', delimiter='\t')
        np.savetxt(resultPath + runLabel + '_Ctrl_edge_weight_alpha' + str(alpha[0]) + '_raw.txt', edge_weight_ctrl, fmt='%s',
                   delimiter='\t')
        np.savetxt(resultPath + runLabel + '_Treat_edge_weight_alpha' + str(alpha[0]) + '_raw.txt', edge_weight_treat,
                   fmt='%s', delimiter='\t')
                   
        sparse_sim_Ctrl_Treat = getSimilarity(PC_data_comb_ctrl, PC_data_comb_treat, k=k, alpha=alpha[1])
        np.savetxt(resultPath + runLabel + '_Comb_similarity_alpha' + str(alpha[1]) + '_raw.txt', sparse_sim_Ctrl_Treat)

        return sparse_sim_Ctrl_Treat



    

    def runEmbedding(self, path, size_1ord=100, size_2ord=100, negative=5, samples=500, alpha=20, threads = 20):
        '''
        Run the LINE (Large-scale Information Network Embedding) algorithm to generate embeddings for both the 
        first-order and second-order proximities for control and treated samples.
        
        Parameters
        ----------
        path : str
            The file path to the executable LINE tool.
        size_1ord : int, optional
            The dimension of the first-order embeddings. Default is 100.
        size_2ord : int, optional
            The dimension of the second-order embeddings. Default is 100.
        negative : int, optional
            The number of negative samples used in training the LINE model, which helps in approximating the negative
            sampling portion of the loss function. Default is 5.
        samples : int, optional
            The total number of samples to draw during the stochastic gradient descent optimization. Default is 500.
        alpha : int, optional
            The alpha parameter for the LINE algorithm, controlling the importance of the proximity information.
            Default is 20.

        Returns
        -------
        None
            Save embedding as txt file in resultPath.  
        '''

        resultPath = self.resultPath
        runLabel = self.runLabel

        time_start = time.time()
        edges = resultPath + runLabel + '_Ctrl_edge_weight_alpha' + str(alpha) + '_raw.txt'
        output = resultPath + runLabel + '_Ctrl_1st_embeddings.txt'
        runLINE(path, edges=edges, output=output, order=1, size=size_1ord, negative=negative, samples=samples, threads = threads)

        edges = resultPath + runLabel + '_Treat_edge_weight_alpha' + str(alpha) + '_raw.txt'
        output = resultPath + runLabel + '_Treat_1st_embeddings.txt'
        runLINE(path, edges=edges, output=output, order=1, size=size_1ord, negative=negative, samples=samples, threads = threads)

        edges = resultPath + runLabel + '_Comb_edge_weight_alpha' + str(alpha) + '_raw.txt'
        output = resultPath + runLabel + '_2nd_embeddings.txt'
        runLINE(path, edges=edges, output=output, order=2, size=size_2ord, negative=negative, samples=samples, threads = threads)

        time_end = time.time()
        print('time cost', time_end - time_start, 's')





    def scConGraphCluster(self,  resolution_ctrl = 0.5, resolution_treat = 0.5):
        '''
        Perform clustering on the LINE embeddings within the scConGraph object. 

        Parameters
        ----------
        resolution_ctrl : float, optional
            The resolution parameter for the clustering algorithm applied to the control dataset. Higher values
            lead to more refined clustering, potentially increasing the number of clusters identified. Default is 0.5.
        resolution_treat : float, optional
            The resolution parameter for the clustering algorithm applied to the treated dataset. Like the control,
            higher values refine the granularity of clustering. Default is 0.5.

        Returns
        -------
        None
            The method updates the scConGraph object in-place, adding new clustering labels  `SCG` to the `.obs` 
            attribute and new UMAP labels `SCG_UMAP` to the `.obsm` attribute.
        '''
        resultPath = self.resultPath
        runLabel = self.runLabel
        
        n_cells = np.array([self.data_ctrl.obs.shape[0], self.data_treat.obs.shape[0]])
        cell_name_ctrl_id = ['C' + str(i) for i in range(n_cells[0])]
        cell_name_treat_id = ['T' + str(i) for i in range(n_cells[1])]

        cell_name_ctrl_treat_id = cell_name_ctrl_id.copy()
        cell_name_ctrl_treat_id.extend(cell_name_treat_id)

        self.data_ctrl.obs['simple_id'] = cell_name_ctrl_id
        self.data_treat.obs['simple_id'] = cell_name_treat_id

        cell_name_ctrl = self.data_ctrl.obs_names.tolist()
        cell_name_treat = self.data_treat.obs_names.tolist()
        cell_name_ctrl_treat = self.data_ctrl.obs_names.tolist()
        cell_name_ctrl_treat.extend(self.data_treat.obs_names.tolist())

        
        embed_ctrl_1st_raw = pd.read_csv(resultPath + runLabel + '_Ctrl_1st_embeddings.txt',
                                       skiprows=1, header=None, index_col=0,
                                       delim_whitespace=True)
        embed_ctrl_1st = embed_ctrl_1st_raw.loc[cell_name_ctrl_id, :]
        embed_ctrl_1st = embed_ctrl_1st.set_axis(cell_name_ctrl, axis=0)

        embed_treat_1st_raw = pd.read_csv(resultPath + runLabel + '_Treat_1st_embeddings.txt',
                                       skiprows=1, header=None, index_col=0,
                                       delim_whitespace=True)
        embed_treat_1st = embed_treat_1st_raw.loc[cell_name_treat_id, :]
        embed_treat_1st = embed_treat_1st.set_axis(cell_name_treat, axis=0)

        embed_ctrl_1st = embed_ctrl_1st.loc[self.data_ctrl.obs_names]
        embed_treat_1st = embed_treat_1st.loc[self.data_treat.obs_names]

        # cate and normalize embedding
        embed_ctrl_treat_2nd_raw = pd.read_csv(resultPath + runLabel + '_2nd_embeddings.txt',
                                         skiprows=1, header=None, index_col=0,
                                         delim_whitespace=True)
        embed_ctrl_treat_2nd = embed_ctrl_treat_2nd_raw.loc[cell_name_ctrl_treat_id,]
        embed_ctrl_treat_2nd = embed_ctrl_treat_2nd.set_axis(cell_name_ctrl_treat, axis=0)

        embed_ctrl_2nd = embed_ctrl_treat_2nd.loc[self.data_ctrl.obs_names]
        embed_treat_2nd = embed_ctrl_treat_2nd.loc[self.data_treat.obs_names]

        embed_ctrl = norm_embedding(embed_1st=embed_ctrl_1st, embed_2nd_1=embed_ctrl_2nd,
                                  mode=0, weight=[1, 1])
        embed_treat = norm_embedding(embed_1st=embed_treat_1st, embed_2nd_1=embed_treat_2nd,
                                  mode=0, weight=[1, 1])

        # Run clustering
        self.data_ctrl = runLINEClustering(self.data_ctrl, embed_ctrl, key='SCG', resolution=resolution_ctrl)
        self.data_treat = runLINEClustering(self.data_treat, embed_treat, key='SCG', resolution=resolution_treat)
        self.data_ctrl.obsm['SCG_UMAP'] = self.data_ctrl.obsm['X_umap']
        self.data_treat.obsm['SCG_UMAP'] = self.data_treat.obsm['X_umap']





    def Alignment(self, sparse_weight_mat, threshold_cell = 0.05, threshold_cls = 0.05):
        '''
        Align the control and treated clusters by calculating transition and origin probability 
        matrices. This method refines the transition probabilities by applying thresholds to filter out 
        low-probability transitions, which are often regarded as noise. 

        Parameters
        ----------
        sparse_weight_mat : ndarray
            A sparse cell-cell affinity matrix that derived from the `calculateAffinity` step. 

        Returns
        -------
        tuple of ndarray
            A tuple containing the transition and origin probability matrices. 
        '''

        # obtain the cls2cls matrix
        resultPath = self.resultPath
        runLabel = self.runLabel
        
        transition_TPM = cls2cls_prob(self.data_ctrl, self.data_treat, sparse_weight_mat, thre_cell = threshold_cell, thre_cls = threshold_cls)
        transition_TPM_raw = cls2cls_prob(self.data_ctrl, self.data_treat, sparse_weight_mat, thre_cell = threshold_cell, thre_cls=0)
        
        sparse_weight_mat_transpose = np.transpose(sparse_weight_mat)
        origin_TPM = cls2cls_prob(self.data_treat, self.data_ctrl, sparse_weight_mat_transpose, thre_cell = threshold_cell, thre_cls = threshold_cls)
        origin_TPM = np.transpose(origin_TPM)
        origin_TPM_raw = cls2cls_prob(self.data_treat, self.data_ctrl, sparse_weight_mat_transpose, thre_cell = threshold_cell, thre_cls=0)
        origin_TPM_raw = np.transpose(origin_TPM_raw)
        
        tmp = rapid_QC(transition_TPM, origin_TPM, threshold = threshold_cls*2)
        transition_TPM = tmp[0]
        origin_TPM = tmp[1]
#         print(origin_TPM)
#         print(transition_TPM)
        
        for i in range(transition_TPM.shape[0]):
            for j in range(transition_TPM.shape[1]):
                if ((transition_TPM[i, j] ==0) and (origin_TPM[i, j] == 0)):
                    origin_TPM_raw[i, j] = 0
                    transition_TPM_raw[i, j] = 0
                    
        transition_TPM_final = normalize(origin_TPM_raw,  norm='l1')
        origin_TPM_final = normalize(origin_TPM_raw,  axis = 0,  norm='l1')
        
        self.transition_TPM = transition_TPM_final
        self.origin_TPM = origin_TPM_final

        np.savetxt(resultPath + runLabel + '_transition_TPM_probability.txt', transition_TPM_final)
        np.savetxt(resultPath + runLabel + '_origin_TPM_probability.txt', origin_TPM_final)


        return transition_TPM_final, origin_TPM_final





    def showUMAP(self, UMAP = 'SCG_UMAP', Cluster = 'SCG_cluster'):
        '''
        Display UMAP plots for the control and treated samples processed within the scConGraph object. 
        Parameters
        ----------
        UMAP : str, optional
            The key in the `.obsm` attribute of the AnnData objects containing the UMAP coordinates. Default is 'SCG_UMAP'.
        Cluster : str, optional
            The clustering label in the `.obs` attribute of the AnnData objects used to color the cells in the UMAP plot.
            Default is 'SCG_cluster'.

        Returns
        -------
        None
            This function does not return a value but displays the UMAP plots directly.
        '''
        
        titles = [Cluster + '_Ctrl', Cluster + '_Treat']
        fig = plt.figure(figsize=(10, 4.5), dpi=300)
        ax1 = fig.add_subplot(1, 2, 1)
        plotUMAP(self.data_ctrl, key=UMAP, group_by= Cluster, title = titles[0], colors = self.ctrl_colors)
        ax2 = fig.add_subplot(1, 2, 2)
        plotUMAP(self.data_treat, key=UMAP, group_by=Cluster, title =  titles[1], colors = self.treat_colors)
        fig.tight_layout()


        

    def showHeatmap(self, transition_TPM = None, origin_TPM = None):
        '''
        Display heatmaps for the transition and origin probability matrices within the scConGraph object.

        Parameters
        ----------
        transition_TPM : ndarray, optional
            The transition probability matrix, where each element (i, j) represents the probability of a
            cell transitioning from cluster i under the control condition to cluster j under the treated 
            condition. If None, the internal stored matrix in the scConGraph object will be used. Default is None.
        origin_TPM : ndarray, optional
            The origin probability matrix, where each element (i, j) represents the probability of a
            cell in cluster j under the treated condition originating from cluster i under the control condition.
            If None, the internal stored matrix in the scConGraph object will be used. Default is None.

        Returns
        -------
        None
            This function does not return a value but directly renders the heatmaps for immediate visualization.
        '''
        if transition_TPM is None:
            transition_TPM = self.transition_TPM
        if origin_TPM is None:
            origin_TPM = self.origin_TPM

        columns = ['T{}'.format(i) for i in range(transition_TPM.shape[1])]
        index = ['C{}'.format(i) for i in range(transition_TPM.shape[0])]
        
        fig = plt.figure(figsize=(10, 4.5), dpi=300)
        ax1 = fig.add_subplot(1, 2, 1)
        sns.heatmap(pd.DataFrame(transition_TPM, index=index, columns=columns), annot=True, linewidths=.5, cmap="YlGnBu")
        plt.xlabel('Treated cluster')
        plt.ylabel('Control cluster')
        plt.title(label = 'Transition TPM', fontsize= 14)  

        ax1 = fig.add_subplot(1, 2, 2)
        sns.heatmap(pd.DataFrame(origin_TPM, index=index, columns=columns), annot=True, linewidths=.5, cmap="YlGnBu")
        plt.xlabel('Treated cluster')
        plt.ylabel('Control cluster')
        plt.title(label = 'Origin TPM', fontsize= 14)  
        fig.tight_layout()





    def rank_genes(self):
        self.data_ctrl.uns['log1p']["base"] = None
        sc.tl.rank_genes_groups(self.data_ctrl, 'SCG_cluster', method='wilcoxon')
        self.data_treat.uns['log1p']["base"] = None
        sc.tl.rank_genes_groups(self.data_treat, 'SCG_cluster', method='wilcoxon')

        sc.pl.rank_genes_groups_heatmap(self.data_ctrl, n_genes=20, use_raw=True, swap_axes=True,
                                        show_gene_labels=True)
        sc.pl.rank_genes_groups_heatmap(self.data_treat, n_genes=20, use_raw=True, swap_axes=True,
                                        show_gene_labels=True)





    def evaluateClustering(self, score_type = 'silhouette_score'):
        '''
        Evaluate the effientness of scConGraph clustering to extract the cells' given-time transcriptomic similarity
        and cross-time cell context information using silhouette score or Davies-Bouldin index.

        Parameters
        ----------
        score_type : str, optional
            The type of clustering evaluation metric to use. Supported options include 'silhouette_score', 
            'calinski_harabasz_score' and 'davies_bouldin_score' for evaluating the compactness and separation of clusters. 
            The default is 'silhouette_score'.

        Returns
        -------
        tuple
            Returns a tuple containing two pandas DataFrames with the calculated scores for each clustering approach
            in the control and treated datasets. These DataFrames facilitate further analysis and comparison of the
            scoring metrics across datasets.
        '''
        fig, axs = plt.subplots(1, 2, figsize=(16, 4.5)) 
        score_df1 = plotScore(self.data_ctrl, title = 'Ctrl', ax = axs[0], score_type = score_type, show_legend=False)
        score_df2 = plotScore(self.data_treat, title = 'Treat', ax = axs[1], score_type = score_type, show_legend=True)

        plt.show()
        return score_df1, score_df2




    def showClusterGraph(self,  graph_comb_QC = True):
        '''
        Generate cluster-level graph where nodes represent clusters and edges represent their transcriptomic similarities.

        Returns
        -------
        None
            This function does not return a value but displays the cluster transition graph directly for visual analysis.
        '''
        transition_TPM = self.transition_TPM
        origin_TPM = self.origin_TPM
        self.data_ctrl.obs['SCG_cluster_name'] = ['C' + str(i) for i in self.data_ctrl.obs['SCG_cluster']]
        self.data_treat.obs['SCG_cluster_name'] = ['T' + str(i) for i in self.data_treat.obs['SCG_cluster']]

        self.data_comb.obs['SCG_cluster_name'] = pd.concat([self.data_ctrl.obs['SCG_cluster_name'],
                                                         self.data_treat.obs['SCG_cluster_name']], axis=0)

        clus_ctrl = np.unique(list(self.data_ctrl.obs['SCG_cluster_name']))
        clus_treat = np.unique(list(self.data_treat.obs['SCG_cluster_name']))

        bulk_expr, clus = getClusterExpr(self.data_comb, self.data_comb, key='SCG_cluster_name')
        bulk_expr = pd.DataFrame(bulk_expr.transpose(), index=self.data_comb.var_names, columns=clus)
        bulk_expr_hrg = bulk_expr.loc[self.data_comb.var.highly_variable, :]
        #print(bulk_expr.shape)
        #print(bulk_expr_hrg.shape)

        k_value = 3
        alpha_value = 5
        sim = getSimilarity(bulk_expr_hrg.transpose(), bulk_expr_hrg.transpose(), k=k_value, alpha=alpha_value)
        sim_df = pd.DataFrame(sim, index=bulk_expr_hrg.columns, columns=bulk_expr_hrg.columns)

        cls2cls_sim = np.array(sim_df.loc[clus_ctrl, clus_treat])
        cls2cls_sim_ctrl = np.array(sim_df.loc[clus_ctrl, clus_ctrl])
        cls2cls_sim_treat = np.array(sim_df.loc[clus_treat, clus_treat])

        G, set1, set2, set_pair, weight_list1, weight_list2, weight_list_pair = getGraph(cls2cls_sim, cls2cls_sim_ctrl,
                                                                                         cls2cls_sim_treat, transition_TPM,
                                                                                         origin_TPM, graph_comb_QC = True)
        plotGraph(G, origin_TPM, set1, set2, set_pair, weight_list1, weight_list2, weight_list_pair)

        plt.show()





    def DrugResponse(self):

        '''
        Analyze drug responses by assessing changes in cell abundances and cell states between control and treated
        clusters.

        Parameters
        ----------
        None

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing detailed information on changes in cell abundances and states between aligned 
            control and treated clusters.

    '''
        transition_TPM = self.transition_TPM
        origin_TPM = self.origin_TPM
        
        n_clus = [transition_TPM.shape[0], transition_TPM.shape[1]]
        
        # Cluster size
        num1 = Counter(self.data_ctrl.obs['SCG_cluster'])
        n_cells_ctrl = [num1[str(i)] for i in range(n_clus[0])]
        num2 = Counter(self.data_treat.obs['SCG_cluster'])
        n_cells_treat = [num2[str(i)] for i in range(n_clus[1])]
        
        # Cluster proportion
        perct_clus_ctrl = [n_cells_ctrl[i] / self.data_ctrl.shape[0] for i in range(len(n_cells_ctrl))]
        perct_clus_treat = [n_cells_treat[i] / self.data_treat.shape[0] for i in range(len(n_cells_treat))]

        # Flow proportion
        transition_TPM_mat = transition_TPM.copy()
        origin_TPM_mat = origin_TPM.copy()
        for i in range(n_clus[0]):
            transition_TPM_mat[i, :] = perct_clus_ctrl[i] * transition_TPM[i, :]

        for j in range(n_clus[1]):
            origin_TPM_mat[:, j] = perct_clus_treat[j] * origin_TPM[:, j]

        #ctrl_clus_perct1 = np.sum(transition_TPM_mat, axis=1)
        #treat_clus_perct2 = np.sum(origin_TPM_mat, axis=0)
        
        # calculate STC
        self.data_ctrl.obs['SCG_cluster_name'] = ['C' + str(i) for i in self.data_ctrl.obs['SCG_cluster']]
        self.data_treat.obs['SCG_cluster_name'] = ['T' + str(i) for i in self.data_treat.obs['SCG_cluster']]

        self.data_comb.obs['SCG_cluster_name'] = pd.concat([self.data_ctrl.obs['SCG_cluster_name'],
                                                         self.data_treat.obs['SCG_cluster_name']], axis=0)
    
        clus_ctrl = np.unique(list(self.data_ctrl.obs['SCG_cluster_name']))
        clus_treat = np.unique(list(self.data_treat.obs['SCG_cluster_name']))
        
        bulk_expr, clus = getClusterExpr(self.data_comb, self.data_comb, key='SCG_cluster_name')
        bulk_expr = pd.DataFrame(bulk_expr.transpose(), index=self.data_comb.var_names, columns=clus)
        bulk_expr_hrg = bulk_expr.loc[self.data_comb.var.highly_variable, :]
        
        k_value = 3
        alpha_value = 5
        sim = getSimilarity(bulk_expr_hrg.transpose(), bulk_expr_hrg.transpose(), k=k_value, alpha=alpha_value)
        sim_df = pd.DataFrame(sim, index=bulk_expr_hrg.columns, columns=bulk_expr_hrg.columns)

        cls2cls_dis = 1 - np.array(sim_df.loc[clus_ctrl, clus_treat])
        cls2cls_dis_ctrl = 1 - np.array(sim_df.loc[clus_ctrl, clus_ctrl])
        cls2cls_dis_treat = 1 - np.array(sim_df.loc[clus_treat, clus_treat])
        G2 = getDiGraph(cls2cls_dis, cls2cls_dis_ctrl, cls2cls_dis_treat, transition_TPM, origin_TPM)

        source = []
        target = []
        flow = []
        prob_ctrl = []
        prob_treat = []
        flow_ctrl_percent = []
        flow_treat_percent = []
        STC = []
        RAC = []

        for i in range(transition_TPM.shape[0]):
            for j in range(transition_TPM.shape[1]):
                if (transition_TPM[i, j] > 0 or origin_TPM[i, j] > 0):
                    source.append('C' + str(i))
                    target.append('T' + str(j))
                    flow.append('C' + str(i) + '->' + 'T' + str(j))
                    prob_ctrl.append(transition_TPM[i, j])
                    prob_treat.append(origin_TPM[i, j])
                    flow_ctrl_percent.append(transition_TPM_mat[i, j])
                    flow_treat_percent.append(origin_TPM_mat[i, j])
                    #path = nx.dijkstra_path(G2, source='C' + str(i), target='T' + str(j))
                    weight = nx.dijkstra_path_length(G2, source='C' + str(i), target='T' + str(j))
                    STC.append(weight)
                    RAC.append(origin_TPM_mat[i, j] / transition_TPM_mat[i, j]  )

        flow_info = pd.DataFrame({'Flow': flow,
                                  'Source': source,
                                  'Target': target,
                                  'Probability of Flow in Control Cluster': prob_ctrl,
                                  'Probability of Flow in Treated Cluster': prob_treat,
                                  'Percent of Flow in Control Sample': flow_ctrl_percent,
                                  'Percent of Flow in Treated Sample': flow_treat_percent,
                                  'STC': STC,
                                  'RAC': RAC
                                  })

        #print(flow_info)
        return flow_info






    def sankey(self, DrugResponseInfo):
        '''
        Visualize the flow of cellular transitions between clusters from control to treated conditions using a Sankey diagram.

        Parameters
        ----------
        DrugResponseInfo : pandas.DataFrame
            A DataFrame containing detailed information on changes in cell abundances and cell states between control and 
            treated clusters. Derived from the `DrugResponse` step. 

        Returns
        -------
        None
            The function renders a Sankey diagram directly and does not return any values. It serves as a visual aid
            to complement quantitative analyses.
        '''
        transition_TPM = self.transition_TPM
        origin_TPM = self.origin_TPM
        
        n_clus = [transition_TPM.shape[0], transition_TPM.shape[1]]
        
        num1 = Counter(self.data_ctrl.obs['SCG_cluster'])
        n_cells_ctrl = [num1[str(i)] for i in range(n_clus[0])]
        num2 = Counter(self.data_treat.obs['SCG_cluster'])
        n_cells_treat = [num2[str(i)] for i in range(n_clus[1])]
        
        perct_clus_ctrl = [n_cells_ctrl[i] / self.data_ctrl.shape[0] for i in range(len(n_cells_ctrl))]
        perct_clus_treat = [n_cells_treat[i] / self.data_treat.shape[0] for i in range(len(n_cells_treat))]
        
        ctrl_fractions = perct_clus_ctrl
        treat_fractions = perct_clus_treat
        fig = plt.figure(figsize=(6, 6))
        gs = GridSpec(DrugResponseInfo.shape[0], 2, width_ratios=[2.1, 1], hspace=0.1, wspace=0.7)
        ax = fig.add_subplot(gs[0:DrugResponseInfo.shape[0], 0])

        for i in range(n_clus[0]):
            bottom = np.sum(ctrl_fractions[(i + 1):])
            rectangle = ax.bar(x=[0], height=ctrl_fractions[i], bottom=bottom, color=self.ctrl_colors[i],
                               edgecolor='black', fill=True, linewidth=0.7, width=0.16)
            text_y =bottom + ctrl_fractions[i] / 2
            ax.text(x=0, y=text_y, s=f'C{i}', color='black',  
            ha='center', va='center', fontsize=10)  

        for i in range(n_clus[1]):
            bottom =np.sum(treat_fractions[(i + 1):])
            rectangle = ax.bar(x=[1], height=treat_fractions[i], bottom=bottom, color=self.treat_colors[i],
                               edgecolor='black', fill=True, linewidth=0.7, width=0.16)
            text_y =bottom + treat_fractions[i] / 2
            ax.text(x=1, y=text_y, s=f'T{i}', color='black',  
            ha='center', va='center', fontsize=10)  


#         legend_labels_ctrl = [f"Ctrl{i}" for i in range(n_clus[0])]
#         legend_handles_ctrl = [Rectangle((0, 0), 1, 1, color=color, ec='black') for color in self.ctrl_colors]

#         legend_labels_treat = [f"Treat{i}" for i in range(n_clus[1])]
#         legend_handles_treat = [Rectangle((0, 0), 1, 1, color=color, ec='black') for color in self.treat_colors]

#         all_legend_labels = legend_labels_ctrl + legend_labels_treat
#         all_legend_handles = legend_handles_ctrl + legend_handles_treat

#         fig.tight_layout()
#         ax.legend(all_legend_handles, all_legend_labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)

        for pos in ['right', 'top', 'bottom', 'left']:
            ax.spines[pos].set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.text(x=0, y=-0.07, s='Ctrl', horizontalalignment='center', verticalalignment='center', fontsize=13)
        ax.text(x=1, y=-0.07, s='Treat', horizontalalignment='center', verticalalignment='center', fontsize=13)

        xs = np.linspace(-5, 5, num=100)
        ys = np.array([sigmoid(x) for x in xs])
        xs = xs / 10 + 0.5
        xs *= 0.83
        xs += 0.085

        y_start_record = [1 - np.sum(ctrl_fractions[0:ii]) for ii in range(len(ctrl_fractions))]
        y_end_record = [1 - np.sum(treat_fractions[0:ii]) for ii in range(len(treat_fractions))]
        

        y_up_start, y_dw_start = 1, 1
        y_up_end, y_dw_end = 1, 1
        axi = 0
        for si in range(n_clus[0]):
            cur_p_df = DrugResponseInfo.loc[DrugResponseInfo['Source'] ==  'C'+str(si), :]
            if cur_p_df.shape[0] > 0:
                for fi in range(cur_p_df.shape[0]):
                    y_up_start = y_start_record[si]
                    y_dw_start = y_up_start - cur_p_df['Percent of Flow in Control Sample'].iloc[fi]
                    y_start_record[si] = y_dw_start

                    ti = cur_p_df['Target'].iloc[fi]
                    ti = int(ti[1:])
                    y_up_end = y_end_record[ti]
                    y_dw_end = y_up_end - cur_p_df['Percent of Flow in Treated Sample'].iloc[fi]
                    y_end_record[ti] = y_dw_end

                    y_up_start -= 0.003
                    y_dw_start += 0.004
                    y_up_end -= 0.003
                    y_dw_end += 0.004

                    ys_up = y_up_start + (y_up_end - y_up_start) * ys
                    ys_dw = y_dw_start + (y_dw_end - y_dw_start) * ys

                    color_s_t = [self.ctrl_colors[si], self.treat_colors[ti]]
                    cmap = LinearSegmentedColormap.from_list('mycmap', [color_s_t[0], color_s_t[1]])
                    grad_colors = cmap(np.linspace(0, 1, len(xs) - 1))
                    grad_colors = [rgb2hex(color) for color in grad_colors]
                    for pi in range(len(xs) - 1):
                        ax.fill_between(xs[pi:(pi + 2)], ys_dw[pi:(pi + 2)], ys_up[pi:(pi + 2)], alpha=0.7,
                                        color=grad_colors[pi], edgecolor=None)


                        
                        
                        
    def sankey_overall(self, DrugResponseInfo, Vol_Ratio):
        '''
        Generate a Sankey diagram that incorporates the Volume Ratio index to visualize the
        flow of cellular transitions between clusters, adjusting for absolute abundance changes. This visualization 
        highlights how treatment affects cell population dynamics, taking into account the reduction in tumor volume 
        or cell count.

        Parameters
        ----------
        DrugResponseInfo : pandas.DataFrame
            A DataFrame containing detailed information on changes in cell abundances and cell states between control and 
            treated clusters. Derived from the `DrugResponse` step. 

        Vol_Ratio : float
            The Volume ratio index, a measure used to quantify the effectiveness of a treatment in inhibiting
            tumor growth. This value adjusts the flow volumes in the Sankey diagram to reflect the reduced cell counts
            or tumor volumes post-treatment.

        Returns
        -------
        None
            The function renders a Sankey diagram directly to the output device (e.g., screen or plot window) and
            does not return any values. It provides a dynamic, intuitive means of exploring data that combines both
            quantitative and qualitative analysis.
        '''
                
        frac = Vol_Ratio
        transition_TPM = self.transition_TPM
        origin_TPM = self.origin_TPM
        
        n_clus = [transition_TPM.shape[0], transition_TPM.shape[1]]
        
        num1 = Counter(self.data_ctrl.obs['SCG_cluster'])
        n_cells_ctrl = [num1[str(i)] for i in range(n_clus[0])]
        num2 = Counter(self.data_treat.obs['SCG_cluster'])
        n_cells_treat = [num2[str(i)] for i in range(n_clus[1])]
        
        perct_clus_ctrl = [n_cells_ctrl[i] / self.data_ctrl.shape[0] for i in range(len(n_cells_ctrl))]
        perct_clus_treat = [n_cells_treat[i] / self.data_treat.shape[0] for i in range(len(n_cells_treat))]
        
        ctrl_fractions = perct_clus_ctrl
        treat_fractions = perct_clus_treat
        fig = plt.figure(figsize=(6, 6))
        gs = GridSpec(DrugResponseInfo.shape[0], 2, width_ratios=[2.1, 1], hspace=0.1, wspace=0.7)
        ax = fig.add_subplot(gs[0:DrugResponseInfo.shape[0], 0])

        for i in range(n_clus[0]):
            bottom = np.sum(ctrl_fractions[(i + 1):])
            rectangle = ax.bar(x=[0], height=ctrl_fractions[i], bottom=bottom, color=self.ctrl_colors[i],
                               edgecolor='black', fill=True, linewidth=0.7, width=0.16)
            text_y =bottom + ctrl_fractions[i] / 2
            ax.text(x=0, y=text_y, s=f'C{i}', color='black',  
            ha='center', va='center', fontsize=10)  

        for i in range(n_clus[1]):
            bottom =np.sum(treat_fractions[(i + 1):])*frac
            rectangle = ax.bar(x=[1], height=treat_fractions[i]*frac, bottom=bottom, color=self.treat_colors[i],
                               edgecolor='black', fill=True, linewidth=0.7, width=0.16)
            text_y =bottom + treat_fractions[i]*frac / 2
            ax.text(x=1, y=text_y, s=f'T{i}', color='black',  
            ha='center', va='center', fontsize=10)  


#         legend_labels_ctrl = [f"Ctrl{i}" for i in range(n_clus[0])]
#         legend_handles_ctrl = [Rectangle((0, 0), 1, 1, color=color, ec='black') for color in self.ctrl_colors]

#         legend_labels_treat = [f"Treat{i}" for i in range(n_clus[1])]
#         legend_handles_treat = [Rectangle((0, 0), 1, 1, color=color, ec='black') for color in self.treat_colors]

#         all_legend_labels = legend_labels_ctrl + legend_labels_treat
#         all_legend_handles = legend_handles_ctrl + legend_handles_treat

#         fig.tight_layout()
#         ax.legend(all_legend_handles, all_legend_labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)

        for pos in ['right', 'top', 'bottom', 'left']:
            ax.spines[pos].set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.text(x=0, y=-0.07, s='Ctrl', horizontalalignment='center', verticalalignment='center', fontsize=13)
        ax.text(x=1, y=-0.07, s='Treat', horizontalalignment='center', verticalalignment='center', fontsize=13)

        xs = np.linspace(-5, 5, num=100)
        ys = np.array([sigmoid(x) for x in xs])
        xs = xs / 10 + 0.5
        xs *= 0.83
        xs += 0.085

        y_start_record = [1 - np.sum(ctrl_fractions[0:ii]) for ii in range(len(ctrl_fractions))]
        y_end_record = [(1 - np.sum(treat_fractions[0:ii]))*frac  for ii in range(len(treat_fractions))]
        
        y_up_start, y_dw_start = 1, 1
        y_up_end, y_dw_end = 1, 1
        axi = 0
        for si in range(n_clus[0]):
            cur_p_df = DrugResponseInfo.loc[DrugResponseInfo['Source'] ==  'C'+str(si), :]
            if cur_p_df.shape[0] > 0:
                for fi in range(cur_p_df.shape[0]):
                    y_up_start = y_start_record[si]
                    y_dw_start = y_up_start - cur_p_df['Percent of Flow in Control Sample'].iloc[fi]
                    y_start_record[si] = y_dw_start

                    ti = cur_p_df['Target'].iloc[fi]
                    ti = int(ti[1:])
                    y_up_end = y_end_record[ti]
                    y_dw_end = y_up_end - frac * cur_p_df['Percent of Flow in Treated Sample'].iloc[fi]
                    y_end_record[ti] = y_dw_end

                    y_up_start -= 0.000003
                    y_dw_start += 0.000004
                    y_up_end -= 0.000003
                    y_dw_end += 0.000004

                    ys_up = y_up_start + (y_up_end - y_up_start) * ys
                    ys_dw = y_dw_start + (y_dw_end - y_dw_start) * ys

                    color_s_t = [self.ctrl_colors[si], self.treat_colors[ti]]
                    cmap = LinearSegmentedColormap.from_list('mycmap', [color_s_t[0], color_s_t[1]])
                    grad_colors = cmap(np.linspace(0, 1, len(xs) - 1))
                    grad_colors = [rgb2hex(color) for color in grad_colors]
                    for pi in range(len(xs) - 1):
                        ax.fill_between(xs[pi:(pi + 2)], ys_dw[pi:(pi + 2)], ys_up[pi:(pi + 2)], alpha=0.7,
                                        color=grad_colors[pi], edgecolor=None)




    def plot_degree_distributions(
        self,
        object_type = 'Ctrl',
        alpha_list=None,
        k_list=None,
        threshold=1e-40,
        bins=100,
        figsize=(15, 13),
        save_plot=False,
        save_path="degree_distributions.png"
    ):
        """
        绘制不同 alpha (alpha_list) 和 k (k_list) 值下的度分布直方图。

        参数：
        --------
        data : np.ndarray
            用于计算相似度矩阵的输入数据 (e.g. principal components)。
        alpha_list : list, optional
            要测试的 alpha 值。默认 [5, 10, 20, 30, 50]
        k_list : list, optional
            要测试的 k 值。默认 [2, 3, 5, 10]
        threshold : float, optional
            判断是否有边的阈值，默认 1e-40
        bins : int, optional
            绘图直方图的 bins 数量，默认 100
        figsize : tuple, optional
            整体图像的大小，默认 (15, 13)
        save_plot : bool, optional
            是否保存生成的图像。默认 False（不保存）。
        save_path : str, optional
            当 save_plot=True 时，保存图像的文件路径或文件名。默认 "degree_distributions.png"

        返回值：
        --------
        fig, axes : (matplotlib.figure.Figure, np.ndarray)
            返回创建的 Figure 和子图 axes 数组，便于后续操作或保存。
        """
        if object_type == 'Ctrl':
            data_use = self.data_ctrl.obsm['X_pca'].copy()
        elif object_type == 'Treat':
            data_use = self.data_treat.obsm['X_pca'].copy()
        elif object_type =='Comb':
            data_use = self.data_comb.obsm['X_pca'].copy()

        if alpha_list is None:
            alpha_list = [5, 10, 20, 30, 50]
        if k_list is None:
            k_list = [2, 3, 5, 10]

        # 创建子图
        fig, axes = plt.subplots(len(alpha_list), len(k_list), figsize=figsize, sharex=True, sharey=True)

        # 如果只有一个 alpha 或一个 k 的情况，需要保证 axes 可正确索引
        if len(alpha_list) == 1 and len(k_list) == 1:
            # 只有一个子图
            axes = np.array([[axes]])
        elif len(alpha_list) == 1:
            # 只有一行
            axes = np.array([axes])
        elif len(k_list) == 1:
            # 只有一列
            axes = axes[:, np.newaxis]
        
        for i, alpha in enumerate(alpha_list):
            for j, k in enumerate(k_list):
                # 计算相似度矩阵
                
                affinity_mat = getSimilarity(data_use, data_use, k, alpha)

                # 计算每个细胞(行)对应的度（大于 threshold 即视为有边）
                degree_arr = np.sum(affinity_mat > threshold, axis=1)
                median_degree = int(np.median(degree_arr))

                # 绘制直方图
                ax = axes[i, j]
                ax.hist(degree_arr, bins=bins, color='blue', alpha=0.7)
                ax.set(
                    title=f'α={alpha}, k={k}, Median={median_degree}',
                    xlabel='Degree',
                    ylabel='Frequency'
                )
                ax.grid(True)

        plt.tight_layout()

        # 如果用户希望保存图像，则执行保存
        if save_plot:
            fig.savefig(save_path, dpi=300)  # 可以根据需求调整 dpi
        
        plt.show()
