import networkx as nx
import matplotlib.pyplot as plt
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import distance
from scipy import sparse
import random
import community
import time
import dgl
import time
import scanpy as sc
from sklearn.preprocessing import normalize
from matplotlib.sankey import Sankey
from matplotlib.gridspec import GridSpec
from sklearn import metrics
from sklearn.mixture import GaussianMixture
from collections import Counter
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, rgb2hex
from matplotlib.patches import Rectangle
import subprocess
import os

def createScConGraphObj(sc_obj_comb, run_label, key, pre_name, pos_name, pre_colors=None, pos_colors=None,
                        cls_prefixes=['C', 'T']):
    data_ctrl, data_treat = splitScObjects(sc_obj_comb, sam_key=key, sam_values=[pre_name, pos_name])
    scg = scConGraph(sc_obj_comb, data_ctrl, data_treat, key, pre_name, pos_name, run_label,pre_colors, pos_colors,
                     cls_prefixes)
    return (scg)

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def splitScObjects(scobj, sam_key, sam_values):
    data_split = [scobj[scobj.obs[sam_key] == v,:] for v in sam_values]
    return(data_split)


def getSimilarity(PC_data1, PC_data2, k, alpha):
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


def runLINE(path, edges, output, order=1, size=128, negative=5, samples=500):
    LINE = path
    cmd = [
        LINE,
        '-train', edges,
        '-output', output,
        '-binary', '0',
        '-size', str(size),
        '-order', str(order),
        '-negative', str(negative),
        '-samples', str(samples),
        '-threads', '20'
    ]
    print('Command line - Run LINE: %s' % ' '.join(cmd))

    subprocess.run(cmd)

def norm_embedding(embed_1st = None, embed_2nd_1 = None, embed_2nd_2 = None, mode = 0, weight = [0.5, 0.5]):
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

def runLineClustering(obj, embed, key = 'LINE', resolution = 0.5):
    obj.obsm[key] = embed
    sc.pp.neighbors(obj, use_rep = key)
    sc.tl.umap(obj)
    sc.tl.leiden(obj, key_added = key+'_cluster', resolution = resolution)
    return obj


def cls2cls_prob(data_pre, data_pos, weight_df,
                 thre_cell=0.05, thre_cls=0.05):
    # probability of cell to cell
    n_clus = np.array([len(data_pre.obs["LINE_cluster"].unique()), len(data_pos.obs["LINE_cluster"].unique())])
    n_cells = np.array([data_pre.obs.shape[0], data_pos.obs.shape[0]])
    cell2cell = normalize(np.array(weight_df), norm='l1')

    # probability of cell to cluster
    cell_cls_pre = data_pre.obs['LINE_cluster']
    cell_cls_pos = data_pos.obs['LINE_cluster']
    cls_cellnum_pos = [np.sum(cell_cls_pos == str(j)) for j in range(n_clus[1])]

    cell2cls = [np.sum(cell2cell[:, np.where(cell_cls_pos == str(j))[0]], axis=1) / cls_cellnum_pos[j] for j in
                range(n_clus[1])]
    cell2cls = np.transpose(np.array(cell2cls))
    cell2cls = normalize(np.array(cell2cls), norm='l1')
    cell2cls[cell2cls < thre_cell] = 0
    cell2cls = normalize(np.array(cell2cls), norm='l1')

    # probability of cluster to cluster
    cls2cls = [np.sum(cell2cls[np.where(cell_cls_pre == str(i))[0], :], axis=0) for i in range(n_clus[0])]
    cls2cls = normalize(np.array(cls2cls), norm='l1')

    cls2cls[cls2cls < thre_cls] = 0
    cls2cls = normalize(cls2cls, norm='l1')
    return (cls2cls)


def rapid_small(align_pre, align_pos):
    n_clus = align_pre.shape
    for i in range(n_clus[0]):
        for j in range(n_clus[1]):
            if align_pre[i, j] + align_pos[i, j] < 0.1:
                align_pre[i, j] = 0
                align_pos[i, j] = 0
                align_pre = normalize(align_pre, norm='l1', axis=1)
                align_pos = normalize(align_pos, norm='l1', axis=0)
    return align_pre, align_pos

def plotUMAP(sc_obj, key = 'umap_cell_embeddings', group_by = 'leiden', title = '', label ='UMAP', size = 15):
    plot_data = pd.DataFrame({'UMAP1':sc_obj.obsm[key][:,0],
                              'UMAP2':sc_obj.obsm[key][:,1],
                              'Cluster':list(sc_obj.obs[ group_by ])})
    plot_data = plot_data.sort_values(by='Cluster', axis=0, ascending=True)
    ax = sns.scatterplot(data=plot_data, x="UMAP1", y="UMAP2", hue="Cluster", linewidth=0, s=size)
    ax.legend(loc = 2, bbox_to_anchor = (1.01,1), ncol=1)
    ax.set(title=title, xlabel = label+'1', ylabel = label+'2')
    return(ax)


def plotScore(data, type, Time):
    size = 3

    cluster_single = data.obs['leiden_single']
    cluster_whole = data.obs['leiden_all']
    cluster_ALINE = data.obs['LINE_cluster']

    umap = ['X_umap_all', 'X_umap', 'LINE_umap']
    if type == 'silhouette_score':
        a = np.array([metrics.silhouette_score(data.obsm[i], cluster_whole) for i in umap])
        b = np.array([metrics.silhouette_score(data.obsm[i], cluster_single) for i in umap])
        c = np.array([metrics.silhouette_score(data.obsm[i], cluster_ALINE) for i in umap])
        print(a, b, c)
    if type == 'davies_bouldin_score':
        a = np.array([metrics.davies_bouldin_score(data.obsm[i], cluster_whole) for i in umap])
        b = np.array([metrics.davies_bouldin_score(data.obsm[i], cluster_single) for i in umap])
        c = np.array([metrics.davies_bouldin_score(data.obsm[i], cluster_ALINE) for i in umap])
        print(a, b, c)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(size)
    x_labels = ["UMAP_whole", "UMAP_single", "UMAP_scConGraph"]
    plt.xticks(x, x_labels)

    total_width, n = 0.8, size
    width = total_width / n
    x = x - (total_width - width) / 2

    ax.bar(x, a, width=width, label="Cluster_whole")
    ax.bar(x + width, b, width=width, label="Cluster_single")
    ax.bar(x + 2 * width, c, width=width, label="scConGraph_cluster")

    ax.legend()
    ax.set_title(Time + ': ' + type)
    plt.show()

def getMetacellExpr(scobj, adata, key):
    expr = adata[scobj.obs.index,:].copy()
    expr_data = expr.X.toarray()
    clqs = np.unique(list(scobj.obs[key]))
    clq_pc_data = [np.average(expr_data[np.where(scobj.obs[key] == x)[0], :], axis=0) for x in clqs]
    clq_pc_data = np.array(clq_pc_data)
    return (clq_pc_data, clqs)

def getGraph(cls2cls_sim, cls2cls_sim_tp0, cls2cls_sim_tp1, align_pre, align_pos):
    import networkx as nx
    G = nx.Graph()
    weight_list=[]
    set1 = []
    set2 = []
    set_pair = []

    weight_list1 = []
    weight_list2= []
    weight_list_pair = []

    for i in range(align_pos.shape[0]):
        for j in range(align_pos.shape[1]):
            if(align_pos[i, j] > 0  or align_pre[i, j] > 0):
                nodes1 = 'C_'+str(i)
                nodes2 = 'T_'+str(j)
                weight = cls2cls_sim[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                set_pair.append((nodes1, nodes2))
                weight_list_pair.append(weight)

    for i in range(align_pos.shape[0]):
        for j in range(align_pos.shape[0]):
            if (i != j ):
                nodes1 = 'C_'+str(i)
                nodes2 = 'C_'+str(j)
                weight = cls2cls_sim_tp0[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                set1.append((nodes1, nodes2))
                weight_list1.append(weight)

    for i in range(align_pos.shape[1]):
        for j in range(align_pos.shape[1]):
            if (i != j ):
                nodes1 = 'T_'+str(i)
                nodes2 = 'T_'+str(j)
                weight = cls2cls_sim_tp1[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                set2.append((nodes1, nodes2))
                weight_list2.append(weight)
    return G, set1, set2, set_pair, weight_list1, weight_list2, weight_list_pair

def plotGraph(G, align_pos, set1, set2, set_pair,  weight_list1, weight_list2, weight_list_pair):
    fig = plt.figure(figsize=(5, 3), dpi=300)
    pos2 = nx.spring_layout(G, iterations=200, seed=2024)

    options = { "node_size": 600, "alpha": 0.9}
    nx.draw_networkx_nodes(G, pos2, nodelist=['C_'+str(i) for i in range(align_pos.shape[0])], node_color="tab:blue", **options)
    nx.draw_networkx_nodes(G, pos2, nodelist=['T_'+str(i) for i in range(align_pos.shape[1])], node_color="tab:red", **options)

    nx.draw_networkx_edges(G,
        pos2,
        edgelist= set_pair,
        width= 2,
        alpha=0.5,
        edge_color = weight_list_pair,
        edge_cmap = plt.cm.YlGn
                          )


    nx.draw_networkx_edges(
        G,
        pos2,
        edgelist=set2,
        width= 2,
        alpha=0.5,
        edge_color = weight_list2,
        edge_cmap = plt.cm.Reds
    )


    nx.draw_networkx_edges(
        G,
        pos2,
        edgelist= set1,
        width=2,
        alpha=0.5,
        edge_color = weight_list1,
        edge_cmap = plt.cm.Blues
    )

    nx.draw_networkx_labels(G, pos2, font_size=4, font_color="whitesmoke")

    plt.tight_layout()
    plt.axis("off")


def GMM_flow(flow_info):
    ks = np.arange(1, 4)
    models = [GaussianMixture(k, random_state=0).fit(np.array(flow_info['Min_weight']).reshape(-1, 1)) for k in ks]

    gmm_model_comp = pd.DataFrame({"ks": ks,
                                   "BIC": [m.bic(np.array(flow_info['Min_weight']).reshape(-1, 1)) for m in models],
                                   "AIC": [m.aic(np.array(flow_info['Min_weight']).reshape(-1, 1)) for m in models]})
    return (models, gmm_model_comp)


def assignFlowType(flow_info, models, k):
    gmm_model = models[k - 1]
    gmm_types = gmm_model.predict(np.array(flow_info['Min_weight']).reshape(-1, 1))
    gmm_means = gmm_model.means_.reshape(1, -1)
    flow_info['GMM_Type'] = \
    np.array(['almost no change', 'change a little', 'change a lot'])[np.argsort(np.argsort(gmm_means))[0]][gmm_types]
    return (flow_info, gmm_model)


def judgePercentVariance(perct):
    cluster_type = []
    for i in range(len(perct)):
        if perct[i] >= 1.2:
            cluster_type.append('Increase')
        elif perct[i] < 1.2 and perct[i] > 0.8:
            cluster_type.append('Invariant')
        elif perct[i] <= 0.8:
            cluster_type.append('Decrease')
    return cluster_type

def getGraph_notdel_dir(cls2cls_sim, cls2cls_sim_tp0, cls2cls_sim_tp1, align_pre, align_pos):

    G = nx.DiGraph()

    for i in range(align_pos.shape[0]):
        for j in range(align_pos.shape[1]):
            if True:
                nodes1 = 'C_'+str(i)
                nodes2 = 'T_'+str(j)
                weight = cls2cls_sim[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)

    for i in range(align_pos.shape[0]):
        for j in range(align_pos.shape[0]):
            if (i != j ):
                nodes1 = 'C_'+str(i)
                nodes2 = 'C_'+str(j)
                weight = cls2cls_sim_tp0[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                G.add_edge(nodes2, nodes1, weight = weight)

    for i in range(align_pos.shape[1]):
        for j in range(align_pos.shape[1]):
            if (i != j ):
                nodes1 = 'T_'+str(i)
                nodes2 = 'T_'+str(j)
                weight = cls2cls_sim_tp1[i, j]
                G.add_edge(nodes1, nodes2, weight = weight)
                G.add_edge(nodes2, nodes1, weight = weight)

    return G


class scConGraph:

    def __init__(self, data_comb, data_ctrl, data_treat, key, pre_name, pos_name, run_label, pre_colors=None,
                 pos_colors=None,
                 cls_prefixes=['C', 'G']):
        self.data_comb = data_comb
        self.data_ctrl = data_ctrl
        self.data_treat = data_treat
        self.key = key
        self.pre_name = pre_name
        self.pos_name = pos_name
        self.run_label = run_label
        self.pre_colors = pre_colors
        self.pos_colors = pos_colors
        if self.pre_colors is None:
            self.pre_colors = ["#FEC643", "#437BFE", "#43FE69", "#FE6943", "#E78AC3",
                               "#43D9FE", "#FFEC1A", "#E5C494", "#A6D854", "#33AEB1",
                               "#EA6F5D", "#FEE8C3", "#3DBA79", "#D6EBF9", "#7F6699",
                               "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                               "#333399", "#679966", "#c12e34", "#66d5a5", "#5599ff"]
        if self.pos_colors is None:
            self.pos_colors = ['#EE854A', '#4878D0', '#6ACC64', '#D65F5F', '#956CB4',
                               '#82C6E2', '#D5BB67', '#8C613C', '#DC7EC0', '#797979']
        self.cls_prefixes = cls_prefixes





    # proProcess(PCA)
    def proProcess(scobj, npcs=40):
        # scobj.raw = scobj
        # scobj = scobj[:, scobj.var.highly_variable]
        #
        # sc.pp.regress_out(scobj, ['nCount_RNA', 'percent_mito', 'percent_ribo'])
        # sc.pp.scale(scobj, max_value=10)

        sc.tl.pca(scobj.data_comb, svd_solver='arpack')
        sc.pp.neighbors(scobj.data_comb, n_neighbors=10, n_pcs=npcs)
        sc.tl.umap(scobj.data_comb)

        sc.tl.pca(scobj.data_ctrl, svd_solver='arpack')
        sc.pp.neighbors(scobj.data_ctrl, n_neighbors=10, n_pcs=npcs)
        sc.tl.umap(scobj.data_ctrl)

        sc.tl.pca(scobj.data_treat, svd_solver='arpack')
        sc.pp.neighbors(scobj.data_treat, n_neighbors=10, n_pcs=npcs)
        sc.tl.umap(scobj.data_treat)

        return (scobj)

    # runCluster
    def runClustering(scobj, res=0.5):
        sc.tl.leiden(scobj.data_ctrl, key_added='leiden_single', resolution=res)
        sc.tl.leiden(scobj.data_treat, key_added='leiden_single', resolution=res)
        sc.tl.leiden(scobj.data_comb, key_added='leiden_all', resolution=res)

        scobj.data_ctrl.obs['leiden_all'] = scobj.data_comb.obs['leiden_all'][
            scobj.data_comb.obs[scobj.key] == scobj.pre_name]
        scobj.data_treat.obs['leiden_all'] = scobj.data_comb.obs['leiden_all'][
            scobj.data_comb.obs[scobj.key] == scobj.pos_name]
        scobj.data_comb.obs['leiden_single'] = pd.concat([scobj.data_ctrl.obs['leiden_single'],
                                                          scobj.data_treat.obs['leiden_single']], axis=0)

        scobj.data_ctrl.obsm['X_umap_all'] = scobj.data_comb.obsm['X_umap'][
                                             scobj.data_comb.obs[scobj.key] == scobj.pre_name, :]
        scobj.data_treat.obsm['X_umap_all'] = scobj.data_comb.obsm['X_umap'][
                                              scobj.data_comb.obs[scobj.key] == scobj.pos_name, :]

    def calculateAffinity(scobj,  resultpath,runLabel,k=3, alpha=20):
        PC_data_B1 = scobj.data_ctrl.obsm['X_pca']
        PC_data_B2 = scobj.data_treat.obsm['X_pca']
        PC_data_all = scobj.data_comb.obsm['X_pca']
        n_clus = np.array(
            [len(scobj.data_ctrl.obs["leiden_single"].unique()), len(scobj.data_treat.obs["leiden_single"].unique())])
        n_cells = np.array([scobj.data_ctrl.shape[0], scobj.data_treat.obs.shape[0]])
        PC_data_all_B1 = PC_data_all[0: n_cells[0], :]

        PC_data_all_B2 = PC_data_all[n_cells[0]: (n_cells[0] + n_cells[1]), :]

        global_sim_B1 = getSimilarity(PC_data_B1, PC_data_B1, k=k, alpha=alpha)
        global_sim_B2 = getSimilarity(PC_data_B2, PC_data_B2, k=k, alpha=alpha)
        pair_sim_B1B2 = getSimilarity(PC_data_all_B1, PC_data_all_B2, k=k, alpha=alpha)

        np.savetxt(resultpath + runLabel + '_CG_similarity_alpha' + str(alpha) + '_raw.txt', pair_sim_B1B2)
        np.savetxt(resultpath + runLabel + '_Ctrl_similarity_alpha' + str(alpha) + '_raw.txt', global_sim_B1)
        np.savetxt(resultpath + runLabel + '_Gemc_similarity_alpha' + str(alpha) + '_raw.txt', global_sim_B2)

        rows, cols = np.where(pair_sim_B1B2 != 0)
        edge_weight_comb = [(f'C{i}', f'G{j}', pair_sim_B1B2[i, j]) for i, j in zip(rows, cols)]

        rows, cols = np.where(np.triu(global_sim_B1, k=1) != 0)
        edge_weight_ctrl = [(f'C{i}', f'C{j}', global_sim_B1[i, j]) for i, j in zip(rows, cols)]

        rows, cols = np.where(np.triu(global_sim_B2, k=1) != 0)
        edge_weight_treat = [(f'G{i}', f'G{j}', global_sim_B2[i, j]) for i, j in zip(rows, cols)]

        np.savetxt(resultpath + runLabel + '_CG_edge_weight_alpha' + str(alpha) + '_raw.txt', edge_weight_comb,
                   fmt='%s', delimiter='\t')
        np.savetxt(resultpath + runLabel + '_C_edge_weight_alpha' + str(alpha) + '_raw.txt', edge_weight_ctrl, fmt='%s',
                   delimiter='\t')
        np.savetxt(resultpath + runLabel + '_G_edge_weight_alpha' + str(alpha) + '_raw.txt', edge_weight_treat,
                   fmt='%s', delimiter='\t')

        return pair_sim_B1B2, global_sim_B1, global_sim_B2

    def ALINE(path, resultpath,runLabel, order=1, size_1ord=100, size_2ord=100, negative=5, samples=500, alpha=20):
        time_start = time.time()
        edges = resultpath + runLabel + '_C_edge_weight_alpha' + str(alpha) + '_raw.txt'
        output = resultpath + runLabel + '_Ctrl_1st_embed' + str(size_1ord) + '.txt'
        runLINE(path, edges=edges, output=output, order=1, size=size_1ord, negative=negative, samples=samples)

        edges = resultpath + runLabel + '_G_edge_weight_alpha' + str(alpha) + '_raw.txt'
        output = resultpath + runLabel + '_Gemc_1st_embed' + str(size_1ord) + '.txt'
        runLINE(path, edges=edges, output=output, order=order, size=size_1ord, negative=negative, samples=samples)

        edges = resultpath + runLabel + '_CG_edge_weight_alpha' + str(alpha) + '_raw.txt'
        output = resultpath + runLabel + '_2nd_embed' + str(size_2ord) + '.txt'
        runLINE(path, edges=edges, output=output, order=order, size=size_2ord, negative=negative, samples=samples)

        time_end = time.time()
        print('time cost', time_end - time_start, 's')

    def ConGraphCluster(scobj, resultpath,runLabel, size_1ord=100, size_2ord=100, negative=5, samples=500, alpha=20,
                        res=0.3):
        n_cells = np.array([scobj.data_ctrl.obs.shape[0], scobj.data_treat.obs.shape[0]])
        cell_name_B1_id = ['C' + str(i) for i in range(n_cells[0])]
        # print(len(cell_name_B1_id))
        cell_name_B2_id = ['G' + str(i) for i in range(n_cells[1])]
        # print(len(cell_name_B2_id))
        cell_name_B1B2_id = cell_name_B1_id.copy()
        cell_name_B1B2_id.extend(cell_name_B2_id)
        # print(len(cell_name_B1B2_id))

        scobj.data_ctrl.obs['simple_id'] = cell_name_B1_id
        scobj.data_treat.obs['simple_id'] = cell_name_B2_id

        cell_name_B1 = scobj.data_ctrl.obs_names.tolist()
        cell_name_B2 = scobj.data_treat.obs_names.tolist()
        cell_name_B1B2 = scobj.data_ctrl.obs_names.tolist()
        cell_name_B1B2.extend(scobj.data_treat.obs_names.tolist())
        # print(len(cell_name_B1B2))

        embed_B1_1st_raw = pd.read_csv(resultpath + runLabel + '_Ctrl_1st_embed' + str(size_1ord) + '.txt',
                                       skiprows=1, header=None, index_col=0,
                                       delim_whitespace=True)
        embed_B1_1st = embed_B1_1st_raw.loc[cell_name_B1_id, :]
        embed_B1_1st = embed_B1_1st.set_axis(cell_name_B1, axis=0)

        embed_B2_1st_raw = pd.read_csv(resultpath + runLabel + '_Gemc_1st_embed' + str(size_1ord) + '.txt',
                                       skiprows=1, header=None, index_col=0,
                                       delim_whitespace=True)
        embed_B2_1st = embed_B2_1st_raw.loc[cell_name_B2_id, :]
        embed_B2_1st = embed_B2_1st.set_axis(cell_name_B2, axis=0)

        embed_B1_1st = embed_B1_1st.loc[scobj.data_ctrl.obs_names]
        embed_B2_1st = embed_B2_1st.loc[scobj.data_treat.obs_names]

        # cat and normalize embedding
        embed_B1B2_2nd_raw = pd.read_csv(resultpath + runLabel + '_2nd_embed' + str(size_2ord) + '.txt',
                                         skiprows=1, header=None, index_col=0,
                                         delim_whitespace=True)
        embed_B1B2_2nd = embed_B1B2_2nd_raw.loc[cell_name_B1B2_id,]
        # print(embed_B1B2_2nd.index)
        embed_B1B2_2nd = embed_B1B2_2nd.set_axis(cell_name_B1B2, axis=0)
        # print(embed_B1B2_2nd.index)

        embed_B1_2nd = embed_B1B2_2nd.loc[scobj.data_ctrl.obs_names]
        embed_B2_2nd = embed_B1B2_2nd.loc[scobj.data_treat.obs_names]

        embed_B1 = norm_embedding(embed_1st=embed_B1_1st, embed_2nd_1=embed_B1_2nd,
                                  mode=0, weight=[1, 1])
        embed_B2 = norm_embedding(embed_1st=embed_B2_1st, embed_2nd_1=embed_B2_2nd,
                                  mode=0, weight=[1, 1])

        # Run clustering
        scobj.data_ctrl = runLineClustering(scobj.data_ctrl, embed_B1, key='LINE', resolution=res)
        scobj.data_treat = runLineClustering(scobj.data_treat, embed_B2, key='LINE', resolution=res)
        scobj.data_ctrl.obsm['LINE_umap'] = scobj.data_ctrl.obsm['X_umap']
        scobj.data_treat.obsm['LINE_umap'] = scobj.data_treat.obsm['X_umap']

    def cls2cls(scobj, weight_qc_mat, resultpath,runLabel, res=0.3):
        # obtain the cls2cls matrix
        align_pre = cls2cls_prob(scobj.data_ctrl, scobj.data_treat, weight_qc_mat, thre_cell=0.05, thre_cls=0.05)
        weight_mat_t = np.transpose(weight_qc_mat)

        align_pos = cls2cls_prob(scobj.data_treat, scobj.data_ctrl, weight_mat_t, thre_cell=0.05, thre_cls=0.05)
        align_pos = np.transpose(align_pos)

        tmp = rapid_small(align_pre, align_pos)
        align_pre = tmp[0]
        align_pos = tmp[1]

        np.savetxt(resultpath + runLabel + '_align_pre_probability_res' + str(res) + '.txt', align_pre)
        np.savetxt(resultpath + runLabel + '_align_pos_probability_res' + str(res) + '.txt', align_pos)

        return align_pre, align_pos

    def showUMAP(scobj):
        fig = plt.figure(figsize=(20, 8), dpi=300)
        ax1 = fig.add_subplot(1, 2, 1)
        plotUMAP(scobj.data_ctrl, key='LINE_umap', group_by='LINE_cluster', title='LINE-Ctrl')
        ax2 = fig.add_subplot(1, 2, 2)
        plotUMAP(scobj.data_treat, key='LINE_umap', group_by='LINE_cluster', title='LINE-Treat')
        fig.tight_layout()

    def showHeatmap(align_pre, align_pos):
        fig = plt.figure(figsize=(20, 8), dpi=300)
        ax1 = fig.add_subplot(1, 2, 1)
        sns.heatmap(pd.DataFrame(align_pre), annot=True, linewidths=.5, cmap="YlGnBu")
        plt.xlabel('Post cluster')
        plt.ylabel('Pre cluster')

        ax1 = fig.add_subplot(1, 2, 2)
        sns.heatmap(pd.DataFrame(align_pos), annot=True, linewidths=.5, cmap="YlGnBu")
        plt.xlabel('Post cluster')
        plt.ylabel('Pre cluster')
        fig.tight_layout()

    def rank_genes(scobj):
        scobj.data_ctrl.uns['log1p']["base"] = None
        sc.tl.rank_genes_groups(scobj.data_ctrl, 'LINE_cluster', method='wilcoxon')
        scobj.data_treat.uns['log1p']["base"] = None
        sc.tl.rank_genes_groups(scobj.data_treat, 'LINE_cluster', method='wilcoxon')

        sc.pl.rank_genes_groups_heatmap(scobj.data_ctrl, n_genes=20, use_raw=True, swap_axes=True,
                                        show_gene_labels=True)
        sc.pl.rank_genes_groups_heatmap(scobj.data_treat, n_genes=20, use_raw=True, swap_axes=True,
                                        show_gene_labels=True)

    def showScore(scobj, type):
        plotScore(scobj.data_ctrl, type, 'Ctrl')
        plotScore(scobj.data_treat, type, 'Treat')

    def clusterlevel(scobj,align_pre,align_pos):
        scobj.data_ctrl.obs['LINE_cluster_name'] = ['C' + str(i) for i in scobj.data_ctrl.obs['LINE_cluster']]
        scobj.data_treat.obs['LINE_cluster_name'] = ['T' + str(i) for i in scobj.data_treat.obs['LINE_cluster']]

        scobj.data_comb.obs['LINE_cluster'] = pd.concat([scobj.data_ctrl.obs['LINE_cluster_name'],
                                                         scobj.data_treat.obs['LINE_cluster_name']], axis=0)

        bulk_expr_B1, clus_B1 = getMetacellExpr(scobj.data_ctrl, scobj.data_comb, key='LINE_cluster')
        bulk_expr_B2, clus_B2 = getMetacellExpr(scobj.data_treat, scobj.data_comb, key='LINE_cluster')
        bulk_expr, clus = getMetacellExpr(scobj.data_comb, scobj.data_comb, key='LINE_cluster')

        clus_B1 = ['C' + str(i) for i in clus_B1]
        bulk_expr_B1 = pd.DataFrame(bulk_expr_B1.transpose(), index=scobj.data_comb.var_names, columns=clus_B1)

        clus_B2 = ['T' + str(i) for i in clus_B2]
        bulk_expr_B2 = pd.DataFrame(bulk_expr_B2.transpose(), index=scobj.data_comb.var_names, columns=clus_B2)

        bulk_expr = pd.DataFrame(bulk_expr.transpose(), index=scobj.data_comb.var_names, columns=clus)

        k_value = 3
        alpha_value = 5
        sim = getSimilarity(bulk_expr.transpose(), bulk_expr.transpose(), k=k_value, alpha=alpha_value)
        sim_df = pd.DataFrame(sim, index=bulk_expr.columns, columns=bulk_expr.columns)

        cls2cls_sim = np.array(sim_df.loc[clus_B1, clus_B2])
        cls2cls_sim_B1 = np.array(sim_df.loc[clus_B1, clus_B1])
        cls2cls_sim_B2 = np.array(sim_df.loc[clus_B2, clus_B2])

        G, set1, set2, set_pair, weight_list1, weight_list2, weight_list_pair = getGraph(cls2cls_sim, cls2cls_sim_B1,
                                                                                         cls2cls_sim_B2, align_pre,
                                                                                         align_pos)
        plotGraph(G, align_pos, set1, set2, set_pair, weight_list1, weight_list2, weight_list_pair)

        plt.show()

    def drugresponse(scobj, weight_qc_mat):

        align_pre = cls2cls_prob(scobj.data_ctrl, scobj.data_treat, weight_qc_mat, thre_cell=0.05, thre_cls=0.05)
        weight_mat_t = np.transpose(weight_qc_mat)
        align_pos = cls2cls_prob(scobj.data_treat, scobj.data_ctrl, weight_mat_t, thre_cell=0.05, thre_cls=0.05)
        align_pos = np.transpose(align_pos)
        tmp = rapid_small(align_pre, align_pos)
        align_pre = tmp[0]
        align_pos = tmp[1]

        align_pre_raw = cls2cls_prob(scobj.data_ctrl, scobj.data_treat, weight_qc_mat, thre_cell=0, thre_cls=0)
        weight_mat_t = np.transpose(weight_qc_mat)
        align_pos_raw = cls2cls_prob(scobj.data_treat, scobj.data_ctrl, weight_mat_t, thre_cell=0, thre_cls=0)
        align_pos_raw = np.transpose(align_pos_raw)

        n_clus = [align_pre.shape[0], align_pre.shape[1]]

        num1 = Counter(scobj.data_ctrl.obs['LINE_cluster'])
        n_cells_B1 = [num1[str(i)] for i in range(n_clus[0])]
        num2 = Counter(scobj.data_treat.obs['LINE_cluster'])
        n_cells_B2 = [num2[str(i)] for i in range(n_clus[1])]

        perct_clus_B1 = [n_cells_B1[i] / scobj.data_ctrl.shape[0] for i in range(len(n_cells_B1))]
        perct_clus_B2 = [n_cells_B2[i] / scobj.data_treat.shape[0] for i in range(len(n_cells_B2))]

        align_pre_mat = align_pre.copy()
        align_pos_mat = align_pos.copy()
        for i in range(n_clus[0]):
            align_pre_mat[i, :] = perct_clus_B1[i] * align_pre[i, :]

        for j in range(n_clus[1]):
            align_pos_mat[:, j] = perct_clus_B2[j] * align_pos[:, j]

        align_pre_mat_raw = align_pre_raw.copy()
        align_pos_mat_raw = align_pos_raw.copy()
        for i in range(n_clus[0]):
            align_pre_mat_raw[i, :] = perct_clus_B1[i] * align_pre_raw[i, :]

        for j in range(n_clus[1]):
            align_pos_mat_raw[:, j] = perct_clus_B2[j] * align_pos_raw[:, j]

        pre_clus_perct1 = np.sum(align_pre_mat, axis=1)
        pre_clus_perct2 = np.sum(align_pos_mat, axis=1)
        pos_clus_perct1 = np.sum(align_pre_mat, axis=0)
        pos_clus_perct2 = np.sum(align_pos_mat, axis=0)

        scobj.data_ctrl.obs['LINE_cluster_name'] = ['C' + str(i) for i in scobj.data_ctrl.obs['LINE_cluster']]
        scobj.data_treat.obs['LINE_cluster_name'] = ['T' + str(i) for i in scobj.data_treat.obs['LINE_cluster']]

        scobj.data_comb.obs['LINE_cluster'] = pd.concat([scobj.data_ctrl.obs['LINE_cluster_name'],
                                                         scobj.data_treat.obs['LINE_cluster_name']], axis=0)
        bulk_expr_B1, clus_B1 = getMetacellExpr(scobj.data_ctrl, scobj.data_comb, key='LINE_cluster')
        bulk_expr_B2, clus_B2 = getMetacellExpr(scobj.data_treat, scobj.data_comb, key='LINE_cluster')
        bulk_expr, clus = getMetacellExpr(scobj.data_comb, scobj.data_comb, key='LINE_cluster')
        clus_B1 = ['C' + str(i) for i in clus_B1]
        bulk_expr_B1 = pd.DataFrame(bulk_expr_B1.transpose(), index=scobj.data_comb.var_names, columns=clus_B1)
        clus_B2 = ['T' + str(i) for i in clus_B2]
        bulk_expr_B2 = pd.DataFrame(bulk_expr_B2.transpose(), index=scobj.data_comb.var_names, columns=clus_B2)
        bulk_expr = pd.DataFrame(bulk_expr.transpose(), index=scobj.data_comb.var_names, columns=clus)

        k_value = 3
        alpha_value = 5
        bulk_expr_hrg = bulk_expr.loc[scobj.data_comb.var_names, :]
        sim2 = getSimilarity(bulk_expr_hrg.transpose(), bulk_expr_hrg.transpose(), k=k_value, alpha=alpha_value)
        sim_df2 = pd.DataFrame(sim2, index=bulk_expr.columns, columns=bulk_expr.columns)

        cls2cls_dis = 1 - np.array(sim_df2.loc[clus_B1, clus_B2])
        cls2cls_dis_B1 = 1 - np.array(sim_df2.loc[clus_B1, clus_B1])
        cls2cls_dis_B2 = 1 - np.array(sim_df2.loc[clus_B2, clus_B2])
        G2 = getGraph_notdel_dir(cls2cls_dis, cls2cls_dis_B1, cls2cls_dis_B2, align_pre, align_pos)

        source = []
        target = []
        flow = []
        prob_pre = []
        prob_pos = []
        flow_pre_percent = []
        flow_pos_percent = []
        STC = []
        RAC = []

        for i in range(align_pre.shape[0]):
            for j in range(align_pre.shape[1]):
                if (align_pre[i, j] > 0 or align_pos[i, j] > 0):
                    source.append('C' + str(i))
                    target.append('T' + str(j))
                    flow.append('C' + str(i) + '->' + 'T' + str(j))
                    prob_pre.append(align_pre[i, j])
                    prob_pos.append(align_pos[i, j])
                    flow_pre_percent.append(align_pre_mat[i, j])
                    flow_pos_percent.append(align_pos_mat[i, j])
                    path = nx.dijkstra_path(G2, source='C_' + str(i), target='T_' + str(j))
                    weight = nx.dijkstra_path_length(G2, source='C_' + str(i), target='T_' + str(j))
                    STC.append(weight)
                    RAC.append(align_pos_mat_raw[i, j] / align_pre_mat_raw[i, j])

        flow_info = pd.DataFrame({'Flow': flow,
                                  'Source': source,
                                  'Target': target,
                                  'Probability of Flow in Control Cluster': prob_pre,
                                  'Probability of Flow in Gemcitabine-treated Cluster': prob_pos,
                                  'Percent of Flow in Control Sample': flow_pre_percent,
                                  'Percent of Flow in Gemcitabine-treated Sample': flow_pos_percent,
                                  'STC': STC,
                                  'RAC': RAC
                                  })

        print(flow_info)
        return flow_info

    def sankey(scobj, weight_qc_mat):
        align_pre = cls2cls_prob(scobj.data_ctrl, scobj.data_treat, weight_qc_mat, thre_cell=0.05, thre_cls=0.05)
        weight_mat_t = np.transpose(weight_qc_mat)
        align_pos = cls2cls_prob(scobj.data_treat, scobj.data_ctrl, weight_mat_t, thre_cell=0.05, thre_cls=0.05)
        align_pos = np.transpose(align_pos)
        tmp = rapid_small(align_pre, align_pos)
        align_pre = tmp[0]
        align_pos = tmp[1]
        n_clus = [align_pre.shape[0], align_pre.shape[1]]
        num1 = Counter(scobj.data_ctrl.obs['LINE_cluster'])
        n_cells_B1 = [num1[str(i)] for i in range(n_clus[0])]
        num2 = Counter(scobj.data_treat.obs['LINE_cluster'])
        n_cells_B2 = [num2[str(i)] for i in range(n_clus[1])]
        perct_clus_B1 = [n_cells_B1[i] / scobj.data_ctrl.shape[0] for i in range(len(n_cells_B1))]
        perct_clus_B2 = [n_cells_B2[i] / scobj.data_treat.shape[0] for i in range(len(n_cells_B2))]
        align_pre_mat = align_pre.copy()
        align_pos_mat = align_pos.copy()
        for i in range(n_clus[0]):
            align_pre_mat[i, :] = perct_clus_B1[i] * align_pre[i, :]

        for j in range(n_clus[1]):
            align_pos_mat[:, j] = perct_clus_B2[j] * align_pos[:, j]
        pre_clus_perct1 = np.sum(align_pre_mat, axis=1)
        pre_clus_perct2 = np.sum(align_pos_mat, axis=1)
        pos_clus_perct1 = np.sum(align_pre_mat, axis=0)
        pos_clus_perct2 = np.sum(align_pos_mat, axis=0)

        scobj.data_ctrl.obs['LINE_cluster_name'] = ['C' + str(i) for i in scobj.data_ctrl.obs['LINE_cluster']]
        scobj.data_treat.obs['LINE_cluster_name'] = ['T' + str(i) for i in scobj.data_treat.obs['LINE_cluster']]

        scobj.data_comb.obs['LINE_cluster'] = pd.concat([scobj.data_ctrl.obs['LINE_cluster_name'],
                                                         scobj.data_treat.obs['LINE_cluster_name']], axis=0)
        bulk_expr_B1, clus_B1 = getMetacellExpr(scobj.data_ctrl, scobj.data_comb, key='LINE_cluster')
        bulk_expr_B2, clus_B2 = getMetacellExpr(scobj.data_treat, scobj.data_comb, key='LINE_cluster')
        bulk_expr, clus = getMetacellExpr(scobj.data_comb, scobj.data_comb, key='LINE_cluster')
        clus_B1 = ['C' + str(i) for i in clus_B1]
        bulk_expr_B1 = pd.DataFrame(bulk_expr_B1.transpose(), index=scobj.data_comb.var_names, columns=clus_B1)
        clus_B2 = ['T' + str(i) for i in clus_B2]
        bulk_expr_B2 = pd.DataFrame(bulk_expr_B2.transpose(), index=scobj.data_comb.var_names, columns=clus_B2)
        bulk_expr = pd.DataFrame(bulk_expr.transpose(), index=scobj.data_comb.var_names, columns=clus)

        k_value = 3
        alpha_value = 5
        bulk_expr_hrg = bulk_expr.loc[scobj.data_comb.var_names, :]
        sim2 = getSimilarity(bulk_expr_hrg.transpose(), bulk_expr_hrg.transpose(), k=k_value, alpha=alpha_value)
        sim_df2 = pd.DataFrame(sim2, index=bulk_expr.columns, columns=bulk_expr.columns)

        cls2cls_dis = 1 - np.array(sim_df2.loc[clus_B1, clus_B2])
        cls2cls_dis_B1 = 1 - np.array(sim_df2.loc[clus_B1, clus_B1])
        cls2cls_dis_B2 = 1 - np.array(sim_df2.loc[clus_B2, clus_B2])
        G2 = getGraph_notdel_dir(cls2cls_dis, cls2cls_dis_B1, cls2cls_dis_B2, align_pre, align_pos)

        source = []
        target = []
        flow = []
        flow_pre_percent = []
        flow_pos_percent = []
        STC = []
        pre_sum = []
        pos_sum = []

        for i in range(align_pre.shape[0]):
            for j in range(align_pre.shape[1]):
                if (align_pre[i, j] > 0 or align_pos[i, j] > 0):
                    source.append(i)
                    target.append(j)
                    flow_pre_percent.append(align_pre_mat[i, j])
                    flow_pos_percent.append(align_pos_mat[i, j])
                    weight = nx.dijkstra_path_length(G2, source='C_' + str(i), target='T_' + str(j))
                    STC.append(weight)

        p_df = pd.DataFrame({
            'c': source,
            't': target,
            'c_p': flow_pre_percent,
            't_p': flow_pos_percent,
            'STC': STC
        })
        pre_fractions = p_df[p_df['c'].isin(range(n_clus[0]))].groupby('c')['c_p'].sum().values
        pos_fractions = p_df[p_df['t'].isin(range(n_clus[1]))].groupby('t')['t_p'].sum().values
        fig = plt.figure(figsize=(4, 4))
        gs = GridSpec(p_df.shape[0], 2, width_ratios=[2.1, 1], hspace=0.1, wspace=0.7)
        ax = fig.add_subplot(gs[0:p_df.shape[0], 0])

        dist_cmap = cm.get_cmap('viridis')
        dist_range = max(p_df['STC']) - min(p_df['STC'])
        smap = cm.ScalarMappable(cmap=dist_cmap)
        smap.set_array(p_df['STC'])

        for i in range(n_clus[0]):
            bottom = pre_fractions[(i + 1):].sum()
            rectangle = ax.bar(x=[0], height=pre_fractions[i], bottom=bottom, color=scobj.pre_colors[i],
                               edgecolor='black', fill=True, linewidth=0.7, width=0.16)
            text_y = rectangle[0].get_height() / 2 + bottom
            # ax.text(x=1.16, y=text_y, s='Ctrl'+str(i), color=scobj.pre_colors[i],
            # horizontalalignment='left', verticalalignment='center', fontsize=13)

        for i in range(n_clus[1]):
            bottom = pos_fractions[(i + 1):].sum()
            rectangle = ax.bar(x=[1], height=pos_fractions[i], bottom=bottom, color=scobj.pos_colors[i],
                               edgecolor='black', fill=True, linewidth=0.7, width=0.16)
            text_y = rectangle[0].get_height() / 2 + bottom
            # ax.text(x=1.26, y=text_y, s='Treat'+str(i), color=scobj.pos_colors[i],
            # horizontalalignment='left', verticalalignment='center', fontsize=13)

        legend_labels_pre = [f"Ctrl{i}" for i in range(n_clus[0])]
        legend_handles_pre = [Rectangle((0, 0), 1, 1, color=color, ec='black') for color in scobj.pre_colors]

        legend_labels_pos = [f"Treat{i}" for i in range(n_clus[1])]
        legend_handles_pos = [Rectangle((0, 0), 1, 1, color=color, ec='black') for color in scobj.pos_colors]

        all_legend_labels = legend_labels_pre + legend_labels_pos
        all_legend_handles = legend_handles_pre + legend_handles_pos

        fig.tight_layout()
        ax.legend(all_legend_handles, all_legend_labels, loc='upper left', bbox_to_anchor=(1, 1), ncol=1)

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

        y_start_record = [1 - pre_fractions[0:ii].sum() for ii in range(len(pre_fractions))]
        y_end_record = [1 - pos_fractions[0:ii].sum() for ii in range(len(pos_fractions))]

        y_up_start, y_dw_start = 1, 1
        y_up_end, y_dw_end = 1, 1
        axi = 0
        for si in range(n_clus[0]):
            cur_p_df = p_df.loc[p_df['c'] == si, :]
            if cur_p_df.shape[0] > 0:
                for fi in range(cur_p_df.shape[0]):
                    y_up_start = y_start_record[si]
                    y_dw_start = y_up_start - cur_p_df['c_p'].iloc[fi]
                    y_start_record[si] = y_dw_start

                    ti = cur_p_df['t'].iloc[fi]
                    y_up_end = y_end_record[ti]
                    y_dw_end = y_up_end - cur_p_df['t_p'].iloc[fi]
                    y_end_record[ti] = y_dw_end

                    y_up_start -= 0.003
                    y_dw_start += 0.004
                    y_up_end -= 0.003
                    y_dw_end += 0.004

                    ys_up = y_up_start + (y_up_end - y_up_start) * ys
                    ys_dw = y_dw_start + (y_dw_end - y_dw_start) * ys

                    color_s_t = [scobj.pre_colors[si], scobj.pos_colors[ti]]
                    cmap = LinearSegmentedColormap.from_list('mycmap', [color_s_t[0], color_s_t[1]])
                    grad_colors = cmap(np.linspace(0, 1, len(xs) - 1))
                    grad_colors = [rgb2hex(color) for color in grad_colors]
                    for pi in range(len(xs) - 1):
                        ax.fill_between(xs[pi:(pi + 2)], ys_dw[pi:(pi + 2)], ys_up[pi:(pi + 2)], alpha=0.7,
                                        color=grad_colors[pi], edgecolor=None)


