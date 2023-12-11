
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 0.1 #set the value globally
import matplotlib.pyplot as plt
plt.rc('font', family='Helvetica')

from matplotlib.colors import LinearSegmentedColormap
from matplotlib import gridspec
#from scipy.spatial import Voronoi, ConvexHull,voronoi_plot_2d, Delaunay
from numpy.linalg import norm

from sklearn.datasets import make_classification
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import LogisticRegression,LogisticRegressionCV, Lasso,Ridge, RidgeCV,LassoCV, LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import RandomizedSearchCV,GridSearchCV,cross_val_predict, cross_val_score,RepeatedKFold,RepeatedStratifiedKFold,StratifiedShuffleSplit
#from sklearn.metrics import make_scorer,accuracy_score, f1_score, classification_report,confusion_matrix,roc_curve, roc_auc_score, precision_score, recall_score, precision_recall_curve
from sklearn.metrics import confusion_matrix,r2_score,mean_absolute_error,mean_squared_error,mean_squared_log_error,mean_absolute_percentage_error,median_absolute_error, max_error, explained_variance_score

from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import Pipeline
#from sklearn.metrics import precision_recall_fscore_support as score
#from imblearn.over_sampling import SMOTE, SMOTEN,ADASYN, KMeansSMOTE, SVMSMOTE
from sklearn.utils import class_weight
from sklearn.metrics import roc_curve, auc,consensus_score

#from sklearn.datasets import make_checkerboard
from sklearn.cluster import SpectralBiclustering
from sklearn.decomposition import NMF
from matplotlib.tri import Triangulation
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import SubplotSpec

from sklearn.metrics.pairwise import cosine_similarity
#bicluster

from gseapy.plot import gseaplot, heatmap
import gseapy
from sklearn.decomposition import PCA as skPCA

from scipy.spatial import cKDTree
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity

#Metrics
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import log_loss
from sklearn.metrics import zero_one_loss
from sklearn.metrics import matthews_corrcoef
from scipy.stats import pearsonr,spearmanr

import pandas as pd
import numpy as np
import seaborn as snn
import os
import random
import warnings
import time
import scanpy as sc
import pickle
import xlsxwriter
from types import SimpleNamespace
import math
import scipy
from sklearn.utils.extmath import svd_flip
import statsmodels.api as sm
#import shap
from scipy.stats import entropy
import sys
from sklearn.model_selection import KFold


fpath=os.path.join(os.path.dirname(__file__),'utils')
sys.path.append(fpath)
#sys.path.insert(1,'./utilities/')
#sys.path.insert(1,'./ionmf/factorization/')
#from ionmf.factorization.onmf import onmf
#from ionmf.factorization.model import iONMF
from pyliger_utilities import nnlsm_blockpivot,iNMF,NMF_obj_eval


def gene_covariation_analysis(Radius=0,niche_prediction_dir='./spatial_ct_ct_interactions/',
refpath='./inputRef/',quepath='./inputQuery/',ref_cluster_tag='cluster',
ref_common_counts='common_counts_sc.h5ad',
ref_original_counts='Original_counts.h5ad',
que_common_counts='common_counts_sp.h5ad',
LRdbFilename='NiCoLRdb.txt',iNMFmode=True,no_of_factors=3,
shap_analysis=False,shap_cluster_cutoff=0.5,
cutoff_to_count_exp_cell_population=0,
lambda_c=list(np.power(2.0, np.arange(-10, 10))),
#lambda_c=list(10 * 0.90 ** np.arange(1,100)),
#lambda_c=[1],
coeff_cutoff_for_rid_reg=0,logistic_coef_cutoff=0):

    #outputdir='./spatial_ct_ct_interactions/',
    ####Random seed used in RepeatedStratifiedKFold
    #####seed=3685134seed=3685134,

    """
    | **This is the primary function called by the user to perform gene covariation analysis in the niche**.
    | **Before calling this function user must call spatial_neighborhood_analysis function from the interaction module.**
    | Please provide Original_counts.h5ad, common_counts_sc.h5ad, sct_singleCell.h5ad files from scRNAseq data.
    | And common_counts_sp.h5ad, sct_spatial.h5ad files for the single-cell resolution of spatial data.
    | Original_counts.h5ad should also have the cluster information in the .obs[] and the umap information in .obsm['X_umap']

    Inputs:

    | The previous niche interaction runs output directory location from the function spatial_neighborhood_analysis
    | (default) niche_prediction_dir='./spatial_ct_ct_interactions/'

    | Reference scRNAseq count matrix in scTransform-like normalization in the common gene space. The filename must be sct_singleCell.h5ad
    | (default) refpath='./inputRef/'

    | Queried spatial count matrix in scTransform-like normalization in the common gene space. The filename must be sct_spatial.h5ad
    | (default) quepath='./inputQuery/'

    | Original count data of scRNAseq in h5ad format
    | (default) ref_original_counts='Original_counts.h5ad'
    | Must have the cluster information in the .obs and the umap information in .obsm['X_umap']
    | .raw layer should have integer count matrix data. It will used to find the Spearman correlation and cosine similarity.  

    | The tag in reference h5ad file where cluster information is stored
    | (default) ref_cluster_tag='cluster'

    | Common genes scRNAseq count data in h5ad format
    | (default) ref_common_counts='common_counts_sc.h5ad'

    | Common gene spatial count data in h5ad format
    | (default) que_common_counts='common_counts_sp.h5ad'

    | The filename of the ligand-receptor database (the first column is Ligand, the second column is Receptor, a third is the resource list)
    | (default) LRdbFilename='NiCoLRdb.txt'

    | This radius parameter should be the same as used in spatial neighborhood analysis to find the niche interactions.
    | (default) Radius=0

    | Number of factors used in NMF for finding the common gene latent dimension space
    | (default) no_of_factors=3

    | By default, the common gene latent dimension space uses an integrated NMF approach to learn a gene-by-factor submatrix from both modalities.
    | If this value is false, then it uses an ordinary NMF approach to learn a gene by factor submatrix only from scRNAseq data.
    | (default) iNMFmode=True

    | The initial range of regularization parameters used in the ridge regression step to find the best-regularized parameter
    | (default) lambda_c=list(np.power(2.0, np.arange(-10, 10)))

    | Do you want to perform the shap analysis
    | shap_analysis=False

    | Shap analysis cutoff parameter
    | shap_cluster_cutoff=0.5

    | This cutoff is used to make the list of significant celltype_factor-celltype_factor niche interactions whose absolute regression coefficient is greater than this
    | coeff_cutoff_for_rid_reg=0

    | This cutoff is used to know the positive niche interactions (cell type -cell type). When it is >0, then cell type - cell type likely to interact with each other
    | logistic_coef_cutoff=0

    | Used to find the percentage of the cell population that expressed a given gene in a given cell type. Value 0 is acceptable with count data.
    | cutoff_to_count_exp_cell_population=0

    Outputs:

    | The output of covariations in the interacted niche
    | (default) niche_prediction_dir=./spatial_ct_ct_interactions/covariations_R*_F*
    """

    #ref_h5ad=refpath+'sct_singleCell.h5ad'
    #que_h5ad=quepath+'sct_spatial.h5ad'

    que_common_counts=quepath+que_common_counts
    ref_common_counts=refpath+ref_common_counts
    ref_original_counts=refpath+ref_original_counts

    original_h5ad=sc.read_h5ad(ref_original_counts)
    df=original_h5ad.obs[ref_cluster_tag]
    annotation_singlecell_celltypename=df.to_numpy()
    cellname=df.index.to_numpy()


    sc_ct_name=[]
    A=list(sorted(np.unique(annotation_singlecell_celltypename)))
    d={}
    for i in range(len(A)):
        sc_ct_name.append([i,A[i]])
        d[A[i]]=i
    sc_ct_name=np.array(sc_ct_name)

    sc_cluster=[]
    for j in range(len(annotation_singlecell_celltypename)):
        sc_cluster.append([cellname[j],d[annotation_singlecell_celltypename[j]]])
    sc_cluster=np.array(sc_cluster)
    annotation_singlecell_barcode_id=sc_cluster[:,0]
    annotation_singlecell_cluster_id=sc_cluster[:,1]
    singlecell_unique_clustername=sc_ct_name[:,1]
    singlecell_unique_clusterid=sc_ct_name[:,0]


    f=open(LRdbFilename,'r')
    LRdb=f.readlines()
    f.close()

    sct_ad_sp=sc.read_h5ad(que_common_counts)
    sct_ad_sc=sc.read_h5ad(ref_common_counts)
    full_ad_sc=original_h5ad.raw
    covariation_outdir=niche_prediction_dir+'covariations_'

    strategy='niche_prediction_linear/'
    gene_set_names=[]

    #print('sc1 annotation_singlecell_cluster_id',len(annotation_singlecell_cluster_id))
    #print('sc2 annotation_singlecell_barcode_id',len(annotation_singlecell_barcode_id))
    #print('sc3 annotation_singlecell_celltypename',len(annotation_singlecell_celltypename))
    #print('sc4 singlecell_unique_clustername', len(singlecell_unique_clustername))


    # load spatial dat
    sp_genename=sct_ad_sp.var_names.to_numpy()
    sc_genename=sct_ad_sc.var_names.to_numpy()
    index_sp,index_sc=find_index(sp_genename,sc_genename)
    #print('common genes between sc and sp',len(index_sp),len(index_sc))
    ad_sp_ori=sct_ad_sp[:,index_sp].copy()
    ad_sc_ori=sct_ad_sc[:,index_sc].copy()


    inputRadius=[Radius]
    for radius in inputRadius:
        celltypeFilename=niche_prediction_dir+'used_CT.txt'
        clusterFilename=niche_prediction_dir+'used_Clusters'+str(radius)+'.csv'

        annotation_spatial_celltypename,annotation_spatial_barcode_id,annotation_spatial_cluster_id,spatialcell_unique_clustername,spatialcell_unique_clusterid=read_spatial_data(clusterFilename,celltypeFilename)


        neighbors=pickle.load( open(niche_prediction_dir+'neighbors_'+str(radius)+'.p', "rb" ) )
        distances=pickle.load( open(niche_prediction_dir+'distances_'+str(radius)+'.p', "rb" ) )

        covariation_dir=covariation_outdir+'R'+str(radius)+'_F'+str(no_of_factors)+'/'
        create_directory(covariation_dir)
        outputname=covariation_dir+'Principal_component_feature_matrix.npz'
        inputdata={}
        inputdata['no_of_pc']=no_of_factors
        inputdata['outputname']=outputname
        inputdata['covariation_dir']=covariation_dir


        fname=niche_prediction_dir+strategy+'/classifier_matrices_'+str(radius)+'.npz'
        data=np.load(fname,allow_pickle=True)
        logistic_coef=data['coef']
        logistic_cmn=data['cmn']
        logistic_cmn_std=data['cmn_std']
        logistic_coef_std=data['coef_std']
        logistic_CTFeatures=data['CTFeatures']
        #f=open(input_spatial+'BiologicalNameOfCT.dat')
        f=open(celltypeFilename)
        nameOfCellType={}
        for line in f:
            l=line[0:-1].split('\t')
            nameOfCellType[int(l[0])]=l[1]

        logistic_predicted_interactions=find_logistic_regression_interacting_score(logistic_cmn,logistic_coef,logistic_CTFeatures,nameOfCellType,logistic_coef_cutoff)



        inputdata['ad_sp']=ad_sp_ori #sct_ad_sp
        inputdata['ad_sc']=ad_sc_ori#sct_ad_sc#
        inputdata['annotation_spatial_cluster_id']=annotation_spatial_cluster_id
        inputdata['annotation_spatial_barcode_id']=annotation_spatial_barcode_id
        inputdata['annotation_spatial_celltypename']=annotation_spatial_celltypename
        inputdata['spatialcell_unique_clustername']=spatialcell_unique_clustername
        inputdata['spatialcell_unique_clusterid']=spatialcell_unique_clusterid

        inputdata['annotation_singlecell_cluster_id']=annotation_singlecell_cluster_id
        inputdata['annotation_singlecell_barcode_id']=annotation_singlecell_barcode_id
        inputdata['annotation_singlecell_celltypename']=annotation_singlecell_celltypename
        inputdata['singlecell_unique_clustername']=singlecell_unique_clustername
        inputdata['singlecell_unique_clusterid']=singlecell_unique_clusterid
        inputdata['neighbors']=neighbors
        inputdata['neigh_distances']=distances
        inputdata['nmf_output']=covariation_dir+'NMF_output/'

        regression_outdir=covariation_dir+'Regression_outputs'+'/'
        create_directory(regression_outdir)
        #inputdata['seed']=seed
        inputdata['lambda_c']=lambda_c
        inputdata['iNMFmode']=iNMFmode
        inputdata['regression_outdir']=regression_outdir
        #inputdata['K_fold']=K_fold
        #inputdata['n_repeats']=n_repeats
        #inputdata['n_jobs']=n_jobs
        inputdata['shap_analysis']=shap_analysis
        inputdata['shap_cluster_cutoff']=shap_cluster_cutoff

        inputdata['logistic_coef_cutoff']=logistic_coef_cutoff
        inputdata['coeff_cutoff_for_rid_reg']=coeff_cutoff_for_rid_reg
        inputdata['gene_set_names']=gene_set_names
        #inputdata['pvalueCutoff']=pvalueCutoff

        inputdata['cutoff_to_count_exp_cell_population']=cutoff_to_count_exp_cell_population
        inputdata['LRdb']=LRdb

        input=SimpleNamespace(**inputdata)

        flag=1
        if os.path.isfile(outputname):
            filesize = os.path.getsize(outputname)
            if filesize>0: #If file is already exist and have size greater than 0 then no need to run again. It will save some time if you want to run it again with different parameters
                flag=0

        if flag==1:
            pc_of_sp_clusterid,PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=compute_PC_space(input,full_ad_sc)
            # full_ad_sc use in only find_PC_of_invidualCluster_in_SC function
            # ideally it should be sctransform way of normalized matrix equivalent to sct_ad_sc but
            # if not then need to do perform scaling HVG etc
            pickle.dump((PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors),open(covariation_dir+'factors_info.p', 'wb'))
            inputdata['pc_of_sp_clusterid']=pc_of_sp_clusterid
            input=SimpleNamespace(**inputdata)
            makePCneighboorhoodFeatureMatrix(input)


        save_reg_coef=model_linear_regression(input,logistic_predicted_interactions)
        inputdata['save_reg_coef']=save_reg_coef
        input=SimpleNamespace(**inputdata)
    return input


def plot_cosine_and_spearman_correlation_to_factors(input,choose_celltypes=[],NOG_Fa=30,saveas='pdf',transparent_mode=False,showit=True,figsize=(15,10)):
    """
    Inputs:

    The main input is the output from gene_covariation_analysis.

    | The cell type that you want to see the covariation regression pattern.
    | (default) choose_celltypes=[]
    | If the list is empty, then the output will show for all the cell types.

    | Number of genes to visualize in each factor
    | (default) NOG_Fa=30

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Dimension of the figure size
    | (default) figsize=(15,10)

    | Background color in the figures
    | (default) transparent_mode=False

    Outputs:

    | The output NMF plots are saved in ./spatial_ct_ct_interactions/covariations_R*_F*/NMF_output

    """

    create_directory(input.nmf_output)

    n=len(input.spatialcell_unique_clusterid)

    perform=[]
    found=[]
    for fi in range(n):
          CC_celltype_name=input.spatialcell_unique_clustername[fi]
          if len(choose_celltypes)==0:
              perform.append(fi)
          else:
              if CC_celltype_name in choose_celltypes:
                  perform.append(fi)
                  found.append(CC_celltype_name)
    if len(choose_celltypes)!=0:
          print("cell types found ",found)

    print("These figures are saved in following path ", input.nmf_output)


    xlabels=[]
    for i in range(input.no_of_pc):
        xlabels.append('NMF'+str(i+1))

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)

    for fi in perform:
        clid=input.spatialcell_unique_clusterid[fi]
        spearman_factors,CC_PCA,CC_gene,CC_meanExpression,CC_popExpression,cosine_factors,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        CC_celltype_name=input.spatialcell_unique_clustername[fi]

        genename_full=CC_gene
        #sc_cosine=find_correlation_bw_genes_and_PC_component_in_singlecell_cosine(H.T,CbyG)
        gname2b,geneNMF2b=top_genes_in_correlation_list_without(genename_full,cosine_factors,NOG_Fa)
        #sc_spearman=find_correlation_bw_genes_and_PC_component_in_singlecell(H.T,CbyG)
        gname3b,geneNMF3b=top_genes_in_correlation_list_without(genename_full,spearman_factors,NOG_Fa)

        selectedGenesAvgExp_cosine=np.zeros( (len(gname2b),1) )
        for i in range(len(gname2b)):
            ind=np.where(genename_full==gname2b[i])
            selectedGenesAvgExp_cosine[i,0]=np.log10(CC_meanExpression[ind[0]])

        selectedGenesAvgExp=np.zeros( (len(gname3b),1) )
        for i in range(len(gname3b)):
            ind=np.where(genename_full==gname3b[i])
            selectedGenesAvgExp[i,0]=np.log10(CC_meanExpression[ind[0]])


        fig=plt.figure(figsize=figsize)
        gs = fig.add_gridspec(ncols=4, nrows=1, wspace=0.5,width_ratios=[0.5, 2,2,0.5])
        ax0=fig.add_subplot(gs[0])
        ax1=fig.add_subplot(gs[1])
        ax2=fig.add_subplot(gs[2])
        ax3=fig.add_subplot(gs[3])
        b=snn.heatmap(selectedGenesAvgExp_cosine,yticklabels=gname2b,ax=ax0)#componentlabel,ax=ax
        b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)
        b.set_title('log(avg exp)')


        b=snn.heatmap(geneNMF2b,yticklabels=gname2b,ax=ax1)#componentlabel,ax=ax
        b.set_xticklabels(xlabels,size = 8,rotation=90)
        b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)
        b.set_title('cosine')
        #b.set_title('spatial'+entropy_SH)
        b=snn.heatmap(geneNMF3b,yticklabels=gname3b,ax=ax2)#componentlabel,ax=ax
        b.set_xticklabels(xlabels,size = 8,rotation=90)
        b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)
        b.set_title('spearman corr '+CC_celltype_name+', lambda = '+ str(alpha))

        b=snn.heatmap(selectedGenesAvgExp,yticklabels=gname3b,ax=ax3)#componentlabel,ax=ax
        #b.set_xticklabels('exp',size = 8,rotation=90)
        b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)
        b.set_title('log(avg exp)')
        #plt.tight_layout()
        print("The figures are saved: ", input.nmf_output+remove_extra_character_from_name(CC_celltype_name)+'.'+saveas)
        fig.savefig(input.nmf_output+remove_extra_character_from_name(CC_celltype_name)+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
        if showit:
            pass
        else:
            plt.close('all')

def plot_feature_matrices(input,showit=True,saveas='pdf',transparent_mode=False,figsize=(10,10)):
    """
    Plots features vectors of the spatial factors from the regression step.
    """
    #maindir1=outputdir+'covariations_'
    #maindir=maindir1+str(radius)+'/'
    #outputname=maindir+'Principal_component_feature_matrix'+str(no_of_factors)+'.npz'
    ylabelname=[]
    for i in range(len(input.spatialcell_unique_clustername)):
        for j in range(input.no_of_pc):
            ylabelname.append(input.spatialcell_unique_clustername[i]+'_'+'Fa'+str(j+1))

    data1=np.load(input.outputname,allow_pickle=True)
    data=data1['weighted_neighborhood_of_factors_in_niche']
    fig,axs=plt.subplots(1,1,figsize=figsize)
    #data=np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
    Feature=data[:,(1+input.no_of_pc):data.shape[1]]
    index=np.argsort(input.annotation_spatial_cluster_id)
    snn.heatmap(np.log10(Feature[index,:]),xticklabels=ylabelname)
    #fig.tight_layout()
    print("The figures are saved: ", input.covariation_dir+'Feature_matrix_PC'+'.'+saveas)
    fig.savefig(input.covariation_dir+'Feature_matrix_PC'+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
    if showit:
        pass
    else:
        plt.close('all')

def plot_significant_regression_covariations_as_circleplot(input,choose_celltypes=[],saveas='pdf',pvalue_cutoff=0.05,mention_pvalue=True,
transparent_mode=False,showit=True,figsize=(6,1.25)):

    """
    Inputs:

    | The main input is the output from gene_covariation_analysis.

    | The cell type that you want to see the covariation regression pattern.
    | (default) choose_celltypes=[]
    | If the list is empty, then the output will show for all the cell types.

    | The pvalue cutoff used to print the -log10(pvalue) on the top of circle.
    | (default) pvalueCutoff=0.05

    | Whether to mention the p-value on top of the circle plot. If false, then it would not show.
    | (default) mention_pvalue=True

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Background color in the figures
    | (default) transparent_mode=False

    | Dimension of the figure size
    | (default) figsize=(6,1.25)

    Outputs:

    | The regression figures are saved in ./spatial_ct_ct_interactions/covariations_R*_F*/Regression_outputs/

    """

    n=len(input.spatialcell_unique_clusterid)

    perform=[]
    found=[]
    for fi in range(n):
          CC_celltype_name=input.spatialcell_unique_clustername[fi]
          if len(choose_celltypes)==0:
              perform.append(fi)
          else:
              if CC_celltype_name in choose_celltypes:
                  perform.append(fi)
                  found.append(CC_celltype_name)
    if len(choose_celltypes)!=0:
          print("cell types found ",found)

    print("The regression figures as pvalue circle plots are saved in following path ", input.regression_outdir+'pvalue_coeff_circleplot_*')


    for i in perform:
        filename=input.spatialcell_unique_clustername[i]
        temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        index=temp[0]

        data=input.save_reg_coef[input.spatialcell_unique_clusterid[i]]
        coef_mu,intercept,alpha,xlabel,score,target,neighborhoodClass,pvalue,pve,rve=data

        '''
        savedata=input.regression_outdir+'coef'+str(input.spatialcell_unique_clusterid[i])+'.npz'
        data=np.load(savedata,allow_pickle=True)
        coef_mu=data['coef_mu']
        intercept=data['intercept']
        pve=data['pve'] # percentage variance explanined
        rve=data['rve'] # residual variance explained
        pvalue=data['pvalue']
        #coef_std=data['coef_std']
        #comp_score=data['comp_score']
        #comp_score_std=data['comp_score_std']
        alpha=data['alpha']
        xlabel=data['xlabel']
        score=data['score']
        '''

        componentlabel=[]
        for j in range(input.no_of_pc):
            componentlabel.append('Fa'+str(j+1))

        percentVE=''
        percentRE=''

        for j in range(len(pve)):
            if j!=0:
                percentVE+=', '
                percentRE+=', '
            percentVE+='%0.3f'%pve[j]
            #percentRE+='%0.1f'%rve[j]

        ylabelname=[]
        for k in range(len(xlabel)):
            for j in range(input.no_of_pc):
                ylabelname.append(xlabel[k]+'_s'+'%0.3f'%score[k]+'_Fa'+str(j+1))


        #tempG=pvalue<0.1
        #m1,m2=tempG.nonzero()
        #coef_mu[m1,m2]=0
        #coef_mu=tempG.astype(int)

        #pvalue=pvalue<0.05
        pvalue[pvalue<10**-10]=10**-10
        pvalue=-np.log10(pvalue)
        pvalue=np.nan_to_num(pvalue)
        pvcut=-np.log10(pvalue_cutoff)

        factor=0.5
        newfigsize=(factor*len(ylabelname),figsize[1])

        fig, ax = plt.subplots(1,1,figsize=newfigsize)
        M=pvalue.shape[1]
        N=pvalue.shape[0]
        c=coef_mu
        x, y = np.meshgrid(np.arange(M), np.arange(N))
        R = pvalue/10.0/2
        maxp=pvalue.max()
        circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
        col = PatchCollection(circles, array=c.flatten(), cmap='jet')#cmap="RdYlGn")
        ax.add_collection(col)
        ax.set(xticks=np.arange(M), yticks=np.arange(N),
               xticklabels=ylabelname, yticklabels=componentlabel)
        ax.set_xticks(np.arange(M+1)-0.5, minor=True)
        ax.set_yticks(np.arange(N+1)-0.5, minor=True)
        ax.set_xticklabels(ylabelname,size = 8,rotation=90)
        ax.set_title(filename+r',$\alpha$='+str(alpha)+',EVS='+percentVE,fontsize=6)
        ax.grid(which='minor')

        if mention_pvalue:
            for e in range(M):
                for f in range(N):
                    if pvalue[f,e]>pvcut:
                        ax.text(e,f,'%0.2f'%pvalue[f,e],fontsize=4)
        fig.colorbar(col)
        #fig.tight_layout()
        savefname=remove_extra_character_from_name(str(input.spatialcell_unique_clusterid[i])+'_'+filename)
        #print('\n\n\n',input.regression_outdir+'pvalue_coeff_circleplot_'+savefname+'.'+saveas)
        #print("The figures are saved: ", input.regression_outdir+'pvalue_significance_coeff_matrix_'+savefname+'.'+saveas)
        fig.savefig(input.regression_outdir+'pvalue_coeff_circleplot_'+savefname+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
        if showit:
            pass
        else:
            plt.close('all')


def plot_significant_regression_covariations_as_heatmap(input,choose_celltypes=[],saveas='pdf',transparent_mode=False,showit=True,figsize=(6,10)):

    """
    Inputs:

    | The main input is the output from gene_covariation_analysis.

    | The cell type that you want to see the covariation regression pattern.
    | (default) choose_celltypes=[]
    | If the list is empty, then the output will show for all the cell types.

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Background color in the figures
    | (default) transparent_mode=False

    | Dimension of the figure size
    | (default) figsize=(6,1.25)

    Outputs:

    | The regression figures are saved in ./spatial_ct_ct_interactions/covariations_R*_F*/Regression_outputs/

    """

    n=len(input.spatialcell_unique_clusterid)

    perform=[]
    found=[]
    for fi in range(n):
          CC_celltype_name=input.spatialcell_unique_clustername[fi]
          if len(choose_celltypes)==0:
              perform.append(fi)
          else:
              if CC_celltype_name in choose_celltypes:
                  perform.append(fi)
                  found.append(CC_celltype_name)
    if len(choose_celltypes)!=0:
          print("cell types found ",found)

    print("The regression figures as pvalue heatmap plots are saved in following path ", input.regression_outdir+'pvalue_coeff_heatmap_*')


    for i in perform:
        filename=input.spatialcell_unique_clustername[i]
        temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        index=temp[0]

        data=input.save_reg_coef[input.spatialcell_unique_clusterid[i]]
        coef_mu,intercept,alpha,xlabel,score,target,neighborhoodClass,pvalue,pve,rve=data

        componentlabel=[]
        for j in range(input.no_of_pc):
            componentlabel.append('Fa'+str(j+1))

        percentVE=''
        percentRE=''

        for j in range(len(pve)):
            if j!=0:
                percentVE+=', '
                percentRE+=', '
            percentVE+='%0.3f'%pve[j]
            #percentRE+='%0.1f'%rve[j]

        xlabelname=[]
        for k in range(len(xlabel)):
            for j in range(input.no_of_pc):
                xlabelname.append(xlabel[k]+'_s'+'%0.3f'%score[k]+'_Fa'+str(j+1))

        pvalue[pvalue<10**-10]=10**-10
        pvalue=-np.log10(pvalue)
        pvalue=np.nan_to_num(pvalue)

        factor=0.5
        newfigsize=(factor*len(xlabelname),figsize[1])

        fig=plt.figure(figsize=newfigsize)
        gs = fig.add_gridspec(ncols=1, nrows=2, height_ratios=[3, 1])
        ax0=fig.add_subplot(gs[0])
        ax1=fig.add_subplot(gs[1])


        a=snn.heatmap(coef_mu,xticklabels=xlabelname,yticklabels=componentlabel,ax=ax0,cbar_kws={"shrink": 1})
        #snn.axes_style(xtick.top=True)
        xlabels= a.get_xticklabels()
        a.set_xticklabels([])
        #a.set_xticklabels(xlabels,size = 8,rotation=90)
        a.set_yticklabels(componentlabel,rotation = 0,size=6)
        #a.set_ylabel('Principal components')
        a.set_title(filename+r',$\alpha$='+str(alpha)+',EVS='+percentVE,fontsize=6)


        b=snn.heatmap(pvalue,annot=True, fmt='.2f',cmap=snn.cm.rocket_r,annot_kws={"size": 3},xticklabels=xlabelname,yticklabels=componentlabel,ax=ax1,cbar_kws={"shrink": 1})
        #b.xaxis.tick_top()
        xlabels= b.get_xticks()
        xlabels= b.get_xticklabels()

        #b.set_xticklabels([])#xlabels,size = 0)
        b.set_xticklabels(xlabels,size = 6,rotation=90)
        b.set_yticklabels(componentlabel,rotation = 0,size=6)




        #fig.tight_layout()
        savefname=remove_extra_character_from_name(str(input.spatialcell_unique_clusterid[i])+'_'+filename)
        #print("The figures are saved: ", input.regression_outdir+'pvalue_significance_coeff_matrix_'+savefname+'.'+saveas)
        fig.savefig(input.regression_outdir+'pvalue_coeff_heatmap_'+savefname+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
        if showit:
            pass
        else:
            plt.close('all')


def save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types(input,pvalueCutoff=0.05, correlation_with_spearman=True,
LR_plot_NMF_Fa_thres=0.2, LR_plot_Exp_thres=0.2,number_of_top_genes_to_print=20):

    """
    Inputs:

    | The main input is the output from gene_covariation_analysis.

    | The cutoff used to find the significant central cell type factor and niche cell type factor interactions
    | (default) pvalueCutoff=0.05

    | The number of top correlating genes to print in regression summary text file.
    | (default) number_of_top_genes_to_print=20

    | If True, visualize genes factor correlation with Spearman; otherwise, If False, then cosine.
    | (default) correlation_with_spearman=True

    | Ligand Receptor analysis uses factor cutoff to find the enriched LR pairs
    | (default) LR_plot_NMF_Fa_thres=0.2

    | Ligand Receptor analysis uses % expressed cell population cutoff to find the enriched LR pairs.
    | (default) LR_plot_Exp_thres=0.2

    Outputs:

    | The Ligand receptors plots are also saved in an Excel sheet for easy access to this information.
    | In the sheet, columns are structured as follows:
    | column 1st denotes the ID of the cell type-cell type interaction,
    | 2nd and 3rd columns identify the interacting cell types A and B,
    | 4th column contains the normalized colocalized scores from the classifier,
    | 5th and 6th represent the NMF factor IDs (metagenes) in cell types A and B,
    | 7th column displays the ridge regression coefficient of the factors’ interaction,
    | 8th column contains the ligand in cell type A,
    | 9th column holds the receptor in cell type B,
    | 10th and 11th exhibit the Pearson correlation of ligand and receptor genes, respectively, in cell types A and B, corresponding to the factors.
    | 12th and 13th showcase the average expression of ligands and receptors in cell types A and B
    | 14th and 15th reveal the population of cells expressing these genes with counts greater than zero in cell types A and B.

    | The regression summary text mentioned the central cell type factor as (CC-Fa(i)), niche cell type factor as (NC-Fa(j))
    | The first-row information is as follows: CC-Fa(i), CC (cell type),niche_score (from classifier), NC-Fa*, NC (cell type), RegCoeff (covariation score), p-value in normal scale, p-value in -log10 scale
    | The second row inforamtion contains: Top 20 (number_of_top_genes_to_print) genes from Fa(i) of central cell type. Genes and their factor ID are mentioned in the pair.
    | The third row inforamtion contains: Top 20 (number_of_top_genes_to_print) genes from Fa(j) of niche cell type. Genes and their factor ID are mentioned in the pair.

    | Our analysis accounts for bidirectional cellular crosstalk interactions of ligands and receptors in cell types A and B.
    | It means, in one way, it assumes a ligand diffused from cell type A, and the receptor is detected in cell type B,
    | while in another way, it assumes a ligand is diffusing from cell type B, and the receptor is detected in cell type A.
    | Both ligand-receptor plot and sheet show the bidirectional cellular crosstalk interaction of ligand and receptor in cell types A and B.
    | Each central cell type is represented in a separate Excel sheet, while the LR enrichment sheet aggregates all interactions across central cell types.

    """

    totalLRpairs,ligand,receptor,either=read_LigRecDb(input.LRdb)
    coeff_cutoff_for_log_reg=input.logistic_coef_cutoff
    coeff_cutoff_for_rid_reg=input.coeff_cutoff_for_rid_reg
    gene_set_names=input.gene_set_names
    LRcutoff=LR_plot_NMF_Fa_thres #Used in excel sheet to show the enrichment of ligand receptor intera

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)

    workbook = xlsxwriter.Workbook(input.covariation_dir+'Lig_and_Rec_enrichment_in_interacting_celltypes.xlsx')
    fout=open(input.covariation_dir+'Regression_summary.txt','w')
    worksheet = workbook.add_worksheet('LR enrichment')
    worksheetrow=0
    main_header=['Id','A','B','localized score','Fa(A)','Fa(B)', 'Coeff' ,'Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)','AvgExp(A)','AvgExp(B)','PopExp(A)','PopExp(B)']
    for ri in range(len(main_header)):
        worksheet.write(worksheetrow,ri,main_header[ri])
    worksheetrow+=1

    d={}
    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        clname=input.spatialcell_unique_clustername[i]
        d[clname]=clid

    print("The Excel sheet is saved: ",input.covariation_dir+'Lig_and_Rec_enrichment_in_interacting_celltypes.xlsx')
    print("The text file is saved:",input.covariation_dir+'Regression_summary.txt')

    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        CC_corr_spearman,CC_PCA,CC_gene,CC_meanExpression,CC_popExpression,CC_corr_cosine,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        CC_celltype_name=input.spatialcell_unique_clustername[i]
        #temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        #index=temp[0]

        data=input.save_reg_coef[input.spatialcell_unique_clusterid[i]]
        coef_mu,intercept,alpha,xlabel,score,target,neighborhoodClass,pvalue,pve,rve=data
        NC_celltype_name=xlabel
        largest=np.max(abs(coef_mu))
        normalized_ridge_coef=coef_mu/largest

        ylabelname=[]
        componentlabel=[]
        for j in range(input.no_of_pc):
            ylabelname.append('CC_'+CC_celltype_name+'_Fa'+str(j+1))
            componentlabel.append('Fa'+str(j+1))

        for k in range(len(NC_celltype_name)):
            if score[k]>coeff_cutoff_for_log_reg:
                if CC_celltype_name!=NC_celltype_name[k]:
                    for j in range(input.no_of_pc):
                        ylabelname.append('NC_'+NC_celltype_name[k]+'_s'+'%0.3f'%score[k]+'_Fa'+str(j+1))

        pc_index_nc=[]
        for k in range(len(NC_celltype_name)):
            for j in range(input.no_of_pc):
                pc_index_nc.append(j)

        CC_celltype_sheetname=remove_extra_character_from_name(CC_celltype_name)
        worksheet_local = workbook.add_worksheet(CC_celltype_sheetname)
        worksheetrow_local=0
        for ri in range(len(main_header)):
            worksheet_local.write(worksheetrow_local,ri,main_header[ri])
        worksheetrow_local+=1

        interaction_id=0
        for k in range(normalized_ridge_coef.shape[0]):
            #k is PC of central cell type
            for j in range(normalized_ridge_coef.shape[1]):
                interaction_id+=1
                index=math.floor(j/input.no_of_pc)
                if (pvalue[k,j]<pvalueCutoff)&(abs(normalized_ridge_coef[k,j])>coeff_cutoff_for_rid_reg):
                #if True:
                    if score[index]>coeff_cutoff_for_log_reg:
                        NC_corr_spearman,NC_PCA,NC_gene,NC_meanExpression,NC_popExpression,NC_corr_cosine,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[d[NC_celltype_name[index]]]
                        if correlation_with_spearman:
                            top_genes_in_CC,top_genes_in_NC,genesWithUP,genesWithDown,Found1,Found2=find_fold_change(CC_corr_spearman,NC_corr_spearman,CC_gene,k,pc_index_nc[j],totalLRpairs,LRcutoff,CC_meanExpression,NC_meanExpression,CC_popExpression,NC_popExpression,number_of_top_genes_to_print)
                        else:
                            top_genes_in_CC,top_genes_in_NC,genesWithUP,genesWithDown,Found1,Found2=find_fold_change(CC_corr_cosine,NC_corr_cosine,CC_gene,k,pc_index_nc[j],totalLRpairs,LRcutoff,CC_meanExpression,NC_meanExpression,CC_popExpression,NC_popExpression,number_of_top_genes_to_print)

                        common_genes=list(set(top_genes_in_CC).intersection(set(top_genes_in_NC)))


                        if CC_celltype_name!=NC_celltype_name[index]:
                                for ele in range(len(Found1)):
                                    header=[str(i)+'-'+str(interaction_id),CC_celltype_name+'(cc)',NC_celltype_name[index]+'(nc)',score[index],k+1,1+pc_index_nc[j],normalized_ridge_coef[k,j] ,'Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)','Receptor(A)','Ligand(B)','GeneCor(Rec)','GeneCor(Lig)']
                                    header[7]=Found1[ele][0][0]
                                    header[8]=Found1[ele][1][0]
                                    header[9]=Found1[ele][0][1]
                                    header[10]=Found1[ele][1][1]
                                    header[11]=Found1[ele][0][2]
                                    header[12]=Found1[ele][1][2]
                                    header[13]=Found1[ele][0][3]
                                    header[14]=Found1[ele][1][3]
                                    for ri in range(15):
                                        worksheet.write(worksheetrow,ri,header[ri])
                                        worksheet_local.write(worksheetrow_local,ri,header[ri])
                                    worksheetrow+=1
                                    worksheetrow_local+=1


                                for ele in range(len(Found2)):
                                    header=[str(i)+'-'+str(interaction_id),NC_celltype_name[index]+'(nc)',CC_celltype_name+'(cc)',score[index],1+pc_index_nc[j],k+1,normalized_ridge_coef[k,j] ,'Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)','Receptor(A)','Ligand(B)','GeneCor(Rec)','GeneCor(Lig)']
                                    header[7]=Found2[ele][0][0]
                                    header[8]=Found2[ele][1][0]
                                    header[9]=Found2[ele][0][1]
                                    header[10]=Found2[ele][1][1]
                                    header[11]=Found2[ele][0][2]
                                    header[12]=Found2[ele][1][2]
                                    header[13]=Found2[ele][0][3]
                                    header[14]=Found2[ele][1][3]
                                    for ri in range(15):
                                        worksheet.write(worksheetrow,ri,header[ri])
                                        worksheet_local.write(worksheetrow_local,ri,header[ri])
                                    worksheetrow+=1
                                    worksheetrow_local+=1

                        fout.write('CC-Fa'+str(k+1)+'\t'+CC_celltype_name+'\t'+'%0.3f'%(score[index])+'\tNC-Fa'+str(1+pc_index_nc[j])+'\t'+NC_celltype_name[index]+'\tRegCoeff=%0.3f'%(normalized_ridge_coef[k,j])+'\t'+'pvalue=%0.2e'%pvalue[k,j]+'\t-log10(pvalue)=%0.2f'%(-np.log10(pvalue[k,j])))#str(interaction_id)
                        fout.write('\n')

                        fout.write('CC'+str(genesWithUP)+'\n')
                        fout.write('NC'+str(genesWithDown)+'\n')
                        fout.write('\n')
        fout.write('\n\n')
    workbook.close()





def find_LR_interactions_in_interacting_cell_types(input,choose_interacting_celltype_pair=[],choose_factors_id=[],pvalueCutoff=0.05,
correlation_with_spearman=True, LR_plot_NMF_Fa_thres=0.2, LR_plot_Exp_thres=0.2,saveas='pdf',transparent_mode=False,showit=True,figsize=(12,10)):
    """
    Inputs:

    | The main input is the output from gene_covariation_analysis.

    | Define the cell type pairs in the list you want to see the ligand-receptor communication in the central (CC) and niche cell types (NC).
    | The first element of the list is the central cell type, and the second element of the list is the niche one.
    | LR will perform for all interacting niches if the niche cell type is not defined.
    | choose_interacting_celltype_pair=[?,?]
    | If the list is Null, then LR plots will be saved for all the significant interacting niche cell types.

    | Define the factors id via the list you want to visualize in the ligand-receptor communication.
    | The first element of the list is the factor ID of the central cell type, and the second element of the list is the factor ID of the niche cell type.
    | choose_factors_id=[?,?]
    | If the list is Null, then LR plots will be saved for all the significant niche cell type factor interactions.


    | If True, visualize genes factor correlation with Spearman; otherwise, If False, then cosine.
    | (default) correlation_with_spearman=True

    | The cutoff used to find the significant central cell type factor and niche cell type factor interactions
    | (default) pvalueCutoff=0.05

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Ligand Receptor analysis uses factor cutoff to find the enriched LR pairs
    | (default) LR_plot_NMF_Fa_thres=0.2

    | Ligand Receptor analysis uses % expressed cell population cutoff to find the enriched LR pairs.
    | (default) LR_plot_Exp_thres=0.2

    | Background color in the figures
    | (default) transparent_mode=False

    | Dimension of the figure size.
    | Figure size on the X-axis direction is the (number of genes) multiplied by factor 12/34.
    | Figure size on the Y-axis direction is the (number of genes) multiplied by factor 10/44
    | All generated figure size are scaled according to the above factors
    | (initital figure size) figsize=(12,10)

    Outputs:

    | The LR interaction figures are saved in "./spatial_ct_ct_interactions/covariations_R*_F*/Plot_ligand_receptor_in_niche*"

    | Our analysis accounts for bidirectional cellular crosstalk interactions of ligands and receptors in cell types A and B.
    | It means in one way, it assumes a ligand diffused from cell type A and the receptor is detected in cell type B,
    | while in another way, it assumes a ligand is diffusing from cell type B, and the receptor is detected in cell type A.
    | Ligand-receptor plots show the bidirectional cellular crosstalk interaction of ligand and receptor in cell types A and B.

    """


    totalLRpairs,ligand,receptor,either=read_LigRecDb(input.LRdb)
    coeff_cutoff_for_log_reg=input.logistic_coef_cutoff
    coeff_cutoff_for_rid_reg=input.coeff_cutoff_for_rid_reg
    gene_set_names=input.gene_set_names
    LRcutoff=LR_plot_NMF_Fa_thres #Used in excel sheet to show the enrichment of ligand receptor intera

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)




    saveLRplots=input.covariation_dir+'Plot_ligand_receptor_in_niche/'
    create_directory(saveLRplots)
    saveLRplotsFirst=input.covariation_dir+'Plot_ligand_receptor_in_niche_cc_vs_nc/'
    create_directory(saveLRplotsFirst)
    saveLRplotsSecond=input.covariation_dir+'Plot_ligand_receptor_in_niche_nc_vs_cc/'
    create_directory(saveLRplotsSecond)

    print("LR figures for both ways are saved in following path ", saveLRplots)
    print("LR figures for CC to NC are saved in following path ", saveLRplotsFirst)
    print("LR figures for NC to CC are saved in following path ", saveLRplotsSecond)


    d={}
    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        clname=input.spatialcell_unique_clustername[i]
        d[clname]=clid


    if len(choose_interacting_celltype_pair)>0:
        choose_CC_celltypes=[choose_interacting_celltype_pair[0]]
    else:
        choose_CC_celltypes=[]

    perform=[]
    for fi in range(n):
        CC_celltype_name=input.spatialcell_unique_clustername[fi]
        if len(choose_CC_celltypes)==0:
            perform.append(fi)
        else:
            if CC_celltype_name in choose_CC_celltypes:
                perform.append(fi)

    for i in perform:
        clid=input.spatialcell_unique_clusterid[i]
        CC_corr_spearman,CC_PCA,CC_gene,CC_meanExpression,CC_popExpression,CC_corr_cosine,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        CC_celltype_name=input.spatialcell_unique_clustername[i]
        #temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        #index=temp[0]

        data=input.save_reg_coef[input.spatialcell_unique_clusterid[i]]
        coef_mu,intercept,alpha,xlabel,score,target,neighborhoodClass,pvalue,pve,rve=data

        NC_celltype_name=xlabel
        largest=np.max(abs(coef_mu))
        normalized_ridge_coef=coef_mu/largest

        ylabelname=[]
        componentlabel=[]
        for j in range(input.no_of_pc):
            ylabelname.append('CC_'+CC_celltype_name+'_Fa'+str(j+1))
            componentlabel.append('Fa'+str(j+1))

        for k in range(len(NC_celltype_name)):
            if score[k]>coeff_cutoff_for_log_reg:
                #in ylabelname first (# of pc) is the central cell type
                #and remaining are (# of pc) from the negihborhood cell type
                if CC_celltype_name!=NC_celltype_name[k]:
                    for j in range(input.no_of_pc):
                        ylabelname.append('NC_'+NC_celltype_name[k]+'_s'+'%0.3f'%score[k]+'_Fa'+str(j+1))

        pc_index_nc=[]
        for k in range(len(NC_celltype_name)):
            for j in range(input.no_of_pc):
                pc_index_nc.append(j)

        #normalized_ridge_coef  noofPC x (noofPC x +ve coff in log reg)


        interaction_id=0
        for k in range(normalized_ridge_coef.shape[0]):
            #k is PC of central cell type
            for j in range(normalized_ridge_coef.shape[1]):
                interaction_id+=1
                index=math.floor(j/input.no_of_pc)
                #index is the id neighboring cell type
                #if abs(normalized_ridge_coef[k,j])>coeff_cutoff_for_rid_reg:
                #pvalueCutoff=1
                if (pvalue[k,j]<pvalueCutoff)&(abs(normalized_ridge_coef[k,j])>coeff_cutoff_for_rid_reg):
                #if True:
                    if score[index]>coeff_cutoff_for_log_reg:
                        NC_corr_spearman,NC_PCA,NC_gene,NC_meanExpression,NC_popExpression,NC_corr_cosine,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[d[NC_celltype_name[index]]]
                        if correlation_with_spearman:
                            top_genes_in_CC,top_genes_in_NC,genesWithUP,genesWithDown,Found1,Found2=find_fold_change(CC_corr_spearman,NC_corr_spearman,CC_gene,k,pc_index_nc[j],totalLRpairs,LRcutoff,CC_meanExpression,NC_meanExpression,CC_popExpression,NC_popExpression,1)
                        else:
                            top_genes_in_CC,top_genes_in_NC,genesWithUP,genesWithDown,Found1,Found2=find_fold_change(CC_corr_cosine,NC_corr_cosine,CC_gene,k,pc_index_nc[j],totalLRpairs,LRcutoff,CC_meanExpression,NC_meanExpression,CC_popExpression,NC_popExpression,1)

                        common_genes=list(set(top_genes_in_CC).intersection(set(top_genes_in_NC)))

                        if len(choose_interacting_celltype_pair)>1:
                            if NC_celltype_name[index]==choose_interacting_celltype_pair[1]:
                                ncflag=1
                            else:
                                ncflag=0
                        else:
                            ncflag=1
                        if ncflag==1:
                            if len(choose_factors_id)==2:
                                flag=0
                                CC_factor=choose_factors_id[0]
                                NC_factor=choose_factors_id[1]
                                if (CC_factor==(k+1))& (NC_factor==(1+pc_index_nc[j])):
                                    flag=1
                            else:
                                flag=1
                            if flag==1:
                                plot_ligand_receptor_in_interacting_celltypes(CC_celltype_name,NC_celltype_name[index],score[index],k+1,1+pc_index_nc[j],normalized_ridge_coef[k,j],pvalue[k,j],Found1,Found2,saveLRplots,LR_plot_Exp_thres,saveas,transparent_mode,showit,figsize,'Both')
                                plot_ligand_receptor_in_interacting_celltypes(CC_celltype_name,NC_celltype_name[index],score[index],k+1,1+pc_index_nc[j],normalized_ridge_coef[k,j],pvalue[k,j],Found1,Found2,saveLRplotsFirst,LR_plot_Exp_thres,saveas,transparent_mode,showit,figsize,'First')
                                plot_ligand_receptor_in_interacting_celltypes(CC_celltype_name,NC_celltype_name[index],score[index],k+1,1+pc_index_nc[j],normalized_ridge_coef[k,j],pvalue[k,j],Found1,Found2,saveLRplotsSecond,LR_plot_Exp_thres,saveas,transparent_mode,showit,figsize,'Second')

    return 0


def make_excel_sheet_for_gene_correlation(input):
    """
    This function create the excel sheet for all compiled genes related to factors across different cell types.

    This organization facilitates a more structured and accessible
    representation of gene factors association in the cell types.

    These sheets were categorized to accommodate various types of information, including average gene
    expression (‘avg gene exp’ sheet), Spearman correlation values
    for different factors within scRNASeq data (‘spearman scRNAseq Fa(i)’ sheets),
    cosine similarity values within scRNASeq data (‘cosine scRNAseq Fa(i)’ sheet),
    as well as a subset of genes unique to spatial data (‘spearman spatial Fa(i)’and
    ‘cosine spatial Fa(i)’ sheets).

    | In these sheet names, ‘i’ corresponds to the factor ID.
    | As columns are traversed within each sheet, factors representing various cell types are presented.
    | These genes in factors are sorted based on their factor ID present in the sheet name.
    | To enhance clarity and insight, we employed a color coding scheme that distinguishes genes as ligands (depicted in blue), receptors (in red), or both ligand and receptor functions (in magenta).
    """



    totalLRpairs,ligand,receptor,either=read_LigRecDb(input.LRdb)

    workbook = xlsxwriter.Workbook(input.covariation_dir+'gene_correlation.xlsx')

    worksheetAvgGeneExp= workbook.add_worksheet('avg gene exp')
    worksheetFullGene_spearman=[]
    worksheetFullGene_cosine=[]
    for i in range(input.no_of_pc):
        worksheetFullGene_spearman.append( workbook.add_worksheet('spearman scRNAseq Fa'+str(i+1)))
        worksheetFullGene_cosine.append( workbook.add_worksheet('cosine scRNAseq Fa'+str(i+1)))

    worksheetSpatialGene_spearman=[]
    worksheetSpatialGene_cosine=[]
    for i in range(input.no_of_pc):
        worksheetSpatialGene_spearman.append( workbook.add_worksheet('spearman spatial Fa'+str(i+1)))
        worksheetSpatialGene_cosine.append( workbook.add_worksheet('cosine spatial Fa'+str(i+1)))


    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))


    #outputFolder=maindir+'geneCorr'+str(input.no_of_pc)+'/'
    #create_directory(outputFolder)


    genenames=sorted(list(input.ad_sp.var_names.to_numpy()))
    n=len(input.spatialcell_unique_clustername)

    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        CC_corr_spearman,CC_PCA,gene,CC_meanExpression,CC_popExpression,CC_corr_cosine,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        worksheetrow=0
        worksheetAvgGeneExp.write(worksheetrow,3*i,input.spatialcell_unique_clustername[i])
        for j in range(input.no_of_pc):
            worksheetFullGene_spearman[j].write(worksheetrow,(input.no_of_pc+2)*i,input.spatialcell_unique_clustername[i])
            worksheetSpatialGene_spearman[j].write(worksheetrow,(input.no_of_pc+2)*i,input.spatialcell_unique_clustername[i])
            worksheetFullGene_cosine[j].write(worksheetrow,(input.no_of_pc+2)*i,input.spatialcell_unique_clustername[i])
            worksheetSpatialGene_cosine[j].write(worksheetrow,(input.no_of_pc+2)*i,input.spatialcell_unique_clustername[i])
        worksheetrow+=1
        fixvalue=worksheetrow


        index=np.argsort(-CC_meanExpression)
        for j in range(len(index)):
            worksheetAvgGeneExp.write(j+2,3*i+1,CC_meanExpression[index[j]])
            worksheetAvgGeneExp.write(j+2,3*i,gene[index[j]])


        red = workbook.add_format({'bold': True,'color': 'red'})
        green = workbook.add_format({'bold': True,'color': 'green'})
        blue = workbook.add_format({'bold': True,'color': 'blue'})
        magenta = workbook.add_format({'bold': True,'color': 'magenta'})


        '''
        fig,(ax)=plt.subplots(1,1,figsize=(8,6))
        ax.plot(CC_corr[:,0],CC_corr[:,1],'.',markersize=1)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(input.spatialcell_unique_clustername[i])
        #fig.tight_layout()
        fig.savefig(outputFolder+'correlation_'+input.spatialcell_unique_clustername[i]+'.png',bbox_inches='tight',transparent=True,dpi=300)
        plt.close('all')
        '''

        headersave_full,headersave_common,sort_full,sort_common=sorting_of_factors_for_showing_the_value_in_excelsheet(CC_corr_spearman,input.no_of_pc,gene,genenames)

        for k in range(input.no_of_pc):
            worksheetrow=fixvalue
            indsort=np.argsort(-np.array(sort_full[k]))
            for rj in range(len(indsort)):
                header=headersave_full[indsort[rj]]
                mygene=header[0]
                genecolor=''
                if mygene.upper() in ligand:
                    genecolor=blue
                elif mygene.upper() in receptor:
                    genecolor=red
                elif mygene.upper() in either:
                    genecolor=magenta

                for ri in range(len(header)):
                    worksheetFullGene_spearman[k].write(worksheetrow,(input.no_of_pc+2)*i+ri,header[ri],genecolor)
                worksheetrow+=1

            worksheetrow=fixvalue
            indsort=np.argsort(-np.array(sort_common[k]))
            for rj in range(len(indsort)):
                header=headersave_common[indsort[rj]]
                for ri in range(len(header)):
                    worksheetSpatialGene_spearman[k].write(worksheetrow,(input.no_of_pc+2)*i+ri,header[ri])
                worksheetrow+=1



        headersave_full,headersave_common,sort_full,sort_common=sorting_of_factors_for_showing_the_value_in_excelsheet(CC_corr_cosine,input.no_of_pc,gene,genenames)
        for k in range(input.no_of_pc):
            worksheetrow=fixvalue
            indsort=np.argsort(-np.array(sort_full[k]))
            for rj in range(len(indsort)):
                header=headersave_full[indsort[rj]]
                mygene=header[0]
                genecolor=''
                if mygene.upper() in ligand:
                    genecolor=blue
                elif mygene.upper() in receptor:
                    genecolor=red
                elif mygene.upper() in either:
                    genecolor=magenta

                for ri in range(len(header)):
                    worksheetFullGene_cosine[k].write(worksheetrow,(input.no_of_pc+2)*i+ri,header[ri],genecolor)
                worksheetrow+=1

            worksheetrow=fixvalue
            indsort=np.argsort(-np.array(sort_common[k]))
            for rj in range(len(indsort)):
                header=headersave_common[indsort[rj]]
                for ri in range(len(header)):
                    worksheetSpatialGene_cosine[k].write(worksheetrow,(input.no_of_pc+2)*i+ri,header[ri])
                worksheetrow+=1

    workbook.close()


def pathway_analysis(input,NOG_pathway=50, choose_factors_id=[],correlation_with_spearman=True,saveas='pdf',savefigure=False,
 positively_correlated=True,rps_rpl_mt_genes_included=True,choose_celltypes=[],circlesize=12,pathwayCutoff=0.5,
 pathwayorganism='Mouse',database=['GO_Biological_Process_2021','BioPlanet_2019','Reactome_2016']):#background_geneName,background_expression

    """
    Inputs:

    The main input is the output from gene_covariation_analysis.

    | Number of genes associated with NMF factors to search in the pathway database. If you don't observe any pathway, then increase the number of genes and try with 100 or 200 or try with different databases.
    | (default) NOG_pathway=50

    | Define the factors id via the list you want to visualize in the pathway.
    | choose_factors_id=[]
    | If the list is Null, then pathway will be saved for all the factors.

    | The cell type that you want to analyze via the pathway Enrichr library.
    | (default) choose_celltypes=[]
    | If the list is empty, the output will show for all the cell types.

    | If True, visualize genes factor correlation with Spearman; otherwise, If False, then cosine.
    | (default) correlation_with_spearman=True

    | If the genes factor association is chosen as Spearman correlation, the belonging genes can be selected as either positively correlated (True) or negatively correlated (False).
    | (default) positively_correlated=True

    | For pathway analysis, decide whether to include rps, rpl, and mt genes. If " True, " they are included; otherwise, " No. "
    | (default) rps_rpl_mt_genes_included=True

    | The gseapy parameter to find the pathway-enriched library from the top gene list from each factor of a given cell type
    | (default) pathwayCutoff=0.5

    | Organisms used in the gseapy package
    | (default) pathwayorganism='Mouse'

    | The database used in the gseapy package for pathway analysis. The default uses all three databases. For other non-mentioned databases, please follow the GSEApy tutorial.
    | (default) database=['GO_Biological_Process_2021','BioPlanet_2019','Reactome_2016']
    | If you are interested in knowing what other database is available, then you can check for these species ‘Human,’ ‘Mouse,’ ‘Yeast,’ ‘Fly,’ ‘Fish,’ and ‘Worm’ in the following way:
    | >>> import gseapy as gp
    | >>> mouse = gp.get_library_name(organism='Mouse')
    | >>> human = gp.get_library_name(organism='Human')

    | The cell type that you want to see the covariation regression pattern
    | (default) choose_celltypes=[]
    | If the list is empty, the output will show for all the cell types.

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | The size of the point in the pathway figure. Increase if the dots are too small.
    | (default) circlesize=12

    Outputs:

    | The pathways figures are saved in "./spatial_ct_ct_interactions/covariations_R*_F*/Pathway_figures/"

    """


    savename=input.covariation_dir+'Pathway_figures/'
    create_directory(savename)

    print("The pathway figures are saved in ", savename)


    #coeff_cutoff_for_log_reg=input.logistic_coef_cutoff
    #coeff_cutoff_for_rid_reg=input.coeff_cutoff_for_rid_reg
    #gene_set_names=input.gene_set_names
    nog=NOG_pathway
    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)

    perform=[]
    found=[]
    for fi in range(n):
        CC_celltype_name=input.spatialcell_unique_clustername[fi]
        if len(choose_celltypes)==0:
            perform.append(fi)
        else:
            if CC_celltype_name in choose_celltypes:
                perform.append(fi)
                found.append(CC_celltype_name)
    if len(choose_celltypes)!=0:
        print("cell types found ",found)

    for fi in perform:
        clid=input.spatialcell_unique_clusterid[fi]
        spearman_factors,CC_PCA,CC_gene,CC_meanExpression,CC_popExpression,cosine_factors,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        CC_celltype_name=input.spatialcell_unique_clustername[fi]

        for j in range(input.no_of_pc):
            if correlation_with_spearman:
                source=spearman_factors[:,j]
            else:
                source=cosine_factors[:,j]
            ind=np.argsort(-source)
            interestofGene=[]
            value=[]
            for k in range(len(source)):
                temp=CC_gene[ind[k]]
                if rps_rpl_mt_genes_included:
                    flag=1
                else:
                    flag=1
                    if temp[0:3]=='Rps':
                        flag=0
                    if temp[i,0][0:3]=='Rpl':
                        flag=0
                    if temp[i,0][0:3]=='mt-':
                        flag=0
                if flag==1:
                    interestofGene.append(CC_gene[ind[k]])
                    value.append(source[ind[k]])

            value=np.array(value)
            interestofGene=np.array(interestofGene)
            if positively_correlated:
                index=np.argsort(-value)
                tname='pos'
            else:
                index=np.argsort(value)
                tname='neg'

            value=list(value[index])
            interestofGene=list(interestofGene[index])

            if len(interestofGene)>nog:
                va1=value[0:nog]
                ga1=interestofGene[0:nog]
                cutoff=va1[-1]
            else:
                ga1=interestofGene
                va1=value

            ccname=remove_extra_character_from_name(CC_celltype_name)
            titlename=tname+' Fa'+str(j+1)+' '+CC_celltype_name+' c'+str(int(100*cutoff))
            sname1=tname+'Fa'+str(j+1)+'_'+ccname+'_c'+str(int(100*cutoff))

            if len(choose_factors_id)>0:
                if (j+1) in choose_factors_id:
                    flag=1
                else:
                    flag=0
            else:
                flag=1

            if flag==1:
                for i in range(len(database)):
                    titlename1=titlename+'['+database[i]+']'+' #G='+str(len(ga1))
                    sname2=sname1+'_'+database[i]
                    finalsavename=savename+sname2+'.'+saveas

                    enr_res1 = gseapy.enrichr(gene_list=ga1,organism=pathwayorganism,gene_sets=database[i], cutoff = pathwayCutoff)
                    #enr_res1 = gseapy.enrichr(gene_list=g1,organism='Mouse',gene_sets=background_model,description='pathway',cutoff = 0.5)
                    finalsavename.replace(' ','_')
                    try:
                        #gseapy.barplot(enr_res1.res2d,title=titlename1,ofname=finalsavename,fontsize=12)#database[i]+titlename
                        if savefigure:
                            gseapy.dotplot(enr_res1.res2d,title=titlename1,ofname=finalsavename,fontsize=12,size=circlesize,cmap = plt.cm.autumn_r)
                        else:
                            gseapy.dotplot(enr_res1.res2d,title=titlename1,fontsize=12,size=circlesize,cmap = plt.cm.autumn_r)
                    except Exception as e: #Exception: Error getting the Enrichr libraries
                        pass


def extract_and_plot_top_genes_from_chosen_factor_in_celltype(input,choose_celltype,choose_factor_id,top_NOG=30,rps_rpl_mt_genes_included=True,
correlation_with_spearman=True,positively_correlated=True,saveas='pdf',cmap='RdBu_r',transparent_mode=False,showit=True,figsize=(5, 6)):

    """
    Inputs:

    The main input is the output from gene_covariation_analysis.

    Define the cell type that you want to see.

    Define the factor id of the cell type that you want to see.

    | Number of genes to visualize
    | (default) top_NOG=30

    | For pathway analysis, decide whether to include rps, rpl, and mt genes. If " True, " they are included; otherwise, " No. "
    | (default) rps_rpl_mt_genes_included=True

    | If True, visualize genes factor correlation with Spearman; otherwise, If False, then cosine.
    | (default) correlation_with_spearman=True

    | If the genes factor association is chosen as Spearman correlation, the belonging genes can be selected as either positively correlated (True) or negatively correlated (False).
    | (default) positively_correlated=True

    | Define the colormap for visualizing factors
    | (default) cmap='RdBu_r'

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas=’pdf’

    | Dimension of the figure size
    | (default) figsize=(15,10)

    | Background color in the figures
    | (default) transparent_mode=False

    Outputs:

    | Return the data frame of the gene, factor, average expression, and proportion of population expressed that gene.
    | The figures are saved in spatial_ct_ct_interactions/covariations_R*_F*/dotplots/Factors*

    """

    savefigdir=input.covariation_dir+ 'dotplots/'
    create_directory(savefigdir)

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)
    perform=[]
    for fi in range(n):
        CC_celltype_name=input.spatialcell_unique_clustername[fi]
        if CC_celltype_name==choose_celltype:
                perform.append(fi)
    if len(perform)==0:
        print("Cell type name do not match")
        flag_correct=0
    else:
        flag_correct=1

    if 1<=choose_factor_id<=input.no_of_pc:
        flag_correct=1
    else:
        print("Factor ID is wrong")
        flag_correct=0

    df=0
    if flag_correct==1:
        for fi in perform:
            clid=input.spatialcell_unique_clusterid[fi]
            spearman_factors,CC_PCA,CC_gene,CC_meanExpression,CC_popExpression,cosine_factors,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
            CC_celltype_name=input.spatialcell_unique_clustername[fi]
            mu=CC_meanExpression
            pop=CC_popExpression

            #for j in range(input.no_of_pc):
            if True:
                if correlation_with_spearman:
                    source=spearman_factors[:,choose_factor_id-1]
                else:
                    source=cosine_factors[:,choose_factor_id-1]
                ind=np.argsort(-source)
                interestofGene=[]
                value_fact=[]
                value_pop=[]
                value_avgexp=[]
                for k in range(len(source)):
                    temp=CC_gene[ind[k]]
                    if rps_rpl_mt_genes_included:
                        flag=1
                    else:
                        flag=1
                        if temp[0:3]=='Rps':
                            flag=0
                        if temp[i,0][0:3]=='Rpl':
                            flag=0
                        if temp[i,0][0:3]=='mt-':
                            flag=0
                    if flag==1:
                        interestofGene.append(CC_gene[ind[k]])
                        value_fact.append(source[ind[k]])
                        value_pop.append(pop[ind[k]])
                        value_avgexp.append(mu[ind[k]])

                value_fact=np.array(value_fact)
                value_pop=np.array(value_pop)
                value_avgexp=np.array(value_avgexp)

                interestofGene=np.array(interestofGene)
                index_pos=np.argsort(-value_fact)
                index_neg=np.argsort(value_fact)

                gp1=list(interestofGene[index_pos])
                gn1=list(interestofGene[index_neg])

                gex=np.zeros((top_NOG,1),dtype=float)

                vp1=list(value_fact[index_pos])[0:top_NOG]
                vn1=list(value_fact[index_neg])[0:top_NOG]
                pos_pop1=list(value_pop[index_pos])[0:top_NOG]
                neg_pop1=list(value_pop[index_neg])[0:top_NOG]
                pos_avg1=list(value_avgexp[index_pos])[0:top_NOG]
                neg_avg1=list(value_avgexp[index_neg])[0:top_NOG]

                if positively_correlated:
                    nvr1=np.hstack((np.reshape(vp1,(len(vp1),1)),gex))
                    nvr2=np.hstack((np.reshape(pos_pop1,(len(pos_pop1),1)),gex))
                    nvr3=np.hstack((np.reshape(pos_avg1,(len(pos_avg1),1)),gex))
                    comgene=gp1[0:top_NOG]
                    title='Pos Fa'+str(choose_factor_id)

                    d = {'Gene': comgene, 'Fa': vp1,'mean_expression':pos_avg1,'proportion_of_population_expressed':pos_pop1}
                    df = pd.DataFrame(data=d)
                    df.set_index('Gene')
                else:
                    nvr1=np.hstack((np.reshape(vn1,(len(vn1),1)),gex))
                    nvr2=np.hstack((np.reshape(neg_pop1,(len(neg_pop1),1)),gex))
                    nvr3=np.hstack((np.reshape(neg_avg1,(len(neg_avg1),1)),gex))
                    comgene=gn1[0:top_NOG]
                    title='Neg Fa'+str(choose_factor_id)
                    d = {'Gene': comgene, 'Fa': vn1,'mean_expression':neg_avg1,'proportion_of_population_expressed':neg_pop1}
                    df = pd.DataFrame(data=d)
                    df.set_index('Gene')


            fig, ax = plt.subplots(1,2,figsize=figsize)
            x,y,z,bigs=findXYZC(nvr1,nvr2)
            p0=ax[0].scatter(x,y,s=bigs,marker='o',c=z,cmap=cmap) #'cm.cmap_name
            x,y,z,bigs=findXYZC(nvr3,nvr2)
            p1=ax[1].scatter(x,y,s=bigs,marker='o',c=z,cmap=cmap)
            plt.colorbar(p0,ax=ax[0],shrink=0.5)
            plt.colorbar(p1,ax=ax[1],shrink=0.5)
            kw = dict(prop="sizes", num=4, alpha=0.6, fmt="% {x:.0f}")
            legend2 = ax[0].legend(*p0.legend_elements(**kw),loc="lower center", bbox_to_anchor=(0.25, -0.25),title="Fraction of cells expressed",frameon=False)#
            ax[0].set_title(title)
            ax[1].set_title('Avg expression')

            for j in range(2):
                ax[j].set_yticks(range(len(nvr1)))
                ax[j].set_yticklabels(comgene,style='italic')
                ax[j].set_xticks([])#range(1))
                ax[j].set_xticklabels([])#xlabels[i],rotation=30)
                ax[j].set_xlim([-0.5,0.5])
                ax[j].set_ylim([-0.5,len(nvr1)+0.5])


            #create_subtitle(fig, grid[0, ::], CC_celltype_name+' Spearman correlation')
            #create_subtitle(fig, grid[1, ::],  CC_celltype_name+' log(avg expression)')
            fig.tight_layout()
            print("The figures are saved: ", savefigdir+'Factors_'+remove_extra_character_from_name(CC_celltype_name)+'.'+saveas)
            plt.savefig(savefigdir+'Factors_'+remove_extra_character_from_name(CC_celltype_name)+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
            if showit:
                pass
            else:
                plt.close('all')

    return df






def plot_top_selected_genes_for_all_factors_from_chosen_celltype(input,choose_celltypes=[],top_NOG=20,rps_rpl_mt_genes_included=True,correlation_with_spearman=True,saveas='pdf',transparent_mode=False,showit=True,figsize=(12, 10)):
    """
    Inputs:

    The main input is the output from gene_covariation_analysis.

    | Number of genes to visualize
    | (default) top_NOG=20

    | If True, visualize genes factor correlation with Spearman; otherwise, If False, then cosine.
    | (default) correlation_with_spearman=True

    | For pathway analysis, decide whether to include rps, rpl, and mt genes. If " True, " they are included; otherwise, " No. "
    | (default) rps_rpl_mt_genes_included=True

    | The cell type that you want to see the covariation regression pattern
    | (default) choose_celltypes=[]
    | If the list is empty, the output will show for all the cell types.

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Background color in the figures
    | (default) transparent_mode=False

    | Dimension of the figure size.
    | (default) figsize=(12,10)

    Outputs:

    | The gene visualization figures are saved in ./spatial_ct_ct_interactions/covariations_R*_F*/dotplots/*

    """
    savefigdir=input.covariation_dir+ 'dotplots/'
    create_directory(savefigdir)

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)
    perform=[]
    found=[]
    for fi in range(n):
        CC_celltype_name=input.spatialcell_unique_clustername[fi]
        if len(choose_celltypes)==0:
            perform.append(fi)
        else:
            if CC_celltype_name in choose_celltypes:
                perform.append(fi)
                found.append(CC_celltype_name)
    if len(choose_celltypes)!=0:
        print("cell types found ",found)

    for fi in perform:
        clid=input.spatialcell_unique_clusterid[fi]
        spearman_factors,CC_PCA,CC_gene,CC_meanExpression,CC_popExpression,cosine_factors,alpha=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        CC_celltype_name=input.spatialcell_unique_clustername[fi]
        mu=np.log(CC_meanExpression)
        pop=CC_popExpression

        nvr1=[]
        nvr2=[]
        nvr3=[]
        comgene=[]
        for j in range(input.no_of_pc):
            if correlation_with_spearman:
                source=spearman_factors[:,j]
            else:
                source=cosine_factors[:,j]
            ind=np.argsort(-source)
            interestofGene=[]
            value_fact=[]
            value_pop=[]
            value_avgexp=[]
            for k in range(len(source)):
                temp=CC_gene[ind[k]]
                if rps_rpl_mt_genes_included:
                    flag=1
                else:
                    flag=1
                    if temp[0:3]=='Rps':
                        flag=0
                    if temp[i,0][0:3]=='Rpl':
                        flag=0
                    if temp[i,0][0:3]=='mt-':
                        flag=0
                if flag==1:
                    interestofGene.append(CC_gene[ind[k]])
                    value_fact.append(source[ind[k]])
                    value_pop.append(pop[ind[k]])
                    value_avgexp.append(mu[ind[k]])

            value_fact=np.array(value_fact)
            value_pop=np.array(value_pop)
            value_avgexp=np.array(value_avgexp)

            interestofGene=np.array(interestofGene)
            index_pos=np.argsort(-value_fact)
            index_neg=np.argsort(value_fact)

            gp1=list(interestofGene[index_pos])
            gn1=list(interestofGene[index_neg])

            comgene.append(gp1[0:top_NOG])
            comgene.append(gn1[0:top_NOG])

            gex=np.zeros((top_NOG,1),dtype=float)

            vp1=list(value_fact[index_pos])[0:top_NOG]
            vn1=list(value_fact[index_neg])[0:top_NOG]
            pos_pop1=list(value_pop[index_pos])[0:top_NOG]
            neg_pop1=list(value_pop[index_neg])[0:top_NOG]
            pos_avg1=list(value_avgexp[index_pos])[0:top_NOG]
            neg_avg1=list(value_avgexp[index_neg])[0:top_NOG]

            nvr1.append(np.hstack((np.reshape(vp1,(len(vp1),1)),gex)))
            nvr1.append(np.hstack((np.reshape(vn1,(len(vn1),1)),gex)))
            nvr2.append(np.hstack((np.reshape(pos_pop1,(len(pos_pop1),1)),gex)))
            nvr2.append(np.hstack((np.reshape(neg_pop1,(len(neg_pop1),1)),gex)))
            nvr3.append(np.hstack((np.reshape(pos_avg1,(len(pos_avg1),1)),gex)))
            nvr3.append(np.hstack((np.reshape(neg_avg1,(len(neg_avg1),1)),gex)))


        fig, ax = plt.subplots(2,6,figsize=figsize)
        title=['Pos Fa1','Neg Fa1','Pos Fa2','Neg Fa2','Pos Fa3','Neg Fa3']

        for i in range(6):
                x,y,z,bigs=findXYZC(nvr1[i],nvr2[i])
                p0=ax[0,i].scatter(x,y,s=bigs,marker='o',c=z,cmap='RdBu_r') #'cm.cmap_name
                x,y,z,bigs=findXYZC(nvr3[i],nvr2[i])
                p1=ax[1,i].scatter(x,y,s=bigs,marker='o',c=z,cmap='RdBu_r')
                plt.colorbar(p0,ax=ax[0,i],shrink=0.5)
                plt.colorbar(p1,ax=ax[1,i],shrink=0.5)
                kw = dict(prop="sizes", num=5, alpha=0.6, fmt="% {x:.0f}")
                legend2 = ax[1,i].legend(*p1.legend_elements(**kw),loc="upper right", title="Fraction of \ncells \nexpressed",bbox_to_anchor=(1.0, 0),frameon=False)

        for i in range(6):
            for j in range(2):
                ax[j,i].set_yticks(range(len(nvr1[i])))
                ax[j,i].set_yticklabels(comgene[i],style='italic')
                ax[j,i].set_xticks([])#range(1))
                ax[j,i].set_xticklabels([])#xlabels[i],rotation=30)
                ax[j,i].set_xlim([-0.5,0.5])
                ax[j,i].set_ylim([-0.5,len(nvr1[i])+0.5])
                ax[0,i].set_title(title[i])

        grid = plt.GridSpec(2, 6)
        create_subtitle(fig, grid[0, ::], CC_celltype_name+' Spearman correlation')
        create_subtitle(fig, grid[1, ::],  CC_celltype_name+' log(avg expression)')

        fig.tight_layout()
        print("The figures are saved: ", savefigdir+remove_extra_character_from_name(CC_celltype_name)+'.'+saveas)
        plt.savefig(savefigdir+remove_extra_character_from_name(CC_celltype_name)+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
        if showit:
            pass
        else:
            plt.close('all')







def create_directory(outputFolder):
    "This function create empty directory."
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)


def find_index(sp_genename,sc_genename):
    "Helper function used in gene_covariation_analysis to find the common gene space submatrix between two modalities."
    index_sc=[]
    index_sp=[]
    d={}
    for j in range(len(sc_genename)):
        name=sc_genename[j]
        d[name]=j

    for i in range(len(sp_genename)):
        name=sp_genename[i]
        try:
            d[name]
            flag=1
        except KeyError:
            flag=0
        if flag==1:
            index_sc.append(d[name])
            index_sp.append(i)
    return index_sp,index_sc

def read_spatial_data(clusterFilename,celltypeFilename):
    """
    Helper function for gene_covariation_analysis to read the clusters information.
    """

    df=pd.read_csv(celltypeFilename,sep='\t',header=None)
    data=df.to_numpy()
    spatialcell_unique_clustername=data[:,1]
    spatialcell_unique_clusterid=data[:,0]
    CTname=spatialcell_unique_clustername
    CTid=spatialcell_unique_clusterid

    df=pd.read_csv(clusterFilename)
    louvainFull=df.to_numpy()


    celltype={}
    cellsinCT={}
    index=[]
    for i in range(len(louvainFull)):
        clu_id=louvainFull[i][1]
        cel_id=louvainFull[i][0]
        if clu_id in CTid:
            index.append(i)
            #celltype[cel_id]=clu_id
            if clu_id not in cellsinCT:
                cellsinCT[clu_id]=[cel_id]
            else:
                cellsinCT[clu_id].append(cel_id)

    louvain=louvainFull[index,:]
    annotation_spatial_barcode_id= louvain[:,0]
    annotation_spatial_cluster_id= louvain[:,1]

    d={}
    for i in range(len(spatialcell_unique_clustername)):
        d[spatialcell_unique_clusterid[i]]=spatialcell_unique_clustername[i]
    annotation_spatial_celltypename=[]
    for i in range(len(annotation_spatial_cluster_id)):
        annotation_spatial_celltypename.append(d[annotation_spatial_cluster_id[i]])
    annotation_spatial_celltypename=np.array(annotation_spatial_celltypename)

    return annotation_spatial_celltypename,annotation_spatial_barcode_id,annotation_spatial_cluster_id,spatialcell_unique_clustername,spatialcell_unique_clusterid


def find_correlation_bw_genes_and_PC_component_in_singlecell(KcomponentCluster,clusterExpression):
    """
    Helper function used in find_PC_of_invidualCluster_in_SC to find the Spearman correlation between common genes scRNAseq factors and scRNAseq expression.
    """
    mat=np.zeros((clusterExpression.shape[1],KcomponentCluster.shape[1]),dtype=float)
    for i in range(clusterExpression.shape[1]):
        v1=clusterExpression[:,i]
        for j in range(KcomponentCluster.shape[1]):
            v2=KcomponentCluster[:,j]
            #corr,_ = pearsonr(v1,v2)
            corr,_ =spearmanr(v1,v2)
            #corr=np.corrcoef(v1,v2)
            mat[i,j]=corr

    # mat shape is (# of genes x # of pc) it is a correlation between (PC and genes) of the single cell cluster
    # KcomponentCluster shape is (# of single cell in a single cell cluster x # of pc)
    # clusterExpression shape is (# of single cell in a single cell cluster x # of genes)
    mat=np.nan_to_num(mat)
    return mat

def find_correlation_bw_genes_and_PC_component_in_singlecell_cosine(KcomponentCluster,clusterExpression):
    """
    Helper function used in find_PC_of_invidualCluster_in_SC
    to find the cosine similarity between common genes scRNAseq factors and scRNAseq expression.
    """
    #same vector =1 perpendicular vector 0
    #print(KcomponentCluster.shape,clusterExpression.shape)
    #mat=np.zeros((clusterExpression.shape[1],KcomponentCluster.shape[1]),dtype=float)
    mat=cosine_similarity(clusterExpression.T,KcomponentCluster.T)
    return mat



def top_genes_in_correlation_list_without(genename,corr_NMFfactors_genes,n_top_words):
        """
        Helper function for sorting the factors values in plot_cosine_and_spearman_correlation_to_factors.
        """
        top_genes_assoc_factors=[]
        for topic_idx, topic in enumerate(corr_NMFfactors_genes.T):
            top_features_ind = topic.argsort()[: -n_top_words - 1 : -1]
            for i in top_features_ind:
                if i not in top_genes_assoc_factors:
                    top_genes_assoc_factors.append(i)
        gname=genename[top_genes_assoc_factors]
        mat=corr_NMFfactors_genes[top_genes_assoc_factors,:]

        return gname,mat


def alignment_score(H,spH,ind_H,ind_spH):
    """
    The helper function is used in find_PC_of_invidualCluster_in_SC to find the alignment score during the integrated NMF step.
    """

    #print(H.shape,spH.shape,len(ind_H),len(ind_spH))
    r1=H[:,ind_H]
    r2=spH[:,ind_spH]
    comb=np.hstack((r1,r2)).T
    n=len(ind_H)
    knn=max([2,np.ceil(0.01*n) ])
    n_jobs=-1
    k_d,k_ind = cKDTree(comb).query(x=comb, k=knn, workers=n_jobs)

    avgc1=0
    for i in range(n):
        neigh=k_ind[i]
        c1=0
        for j in range(len(neigh)):
            if neigh[j]<n:
                c1=c1+1
        avgc1=avgc1+c1
    avgc1=avgc1/n
    #doi:10.1038/nbt.4096
    score=1 - ((avgc1 - (knn/n) ) / (knn - (knn/n) ))

    return score

def multiplicative_method(W,H,A,max_iter):
    """
    Helper function used in find_PC_of_invidualCluster_in_SC to perform the
    conventional NMF.
    """
    norms = []
    e = 1.0e-10
    for n in range(max_iter):
        # Update H
        W_TA = W.T@A
        W_TWH = W.T@W@H+e
        for i in range(np.size(H, 0)):
            for j in range(np.size(H, 1)):
                H[i, j] = H[i, j] * W_TA[i, j] / W_TWH[i, j]
        # Update W
        #AH_T = A@H.T
        #WHH_T =  W@H@H.T+ e
        #for i in range(np.size(W, 0)):
        #    for j in range(np.size(W, 1)):
        #        W[i, j] = W[i, j] * AH_T[i, j] / WHH_T[i, j]

        norm = np.linalg.norm(A - W@H, 'fro')
        norms.append(norm)
    return W ,H ,norms


def remove_extra_character_from_name(name):
    """
    This function remove the special characters from the cell type names so it should not throw error while saving the figures.
    """
    name=name.replace('/','_')
    name=name.replace(' ','_')
    name=name.replace('"','')
    name=name.replace("'",'')
    name=name.replace(')','')
    name=name.replace('(','')
    name=name.replace('+','p')
    name=name.replace('-','n')
    name=name.replace('.','')
    return name


def find_PC_of_invidualCluster_in_SC(scbarcode,iNMFmode,scadata,no_of_pc,spbarcode,spadata, sct_ad_sc_full,celltype_name,cutoff_to_count_exp_cell_population):

    "Helper function used in compute_PC_space."
    cellname=sct_ad_sc_full.obs_names.to_numpy()
    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i
    index=[]
    for i in range(len(scbarcode)):
        index.append(d[scbarcode[i]])
    full_genes_sc=sct_ad_sc_full[index,:].copy()

    #common gene single cell
    cellname=scadata.obs_names.to_numpy()
    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i
    index=[]
    for i in range(len(scbarcode)):
        index.append(d[scbarcode[i]])
    sct_ad_sc=scadata[index,:].copy()
    sc_cellname=sct_ad_sc.obs_names.to_numpy()

    #common gene spatial
    cellname=spadata.obs_names.to_numpy()
    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i
    index=[]
    for i in range(len(spbarcode)):
        index.append(d[spbarcode[i]])

    sct_ad_sp=spadata[index,:].copy()
    sp_cellname=sct_ad_sp.obs_names.to_numpy()


    if scipy.__version__=='1.7.3':
        matrixtype="<class 'scipy.sparse.csr.csr_matrix'>"
    else:
        matrixtype="<class 'scipy.sparse._csr.csr_matrix'>"

    tp_sc=str(type(full_genes_sc.X))
    if tp_sc==matrixtype:
        CbyG=full_genes_sc.X.toarray()
    else:
        CbyG=full_genes_sc.X


    tp_sc=str(type(sct_ad_sc.X))
    if tp_sc==matrixtype:
        msc=sct_ad_sc.X.toarray()
    else:
        msc=sct_ad_sc.X

    tp_sp=str(type(sct_ad_sp.X))
    if tp_sp==matrixtype:
        msp=sct_ad_sp.X.toarray()
    else:
        msp=sct_ad_sp.X

    #replace nan to zero
    #msp=np.nan_to_num(msp)
    #msc=np.nan_to_num(msc)

    #msc=msc/np.sum(msc)
    #msp=msp/np.sum(msp)
    #CbyG=CbyG/np.sum(CbyG)
    genename_joint=sct_ad_sc.var_names.to_numpy()
    genename_spatial=sct_ad_sp.var_names.to_numpy()


    #Gene based normalization
    #msc=np.log10(1+msc)
    #msp=np.log10(1+msp)
    #CbyG=np.log10(1+CbyG)
    std1=np.std(msc,axis=0)
    std2=np.std(msp,axis=0)
    ind=np.where((std1>0)&(std2>0))
    index=ind[0]
    n=len(index)
    v1=np.zeros((msc.shape[0],n),dtype=float)
    v2=np.zeros((msp.shape[0],n),dtype=float)
    for i in range(n):
        v1[:,i]=msc[:,index[i]]/std1[index[i]]
        v2[:,i]=msp[:,index[i]]/std2[index[i]]
    #sum1=np.std(v1,axis=0)
    #sum2=np.std(v2,axis=0)


    datasets=[v1,v2]
    n1=msc.shape[0]
    n2=msp.shape[0]
    threshold=0.001
    old_score=1

    if iNMFmode==True:
        for alpha in range(0,51,2):
            arr1=[*range(n1)]
            arr2=[*range(n2)]
            if n1>n2:
                np.random.shuffle(arr1)
                arr1=arr1[0:n2]
            else:
                np.random.shuffle(arr2)
                arr2=arr2[0:n1]

            H,spH,W,V,spV = iNMF(datasets,no_of_pc,value_lambda=alpha,print_obj=False)
            spW=W
            score=alignment_score(H,spH,arr1,arr2)
            #print(score,alpha,n1,n2)
            if abs(score-old_score)<threshold:
                break
            old_score=score
    else:
        alpha=0
        model = NMF(n_components=no_of_pc, init = "nndsvda", random_state=1,beta_loss="kullback-leibler",solver="mu",max_iter=1000,alpha_W=0.0,alpha_H=0.0,l1_ratio=0)
        W = model.fit_transform(v1.T)
        H = model.components_
        spW=W
        spH=np.ones((no_of_pc,v2.shape[0]),dtype=float)
        spW ,spH ,norms=multiplicative_method(spW,spH,v2.T,200)


    entropy_H=''
    entropy_SH=''
    entvalue=[]

    for i in range(no_of_pc):
        value=entropy(H[i,:],base=2)  /  np.log2(len(H[i]))
        entvalue.append(value)

    entvalue=np.array(entvalue)
    index=np.argsort(-entvalue)

    H=H[index]
    spH=spH[index]

    for i in range(no_of_pc):
        entropy_H+=',%0.2f'%(entropy(H[i,:],base=2)  /  np.log2(len(H[i])))
        entropy_SH+=',%0.2f'%(entropy(spH[i,:],base=2) / np.log2(len(spH[i]))  )


    #value1=np.sqrt(np.sum((v1.T-np.matmul(W+V,H))**2))
    #value2=np.sqrt(np.sum((v2.T-np.matmul(spW+spV,spH))**2))

    sc_cosine=find_correlation_bw_genes_and_PC_component_in_singlecell_cosine(H.T,CbyG)
    sc_spearman=find_correlation_bw_genes_and_PC_component_in_singlecell(H.T,CbyG)


    sc_cluster_mean_exp=np.mean(CbyG,axis=0)
    sc_cluster_exp_more_than_threshold=CbyG>cutoff_to_count_exp_cell_population
    sc_cluster_exp_more_than_threshold=np.sum(sc_cluster_exp_more_than_threshold,axis=0)
    sc_cluster_exp_more_than_threshold=sc_cluster_exp_more_than_threshold/CbyG.shape[0]

    transfer_sp_com=spH.T
    transfer_sc_com=[]


    sc_barcode=sct_ad_sc.obs_names.to_numpy()
    sp_barcode=sct_ad_sp.obs_names.to_numpy()
    sc_genenames=full_genes_sc.var_names.to_numpy()


    #maximum norm or infinity norm normalization
    for i in range(transfer_sp_com.shape[1]):
        #transfer_sp_com[:,i]=transfer_sp_com[:,i]/max(abs(transfer_sp_com[i:,]))
        #l2norm=np.linalg.norm(transfer_sp_com[:,i],ord=2)
        l2norm=np.std(transfer_sp_com[:,i])
        #l1norm=np.linalg.norm(transfer_sp_com[:,i],ord=1)
        #transfer_sp_com[:,i]=transfer_sp_com[:,i]/l1norm
        transfer_sp_com[:,i]=transfer_sp_com[:,i]/l2norm

    return transfer_sp_com, transfer_sc_com, sp_barcode,sc_barcode, sc_spearman,sc_cosine,sc_genenames, H, spH,sc_cluster_mean_exp,sc_cluster_exp_more_than_threshold,alpha



def makePCneighboorhoodFeatureMatrix(input):
    """
    Helper function in gene_covariation_analysis to find the weighted neighborhood
    average of cell types from the spatial factors.
    """

    n=len(input.spatialcell_unique_clusterid)
    M=np.zeros((len(input.neighbors),n*input.no_of_pc),dtype=float)

    dist_neighbors=input.neigh_distances
    avgdistArray=0
    for i in range(len(dist_neighbors)):
        avgdistArray=avgdistArray+np.mean(dist_neighbors[i])
    avgdist=avgdistArray/len(dist_neighbors)


    for j in range(len(input.neighbors)):
        CC_barcode_id=input.annotation_spatial_barcode_id[j]
        CC_cluster_id=input.annotation_spatial_cluster_id[j]
        PC_component_of_CC=input.pc_of_sp_clusterid[CC_barcode_id]
        PC_component_of_CC=PC_component_of_CC.reshape((1,input.no_of_pc))
        if j==0:
            target=PC_component_of_CC
        else:
            target=np.vstack((target,PC_component_of_CC))

        neigh_dist=np.array(dist_neighbors[j])
        #weightdist=weightdist/avgdist
        neigh_dist=1/neigh_dist
        sum_weight_dist=np.sum(neigh_dist)
        weighted_avg_dist=neigh_dist/sum_weight_dist
        temp={}
        for k in range(len(input.neighbors[j])):
            id=input.neighbors[j][k]
            NC_barcode_id=input.annotation_spatial_barcode_id[id]
            NC_cluster_id=input.annotation_spatial_cluster_id[id]
            PC_component_of_NC=input.pc_of_sp_clusterid[NC_barcode_id]
            PC_component_of_NC=PC_component_of_NC.reshape((1,input.no_of_pc))
            factor=weighted_avg_dist[k]
            if NC_cluster_id not in temp:
                temp[NC_cluster_id]=factor*PC_component_of_NC
            else:
                temp[NC_cluster_id]=np.concatenate((temp[NC_cluster_id],factor*PC_component_of_NC))

        for key in input.spatialcell_unique_clusterid:
            start_index=input.no_of_pc*key
            end_index=start_index+input.no_of_pc
            if key in temp:
                M[j,start_index:end_index]=np.sum(temp[key],axis=0)


    #cluster=input.annotation_spatial_cluster_id
    #cluster=cluster.reshape((len(cluster),1))
    #df=pd.DataFrame(np.hstack((cluster,M)))
    data=np.hstack((target,M))
    #df=pd.DataFrame(np.hstack((target,M)))
    #df.to_csv(input.outputname,index=True,header=None)
    np.savez(input.outputname,weighted_neighborhood_of_factors_in_niche=data)





def compute_PC_space(input,sct_ad_sc_full):
    """
    Helper function used gene_covariation_analysis. The spatial clusters' names must match the scRNA-seq clusters' names. If it is not, then it will throw an error.
    """
    a=set(input.singlecell_unique_clustername)
    b=set(input.spatialcell_unique_clustername)
    common=a.intersection(b)


    print("\n\n Spatial and scRNA-seq number of clusters, respectively ",len(b),len(a))
    print('Common cell types between spatial and scRNA-seq data  ',len(common),common)
    print('\nThe spatial cluster name does not match the scRNA-seq cluster name ', b-common)
    print("If the above answer is Null, then everything is okay. However, If any spatial cell type does not exist in scRNA-seq data, please correct this manually; otherwise, NiCo will not run. ")
    print("\n\n")

    flag=1
    if len(b-common)>0:
        flag=0

    if flag==1:
        n=len(input.spatialcell_unique_clustername)
        pc_of_sp_clusterid={}
        save_scFactors={}
        save_spFactors={}
        PCA_of_sc_cluster_accordingto_spatial_clusterid={}
        for i in range(n):
            clidsp=input.spatialcell_unique_clusterid[i]
            index=np.where(input.annotation_spatial_cluster_id==clidsp)
            spbarcode=input.annotation_spatial_barcode_id[index[0]]
            scbarcode=[]
            for j in range(len(input.singlecell_unique_clustername)):
                if input.singlecell_unique_clustername[j]==input.spatialcell_unique_clustername[i]:
                    clid=input.singlecell_unique_clusterid[j]
                    index=np.where(input.annotation_singlecell_cluster_id==clid)
                    scbarcode=input.annotation_singlecell_barcode_id[index[0]]
                    break

            pc_sp,pc_sc,sp_barcode,sc_barcode,sc_spearman,sc_cosine,sc_genenames,H, spH,sc_cluster_mean_exp,sc_cluster_exp_more_than_threshold,alpha=find_PC_of_invidualCluster_in_SC(scbarcode,input.iNMFmode,input.ad_sc,input.no_of_pc,spbarcode,input.ad_sp, sct_ad_sc_full,input.spatialcell_unique_clustername[i],input.cutoff_to_count_exp_cell_population)

            PCA_of_sc_cluster_accordingto_spatial_clusterid[clidsp]=[sc_spearman,pc_sp,sc_genenames,sc_cluster_mean_exp,sc_cluster_exp_more_than_threshold,sc_cosine,alpha]
            for k in range(len(sp_barcode)):
                pc_of_sp_clusterid[sp_barcode[k]]=pc_sp[k]
                save_spFactors[sp_barcode[k]]=spH[:,k]

            for k in range(len(sc_barcode)):
                save_scFactors[sc_barcode[k]]=H[:,k]

    return pc_of_sp_clusterid,PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors



def model_linear_regression(input,logistic_predicted_interactions):
    """
    Helper function for gene_covariation_analysis that prepare data Y(central cell factors) and X (neighborhood avg spatial cell factors) for each celltype to perform regression.
    """
    shap_cluster_cutoff=input.shap_cluster_cutoff
    data1=np.load(input.outputname,allow_pickle=True)
    data1=data1['weighted_neighborhood_of_factors_in_niche']
    #print(data1.shape)
    #print(data1[0:5])
    #data1 = np.genfromtxt(open(input.outputname, "rb"), delimiter=',', skip_header=0)
    #ind=~np.isnan(data1).any(axis=1)
    #data=data1[ind,:]
    data=np.nan_to_num(data1)

    featureVector=range(input.no_of_pc,data.shape[1]) # #just neighborhood
    AllneighborhoodClass= data[:,featureVector]
    Alltarget= data[:,0:input.no_of_pc]

    count=0
    save_coef={}
    for i in range(len(input.spatialcell_unique_clusterid)):
        temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        index=temp[0]
        neighborhoodClass=AllneighborhoodClass[index,:]
        target=Alltarget[index,:]
        positive_interacted_CT= logistic_predicted_interactions[input.spatialcell_unique_clustername[i]]
        newindex=[]
        xlabel=[]
        score=[]
        for j in range(len(input.spatialcell_unique_clustername)):
            start=j*input.no_of_pc
            end=start+input.no_of_pc
            for k in range(len(positive_interacted_CT)):
                if positive_interacted_CT[k][0]==input.spatialcell_unique_clustername[j]:
                    xlabel.append(positive_interacted_CT[k][0])
                    score.append(positive_interacted_CT[k][1])
                    for kk in range(start,end):
                        newindex.append(kk)

        neighborhoodClass=neighborhoodClass[:,newindex]
        xlabel=np.array(xlabel)
        score=np.array(score)

        ylabelname=[]
        for k in range(len(xlabel)):
            for j in range(input.no_of_pc):
                ylabelname.append(xlabel[k]+'_s'+'%0.3f'%score[k]+'_Fa'+str(j+1))

        count+=neighborhoodClass.shape[0]
        saveoutname=str(input.spatialcell_unique_clusterid[i])+'_'+input.spatialcell_unique_clustername[i]
        coef,intercept,alpha,percent_variance_explained,residual_variance_explained,pv=run_ridge_regression(input,saveoutname,ylabelname,target,neighborhoodClass,shap_cluster_cutoff)
        #coef_mu,comp_score,coef_std,comp_score_std,alpha=run_ridge_regression(input.seed ,input.lambda_c,input.K_fold,input.n_repeats,target,neighborhoodClass)
        #savedata=savedir+'coef'+str(input.spatialcell_unique_clusterid[i])+'.npz'

        save_coef[input.spatialcell_unique_clusterid[i]]=[coef,intercept,alpha,xlabel,score,target,neighborhoodClass,pv,percent_variance_explained,residual_variance_explained]
        #np.savez(savedata,coef_mu=coef,intercept=intercept,alpha=alpha,xlabel=xlabel,score=score,Yreg=target,Xreg=neighborhoodClass,pvalue=pv,pve=percent_variance_explained,rve=residual_variance_explained)
        #np.savez(savedata,coef_mu=coef_mu,coef_std=coef_std,comp_score=comp_score,comp_score_std=comp_score_std,alpha=alpha,xlabel=xlabel,score=score)

    #print(count)
    return save_coef





def run_ridge_regression(input,saveoutname,ylabelname,target,neighborhoodClass,shap_cluster_cutoff):

    """
    Helper function for model_linear_regression to perform the ridge regression per cell type.
    """

    train_index=range(target.shape[0])
    test_index=[]

    x_std=np.std(neighborhoodClass,axis=0)
    y_std=np.std(target,axis=0)
    for i in range(neighborhoodClass.shape[1]):
        if x_std[i]==0:
            x_std[i]=1
        neighborhoodClass[:,i]=neighborhoodClass[:,i]/x_std[i]
    for i in range(target.shape[1]):
        if y_std[i]==0:
            y_std[i]=1
        target[:,i]=target[:,i]/y_std[i]

    #add=np.hstack((target,neighborhoodClass))
    #ind1=~np.isnan(add).any(axis=0) #1 means rows and 0 means columns
    #ind2=~np.isnan(add).any(axis=1) #1 means rows and 0 means columns
    #data=data1[ind,:]
    #print(target.shape,neighborhoodClass.shape,len(ind1),len(ind2))


    #print(neighborhoodClass.shape)
    #print(target.shape,x_std.shape,y_std.shape)
    x_train,x_test=neighborhoodClass[train_index],neighborhoodClass[test_index]
    y_train,y_test=target[train_index],target[test_index]

    #create_directory(savedir+'plot_Y_and_X/')
    if input.shap_analysis:
        dir1=input.regression_outdir+'Shapley_Interventional/'
        dir2=input.regression_outdir+'Shapley_FullConventional/'
        create_directory(dir1)
        create_directory(dir2)

    LRI=[]
    LRC=[]
    yhat=[]
    lambda_c=[]
    Xdata=x_train

    #kf = KFold(10)
    #print(kf)
    for i in range(y_train.shape[1]):
        linear_model = RidgeCV(alphas=input.lambda_c)#,cv=kf,scoring = 'neg_mean_squared_error')
        #pipe=Pipeline([ ('StandardScaler',StandardScaler(with_mean=True)),('ridge_regression',linear_model)])
        pipe=Pipeline([('ridge_regression',linear_model)])
        pipe.fit(Xdata,y_train[:,i])
        yyhat=pipe.predict(Xdata)
        yhat.append(yyhat)
        LR= pipe.named_steps['ridge_regression']
        coef=LR.coef_
        intercept=LR.intercept_
        LRI.append(intercept)
        LRC.append(coef)
        lambda_c.append('%0.2f'%LR.alpha_)

    LRI=np.array(LRI)
    yhat=np.array(yhat).T
    LRC=np.array(LRC)


    #mu=np.mean(y_train,axis=0)
    #total_ss= np.sum((y_train-mu)**2,axis=0)
    #residual_ss=np.sum((y_train-yhat)**2,axis=0)
    #explained_ss= np.sum((yhat-mu)**2,axis=0)
    #percent_variance_explained=100*explained_ss/total_ss
    #residual_variance_explained=100*residual_ss/total_ss

    pv=np.ones(LRC.shape,dtype=float)
    EVS=[]
    rss=[]
    for i in range(y_train.shape[1]):
                #EVS.append(explained_variance_score(save_y_test[:,i], save_y_pred[:,i]))
                EVS.append(explained_variance_score(y_train[:,i], yhat[:,i]))
                rss.append(np.sum((y_train[:,i]-yhat[:,i])**2,axis=0))
                params = np.append(LRI[i],LRC[i,:])
                newX = np.append(np.ones((len(Xdata),1)), Xdata, axis=1)
                MSE = (sum((y_train[:,i]-yhat[:,i])**2))/(len(newX)-len(newX[0]))

                detM=np.linalg.det(np.dot(newX.T,newX))
                if detM>0:
                    flag=0
                    try:
                        var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
                    except np.linalg.LinAlgError as e:
                        if 'Singular matrix' in str(e):
                            var_b=1# your error handling block
                            flag=1
                        else:
                            raise
                    sd_b = np.sqrt(var_b)
                    ts_b = params/ sd_b
                    df = x_train.shape[0] - x_train.shape[1]
                    p_values1 =np.array([[2*(1-scipy.stats.t.cdf(np.abs(ii),df-1)) for ii in ts_b]])
                    pv[i]=p_values1[:,1:]
                    #if flag==1:
                    #    print(i,saveoutname,MSE,"var_b",var_b,"pvalue",pv[i])



        #print("LRC",LRC.shape,LRI.shape)
        #x_train2 = sm.add_constant(Xdata)
        #est1=sm.OLS(y_train[:,0],x_train2).fit()
        #print("summary1",est1.summary())
        #est2=sm.OLS(y_train[:,1],x_train2).fit()
        #print("summary2",est2.summary())
        #est3=sm.OLS(y_train[:,2],x_train2).fit()
        #print("summary3",est3.summary())


    if input.shap_analysis:
            #explainer = shap.LinearExplainer(LR, x_train)
            explainer = shap.explainers.Linear(LR, x_train,feature_names=ylabelname,feature_perturbation="interventional")
            #explainer = shap.Explainer(LR, x_train,feature_names=ylabelname)
            #shap_values = explainer.shap_values(x_train)
            shap_values = explainer(x_train)

            for i in range(y_train.shape[1]):
                #shap.waterfall_plot(explainer.expected_value, shap_values[sample_ind], X.iloc[sample_ind], max_display=14)
                clust = shap.utils.hclust(x_train, y_train[:,i], linkage="single")
                shap.plots.bar(shap_values, clustering=clust, clustering_cutoff=shap_cluster_cutoff, show=False)
                plt.title("True to the model "+saveoutname+'_'+'Fa'+str(i+1)+", EVS = " +'%0.4f'%EVS[i])
                plt.savefig(dir1+saveoutname+'_Fa'+str(i+1)+'.png',dpi=300, bbox_inches = "tight")
                plt.close('all')

                explainer = shap.explainers.Linear(LR, x_train,feature_names=ylabelname,feature_perturbation="correlation_dependent")
                shap_values = explainer(x_train)
                shap.plots.bar(shap_values, clustering=clust, clustering_cutoff=shap_cluster_cutoff, show=False)
                plt.title("True to the data "+saveoutname+'_'+'Fa'+str(i+1)+", EVS = " +'%0.4f'%EVS[i])
                plt.savefig(dir2+saveoutname+'_Fa'+str(i+1)+'.png',dpi=300, bbox_inches = "tight")
                plt.close('all')

    coef=LRC
    intercept=LRI
    residual_variance_explained=0

    return coef,intercept,lambda_c,EVS,residual_variance_explained,pv



def find_logistic_regression_interacting_score(cmn,coef,CTFeatures,nameOfCellType,logistic_coef_cutoff):
    """
    Helper function used in gene_covariation_analysis
    to find niche interaction scores from logistic regression classifier.
    """
    a=np.diag(cmn)
    #b=np.diag(input.cmn_std)
    goodPredictedCellType=np.argsort(-a)
    largest=np.max(abs(coef))
    normalized_coef=coef/largest
    InteractingCTs=[]
    for k in range(len(a)):
            meanCoefficients=normalized_coef[goodPredictedCellType[k]]
            #stdCoefficients=input.coef_std[goodPredictedCellType[k]]
            highestIndex=np.argsort(-abs(meanCoefficients))
            n=len(highestIndex)
            coeff_of_CT=[]
            name_of_the_coeff=[]
            std_of_coeff=[]
            predictedCT=nameOfCellType[goodPredictedCellType[k]]
            positiveprediction=[]
            negativeprediction=[]
            score=[]
            for i in range(n):
                l=CTFeatures[highestIndex[i]].split()
                temp=''
                for j in range(len(l)):
                    temp+=nameOfCellType[int(l[j][1:])]
                    if j!=(len(l)-1):
                        temp+='--'
                if meanCoefficients[ highestIndex[i]]>logistic_coef_cutoff:
                    positiveprediction.append(temp)
                    score.append(meanCoefficients[ highestIndex[i]])
                else:
                    negativeprediction.append(temp)
            InteractingCTs.append([predictedCT,positiveprediction, score   ])

    logistic_predicted_interactions={}
    for i in range(len(InteractingCTs)):
        cCT=InteractingCTs[i][0]
        nCT=InteractingCTs[i][1]
        Interacting_score=InteractingCTs[i][2]
        for j in range(len(nCT)):
            if cCT not in logistic_predicted_interactions:
                logistic_predicted_interactions[cCT]=[[nCT[j],Interacting_score[j]]]
            else:
                logistic_predicted_interactions[cCT].append([nCT[j],Interacting_score[j]])

    return logistic_predicted_interactions


def findXYZC(c,s):
    "Helper function used in plot_top_selected_genes_as_dotplot."
    x=[]
    y=[]
    z=[]
    bigs=[]
    for i in range(len(c)):
        for j in range(len(c[i])):
            x.append(j)
            y.append(i)
            z.append(c[i,j])
            bigs.append(100*s[i,j])
    return x,y,z,bigs



def create_subtitle(fig: plt.Figure, grid: SubplotSpec, title: str):
    """
    Helper function used in plot_top_selected_genes_as_dotplot.
    Sign sets of subplots with title.
    """
    row = fig.add_subplot(grid)
    # the '\n' is important
    row.set_title(f'{title}\n', fontweight='semibold')
    # hide subplot
    row.set_frame_on(False)
    row.axis('off')





def find_fold_change(PCA,NH_PCA,gene,CCPC,NCPC,totalLRpairs,LRcutoff,CC_meanExpression,NC_meanExpression,CC_popExpression,NC_popExpression,number_of_top_genes_to_print):
    "Helper function used to find the ligand-receptor genes for find_LR_interactions_in_interacting_cell_types."
    listofallLR={}
    uniqueLRpairs={}
    for i in range(len(totalLRpairs)):
        l=totalLRpairs[i][0]
        r=totalLRpairs[i][1]
        listofallLR[l]=1
        listofallLR[r]=1
        name=l+'--'+r
        if name not in uniqueLRpairs:
            uniqueLRpairs[name]=1

    first=PCA[:,CCPC]
    second=NH_PCA[:,NCPC]
    ind1=np.argsort(-abs(first))
    ind2=np.argsort(-abs(second))

    cc_genes=[]
    cc_genes2=[]
    cc_genes5=[]

    nc_genes=[]
    nc_genes2=[]
    nc_genes5=[]


    for i in range(number_of_top_genes_to_print):
            cc_genes5.append([gene[ind1[i]],'%0.2f'%first[ind1[i]]])

    for i in range(number_of_top_genes_to_print):
            nc_genes5.append([gene[ ind2[i] ],'%0.2f'%second[ ind2[i] ]])


    for i in range(len(ind1)):
        if (first[ind1[i]])>LRcutoff:
        #if (first[ind1[i]])<-0.4:
            cc_genes.append(gene[ind1[i]])
            if gene[ind1[i]].upper() in listofallLR:
                cc_genes2.append([gene[ind1[i]],'%0.2f'%first[ind1[i]] ,CC_meanExpression[ind1[i]],CC_popExpression[ind1[i]]       ])



    for i in range(len(ind2)):
        if (second[ind2[i]])>LRcutoff:
        #if (second[ind2[i]])<-0.4:
            nc_genes.append(gene[ind2[i]])
            if gene[ind2[i]].upper() in listofallLR:
                nc_genes2.append([gene[ ind2[i] ],'%0.2f'%second[ ind2[i] ], NC_meanExpression[ind2[i]],NC_popExpression[ind2[i]]   ])

    Found1=[]
    Found2=[]
    for i in range(len(cc_genes2)):
        cc=cc_genes2[i][0].upper()
        for j in range(len(nc_genes2)):
            nc=nc_genes2[j][0].upper()
            name1=cc+'--'+nc # lig in CC and rec in NC
            name2=nc+'--'+cc # lig in NC and rec in CC
            if name1 in uniqueLRpairs:
                Found1.append([cc_genes2[i],nc_genes2[j] ])  # lig in CC and rec in NC
            if name2 in uniqueLRpairs:
                Found2.append([nc_genes2[j],cc_genes2[i] ])  # lig in NC and rec in CC

    return cc_genes, nc_genes,cc_genes5,nc_genes5,Found1,Found2


def sorting_of_factors_for_showing_the_value_in_excelsheet(CC_corr,no_of_pc,gene,genenames):
    """
    Helper function for make_excel_sheet_for_gene_correlation to sorting the factors value.
    """

    headersave_full=[]
    headersave_common=[]
    sort_full=[]
    sort_common=[]
    for k in range(no_of_pc):
            sort_full.append([])
            sort_common.append([])
    for j in range(len(CC_corr)):
            ind=~np.isnan(CC_corr[j]).any(axis=0)
            if ind==True:
                #ax.text(CC_corr[j,0],CC_corr[j,1],gene[j],fontsize=5)
                header=[gene[j]]
                for k in range(no_of_pc):
                    sort_full[k].append(CC_corr[j,k]) #without absolute
                    header.append(CC_corr[j,k])
                headersave_full.append(header)
                if gene[j] in genenames:
                    headersave_common.append(header)
                    for k in range(no_of_pc):
                        sort_common[k].append(CC_corr[j,k]) #without absolute
    return headersave_full,headersave_common,sort_full,sort_common




def triangulation_for_triheatmap(M, N):
    """
    Helper function for plotting ligand receptor map in plot_ligand_receptor_in_interacting_celltypes.
    """

    xv, yv = np.meshgrid(np.arange(-0.5, M), np.arange(-0.5, N))  # vertices of the little squares
    xc, yc = np.meshgrid(np.arange(0, M), np.arange(0, N))  # centers of the little squares
    x = np.concatenate([xv.ravel(), xc.ravel()])
    y = np.concatenate([yv.ravel(), yc.ravel()])
    cstart = (M + 1) * (N + 1)  # indices of the centers

    trianglesN = [(i + j * (M + 1), i + 1 + j * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesE = [(i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesS = [(i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesW = [(i + (j + 1) * (M + 1), i + j * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    return [Triangulation(x, y, triangles) for triangles in [trianglesN, trianglesE, trianglesS, trianglesW]]


def  plot_ligand_receptor_in_interacting_celltypes(CC_celltype_name,NC_celltype_name,logRegScore,pc1,pc2,ridgeRegScore,pvalue,Found1,Found2,saveLRplots,LR_plot_Exp_thres,saveas,transparent_mode,showit,figsize,flag):
    """
    Helper function used in find_LR_interactions_in_interacting_cell_types to plot the rectangle p-value figures.

    The Y-axis shows the central cell type factors, and the X-axis shows the neighborhood cell type factors.

    The circle size denotes the p-values. The large circle size is the most significant, and the smaller circle size is the least significant.
    """
    xfact=figsize[0]/34.0
    yfact=figsize[1]/44.0
    LRFigSize=np.zeros(np.array(figsize).shape)
    ligand=[]
    receptor=[]
    fact_lig=[]
    fact_rec=[]
    popExp_lig=[]
    popExp_rec=[]
    A=[]
    B=[]
    if flag=='First':
        for ele in range(len(Found1)):
            ligExpCellPop=Found1[ele][0][3]
            recExpCellPop=Found1[ele][1][3]
            if ((ligExpCellPop>LR_plot_Exp_thres)&(recExpCellPop>LR_plot_Exp_thres)):
                ligand.append(Found1[ele][0][0])
                receptor.append(Found1[ele][1][0])
                fact_lig.append(float(Found1[ele][0][1]))
                fact_rec.append(float(Found1[ele][1][1]))
                popExp_lig.append(Found1[ele][0][3])
                popExp_rec.append(Found1[ele][1][3])
                A.append(CC_celltype_name+'(cc)_Fa'+str(pc1)+'_'+Found1[ele][0][0])
                B.append(NC_celltype_name+'(nc)_Fa'+str(pc2)+'_'+Found1[ele][1][0])

    if flag=='Second':
        for ele in range(len(Found2)):
            ligExpCellPop=Found2[ele][0][3]
            recExpCellPop=Found2[ele][1][3]
            if ((ligExpCellPop>LR_plot_Exp_thres)&(recExpCellPop>LR_plot_Exp_thres)):
                ligand.append(Found2[ele][0][0])
                receptor.append(Found2[ele][1][0])
                fact_lig.append(float(Found2[ele][0][1]))
                fact_rec.append(float(Found2[ele][1][1]))
                popExp_lig.append(Found2[ele][0][3])
                popExp_rec.append(Found2[ele][1][3])
                A.append(NC_celltype_name+'(nc)_Fa'+str(pc2)+'_'+Found2[ele][0][0])
                B.append(CC_celltype_name+'(cc)_Fa'+str(pc1)+'_'+Found2[ele][1][0])
    if flag=='Both':
        for ele in range(len(Found1)):
            #header=['Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)','Receptor(A)','Ligand(B)','GeneCor(Rec)','GeneCor(Lig)']
            #header[11]=Found1[ele][0][2]
            #header[12]=Found1[ele][1][2]
            ligExpCellPop=Found1[ele][0][3]
            recExpCellPop=Found1[ele][1][3]
            if ((ligExpCellPop>LR_plot_Exp_thres)&(recExpCellPop>LR_plot_Exp_thres)):
                ligand.append(Found1[ele][0][0])
                receptor.append(Found1[ele][1][0])
                fact_lig.append(float(Found1[ele][0][1]))
                fact_rec.append(float(Found1[ele][1][1]))
                popExp_lig.append(Found1[ele][0][3])
                popExp_rec.append(Found1[ele][1][3])
                A.append(CC_celltype_name+'(cc)_Fa'+str(pc1)+'_'+Found1[ele][0][0])
                B.append(NC_celltype_name+'(nc)_Fa'+str(pc2)+'_'+Found1[ele][1][0])

        for ele in range(len(Found2)):
            ligExpCellPop=Found2[ele][0][3]
            recExpCellPop=Found2[ele][1][3]
            if ((ligExpCellPop>LR_plot_Exp_thres)&(recExpCellPop>LR_plot_Exp_thres)):
                ligand.append(Found2[ele][0][0])
                receptor.append(Found2[ele][1][0])
                fact_lig.append(float(Found2[ele][0][1]))
                fact_rec.append(float(Found2[ele][1][1]))
                popExp_lig.append(Found2[ele][0][3])
                popExp_rec.append(Found2[ele][1][3])
                A.append(NC_celltype_name+'(nc)_Fa'+str(pc2)+'_'+Found2[ele][0][0])
                B.append(CC_celltype_name+'(cc)_Fa'+str(pc1)+'_'+Found2[ele][1][0])


    #print(CC_celltype_name,NC_celltype_name,flag,logRegScore,len(Found1),len(Found2),len(ligand),len(A),len(B),LRFigSize1) #,len(nA),len(nB)
    if (len(A)>0)&(len(B)>0):
            nA=np.sort(np.unique(A))
            nB=np.sort(np.unique(B))
            fact_lig=np.array(fact_lig)
            fact_rec=np.array(fact_rec)
            popExp_rec=np.array(popExp_rec)
            popExp_lig=np.array(popExp_lig)

            p1=np.max(fact_lig)
            p2=np.max(fact_rec)
            p3=np.max(popExp_rec)
            p4=np.max(popExp_lig)

            q1=np.min(fact_lig)
            q2=np.min(fact_rec)
            q3=np.min(popExp_rec)
            q4=np.min(popExp_lig)

            fmin=min(q1,q2)
            fmax=max(p1,p2)
            pmin=min(q3,q4)
            pmax=max(p3,p4)

            df = pd.DataFrame({'cols': B,
                               'rows': A,
                               'north': fact_lig,
                               'south': fact_rec,
                               'east': popExp_rec,
                               'west':  popExp_lig})

            df['rows'] = pd.Categorical(df['rows'],categories=nA)  # fix an ordering,
            df_piv = df.pivot_table(index='rows', columns='cols')
            M = len(df_piv.columns) // 4
            N = len(df_piv)

            if (len(range(M))<13)|(len(range(N))<13):
                LRFigSize=figsize
            else:
                LRFigSize[0]=len(range(M))*xfact
                LRFigSize[1]=len(range(N))*yfact

            #print('\n\n',flag,range(M),range(N),CC_celltype_name,[xfact,yfact],LRFigSize)

            values = [df_piv[dir] for dir in ['north', 'east', 'south', 'west']]  # these are the 4 column names in df

            triangul = triangulation_for_triheatmap(M, N)
            #cmaps = ['RdYlBu'] * 4
            #cmaps =['cool','copper','cool','copper']
            cmaps =['Blues','copper_r','Blues','copper_r']
            #norms = [plt.Normalize(0, 1) for _ in range(4)]
            norms = [plt.Normalize(fmin, fmax),plt.Normalize(pmin, pmax),plt.Normalize(fmin, fmax),plt.Normalize(pmin, pmax)]

            fig, ax = plt.subplots(figsize=LRFigSize)

            imgs = [ax.tripcolor(t, np.ravel(val), norm=norm,cmap=cmap,ec='white')  #norm=[]
                    for t, val, cmap, norm in zip(triangul, values, cmaps, norms)]

            #ax.tick_params(length=0)
            #ax.set_title('localizationCoef='+'%0.3f'%np.unique(localized)+',regressionCoef='+'%0.3f'%np.unique(regCoff))
            ax.set_xticks(range(M))
            ax.set_xticklabels(df_piv['north'].columns,rotation=90,style='italic')
            ax.set_yticks(range(N))
            ax.set_yticklabels(df_piv.index,style='italic')
            ax.invert_yaxis()
            ax.margins(x=0, y=0)
            pvalue='%0.2f'%(-np.log10(pvalue))
            ax.set_title(CC_celltype_name+'_Fa'+str(pc1)+', '+NC_celltype_name+'_Fa'+str(pc2)+', SS=%0.3f'%logRegScore+', RRS=%0.2f'%ridgeRegScore +', pv='+pvalue  )
            #ax.set_aspect('equal', 'box')  # square cells
            plt.colorbar(imgs[0], ax=ax,label='correlation value in the factors')
            plt.colorbar(imgs[1], ax=ax,label='% expressed cell population')

            #plt.tight_layout()

            savefname=remove_extra_character_from_name(CC_celltype_name)+'_Fa'+str(pc1)+'_'+remove_extra_character_from_name(NC_celltype_name)+'_Fa'+str(pc2)

            fig.savefig(saveLRplots+savefname+'.'+saveas,bbox_inches='tight', transparent=transparent_mode,dpi=300)
            if showit:
                pass
            else:
                plt.close('all')
    return 0



def visualize_factors_in_scRNAseq_umap(input,choose_interacting_celltype_pair,visualize_factors_id,
refpath='./inputRef/',ref_original_counts='Original_counts.h5ad',
msna=0.1,ms=5,cmap=plt.cm.get_cmap('Spectral'),saveas='pdf',transparent_mode=False,showit=True,figsize=(8,3.5)):


    """
    Inputs:

    | The primary input is the output from gene_covariation_analysis.

    | Define the cell type pairs via the list you want to visualize. The first is the central cell type, and the second is the niche cell type.
    | choose_interacting_celltype_pair=['?']

    | Define the factors id via the list you want to visualize in the umap. The first is the factor ID of the central cell type, and the second is the factor ID of the niche cell type.
    | visualize_factors_id=['?']

    | Reference scRNAseq count matrix in scTransform-like normalization in the common gene space. The filename must be sct_singleCell.h5ad
    | (default) refpath='./inputRef/'

    | Original count data of scRNAseq in h5ad format
    | (default) ref_original_counts='Original_counts.h5ad'
    | Must have the umap information in .obsm['X_umap']

    | The marker size of chosen and NA cell types
    | (default) ms=5  (chosen)
    | (default) msna=0.1 (NA)

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Define the colormap for visualizing factors
    | (default) cmap=plt.cm.get_cmap('Spectral')

    | Background color in the figures
    | (default) transparent_mode=False

    | Dimension of the figure size.
    | (default) figsize=(8,3.5)

    Outputs:

    | The factor visualization in scRNAseq figure is saved in "./spatial_ct_ct_interactions/covariations_R*_F*/scRNAseq_factors_in_umap."

    """


    ref_original_counts=refpath+ref_original_counts
    original_h5ad=sc.read_h5ad(ref_original_counts)
    umap=original_h5ad.obsm['X_umap']
    barcode=original_h5ad.obs_names.to_numpy()
    barcode=np.reshape(barcode,(len(barcode),1))
    umap_not_order=np.hstack((barcode,umap))

    cellname=np.reshape(input.annotation_singlecell_barcode_id,(len(input.annotation_singlecell_barcode_id),1) )
    umap_data=sort_index_in_right_order(cellname,umap_not_order)

    #sc_ct_name=np.array([input.singlecell_unique_clusterid,input.singlecell_unique_clustername]))
    #sc_cluster=np.array([input.annotation_singlecell_barcode_id,input.annotation_singlecell_cluster_id])

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))


    t=[]
    for i in range(input.no_of_pc):
        t.append(0)
    H_sc=[]
    for i in range(len(input.annotation_singlecell_barcode_id)):
        cellid=input.annotation_singlecell_barcode_id[i]
        if cellid in save_scFactors:
            H_sc.append(save_scFactors[cellid])
        else:
            H_sc.append(np.array(t))
    H_sc=np.array(H_sc)


    fig,(ax)=plt.subplots(1,2,figsize=figsize)
    v1=H_sc[:,(visualize_factors_id[0]-1)]
    v2=H_sc[:,(visualize_factors_id[1]-1)]

    CTname=input.singlecell_unique_clustername
    cellsinCT={}
    for i in range(len(input.annotation_singlecell_barcode_id)):
        clu_id=input.annotation_singlecell_celltypename[i]
        if clu_id not in cellsinCT:
            cellsinCT[clu_id]=[i]
        else:
            cellsinCT[clu_id].append(i)


    plot_all_ct(CTname,umap_data[:,[1,2]],cellsinCT,ax[0],[choose_interacting_celltype_pair[0]],v1,cmap,ms,msna)
    plot_all_ct(CTname,umap_data[:,[1,2]],cellsinCT,ax[1],[choose_interacting_celltype_pair[1]],v2,cmap,ms,msna)
    #ax1[1].legend(loc='upper right',bbox_to_anchor=(1.50, 1),ncol=1, frameon=False,borderaxespad=0.,prop={"size":10},fancybox=True, shadow=True)

    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    #plt.gca().axes.get_yaxis().set_visible(False)
    ax[0].set_axis_off()
    ax[1].set_axis_off()
    ax[0].set_title(choose_interacting_celltype_pair[0]+':Fa'+str(visualize_factors_id[0]))
    ax[1].set_title(choose_interacting_celltype_pair[1]+':Fa'+str(visualize_factors_id[1]))

    fig.tight_layout()
    print("The figures are saved: ", input.covariation_dir+'scRNAseq_factors_in_umap.'+saveas)
    fig.savefig(input.covariation_dir+'scRNAseq_factors_in_umap.'+saveas,bbox_inches='tight',transparent=False,dpi=300)
    if showit:
        pass
    else:
        plt.close('all')
    return 0


def plot_all_ct(CTname,PP,cellsinCT,ax,mycelltype,Fa,cmap,ms,msna):
    "Helper function used in visualizing factors value in the umap."
    #cmap=plt.cm.get_cmap('Spectral')
    #cmap=plt.cm.get_cmap('jet')
    cumsum=np.linspace(0,1,len(CTname))

    naindex=[]
    for i in range(len(CTname)):
        index=cellsinCT[CTname[i]]
        labelname=str(i)+'-'+CTname[i]+'-'+str(len(index))
        rgba=cmap(cumsum[i])
        if CTname[i] in mycelltype:
            #p1=ax.plot(PP[index,0],PP[index,1],'o',label=labelname,color=rgba,markersize=1)
            p1=ax.scatter(PP[index,0],PP[index,1],s=ms,c=Fa[index],marker='o',cmap=cmap)
        else:
            naindex=naindex+index

    ax.plot(PP[naindex,0],PP[naindex,1],'.',color="0.5",label='NA',markersize=msna)
    plt.colorbar(p1,ax=ax,shrink=0.5)




def visualize_factors_in_spatial_umap(input,choose_interacting_celltype_pair,visualize_factors_id,
quepath='./inputQuery/',msna=0.1,ms=5,cmap=plt.cm.get_cmap('Spectral'),saveas='pdf',transparent_mode=False,showit=True,figsize=(8,3.5)):


    """
    Inputs:

    | The primary input is the output from gene_covariation_analysis.

    | Define the cell type pairs via the list you want to visualize. The first is the central cell type, and the second is the niche cell type.
    | choose_interacting_celltype_pair=['?']

    | Define the factors id via the list you want to visualize in the umap. The first is the factor ID of the central cell type, and the second is the factor ID of the niche cell type.
    | visualize_factors_id=['?']

    | Queried spatial count matrix in scTransform-like normalization in the common gene space. The filename should be sct_spatial.h5ad
    | (default) quepath='./inputQuery/'
    | Must have the umap information in .obsm['X_umap']

    | The marker size of chosen and NA cell types
    | (default) ms=5  (chosen)
    | (default) msna=0.1 (NA)

    | Define the colormap for visualizing factors
    | (default) cmap=plt.cm.get_cmap('Spectral')

    | Save the figures in PDF or PNG format (dpi for PNG format is 300)
    | (default) saveas='pdf'

    | Background color in the figures
    | (default) transparent_mode=False

    | Dimension of the figure size.
    | (default) figsize=(8,3.5)

    Output:

    The output figure will be saved in spatial_ct_ct_interactions/covariations_R*_F*/spatial_factors_in_umap*

    """

    que_h5ad=quepath+'sct_spatial.h5ad'
    sct_ad_sp=sc.read_h5ad(que_h5ad)
    umap=sct_ad_sp.obsm['X_umap']
    barcode=sct_ad_sp.obs_names.to_numpy()
    barcode=np.reshape(barcode,(len(barcode),1))
    umap_not_order=np.hstack((barcode,umap))
    cellname=np.reshape(input.annotation_spatial_barcode_id,(len(input.annotation_spatial_barcode_id),1) )
    umap_data=sort_index_in_right_order(cellname,umap_not_order)

    PCA_of_sc_cluster_accordingto_spatial_clusterid,save_scFactors,save_spFactors=pickle.load(open(input.covariation_dir+'factors_info.p', 'rb'))

    t=[]
    for i in range(input.no_of_pc):
        t.append(0)

    H_sp=[]
    for i in range(len(input.annotation_spatial_cluster_id)):
        cellid=input.annotation_spatial_barcode_id[i]
        if cellid in save_spFactors:
            H_sp.append(save_spFactors[cellid])
            #print(save_scFactors[cellid])
        else:
            H_sp.append(np.array(t))
    H_sp=np.array(H_sp)



    fig,(ax)=plt.subplots(1,2,figsize=figsize)
    v1=H_sp[:,(visualize_factors_id[0]-1)]
    v2=H_sp[:,(visualize_factors_id[1]-1)]

    CTname=input.spatialcell_unique_clustername
    cellsinCT={}
    for i in range(len(input.annotation_spatial_cluster_id)):
        clu_id=input.annotation_spatial_celltypename[i]
        #cel_id=sc_cluster[i][0]
        if clu_id not in cellsinCT:
            cellsinCT[clu_id]=[i]
        else:
            cellsinCT[clu_id].append(i)


    plot_all_ct(CTname,umap_data[:,[1,2]],cellsinCT,ax[0],[choose_interacting_celltype_pair[0]],v1,cmap,ms,msna)
    plot_all_ct(CTname,umap_data[:,[1,2]],cellsinCT,ax[1],[choose_interacting_celltype_pair[1]],v2,cmap,ms,msna)
    #ax1[1].legend(loc='upper right',bbox_to_anchor=(1.50, 1),ncol=1, frameon=False,borderaxespad=0.,prop={"size":10},fancybox=True, shadow=True)

    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    #plt.gca().axes.get_yaxis().set_visible(False)
    ax[0].set_axis_off()
    ax[1].set_axis_off()
    ax[0].set_title(choose_interacting_celltype_pair[0]+':Fa'+str(visualize_factors_id[0]))
    ax[1].set_title(choose_interacting_celltype_pair[1]+':Fa'+str(visualize_factors_id[1]))

    fig.tight_layout()
    print("The figures are saved: ", input.covariation_dir+'spatial_factors_in_umap.'+saveas)
    fig.savefig(input.covariation_dir+'spatial_factors_in_umap.'+saveas,bbox_inches='tight',transparent=False,dpi=300)
    if showit:
        pass
    else:
        plt.close('all')
    return 0


def read_LigRecDb(contdb):
    """
    Read the ligand-receptor database file. The first column is the ligand pair, and the second is the receptor pair.
    """
    #f=open('sort_3_db_L_R_high_confident.dat','r')
    totalLRpairs=[]
    ligand={}
    receptor={}
    either={}
    for j in range(len(contdb)):
        l=contdb[j][0:-1].split()
        ligand[l[0].upper()]=1
        receptor[l[1].upper()]=1
        if [l[0], l[1] ] not in totalLRpairs:
            totalLRpairs.append( [l[0].upper(), l[1].upper() ])

    for key in ligand:
        if key in receptor:
            either[key]=1
    for key in either:
        ligand.pop(key, None)
        receptor.pop(key, None)

    return totalLRpairs,ligand,receptor,either

def sort_index_in_right_order(correct,wrong):
    "Helper function used in visualize cell type annotations."
    d={}
    for i in range(len(wrong)):
        d[wrong[i,0]]=i
    index=[]
    for i in range(len(correct)):
        index.append(d[correct[i,0]])
    right=wrong[index]
    return right