import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

results = "../standard/clusters"
resultsLT = "../standard/clustersLT"
resultsMM = "../standard/clustersLFC"
resultsCDF = "../standard/clustersCDF"
genes_csv = "gene_list.csv"
genes_filtered = "gene_list_filtered.csv"

genes = pd.read_csv(genes_filtered, header=None).squeeze().tolist()

def extract_coeff(folder) -> pd.DataFrame:
    cluster_dfs = {}
    
    for filename in sorted(os.listdir(folder)):
        path = os.path.join(folder, filename)
        rows = []
        
        with open(path) as file:
            for line in file:
                if line.startswith("row-"):
                    values = [float(x) for x in line[4:].strip().split(",")[1:4]]
                    error = line[4:].strip().split(",")[4]
                    values = [(x - float(error)) for x in values]
                    rows.append(values)
        
        num_patterns = len(rows[0]) if rows else 0
        pattern_cols = [f"pattern{i+1}" for i in range(num_patterns)]
        
        cluster_name = os.path.splitext(filename)[0].split('_')[0]
        cluster_dfs[cluster_name] = pd.DataFrame(rows, columns=pattern_cols, index=genes)
        cluster_dfs[cluster_name].index.name = "gene"
    
        keep = pd.read_csv(genes_filtered, header=None)[0].tolist()
        cluster_dfs[cluster_name] = cluster_dfs[cluster_name].loc[cluster_dfs[cluster_name].index.intersection(keep)]

    return cluster_dfs

#cluster_dfs = extract_coeff(results)
#clusterLT_dfs = extract_coeff(resultsLT)
#clusterMM_dfs = extract_coeff(resultsMM)
clusterCDF_dfs = extract_coeff(resultsCDF)

def extract_patterns(folder):
    patterns_df = {}
    
    for filename in sorted(os.listdir(folder)):
        path = os.path.join(folder, filename)
        cluster_name = os.path.splitext(filename)[0].split('_')[0]
        
        if cluster_name not in patterns_df:
            patterns_df[cluster_name] = {}
        
        with open(path) as file:
            for line in file:
                if line.startswith('#"'):
                    parts = line.strip().split(',')
                    parts = parts[1].split()
                    vals = [float(x) for x in parts]
                    pattern_key = f"pattern{len(patterns_df[cluster_name]) + 1}"
                    patterns_df[cluster_name][pattern_key] = vals
    
    return patterns_df

#patterns = extract_patterns(results)
#patternsLT = extract_patterns(resultsLT)
#patternsMM = extract_patterns(resultsMM)
patternsCDF = extract_patterns(resultsCDF)

def top_genes_per_pattern(cluster_df, n=20):
    top_ranked_genes = {}

    for cluster_name, df in cluster_df.items():
        top_ranked_genes[cluster_name] = {}
        for pattern in df.columns:
            top_ranked_genes[cluster_name][pattern] = (df[pattern].sort_values(ascending = False).head(n))
    return top_ranked_genes


#top_genes = top_genes_per_pattern(cluster_dfs)
#top_genesLT = top_genes_per_pattern(clusterLT_dfs, 50)
#top_genesMM = top_genes_per_pattern(clusterMM_dfs, 100)
top_genesCDF = top_genes_per_pattern(clusterCDF_dfs, 50)
#print(top_genesLT)
# for cluster in top_genesLT:
#     for pattern in top_genesLT[cluster]:
#         if pattern == "pattern1":
#             print(f"--- Cluster: {cluster} | {pattern} ---")
#             for gene_name in top_genesLT[cluster][pattern].index:
#                 print(gene_name)
#print(top_genesLT)
#print(top_genes["clusterOne"])

#plot patterns identified in each cluster

constraints = {
    "clusterOne": [1.45,12.0,5.21,6.51],
    "clusterTwo": [0.371,0.803,1.22,1.78],
    "clusterThree": [0.0405,0,0.0369,0.0353],
    "clusterFour": [0.00913,0.638,0.0178,0.0243]
}

constraintsLT = {
    "clusterOneLT": [.894475,2.568650,1.826,2.016452],
    "clusterTwoLT": [0.036407,0.589570,0.798347,1.023510],
    "clusterThreeLT": [0.039660,0.0,0.036266,0.034714],
    "clusterFourLT": [0.009090,0.493705,0.017601,0.024007]
}

constraintsLFC = {
    "clusterOneLFC": [0.5000000, 0.8048910, 0.6845230, 0.7166605],
    "clusterTwoLFC": [0.5000000, 0.6113981, 0.6717390, 0.7262386],
    "clusterThreeLFC": [0.5000000, 0, 0.4865699, 0.4801746],
    "clusterFourLFC": [0.500000, 1.112680, 0.596319, 0.641227]
}

constraintsCDF = {
    "clusterOneCDF": [0.5024788,0.9997918,0.9837637,0.9939621],
    "clusterTwoCDF": [0.4972853,0.7214298,0.8182699,0.8848034],
    "clusterThreeCDF": [0.5087043,0.0000000,0.4846397,0.4731947],
    "clusterFourCDF": [0.5024788, 1.0000000, 0.8684869, 0.9493119]
}

def plot_patterns(pattern_df, lt=False):
    x = [0, 3, 6, 9]
    constraint = {}
    constraint = constraintsLT if lt else constraints
    for cluster in pattern_df:
        print(pattern_df[cluster])
        for pattern in pattern_df[cluster]:
            print(pattern_df[cluster][pattern])
            plt.plot(x, pattern_df[cluster][pattern], label = f"{pattern}")
        print('\n')
        plt.plot(x, constraint[cluster], color="crimson", linewidth=2.5, linestyle="--", label=f"{cluster} constraint")
            
        plt.title(f"{cluster} patterns")
        plt.xlabel("timepoint")
        plt.ylabel("coefficients")
        plt.legend()
        plt.xticks(x)
        plt.show()

#plot_patterns(patterns)
#plot_patterns(patternsLT, True)


def plot_genes_per_pattern(patterns, constraints, top_genes, suffix="", n=100):
    x = [0, 3, 6, 9]
    gene_folder = f"gene_counts{suffix}/"

    for cluster_name, pattern_dict in patterns.items():
        print(cluster_name)
        gene_file = os.path.join(gene_folder, f"{cluster_name}_annotated.csv")
        gene_counts = pd.read_csv(gene_file, index_col=0)
        gene_counts.columns = x
        

        for pattern_name, pattern_vals in pattern_dict.items():
            if pattern_name == "pattern1":
                top = top_genes[cluster_name][pattern_name].head(n)
                print(top)
                top_gene_names = top.index.tolist()

                fig, ax = plt.subplots(figsize=(10, 6))

                # for gene in top_gene_names:
                #     if gene in gene_counts.index:
                #         ax.plot(x, gene_counts.loc[gene], color="steelblue",
                #                 alpha=0.4, linewidth=1, label=gene)

                ax.plot(x, pattern_vals, color="blue", linewidth=2.5, linestyle="--", label=f"{pattern_name}")
                ax.plot(x, constraints[cluster_name], color="crimson", linewidth=2.5, linestyle="--", label="constraint")
                
                ax.set_title(f"{cluster_name} — {pattern_name} | Top {n} genes")
                ax.set_xlabel("Timepoint")
                ax.set_ylabel("Gene counts")
                ax.set_xticks(x)

                # handles, labels = ax.get_legend_handles_labels()
                # pattern_handle = [(h, l) for h, l in zip(handles, labels) if "scaled" in l]
                # gene_handles = [(h, l) for h, l in zip(handles, labels) if "scaled" not in l]

                # ax.legend(
                #     [ph[0] for ph in pattern_handle] + [gene_handles[0][0]],
                #     [ph[1] for ph in pattern_handle] + [f"Top {n} genes (n={len(gene_handles)})"],
                #     loc="upper right"
                # )

                plt.tight_layout()
                plt.show()

#plot_genes_per_pattern(patterns, top_genes)
#plot_genes_per_pattern(patternsLT, top_genesLT, lt=True)
#plot_genes_per_pattern(patternsMM, constraintsLFC, top_genesMM, "MM")
# plot_genes_per_pattern(patternsCDF, constraintsCDF, top_genesCDF, "CDF")


# pull from gsea results to get lowest p-val -> identify genes in pathway -> plot average expression
def get_gmt_dict(gmt_file):
    """Parses GMT file into {Pathway: [Gene1, Gene2...]}"""
    pathway_dict = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 2:
                pathway_dict[parts[0]] = parts[2:]
    return pathway_dict

def genes_in_pathway(pathway_name, gmt_dict, expression_df):
    """Finds genes in the pathway that are also in your expression data rows"""
    pathway_list = gmt_dict.get(pathway_name, [])
    # Find intersection of pathway genes and dataset index
    overlap = [gene for gene in pathway_list if gene in expression_df.index]
    return overlap

def top_pathways(file, gmt_file, n_top=5):
    data = pd.read_csv(file, skiprows=1, header=None, sep=' ')
    
    # Column 0: Pathway Name, Column 2: P-value
    pathway_col = data.columns[1]
    pval_col = data.columns[2]
    
    # Calculate gene counts from the GMT
    pathway_genes = get_gmt_dict(gmt_file)
    data['gene_count'] = data[pathway_col].apply(lambda x: len(pathway_genes.get(x, [])))
    
    # Sort by p-value first (lowest is best), then gene count (highest is best)
    sorted_data = data.sort_values(
        by=[pval_col,], 
        ascending=[True, False]
    )
    
    return sorted_data.head(n_top)
        



def plot_pathways(top_df, gmt_dict, expression_df, deg_df=None, pathway_subset=None, plot_mode='expression'):
    # plot_mode: 'expression' for raw average, 'fold_change' for absolute fold change
    
    # 1. Clean the Expression Data Index
    if expression_df.index.dtype == 'int64':
        expression_df = expression_df.set_index(expression_df.columns[0])
    
    expression_df.index = expression_df.index.str.replace('"', '').str.strip()
    
    pathway_col = top_df.columns[1]
    plt.figure(figsize=(12, 7))
    plot_count = 0

    if deg_df:
        with open(deg_df, 'r') as f:
            deg_set = {line.strip().replace('"', '') for line in f if line.strip()}
        print(len(deg_set))
    x = [0, 3, 6, 9]

    def get_plot_values(pathway_means):
        if plot_mode == 'fold_change':
            baseline = pathway_means.iloc[0]
            if baseline == 0:
                baseline = 1e-9  # avoid division by zero
            return pathway_means.apply(lambda v: abs(np.log2((v + 1e-9) / baseline)))
        return pathway_means

    if pathway_subset:
        top_df = top_df[top_df[pathway_col].isin(pathway_subset)]
        for pathway in top_df[pathway_col]:
            raw_gmt_genes = gmt_dict.get(pathway, [])
            
            if deg_df:
                genes = [g for g in raw_gmt_genes if g in deg_set]
            else:
                genes = [g for g in raw_gmt_genes if g in expression_df.index.tolist()]
            
            print(f"Pathway: {pathway[:20]}... | Overlap: {len(genes)}")

            if genes:
                print(genes)
                pathway_means = expression_df.loc[genes].mean(axis=0)
                plot_values = get_plot_values(pathway_means)
                plt.plot(x, plot_values.values, label=pathway, marker='o', alpha=0.8)
                plot_count += 1
    else:
        for pathway in top_df[pathway_col]:
            raw_gmt_genes = gmt_dict.get(pathway, [])
            
            if deg_df:
                genes = [g for g in raw_gmt_genes if g in deg_set]
            else:
                genes = [g for g in raw_gmt_genes if g in expression_df.index.tolist()]
            
            print(f"Pathway: {pathway[:20]}... | Overlap: {len(genes)}")

            if genes:
                print(genes)
                pathway_means = expression_df.loc[genes].mean(axis=0)
                plot_values = get_plot_values(pathway_means)
                plt.plot(x, plot_values.values, label=pathway, marker='o', alpha=0.8)
                plot_count += 1
    
    if plot_count > 0:
        ylabel = "Absolute Log Fold Change" if plot_mode == 'fold_change' else "Average Expression"
        title = "Absolute Log Fold Change for Each Top Pathway" if plot_mode == 'fold_change' else "Averaged Expressions for Each Top Pathway"
        plt.title(title)
        plt.ylabel(ylabel)
        plt.xlabel("timepoint")
        plt.xticks(x, rotation=45)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='x-small')
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.show()
    else:
        print("Still no overlap.")

def plot_tcr_fold_change(tcr_values, label="TCR Pattern", x=[0, 3, 6, 9]):
    baseline = tcr_values[0] if tcr_values[0] != 0 else 1e-9
    fold_change = [abs(np.log2((v + 1e-9) / baseline)) for v in tcr_values]
    
    plt.figure(figsize=(10, 5))
    plt.plot(x, fold_change, marker='o', alpha=0.8, label=label)
    plt.title(f"Absolute Log Fold Change - {label}")
    plt.ylabel("Absolute Log Fold Change")
    plt.xlabel("Timepoint")
    plt.xticks(x)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def plot_selected_vs_unselected(top_df, gmt_dict, expression_df, selected_pathways, deg_df=None, plot_mode='expression', x=[0, 3, 6, 9]):
    pathway_col = top_df.columns[1]
    
    if expression_df.index.dtype == 'int64':
        expression_df = expression_df.set_index(expression_df.columns[0])
    expression_df.index = expression_df.index.str.replace('"', '').str.strip()
    
    if deg_df:
        with open(deg_df, 'r') as f:
            deg_set = {line.strip().replace('"', '') for line in f if line.strip()}

    def get_plot_values(pathway_means):
        if plot_mode == 'fold_change':
            baseline = pathway_means.iloc[0] if pathway_means.iloc[0] != 0 else 1e-9
            return pathway_means.apply(lambda v: abs(np.log2((v + 1e-9) / baseline)))
        return pathway_means

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)
    
    for pathway in top_df[pathway_col]:
        raw_gmt_genes = gmt_dict.get(pathway, [])
        genes = [g for g in raw_gmt_genes if g in (deg_set if deg_df else expression_df.index.tolist())]
        if not genes:
            continue
        
        pathway_means = expression_df.loc[genes].mean(axis=0)
        plot_values = get_plot_values(pathway_means)
        
        ax = axes[0] if pathway in selected_pathways else axes[1]
        ax.plot(x, plot_values.values, marker='o', alpha=0.8, label=pathway)

    axes[0].set_title(f"Selected Pathways (n={len(selected_pathways)})")
    axes[1].set_title("Unselected Pathways")
    for ax in axes:
        ax.set_xlabel("Timepoint")
        ax.set_xticks(x)
        ax.grid(True, linestyle='--', alpha=0.5)
    axes[0].set_ylabel("Absolute Log2 Fold Change" if plot_mode == 'fold_change' else "Average Expression")
    axes[0].legend(fontsize='x-small', loc='upper left')
    axes[1].legend(fontsize='x-small', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.show()

def plot_raw_and_afc(top_df, gmt_dict, expression_df, selected_pathways=None, deg_df=None, deg_only=False, x=[0, 3, 6, 9]):
    pathway_col = top_df.columns[1]
    
    if expression_df.index.dtype == 'int64':
        expression_df = expression_df.set_index(expression_df.columns[0])
    expression_df.index = expression_df.index.str.replace('"', '').str.strip()
    
    deg_set = None
    if deg_df:
        with open(deg_df, 'r') as f:
            deg_set = {line.strip().replace('"', '') for line in f if line.strip()}

    def get_genes(raw_gmt_genes):
        if deg_only and deg_set:
            return [g for g in raw_gmt_genes if g in deg_set]
        return [g for g in raw_gmt_genes if g in expression_df.index.tolist()]

    def get_afc(pathway_means):
        baseline = pathway_means.iloc[0] if pathway_means.iloc[0] != 0 else 1e-9
        return pathway_means.apply(lambda v: abs(np.log2((v + 1e-9) / baseline)))

    fig, axes = plt.subplots(1, 2, figsize=(20, 7))
    axes[0].set_title("Raw Average Expression")
    axes[1].set_title("Absolute Log2 Fold Change")

    top_df_filtered = top_df[top_df[pathway_col].isin(selected_pathways)] if selected_pathways else top_df  # only change
    
    colors = plt.cm.tab10.colors
    handles, labels = [], []

    for i, pathway in enumerate(top_df_filtered[pathway_col]):
        raw_gmt_genes = gmt_dict.get(pathway, [])
        genes = get_genes(raw_gmt_genes)
        
        if not genes:
            print(f"No genes found for {pathway}")
            continue

        pathway_means = expression_df.loc[genes].mean(axis=0)
        afc_values = get_afc(pathway_means)
        color = colors[i % len(colors)]
        short_label = pathway[:40]

        line, = axes[0].plot(x, pathway_means.values, marker='o', alpha=0.8, color=color, label=short_label)
        axes[1].plot(x, afc_values.values, marker='o', alpha=0.8, color=color, label=short_label)
        
        handles.append(line)
        labels.append(short_label)

    for ax in axes:
        ax.set_xlabel("Timepoint")
        ax.set_xticks(x)
        ax.grid(True, linestyle='--', alpha=0.5)
    
    axes[0].set_ylabel("Average Expression")
    axes[1].set_ylabel("Absolute Log2 Fold Change")

    deg_label = "DE Genes Only" if deg_only else "All Genes"
    fig.suptitle(f"Raw vs AFC — {deg_label}", fontsize=13, fontweight='bold')
    axes[1].legend(handles, labels, loc='upper left', ncol=1, fontsize='x-small', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)  # shrink plot area, more room for legend
    plt.show()

c2 = get_gmt_dict("../gsea/msigdb/c2.all.v2026.1.Hs.symbols.gmt")
c5 = get_gmt_dict("../gsea/msigdb/c5.all.v2026.1.Hs.symbols.gmt")
c7 = get_gmt_dict("../gsea/msigdb/c7.all.v2026.1.Hs.symbols.gmt")
hallmark = get_gmt_dict("../gsea/msigdb/h.all.v2026.1.Hs.symbols.gmt")

tcr1 = [1.45, 12.0, 5.21, 6.51]
tcr2 = [0.0371,0.803,1.22,1.78]
tcr3 = [0.0405, 0, 0.0369, 0.0353]
tcr4 = [0.00913, 0.638, 0.0178, 0.0243]

# plot_tcr_fold_change(tcr1, label="Cluster 1 Constraint")
# plot_tcr_fold_change(tcr2, label="Cluster 2 Constraint")
# plot_tcr_fold_change(tcr3, label="Cluster 3 Constraint")
# plot_tcr_fold_change(tcr4, label="Cluster 4 Constraint")


#<---- CLUSTER FOUR CDF ETHAN ---->

clusterFour = pd.read_csv("gene_counts/clusterFour_annotated.csv")

top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterFourDE_GSEA_C2.csv", header=None, sep='\s+').head(10)
c2_arr = top_c2[top_c2.columns[1]].tolist()
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterFourDE_GSEA_C5.csv", header=None, sep='\s+').head(10)
c5_arr = top_c5[top_c5.columns[1]].tolist()
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterFourDE_GSEA_C7.csv", header=None, sep='\s+').head(10)
c7_arr = top_c7[top_c7.columns[1]].tolist()
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterFourDE_GSEA_Hallmark.csv", header=None, sep='\s+').head(10)
hallmark_arr = top_hallmark[top_hallmark.columns[1]].tolist()

#<-- Differentially Expressed Genes -->
# plot_pathways(top_c2, c2, clusterFour, "DE_genes/clusterFourDE3tp.csv")
# plot_pathways(top_c5, c5, clusterFour, "DE_genes/clusterFourDE3tp.csv")
# plot_pathways(top_c7, c7, clusterFour, "DE_genes/clusterFourDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterFour, "DE_genes/clusterFourDE3tp.csv")

#<-- All Genes -->
# plot_pathways(top_c2, c2, clusterFour)
# plot_pathways(top_c5, c5, clusterFour)
# plot_pathways(top_c7, c7, clusterFour)
# plot_pathways(top_hallmark, hallmark, clusterFour)

#<-- ALL GENES ABS FOLD CHANGE -->
# plot_pathways(top_c2, c2, clusterFour, plot_mode="fold_change")
# plot_pathways(top_c5, c5, clusterFour, plot_mode="fold_change")
# plot_pathways(top_c7, c7, clusterFour, plot_mode="fold_change")
# plot_pathways(top_hallmark, hallmark, clusterFour, plot_mode="fold_change")

#<-- ALL GENES RAW + ABS FOLD CHANGE -->
# plot_raw_and_afc(top_c2, c2, clusterFour)
# plot_raw_and_afc(top_c5, c5, clusterFour)
# plot_raw_and_afc(top_c7, c7, clusterFour)
# plot_raw_and_afc(top_hallmark, hallmark, clusterFour)

#<-- DE GENES RAW + DE GENES ABS FC -->
# plot_raw_and_afc(top_c2, c2, clusterFour, deg_df="DE_genes/clusterFourDE3tp.csv", deg_only=True)
# plot_raw_and_afc(top_c5, c5, clusterFour, deg_df="DE_genes/clusterFourDE3tp.csv", deg_only=True)
# plot_raw_and_afc(top_c7, c7, clusterFour, deg_df="DE_genes/clusterFourDE3tp.csv", deg_only=True)
# plot_raw_and_afc(top_hallmark, hallmark, clusterFour, deg_df="DE_genes/clusterFourDE3tp.csv", deg_only=True)

#<-- Selected Pathways vs Non Selected -->
# selected = ["GOBP_ADENYLATE_CYCLASE_MODULATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY", "GOMF_MIRNA_BINDING", "GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION", "GOMF_REGULATORY_RNA_BINDING", "GOBP_MIRNA_CATABOLIC_PROCESS"]
# plot_selected_vs_unselected(top_c5, c5, clusterFour, selected, plot_mode='fold_change')

# selected = ["GSE36476_CTRL_VS_TSST_ACT_72H_MEMORY_CD4_TCELL_YOUNG_UP", "GSE17721_CPG_VS_GARDIQUIMOD_0.5H_BMDC_UP", "GSE37532_TREG_VS_TCONV_PPARG_KO_CD4_TCELL_FROM_VISCERAL_ADIPOSE_TISSUE_DN", "GSE17186_CD21LOW_VS_CD21HIGH_TRANSITIONAL_BCELL_DN", "GSE1460_DP_VS_CD4_THYMOCYTE_UP"]
# plot_selected_vs_unselected(top_c7, c7, clusterFour, selected, plot_mode='fold_change')

# selected = ["HALLMARK_APOPTOSIS", "HALLMARK_UV_RESPONSE_DN", "HALLMARK_MTORC1_SIGNALING"]
# plot_selected_vs_unselected(top_hallmark, hallmark, clusterFour, selected, plot_mode='fold_change')

#<----- CLUSTER ONE MM Ethan ----- >

clusterOne = pd.read_csv("gene_counts/clusterOne_annotated.csv")
top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterOneDE_GSEA_C2.csv", header=None, sep='\s+').head(10)
c2_arr = top_c2[top_c2.columns[1]].tolist()
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterOneDE_GSEA_C5.csv", header=None, sep='\s+').head(10)
c5_arr = top_c5[top_c5.columns[1]].tolist()
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterOneDE_GSEA_C7.csv", header=None, sep='\s+').head(10)
c7_arr = top_c7[top_c7.columns[1]].tolist()
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterOneDE_GSEA_Hallmark.csv", header=None, sep='\s+').head(10)
hallmark_arr = top_hallmark[top_hallmark.columns[1]].tolist()

#<-- Differentially Expressed Genes -->
# plot_pathways(top_c2, c2, clusterOne, "DE_genes/clusterOneDE3tp.csv")
# plot_pathways(top_c5, c5, clusterOne, "DE_genes/clusterOneDE3tp.csv")
# plot_pathways(top_c7, c7, clusterOne, "DE_genes/clusterOneDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterOne, "DE_genes/clusterOneDE3tp.csv")

#<-- All Genes -->
# plot_pathways(top_c2, c2, clusterOne)
# plot_pathways(top_c5, c5, clusterOne)
# plot_pathways(top_c7, c7, clusterOne)
# plot_pathways(top_hallmark, hallmark, clusterOne)

#<-- ALL GENES ABS FOLD CHANGE -->
# plot_pathways(top_c2, c2, clusterOne, plot_mode="fold_change")
# plot_pathways(top_c5, c5, clusterOne, plot_mode="fold_change")
# plot_pathways(top_c7, c7, clusterOne, plot_mode="fold_change")
# plot_pathways(top_hallmark, hallmark, clusterOne, plot_mode="fold_change")

#<-- ALL GENES RAW + ABS FOLD CHANGE -->
# plot_raw_and_afc(top_c2, c2, clusterOne)
# plot_raw_and_afc(top_c5, c5, clusterOne)
# plot_raw_and_afc(top_c7, c7, clusterOne)
# plot_raw_and_afc(top_hallmark, hallmark, clusterOne)

#<-- DE GENES RAW + DE GENES ABS FC -->
plot_raw_and_afc(top_c2, c2, clusterOne, deg_df="DE_genes/clusterOneDE3tp.csv", deg_only=True)
plot_raw_and_afc(top_c5, c5, clusterOne, deg_df="DE_genes/clusterOneDE3tp.csv", deg_only=True)
plot_raw_and_afc(top_c7, c7, clusterOne, deg_df="DE_genes/clusterOneDE3tp.csv", deg_only=True)
plot_raw_and_afc(top_hallmark, hallmark, clusterOne, deg_df="DE_genes/clusterOneDE3tp.csv", deg_only=True)

#<-- Selected Pathways vs Non Selected -->
# selected =["BIOCARTA_EIF2_PATHWAY"]
# plot_selected_vs_unselected(top_c2, c2, clusterOne, selected, plot_mode='fold_change')
# selected = ["GOBP_PRESYNAPTIC_MODULATION_OF_CHEMICAL_SYNAPTIC_TRANSMISSION", "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY", "HP_SLOWLY_PROGRESSIVE", "HP_APLASIA_HYPOPLASIA_OF_THE_EAR"]
# plot_selected_vs_unselected(top_c5, c5, clusterOne, selected, plot_mode='fold_change')
# selected = ["GSE29618_MONOCYTE_VS_PDC_DN", "GSE6259_DEC205_POS_DC_VS_CD4_TCELL_UP", "GSE360_L_DONOVANI_VS_L_MAJOR_MAC_DN", "GSE25123_CTRL_VS_ROSIGLITAZONE_STIM_MACROPHAGE_DN"]
# plot_selected_vs_unselected(top_c7, c7, clusterOne, selected, plot_mode='fold_change')
# selected = ["HALLMARK_KRAS_SIGNALING_UP"]
# plot_selected_vs_unselected(top_hallmark, hallmark, clusterOne, selected, plot_mode="fold_change")

# #<----- CLUSTER THREE MM Ethan ----- >

clusterThree = pd.read_csv("gene_counts/clusterThree_annotated.csv")
top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterThreeDE_GSEA_C2.csv", header=None, sep='\s+').head(10)
c2_arr = top_c2[top_c2.columns[1]].tolist()
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterThreeDE_GSEA_C5.csv", header=None, sep='\s+').head(10)
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterThreeDE_GSEA_C7.csv", header=None, sep='\s+').head(10)
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterThreeDE_GSEA_Hallmark.csv", header=None, sep='\s+').head(10)
hallmark_arr = top_hallmark[top_hallmark.columns[1]].tolist()


# plot_pathways(top_c2, c2, clusterThree, "DE_genes/clusterThreeDE3tp.csv")
# plot_pathways(top_c5, c5, clusterThree, "DE_genes/clusterThreeDE3tp.csv")
# plot_pathways(top_c7, c7, clusterThree, "DE_genes/clusterThreeDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterThree, "DE_genes/clusterThreeDE3tp.csv")

# plot_pathways(top_c2, c2, clusterThree)
# plot_pathways(top_c5, c5, clusterThree)
# plot_pathways(top_c7, c7, clusterThree)
# plot_pathways(top_hallmark, hallmark, clusterThree)

#<-- ALL GENES ABS FOLD CHANGE -->
# plot_pathways(top_c2, c2, clusterThree, plot_mode="fold_change")
# plot_pathways(top_c5, c5, clusterThree, plot_mode="fold_change")
# plot_pathways(top_c7, c7, clusterThree, plot_mode="fold_change")
# plot_pathways(top_hallmark, hallmark, clusterThree, plot_mode="fold_change")

#<-- DE GENES RAW + DE GENES ABS FC -->


#<-- Selected Pathways vs Non Selected -->
# selected = ["SANA_RESPONSE_TO_IFNG_DN"]
# plot_selected_vs_unselected(top_c2, c2, clusterThree, selected)

# selected = ["HALLMARK_HEDGEHOG_SIGNALING"]
# plot_selected_vs_unselected(top_hallmark, hallmark, clusterThree, selected)

#<--------- MATTS ANALYSIS  ---------->
#<--- CLUSTER ONE --->

clusterOne = pd.read_csv("gene_counts/clusterOne_annotated.csv")
top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterOneDE_GSEA_C2_m.csv", header=None, sep='\s+').head(10)
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterOneDE_GSEA_C5_m.csv", header=None, sep='\s+').head(10)
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterOneDE_GSEA_C7_m.csv", header=None, sep='\s+').head(10)
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterOneDE_GSEA_Hallmark_m.csv", header=None, sep='\s+').head(10)
hallmark_arr = top_hallmark[top_hallmark.columns[1]].tolist()

# plot_pathways(top_c2, c2, clusterOne, "DE_genes_m/clusterOneDE3tp.csv")
# plot_pathways(top_c5, c5, clusterOne, "DE_genes_m/clusterOneDE3tp.csv")
# plot_pathways(top_c7, c7, clusterOne, "DE_genes_m/clusterOneDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterOne, "DE_genes_m/clusterOneDE3tp.csv")

# plot_pathways(top_c2, c2, clusterOne)
# plot_pathways(top_c5, c5, clusterOne)
# plot_pathways(top_c7, c7, clusterOne)
# plot_pathways(top_hallmark, hallmark, clusterOne)

#<-- selected pathways vs non-selected -->
# plot_pathways(top_hallmark, hallmark, clusterOne, pathway_subset=["HALLMARK_PANCREAS_BETA_CELLS"])
# hallmark_arr.remove("HALLMARK_PANCREAS_BETA_CELLS")
# plot_pathways(top_hallmark, hallmark, clusterOne, pathway_subset=hallmark_arr)

#<--- CLUSTER TWO --->

clusterTwo = pd.read_csv("gene_counts/clusterTwo_annotated.csv")
top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterTwoDE_GSEA_C2_m.csv", header=None, sep='\s+').head(10)
c2_arr = top_c2[top_c2.columns[1]].tolist()
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterTwoDE_GSEA_C5_m.csv", header=None, sep='\s+').head(10)
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterTwoDE_GSEA_C7_m.csv", header=None, sep='\s+').head(10)
c7_arr = top_c7[top_c7.columns[1]].tolist()
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterTwoDE_GSEA_Hallmark_m.csv", header=None, sep='\s+').head(10)

# plot_pathways(top_c2, c2, clusterTwo, "DE_genes_m/clusterTwoDE3tp.csv")
# plot_pathways(top_c5, c5, clusterTwo, "DE_genes_m/clusterTwoDE3tp.csv")
# plot_pathways(top_c7, c7, clusterTwo, "DE_genes_m/clusterTwoDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterTwo, "DE_genes_m/clusterTwoDE3tp.csv")

# plot_pathways(top_c2, c2, clusterTwo)
# plot_pathways(top_c5, c5, clusterTwo)
# plot_pathways(top_c7, c7, clusterTwo)
# plot_pathways(top_hallmark, hallmark, clusterTwo)

#<--selected pathways vs non-selected -->
# plot_pathways(top_c2, c2, clusterTwo, pathway_subset=["JIANG_MELANOMA_TRM8_CD4", "JIANG_MELANOMA_TRM9_CD8"])
# c2_arr = [x for x in c2_arr if x not in ["JIANG_MELANOMA_TRM8_CD4", "JIANG_MELANOMA_TRM9_CD8"]]
# plot_pathways(top_c2, c2, clusterTwo, pathway_subset=c2_arr)

# plot_pathways(top_c7, c7, clusterTwo, pathway_subset=["GSE3565_CTRL_VS_LPS_INJECTED_DUSP1_KO_SPLENOCYTES_UP", "GSE40685_NAIVE_CD4_TCELL_VS_FOXP3_KO_TREG_PRECURSOR_UP"])
# c7_arr = [x for x in c7_arr if x not in ["GSE3565_CTRL_VS_LPS_INJECTED_DUSP1_KO_SPLENOCYTES_UPP", "GSE40685_NAIVE_CD4_TCELL_VS_FOXP3_KO_TREG_PRECURSOR_UP"]]
# plot_pathways(top_c7, c7, clusterTwo, pathway_subset=c7_arr)

#<--- CLUSTER THREE --->

clusterThree = pd.read_csv("gene_counts/clusterThree_annotated.csv")
top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterThreeDE_GSEA_C2_m.csv", header=None, sep='\s+').head(10)
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterThreeDE_GSEA_C5_m.csv", header=None, sep='\s+').head(10)
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterThreeDE_GSEA_C7_m.csv", header=None, sep='\s+').head(10)
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterThreeDE_GSEA_Hallmark_m.csv", header=None, sep='\s+').head(10)


# plot_pathways(top_c2, c2, clusterThree, "DE_genes_m/clusterThreeDE3tp.csv")
# plot_pathways(top_c5, c5, clusterThree, "DE_genes_m/clusterThreeDE3tp.csv")
# plot_pathways(top_c7, c7, clusterThree, "DE_genes_m/clusterThreeDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterThree, "DE_genes_m/clusterThreeDE3tp.csv")

# plot_pathways(top_c2, c2, clusterThree)
# plot_pathways(top_c5, c5, clusterThree)
# plot_pathways(top_c7, c7, clusterThree)
# plot_pathways(top_hallmark, hallmark, clusterThree)


#<--- Cluster4 --->
clusterFour = pd.read_csv("gene_counts/clusterFour_annotated.csv")

top_c2 = pd.read_csv("../gsea/pathway_results/c2_3tp/clusterFourDE_GSEA_C2_m.csv", header=None, sep='\s+').head(10)
top_c5 = pd.read_csv("../gsea/pathway_results/c5_3tp/clusterFourDE_GSEA_C5_m.csv", header=None, sep='\s+').head(10)
top_c7 = pd.read_csv("../gsea/pathway_results/c7_3tp/clusterFourDE_GSEA_C7_m.csv", header=None, sep='\s+').head(10)
top_hallmark = pd.read_csv("../gsea/pathway_results/hallmark_3tp/clusterFourDE_GSEA_Hallmark_m.csv", header=None, sep='\s+').head(10)

# plot_pathways(top_c2, c2, clusterFour, "DE_genes_m/clusterFourDE3tp.csv")
# plot_pathways(top_c5, c5, clusterFour, "DE_genes_m/clusterFourDE3tp.csv")
# plot_pathways(top_c7, c7, clusterFour, "DE_genes_m/clusterFourDE3tp.csv")
# plot_pathways(top_hallmark, hallmark, clusterFour, "DE_genes/clusterFourDE3tp.csv")


# plot_pathways(top_c2, c2, clusterFour)
# plot_pathways(top_c5, c5, clusterFour)
# plot_pathways(top_c7, c7, clusterFour)
# plot_pathways(top_hallmark, hallmark, clusterFour)
