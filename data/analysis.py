import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

results = "../standard/clusters"
resultsLT = "../standard/clustersLT"
genes_csv = "gene_list.csv"

genes = pd.read_csv(genes_csv, header=None).squeeze().tolist()

def extract_coeff(folder) -> pd.DataFrame:
    cluster_dfs = {}
    
    for filename in sorted(os.listdir(folder)):
        path = os.path.join(folder, filename)
        rows = []
        
        with open(path) as file:
            for line in file:
                if line.startswith("row-"):
                    values = [float(x) for x in line[4:].strip().split(",")[1:4]]
                    rows.append(values)
        
        num_patterns = len(rows[0]) if rows else 0
        pattern_cols = [f"pattern{i+1}" for i in range(num_patterns)]
        
        cluster_name = os.path.splitext(filename)[0].split('_')[0]
        cluster_dfs[cluster_name] = pd.DataFrame(rows, columns=pattern_cols, index=genes)
        cluster_dfs[cluster_name].index.name = "gene"
    
    return cluster_dfs

cluster_dfs = extract_coeff(results)
clusterLT_dfs = extract_coeff(resultsLT)
#print(cluster_dfs)

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

patterns = extract_patterns(results)
patternsLT = extract_patterns(resultsLT)
#print(patterns)
#print(patterns)

def top_genes_per_pattern(cluster_df, n=20):
    top_ranked_genes = {}

    for cluster_name, df in cluster_df.items():
        top_ranked_genes[cluster_name] = {}
        for pattern in df.columns:
            top_ranked_genes[cluster_name][pattern] = (df[pattern].sort_values(ascending = False).head(n))
    return top_ranked_genes


top_genes = top_genes_per_pattern(cluster_dfs)
top_genesLT = top_genes_per_pattern(clusterLT_dfs)

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


def plot_genes_per_pattern(patterns, top_genes, lt=False, n=20):
    x = [0, 3, 6, 9]
    suffix = "LT" if lt else ""
    gene_folder = f"gene_counts{suffix}/"

    for cluster_name, pattern_dict in patterns.items():
        print(cluster_name)
        gene_file = os.path.join(gene_folder, f"{cluster_name}_annotated.csv")
        gene_counts = pd.read_csv(gene_file, index_col=0)
        gene_counts.columns = x
        

        for pattern_name, pattern_vals in pattern_dict.items():
            top = top_genes[cluster_name][pattern_name].head(n)
            print(top)
            top_gene_names = top.index.tolist()

            fig, ax = plt.subplots(figsize=(10, 6))

            for gene in top_gene_names:
                if gene in gene_counts.index:
                    ax.plot(x, gene_counts.loc[gene], color="steelblue",
                            alpha=0.4, linewidth=1, label=gene)

            #pattern_vals = [5 * x for x in pattern_vals]
            ax.plot(x, pattern_vals, color="crimson", linewidth=2.5, linestyle="--", label=f"{pattern_name}")
            
            ax.set_title(f"{cluster_name} — {pattern_name} | Top {n} genes {'(LT)' if lt else ''}")
            ax.set_xlabel("Timepoint")
            ax.set_ylabel("Gene counts")
            ax.set_xticks(x)

            handles, labels = ax.get_legend_handles_labels()
            pattern_handle = [(h, l) for h, l in zip(handles, labels) if "scaled" in l]
            gene_handles = [(h, l) for h, l in zip(handles, labels) if "scaled" not in l]

            ax.legend(
                [ph[0] for ph in pattern_handle] + [gene_handles[0][0]],
                [ph[1] for ph in pattern_handle] + [f"Top {n} genes (n={len(gene_handles)})"],
                loc="upper right"
            )

            plt.tight_layout()
            plt.show()

#plot_genes_per_pattern(patterns, top_genes)
#plot_genes_per_pattern(patternsLT, top_genesLT, lt=True)