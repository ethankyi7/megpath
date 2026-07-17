import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import io


# with open("denseLT_error.txt") as dense_error:
#     error = []
#     for line in dense_error:
#         data = line.split(':')
#         error.append(float(data[1].strip()))
    
#     x = [i for i in range(1, 11)]
#     plt.plot(x, error)
#     plt.title("Log Transformed Data: Total Error vs # of Patterns")
#     plt.ylabel("total error")
#     plt.xlabel("# of patterns")

#     for i,j in zip(x, error):
#         plt.annotate(f'{j:.2f}', (i,j), textcoords="offset points", xytext=(15,10), ha='center')

#     plt.show()

# with open("dense_error.txt") as dense_error:
#     error = []
#     for line in dense_error:
#         data = line.split(':')
#         error.append(float(data[1].strip()))
    
#     x = [i for i in range(1, 11)]
#     plt.plot(x, error)
#     plt.title("Original Data: Total Error vs # of Patterns")
#     plt.ylabel("total error")
#     plt.xlabel("# of patterns")

#     for i,j in zip(x, error):
#         plt.annotate(f'{j:.2f}', (i,j), textcoords="offset points", xytext=(15,10), ha='center')

#     plt.show()
    

# with open("tcr_error.txt") as tcr_error:
#     error = []
#     for line in tcr_error:
#         data = line.split(':')
#         error.append(float(data[1].strip()))
    
#     x = [i for i in range(1, 11)]
#     plt.plot(x, error)
#     plt.title("Log Transformed TCR Data: Total Error vs # of Patterns")
#     plt.ylabel("total error")
#     plt.xlabel("# of patterns")

#     for i,j in zip(x, error):
#         plt.annotate(f'{j:.2f}', (i,j), textcoords="offset points", xytext=(15,10), ha='center')

#     plt.show()


# patterns = []
# weight_data_lines = []

# with open("../standard/tcr_resultsLT/tcr_results_3.csv", 'r') as tcr_results:
#     for line in tcr_results:
#         if line.startswith('#""'):
#             vals = [float(x) for x in line.strip().split(',')[1].split()]
#             patterns.append(vals)
#         elif not line.startswith('#'):
#             weight_data_lines.append(line)
# print(weight_data_lines[0])

# patterns = np.array(patterns)
# w_df = pd.read_csv(io.StringIO("".join(weight_data_lines)), header=None)
# weights = w_df.iloc[:, 1:7].values
# print(weights)

# data_lt = pd.read_csv("tcr_dataLT.csv", header=None).values
    
# dominant_pattern = np.argmax(weights, axis=1)
    
# time_points = ['T1', 'T2', 'T3', 'T4']
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
# colors = plt.cm.tab10(np.linspace(0, 1, len(patterns)))
    

# for i in range(len(patterns)):
#     p = patterns[i]
#     p_norm = p / (np.max(p) if np.max(p) > 0 else 1)
#     ax1.plot(time_points, p_norm, marker='o', label=f'Pattern {i+1}', linewidth=2, color=colors[i])
    
# ax1.set_title('A. NMF Temporal Basis Shapes (Normalized)')
# ax1.set_xlabel('Time Point')
# ax1.set_ylabel('Relative Contribution')
# ax1.legend()
# ax1.grid(True, linestyle='--', alpha=0.6)
    
# for i in range(len(patterns)):
#     cluster_indices = np.where(dominant_pattern == i)[0]
#     if len(cluster_indices) > 0:
#         cluster_mean = np.mean(data_lt[cluster_indices], axis=0)
#         ax2.plot(time_points, cluster_mean, marker='s', label=f'Cluster {i+1} (n={len(cluster_indices)})', color=colors[i], linewidth=2.5)
    
# ax2.set_title('B. Mean Expansion Profiles (Log Transformed Data)')
# ax2.set_xlabel('Time Point')
# ax2.set_ylabel('Clonal Frequency')
# ax2.legend()
# ax2.grid(True, linestyle='--', alpha=0.6)
    
# plt.tight_layout()
# plt.savefig('tcr_nmf_comparison.png', dpi=300)
# plt.show()


# clonotype_df = pd.read_csv("kmeans_clusters.csv")
# clonotype_df = clonotype_df.drop_duplicates(subset=['TCRB.clonotype'])
# clonotype_map = clonotype_df[["TCRB.clonotype", "cluster"]]

# print((clonotype_map['cluster'] == 1).sum())
# print((clonotype_map['cluster'] == 2).sum())
# print((clonotype_map['cluster'] == 3).sum())
# print((clonotype_map['cluster'] == 4).sum())

# print(clonotype_map)
# clonotype_map.to_csv("clonotype_map.csv", index=None)

patternOne = [0.5024788,0.9997918,0.9837637,0.9939621]
nmfPatternOne = [0.999794, 0.983766, 0.993965, 0.502481]

patternTwo = [0.4972853,0.7214298,0.8182699,0.8848034]
nmfPatternTwo = [0.721921, 0.818761, 0.885295, 0.497777]

patternThree = [0.5087043,0.0000000,0.4846397,0.4731947]
nmfPatternThree = [0.0713356, 0.555975, 0.54453, 0.58004]

patternFour = [0.5024788, 1.0000000, 0.8684869, 0.9493119]
nmfPatternFour = [0.243095, 0.756694, 0.336577, 0.374999]

def plot_data(pattern, data, cluster):
    x = [0, 3, 6, 9]
    plt.plot(x, pattern, label="input pattern")
    plt.plot(x, data, label="identified pattern")
    plt.title(f'Cluster{cluster}: Coefficients vs timepoint')
    plt.ylabel("coefficients")
    plt.xlabel("timepoint")
    plt.legend()
    plt.xticks(x)
    for i,j in zip(x, data):
        plt.annotate(f'{j:.2f}', (i,j), textcoords="offset points", xytext=(15,10), ha='center')
    for i,j in zip(x, pattern):
        plt.annotate(f'{j:.2f}', (i,j), textcoords="offset points", xytext=(15,10), ha='center')
    plt.show()

# plot_data(patternOne, nmfPatternOne, 1)
# plot_data(patternTwo, nmfPatternTwo, 2)
# plot_data(patternThree, nmfPatternThree, 3)
# plot_data(patternFour, nmfPatternFour, 4)

lfc1 = [0.231723, 0.438537, 0.434193, 0.417186]
lfc_pat1 = [0.0000000, 1.0000000, 0.3563981, 0.4796209]
plot_data(lfc_pat1, lfc1, 1)

lfc2 = [0.527994, 0.26754, 0.376936, 0.378721]
lfc_pat2 = [0.0, 0.4394, 0.6787, 1.0]
plot_data(lfc_pat2, lfc2, 2)

lfc3 = [0.460237, 0.267471, 0.461853, 0.462834]
lfc_pat3 = [1.0, 0.0, 0.9111, 0.8716]
plot_data(lfc_pat3, lfc3, 3)

lfc4 = [0.418861, 0.416119, 0.450006, 0.217711]
lfc_pat4 = [0.0, 1.0, 0.0138, 0.0241]
plot_data(lfc_pat4, lfc4, 4)

c1 = [1.45, 12.0, 5.21, 6.51]
c2 = [0.0371,0.803,1.22,1.78]
c3 = [0.0405, 0, 0.0369, 0.0353]
c4 = [0.00913, 0.638, 0.0178, 0.0243]
x = [0, 3, 6, 9]
# plt.plot(x, [10, 20, 15, 25])
# # plt.plot(x, [10, 20, 15, 4])
# plt.xticks(x)
# plt.show()