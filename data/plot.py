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


patterns = []
weight_data_lines = []

with open("../standard/tcr_resultsLT/tcr_results_3.csv", 'r') as tcr_results:
    for line in tcr_results:
        if line.startswith('#""'):
            vals = [float(x) for x in line.strip().split(',')[1].split()]
            patterns.append(vals)
        elif not line.startswith('#'):
            weight_data_lines.append(line)
print(weight_data_lines[0])

patterns = np.array(patterns)
w_df = pd.read_csv(io.StringIO("".join(weight_data_lines)), header=None)
weights = w_df.iloc[:, 1:7].values
print(weights)

data_lt = pd.read_csv("tcr_dataLT.csv", header=None).values
    
dominant_pattern = np.argmax(weights, axis=1)
    
time_points = ['T1', 'T2', 'T3', 'T4']
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
colors = plt.cm.tab10(np.linspace(0, 1, len(patterns)))
    

for i in range(len(patterns)):
    p = patterns[i]
    p_norm = p / (np.max(p) if np.max(p) > 0 else 1)
    ax1.plot(time_points, p_norm, marker='o', label=f'Pattern {i+1}', linewidth=2, color=colors[i])
    
ax1.set_title('A. NMF Temporal Basis Shapes (Normalized)')
ax1.set_xlabel('Time Point')
ax1.set_ylabel('Relative Contribution')
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.6)
    
for i in range(len(patterns)):
    cluster_indices = np.where(dominant_pattern == i)[0]
    if len(cluster_indices) > 0:
        cluster_mean = np.mean(data_lt[cluster_indices], axis=0)
        ax2.plot(time_points, cluster_mean, marker='s', label=f'Cluster {i+1} (n={len(cluster_indices)})', color=colors[i], linewidth=2.5)
    
ax2.set_title('B. Mean Expansion Profiles (Log Transformed Data)')
ax2.set_xlabel('Time Point')
ax2.set_ylabel('Clonal Frequency')
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.6)
    
plt.tight_layout()
plt.savefig('tcr_nmf_comparison.png', dpi=300)
plt.show()