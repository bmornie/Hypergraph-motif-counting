import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib as mpl


mpl.rcParams["pdf.fonttype"] = 42


datasets = ["hospital", "primary-school", "high-school", "sfhh", "eu", "ask-ubuntu", "history", "dblp", "hypertext", "enron", "math-sx", "geology"]

domains = ["social", "social", "social", "social", "email", "online tags", "co-authorship", "co-authorship", "social", "email", "online tags", "co-authorship"]


groups = {d:[] for d in domains}
for i, d in enumerate(domains):
    groups[d].append(i)


abundance_sp = []
abundance_conf = []


for name in datasets:
    abundance = []
    with open(f"{name}_abundance_sp_4.txt", "r") as f:
        for line in f:
            abundance.append(float(line.split("\t")[-1]))
    abundance_sp.append(np.array(abundance))


for name in datasets:
    abundance = []
    with open(f"{name}_abundance_conf_4.txt", "r") as f:
        for line in f:
            abundance.append(float(line.split("\t")[-1]))
    abundance_conf.append(np.array(abundance))


profile_sp = [abundance/np.sqrt(np.sum(abundance**2)) for abundance in abundance_sp]
profile_conf = [abundance/np.sqrt(np.sum(abundance**2)) for abundance in abundance_conf]

correlation_sp = np.corrcoef(profile_sp)
correlation_conf = np.corrcoef(profile_conf)


group_order = []

for d in ["social", "email", "online tags", "co-authorship"]:
    idx = groups[d]
    
    if len(idx) > 1:
        
        # cluster
        sub_profiles = correlation_conf[idx]
        condensed = pdist(sub_profiles, metric='correlation')
        Z = linkage(condensed, method='average')
        order = leaves_list(Z)
        
        # map back to original indices
        ordered_idx = [idx[i] for i in order]
    else:
        ordered_idx = idx
    
    group_order.extend(ordered_idx)


correlation_sp = correlation_sp[group_order][:,group_order]
correlation_conf = correlation_conf[group_order][:,group_order]
domains_ordered = [domains[i] for i in group_order]

labels = ["Hospital", "Primary School", "High School", "SFHH", "EU", "Ask Ubuntu", "History", "dblp", "Hypertext", "Enron", "Math SX", "Geology"]
labels_ordered = [labels[i] for i in group_order]

fontsize = 16


def truncate_colormap(cmap_name, minval=0.15, maxval=0.85, n=256):
    cmap = plt.get_cmap(cmap_name)
    new_colors = cmap(np.linspace(minval, maxval, n))
    return colors.LinearSegmentedColormap.from_list('trunc', new_colors)

def plot_side_by_side(corr1, corr2, labels, domains, title1, title2):
    
    unique_domains = sorted(set(domains))
    
    cmap = plt.get_cmap('tab10')

    domain_to_color = {
    d: cmap(i % cmap.N)
    for i, d in enumerate(unique_domains)
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True, gridspec_kw={'wspace': 0.15})

    vmin, vmax = -1, 1
    cmap = truncate_colormap('RdBu_r', 0.1, 0.9)

    ims = []
    for idx, (ax, corr, title) in enumerate(zip(axes, [corr1, corr2], [title1, title2])):
        im = ax.imshow(corr, cmap=cmap, vmin=vmin, vmax=vmax)
        
        # remove border
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        ims.append(im)

        ax.set_title(title, fontsize=fontsize+2)

        # remove x-axis labels
        ax.set_xticks([])

        # y-axis labels only on the left plot
        if idx == 0:
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=fontsize)
            for label, d in zip(ax.get_yticklabels(), domains):
                label.set_color(domain_to_color[d])
        else:
            ax.set_yticks([])


    # shared colorbar
    cbar = fig.colorbar(ims[0], ax=axes, location='right', shrink=0.9)
    cbar.ax.tick_params(labelsize=fontsize-1)
    
    plt.savefig("correlation.pdf", bbox_inches="tight")
    plt.show()

plot_side_by_side(correlation_sp, correlation_conf, labels_ordered, domains_ordered, "Size-preserving model", "Configuration model")

