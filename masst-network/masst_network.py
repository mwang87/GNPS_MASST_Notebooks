import requests
import pandas as pd
import networkx as nx
import requests_cache
requests_cache.install_cache('demo_cache')

MASST_TASK = "c6b2797224f34d819d20dd7af622bc6b"

# Loading MASST information

spectra_matches_df = pd.read_csv("https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task={}&file=all_dataset_spectra_matches/&block=main".format(MASST_TASK), sep="\t")
spectra_matches_df.head()

# Loading all Datasets Information
dataset_matches = list(set(spectra_matches_df["dataset_id"]))
all_datasets = requests.get("https://massive.ucsd.edu/ProteoSAFe/datasets_json.jsp#%7B%22query%22%3A%7B%7D%2C%22table_sort_history%22%3A%22createdMillis_dsc%22%7D").json()["datasets"]


all_node_usi_list = []

# Source MASST USI
all_node_usi_list.append("mzspec:GNPS:TASK-c6b2797224f34d819d20dd7af622bc6b-spectra/:scan:1")

# Getting all the MASST data
for dataset in dataset_matches:
    print(dataset)
    filtered_dataset = [current_dataset for current_dataset in all_datasets if current_dataset["dataset"] == dataset]
    dataset_task = filtered_dataset[0]["task"]
    continuous_id = requests.get("http://gnps.ucsd.edu/ProteoSAFe/ContinuousIDServlet?task={}".format(dataset_task)).json()
    
    network_url = "https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={}&view=clusters_network_pairs".format(continuous_id["jobs"][0]["task"])
    data = requests.get(network_url).json()['blockData']
    network_df = pd.DataFrame(data)
    
    dataset_spectra_matches = spectra_matches_df[spectra_matches_df["dataset_id"] == dataset]
    clusters_matched = list(set(dataset_spectra_matches["cluster_scan"]))
    print(clusters_matched)
    
    network_df["Node1"] = network_df["Node1"].astype(int)
    filtered_edges = network_df[network_df["Node1"].isin(clusters_matched)]
    print(len(filtered_edges))
    
    for edge in filtered_edges.to_dict(orient="records"):
        cluster = edge["Node2"]
        usi = "mzspec:GNPS:TASK-{}-speccontinuous/speccontinuous-00000.mgf:scan:{}".format(continuous_id["jobs"][0]["task"], cluster)
        all_node_usi_list.append(usi)

print(len(all_node_usi_list), "Total Spectra")

# Now we will load up all the spectra and do stuff with it
from ming_spectrum_library import Spectrum
import spectrum_alignment
all_spectra_list = []
for usi in all_node_usi_list:
    url = "https://metabolomics-usi.ucsd.edu/json/?usi={}".format(usi)
    spectrum_json = requests.get(url).json()

    spectrum = Spectrum("", usi, usi, spectrum_json["peaks"], spectrum_json["precursor_mz"], 1, 2)
    all_spectra_list.append(spectrum)

min_score = 0.7
min_matched_peaks = 4

# Let's create a network now
G = nx.Graph()
from tqdm import tqdm
for i, spectrum1 in tqdm(enumerate(all_spectra_list)):
    for j, spectrum2 in enumerate(all_spectra_list):
        if i <= j:
            continue
        
        # Doing a network here
        total_score, reported_alignments = spectrum_alignment.score_alignment(spectrum1.peaks, spectrum2.peaks, spectrum1.mz, spectrum2.mz, 0.5, max_charge_consideration=1)
        if total_score < min_score:
            continue
        
        if len(reported_alignments) < min_matched_peaks:
            continue

        G.add_edge(spectrum1.scan, spectrum2.scan)

import matplotlib.pyplot as plt
nx.draw(G, with_labels=True, font_weight='bold')        
plt.savefig("Graph.png", format="PNG")
