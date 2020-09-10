# VISION Python API

A library that hooks into a VISION session and allows for direct interfacing with the underlying VISION API using Requests.  
Provides a convenience interface with the VISION data and allows for VISION-like plotting.  
Seldom it is necessary to use the api directly. The plot and data module should provide most functionalities

## Authors

AJ Zwijnenburg

## Requirements

Python >= 3.8.1  
requests >= 2.22.0  
pandas >= 1.1.1  
plotnine >= 0.7.1

## Installation

Copy the folder 'vision' with its contents into the project folder of the project and import directly.

## Usage

API example:  

```python
from vision.api import API

api = API("[VISION session link here]")

# Signatures information
signatures_names = api.signatures.names
signatures_score = api.signatures.score(signature_name[0])
signature_info = api.signatures.signature(signature_name[0])

# Proteins information
proteins_names = api.proteins.names
proteins_value = api.proteins.value(protein_names[0])

# Clusters information
clusters_names = api.clusters.names
clusters_levels = api.clusters.levels
clusters_proteins = api.clusters.proteins(clusters_names[0])
clusters_signatures = api.clusters.signatures(clusters_names[0])
clusters_metadata = api.clusters.metadata(clusters_names[0])

# Latent Class Analysis information
lca_proteins = api.lca.proteins
lca_signatures = api.lca.signatures
lca_metadata = api.lca.metadata

# Projections information
projections_names = api.projections.names
projections_data = api.projections.projection(projections_names[0])

# Expression information
gene_names = api.expression.names
gene_expression = api.expression.expression(gene_names[0])

# Cell metadata
cell_metadata = api.cell.cell("0")
```

Data example:

```python
from vision.data import Data

data = Data("[VISION session link here]")

# Data object contains (most) data of the vision interface
# Use the cells attribute to recover all cell ids
cell_ids = data.cells

# The expression data is stored in the expression attribute
gene_names = data.expression.names
gene_count = data.expression[gene_names[0]]
# values can also always be requested using a list of gene names
gene_counts = data.expression[gene_names[0:3]]

# The protein data is stored in the protein attribute
protein_names = data.protein.names
protein_count = data.protein[protein_names[0]]

# The signature data is stored in the protein attribute
signature_names = data.signature.names
signature_score = data.signature[signature_names[0]]

# The metadata is stored in the meta attribute
meta_names = data.meta.names
meta_levels = data.meta.levels # only for discrete metadata
meta_values = data.meta[cell_ids[0]]

# The projection data is stored in the projection attribute
projection_names = data.projection.names
projection_dimensions = data.projection.dimensions
projection_values = data.projection[projection_names[0]]
```

Plot example:

```python
from vision.plot import Plot

# The plot class handles all data acquisition automatically
plot = Plot("[VISION session link here]")
# You only need to define the type of plot and which parameters to use
scatter = Plot.scatter("projection_name", x="x", y="y", color="gene")

# For discrete parameters a custom color_map can be added
scatter = Plot.scatter("projection_name", x="x", y="y", color="discrete", color_map={"1":"red", "2":"blue"})
```

## Version Info

v1.0 - Implemented the API, Data containers, and scatter plot

## What Can Still Be Done

General - Some cleanup of the API, all caching should be handled by the Data  
API - Implementation of the Yanay, Tree, Analysis, Cells APIs  
Data - Implementation of Clusters information (and the above)  
Plot - Implementation of VISIONs barplots / heatmaps

## License

[MIT](https://choosealicense.com/licenses/mit/)
