# VISION Python API

A library that hooks into a [VISION](https://yoseflab.github.io/software/vision/) session and allows for direct interfacing with the underlying VISION API using Requests. Provides an interface of pandas object with the data from VISION and allows for VISION-like plotting. To reduce data-usage all received data is buffered. Provides access to VISION's signature scores, gene & protein expression, clustering, projections, and latent class analysis.  

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

The data module provides a pythony interface with the data from VISION.

```python
from vision.data import Data

data = Data("[VISION session link here]")

# Data object contains (most) data of the vision interface
# Use the cells attribute to recover all cell ids
cell_ids = data.cells

# The expression data is stored in the expression attribute
# First get the gene name
gene_names = data.expression.names
# Next request the expression counts for a specific gene
gene_count = data.expression[gene_names[0]]
# Data can always be requested using a list of names too
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

# For convenience all data can be requested directly
# This will return the first entree with the specified name
# the search order is: expression->meta->protein->signature
values = data["gene | meta | protein | signature name"]
```

The plot module reimplements the basic VISION plots. This allows for automated plotting and presentation generation.

```python
from vision.plot import Plot

# The plot class handles all data acquisition automatically
plot = Plot("[VISION session link here]")

# The vision.data can be reached using the .data attribute
plot.data

# To plot a scatter plot of a projection, use:
scatter = plot.projection("projection_name", x="x", y="y")

# c(olor) parameter can be used to overlay the (normalized) gene count
scatter = plot.projection("projection_name", x="x", y="y", c="gene")

# Or if you rather plot two genes against eachother, use:
scatter = plot.comparison(x="TBX21", y="GZMB")

# To plot only a subset of the data, you can add a mask
plot.mask = plot.data["discrete_parameter"] == "1"

# For discrete parameters a custom color_map can be added
scatter = plot.projection("projection_name", x="x", y="y", c="discrete c", c_map={"1":"red", "2":"blue"})

# Bar plots can also be created, and can be log10(x+1) scaled
bar = plot.bar(x="x", log=False)

# Stacked bar plots can be made for categorical categories
bar_stacked = plot.bar_stacked(x="x", y="categorical y", y_map={"1":"red", "2":"blue"})
```

The api provides a low-level interface with the VISION api. Use the data module for a user-friendly approach.

```python
from vision.api import API

api = API("[VISION session link here]")

# Signatures information
signatures_names = api.signatures.names()
signatures_score = api.signatures.score(signature_names[0])
signature_info = api.signatures.signature(signature_names[0])
metadata_names = api.signatures.metadata_names()
metadata_value = api.signatures.metadata(metadata_names[0])

# Expression information
gene_names = api.expression.names()
gene_expression = api.expression.expression(gene_names[0])

# Proteins information
proteins_names = api.proteins.names()
proteins_value = api.proteins.value(protein_names[0])

# Clusters information
clusters_names = api.clusters.names()
clusters_levels = api.clusters.levels()
clusters_proteins = api.clusters.proteins(clusters_names[0])
clusters_signatures = api.clusters.signatures(clusters_names[0])
clusters_metadata = api.clusters.metadata(clusters_names[0])

# Projections information
projections_names = api.projections.names()
projections_dimensions = api.projections.dimensions()
projections_data = api.projections.projection(
    projections_names[0],
    projections_dimensions[projections_names[0]][0]
)

# Latent Class Analysis information
lca_proteins = api.lca.proteins()
lca_signatures = api.lca.signatures()
lca_metadata = api.lca.metadata()

# Per cell metadata
cell_metadata = api.cell.cell("0")
```

## Version Info

v1.5 - Updated README  
v1.4 - Implemented bar_stacked plots  
v1.3 - Implemented violin plots  
v1.2 - Implemented bar plots  
v1.1 - Removed some API inconsistencies  
v1.0 - Implemented the API, Data containers, and scatter plot

## What Can Still Be Done

API - Implementation of the Yanay, Tree, Analysis, Cells APIs  
Data - Handle for Clustering and LCA results  
Plot - Implementation of VISIONs heatmaps

## License

[MIT](https://choosealicense.com/licenses/mit/)
