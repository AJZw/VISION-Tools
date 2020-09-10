##############################################################################     ##    ######
#    A.J. Zwijnenburg                   2020-09-02           v1.0                 #  #      ##
#    Copyright (C) 2020 - AJ Zwijnenburg          MIT license                    ######   ##
##############################################################################  ##    ## ######

## Copyright notice ##########################################################
# Copyright 2020 AJ Zwijnenburg
#
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in  
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE # WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
# THE SOFTWARE.
##############################################################################

"""
An API interface with the Yosef Lab Vision server
(Thx to api.js and its nice and clear writing, made it nice and easy to write)

------------------------------------------------------------------------------

Implemented Vision functions:
- Signature API                     -> _Signatures
- Protein API                       -> _Proteins
- Clusters API                      -> _Clusters
- Pearson Correlations (LCA) API    -> _LCA
- Projections API                   -> _Projections
- Expression API                    -> _Expression
- Cell API                          -> _Cell
- Session Info API                  -> _SessionInfo

Not Implemented functions:
- Yanay API
- Tree API                          -> _Tree
- Analysis API
- Cells API                         -> _Cells

------------------------------------------------------------------------------

:class: SessionError
Error to be raised if the session cannot be requested

:class: DataType
Enumeration to specify the data returned/stored in a Cluster or LatentClass object
These classes are used to store the _Clusters and _LCA results

:class: Signature
A storage class for the data of a single signature. The class store the signatures
metadata, included genes, and gene weights.

:class: Cluster
A storage class for _Clusters results. Stores the zscores and pvalues:
.zscores    - zscore result                         -> pd.DataFrame[zscores] indexed on the compared data labels
.pvalues    - pvalue result                         -> pd.DataFrame[pvalues] indexed on the compared data labels

:class: LatentClass
A storage class for _LCA results. Stores the zscores and pvalues:
.zscores    - zscore result                         -> pd.DataFrame[zscores] indexed on the compared data labels
.pvalues    - pvalue result                         -> pd.DataFrame[pvalues] indexed on the compared data labels

:class: API
The main interface with VISION
.signatures - access the Signature API
.proteins   - access the Protein API
.clusters   - access the Clusters API 
.lca        - access the Pearson Correlations (LCA) API
.projections - access the Projections API
.expression - access the Expression API
.cell       - access the Cell API
.session_info - access the Session Info API

:class: _SessionInfo
The interface with VISIONS Session Info API, provides:
.name       - name of the session                   -> name
.ncells     - number of cells in the session        -> cell n
.metadata   - the session-wide metadata names       -> List[metadata_name]

:class: _Signatures
The interface with VISION Signature API, provides:
.names        - the signature names                 -> List[signature_name]
.signature()  - signature metadata                  -> Signature
.score()      - signature score                     -> pd.Series[score] indexed on cell_id
.expression() - unknown score /metadata/gene        -> pd.DataFrame[score] indexed on gene vs metadata_label

:class: _Proteins
The interface with VISION Proteins API, provides:
.names      - protein names                         -> List[protein_name]
.value()    - protein value/cell_id                 -> pd.Series[value] indexed on cell_id

:class: _Clusters
The interface with VISION Clusters API, provides:
.names      - cluster names                         -> List[cluster_name]
.levels     - levels of each cluster name           -> Dict[cluster_name, List[cluster_levels]]
.proteins() - Cluster info protein/metadata         -> Cluster
.signatures() - Cluster info signature/metadata     -> Cluster
.metadata() - Cluster info metadata/metadata        -> Cluster

:class: _LCA
The interface with VISION LCA API, provides:
.proteins() - Latent Class Analysis of proteins     -> LatentClass
.signatures() - Latent Class Analysis of signatures -> LatentClass
.metadata() - Latent Class Analysis of metadata     -> LatentClass

:class: _Projections
The interface with VISION Projections API, provides:
.names      - The projection names                  -> List[projection_name]
.levels     - The projection levels                 -> Dict[projection_name, List[projection_levels]]
.projection() - The value of all events in all projection dimensions -> pd.DataFrame[value] indexed on dimension vs cell_id

:class: _Expression
The interface with VISION Expression API, provides:
.names      - The gene names                        -> List[gene_name]
.expression() - The expression of a gene in all events -> pd.Series[value] indexed on cell_id

:class: _Cell
The interface with VISION Cell API, provides:
.cell()     - The cells metadata                    -> Dict[metadata_name, metadata_value]

:class: _Tree
Incomplete interface with VISION Tree API

:class: _Cells
Incomplete interface with VISION Cells API

------------------------------------------------------------------------------

"""

from __future__ import annotations
import requests
import json
import pandas as pd
from enum import Enum

class SessionError(Exception):
    """
    Error to be raised if the session cannot be requested
    """
    pass

class DataType(Enum):
    METADATA = 1
    PROTEINS = 2
    SIGNATURES = 3

class Signature():
    """
    A class representing a single signature
        :param data: (optional) parses data
    """
    def __init__(self, data: bytes = None) -> None:
        self.name: str = None
        self.source: str = None
        self.description: str = None
        self.genes: Dict[str, float] = [] # gene + importance

        if data:
            self.load_from_vision(data)

    def load_from_vision(self, data: dict) -> None:
        """
        Parses the vision json response into a signature 
        """
        self.name = data["name"]
        self.source = data["source"]
        self.description = data["metaData"]
        self.genes = data["geneImportance"]

    def __repr__(self) -> None:
        return f"Signature({self.name})"

class Cluster():
    """
    A class representing a cluster result
        :param data: (optional) parses data
    """
    def __init__(self, data_type: DataType, cluster_name: str, data: bytes = None) -> None:
        self.type: DataType = data_type
        self.cluster_name: str = cluster_name
        self.zscores: Dict[str, Dict[str, float]] = pd.DataFrame()
        self.pvalues: Dict[str, Dict[str, float]] = pd.DataFrame()

        if data:
            self.load_from_vision(data)

    def load_from_vision(self, data: dict) -> None:
        """
        Parses the vision json response into a cluster
        """
        temp = dict()
        for i, sample in enumerate(data["sig_labels"]):
            temp[sample] = pd.Series(data["zscores"][i], index=data["proj_labels"])
        self.zscores = pd.DataFrame(temp)
        
        temp = dict()
        for i, sample in enumerate(data["sig_labels"]):
            temp[sample] = pd.Series(data["pvals"][i], index=data["proj_labels"])
        self.pvalues = pd.DataFrame(temp)
    
    def __repr__(self) -> str:
        return f"Cluster({self.type.name.lower()} x {self.cluster_name})"

class LatentClass():
    """
    A class representing a result of laten class analysis (LCA)
        :param data: (optional) parses data
    """
    def __init__(self, data_type: DataType, data: bytes = None) -> None:
        self.type: DataType = data_type
        self.zscores: Dict[str, Dict[str, float]] = dict()
        self.pvalues: Dict[str, Dict[str, float]] = dict()

        if data:
            self.load_from_vision(data)

    def load_from_vision(self, data: dict) -> None:
        """
        Parses the vision json response into a cluster
        """
        temp = dict()
        for i, sample in enumerate(data["sig_labels"]):
            temp[sample] = pd.Series(data["zscores"][i], index=data["proj_labels"])
        self.zscores = pd.DataFrame(temp)
        
        temp = dict()
        for i, sample in enumerate(data["sig_labels"]):
            temp[sample] = pd.Series(data["pvals"][i], index=data["proj_labels"])
        self.pvalues = pd.DataFrame(temp)
    
    def __repr__(self) -> str:
        return f"LatentClass({self.type.name.lower()})"

class API():
    """
    Main API interface with the Vision session
        :param session: the link to the vision session
        :raises SessionError: if the session cannot be reached
    """
    def __init__(self, session: str) -> None:
        self.session: str = None

        # Make sure the session contains the correct protocol scheme
        if len(session.split("http://")) == 1:
            self.session = "http://" + session
        else:
            self.session = session

        # Check if website exists
        if not self.valid_session():
            raise SessionError(f"Session '{self.session}' is unavailable")

        # Initiate rest of class
        self.session_info: _SessionInfo = _SessionInfo(self.session)
        self._signatures: _Signatures = _Signatures(self.session)
        self._proteins: _Proteins = _Proteins(self.session)
        self._clusters: _Clusters = _Clusters(self.session)
        self._de: None = None
        self._lca: _LCA = _LCA(self.session)
        self._projections: _Projections = _Projections(self.session)
        self._tree: _Tree = _Tree(self.session)
        self._expression: _Expression = _Expression(self.session)
        self._analysis: None = None
        self._cell: _Cell = _Cell(self.session)
        self._cells: _Cells = _Cells(self.session)
        
    def valid_session(self) -> bool:
        """
        Checks if the sessions exists
            :returns: true is session is available
        """
        try:
            request = requests.head(self.session)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {self.session}") from error

        if request.status_code < 400:
            return True
        else:
            return False

    @property
    def signatures(self) -> _Signatures:
        """
        Getter for signatures
            :raises SessionError: when the current session has no signature info
        """
        if not self.session_info.has_sigs:
            raise SessionError(f"Session '{self.session}' contains no signature information")
        return self._signatures

    @property
    def proteins(self) -> _Proteins:
        """
        Getter for protein
            :raises SessionError: when the current session has no protein info
        """
        if not self.session_info.has_proteins:
            raise SessionError(f"Session '{self.session}' contains no protein information")
        return self._proteins

    @property
    def clusters(self) -> _Clusters:
        """
        Getter for clusters data
        """
        return self._clusters

    @property
    def de(self) -> None:
        """
        Hook for differential expression API. For now not implemented
        """
        raise NotImplementedError("Differential Expression (Yanay) API not implemented")

    @property
    def lca(self) -> _LCA:
        """
        Getter for pearson correlations (LCA)
            :raises SessionError: when the current session has no lca info
        """
        if not self.session_info.has_lca:
            raise SessionError(f"Session '{self.session}' contains no Pearson correlation/lca information")
        return self._lca

    @property
    def projections(self) -> _Projections:
        """
        Getter for projections data
        """
        return self._projections

    @property
    def tree(self) -> _Tree:
        """
        Getter for tree
            :raises SessionError: when the current session has no tree info
            :raises NotImplementedError: when it does
        """
        if not self.session_info.has_tree:
            raise SessionError(f"Session '{self.session}' contains no tree information")
        raise NotImplementedError(f"API for tree has not been implemented")
        return self._tree

    @property
    def expression(self) -> _Expression:
        """
        Getter for expression data
        """
        return self._expression

    @property
    def analysis(self) -> None:
        """
        Hook for analysis API. For now not implemented
        """
        raise NotImplementedError("Analysis API not implemented")

    @property
    def cell(self) -> _Cell:
        """
        Getter for cell metadata
        """
        return self._cell

    @property
    def cells(self) -> _Cells:
        """
        Getter for cell metadata
        """
        raise NotImplementedError(f"API for cells has not been implemented")
        return self._cells

    def __repr__(self) -> str:
        return f"API({self.session})"

class _SessionInfo():
    """
    Main interface with the VISION Session Info API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session
        # links to the internal Vision API
        self._session_info: str = "/SessionInfo"

        # parameters
        self.name: str = None
        self.ncells: int = None

        self.has_tree: bool = None
        self.pooled: bool = None
        self.has_sigs: bool = None
        self.has_proteins: bool = None
        self.has_lca: bool = None

        self.metadata: List[str] = []

        self.load_from_vision()

    def load_from_vision(self) -> None:
        """
        Requests the session info from vision
        """
        path = self.session + self._session_info
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        self.name = parse["name"]
        self.ncells = parse["ncells"]
        self.has_tree = parse["has_tree"]
        self.pooled = parse["pooled"]
        self.has_sigs = parse["has_sigs"]
        self.has_proteins = parse["has_proteins"]
        self.has_lca = parse["has_lca"]
        self.metadata = parse["meta_sigs"]
    
    def __repr__(self) -> str:
        return f"SessionInfo({self.session})"

class _Signatures():
    """
    Main interface with the VISION signature API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session
        # links to the internal Vision API
        # Signature API
        self._signature_info: str = "/Signature/Info/"               # + sig_name
        self._signature_scores: str = "/Signature/Scores/"           # + sig_name
        self._signature_expression: str = "/Signature/Expression/"   # + sig_name + "/" + cluster_var
        self._signature_cluster_normal: str = "/FilterGroup/SigClusters/Normal"

        # Following two entrees provide meta information, not signature specific, so ignore here
        self._signature_meta: str = "/Signature/Meta/"             # + meta_name 
        self._signature_cluster_meta: str = "/FilterGroup/SigClusters/Meta"

        # Data caches
        self._names: List[str] = []
        self._meta: List[str] = []

    @property
    def names(self) -> List[str]:
        """
        Getter for the signature names. First time runs a request, afterwards caches the result
        """
        if not self._names:
            path = self.session + self._signature_cluster_normal
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            self._names = list(parse.keys())

        return self._names

    @property
    def _metadata_names(self) -> List[str]:
        """
        Getter for the meta names. First time runs a request, afterwards caches the result
        """
        if not self._meta:
            path = self.session + self._signature_cluster_meta
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error
    
            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")
    
            self._meta = list(parse.keys())
    
        return self._meta

    def _cells_metadata(self, metadata_name: str) -> pd.Series:
        """
        Getter for the metadata label of all events/cells, sends a request for the specific info.
            :metadata_name: the name of the metadata as in self.meta
            :returns: a pd.Series of metadata_labels indexed on cell_id
        """
        path = self.session + self._signature_meta + metadata_name
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return pd.Series(parse["values"], index=parse["cells"])

    def signature(self, signature_name: str) -> Signature:
        """
        Getter for the metadata of a specific signature, sends a request for the specific info.
            :signature_name: the name of the signature as in self.names
        """
        path = self.session + self._signature_info + signature_name
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        output = Signature(parse)

        return output

    def score(self, signature_name: str) -> pd.Series:
        """
        Getter for the score/event of a specific signature, sends a request for the specific info.
            :signature_name: the name of the signature as in self.names
            :returns: a pd.Series of signature_scores indexed on cell_id
        """
        path = self.session + self._signature_scores + signature_name
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return pd.Series(parse["values"], index=parse["cells"])

    def expression(self, signature_name: str, meta_name: str) -> pd.DataFrame:
        """
        Not 100% sure what this expression score represents...
        :returns: a pd.DataFrame of score in a matrix of signature genes vs metadata_labels
        """
        path = self.session + self._signature_expression + signature_name + '/' + meta_name
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        output = dict()
        for i, sample in enumerate(parse["gene_labels"]):
            output[sample] = pd.Series(parse["data"][i], index=parse["sample_labels"])

        return pd.DataFrame(output)

    def __repr__(self) -> str:
        return f"Signatures({self.session})"

class _Proteins():
    """
    Main interface with the VISION Protein API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._protein_clusters: str = "/FilterGroup/SigClusters/Proteins"
        self._protein_values_a: str = "/Proteins/"
        self._protein_values_b: str = "/Values"

        self._names: List[str] = []

    @property
    def names(self) -> List[str]:
        """
        Getter for the protein names. First time runs a request, afterwards caches the result
        """
        if not self._names:
            path = self.session + self._protein_clusters
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            self._names = list(parse.keys())

        return self._names

    def value(self, protein_name: str) -> pd.Series:
        """
        Getter for the value/event of a specific protein, sends a request for the specific info.
            :protein_name: the name of the protein as in self.names
        """
        path = self.session + self._protein_values_a + protein_name + self._protein_values_b
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return pd.Series(parse["values"], index=parse["cells"])

    def __repr__(self) -> str:
        return f"Proteins({self.session})"

class _Clusters():
    """
    Main interface with the VISION Clusters API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._clusters_base: str = "/Clusters/" # + cluster_variable
        self._clusters_sigproj_matrix_meta: str = "/SigProjMatrix/Meta"     # meta to meta -> zscores/pvals
        self._clusters_sigproj_matrix_normal: str = "/SigProjMatrix/Normal" # meta to signature -> zscores/pvals
        self._clusters_protein_matrix: str = "/ProteinMatrix"               # meta to protein -> zscores/pvals
        self._clusters_cells: str = "/Cells"                                # meta_level/cell

        self._clusters_meta_levels: str = "/Clusters/MetaLevels"            # meta_levels
        self._clusters_list: str = "/Clusters/list"                         # meta_names

        # data caches
        self._meta: List[str] = []
        self._meta_levels: Dict[str, List[str]] = []

    @property
    def names(self) -> List[str]:
        """
        Getter for the cluster_names. First time runs a request, afterwards caches the result
        This is equivalent to the catagorize metadata_names, for the labels see Clusters.levels attribute
        """
        if not self._meta:
            path = self.session + self._clusters_list
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error
    
            parse = json.loads(request.content)

            self._meta = list(parse)
    
        return self._meta

    @property
    def levels(self) -> Dict[str, List[str]]:
        """
        Getter for the cluster levels of all cluster_names, sends a request for the specific info.
        This is equivalent to all labels from all catagorized metadata_levels. 
            :returns: a dictionary of all metadata names with a list of all levels of the specified entree
        """
        if not self._meta_levels:
            path = self.session + self._clusters_meta_levels
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            temp = dict()
            for meta in parse:
                temp[meta] = parse[meta]
            self._meta_levels = temp
        
        return self._meta_levels

    def __cells_metadata(self, cluster_name: str) -> pd.Series:
        """
        Getter for the cluster/metadata label of all events/cells, sends a request for the specific info.
            :cluster_name: the name of the cluster as in self.names
            :returns: a pd.Series of cluster_label indexed on cell_id
        """
        path = self.session + self._clusters_base + cluster_name + self._clusters_cells
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return pd.Series(parse["values"], index=parse["cells"])

    def proteins(self, cluster_name: str) -> Cluster:
        """
        Returns the zscores and pvals of the cluster_name clustering results over the protein data
            :cluster_name: the name of the cluster as in self.names
        """
        path = self.session + self._clusters_base + cluster_name + self._clusters_protein_matrix
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return Cluster(DataType.PROTEINS, cluster_name, parse)

    def signatures(self, cluster_name: str) -> Cluster:
        """
        Returns the zscores and pvals of the cluster_name clustering over the signatures data
            :cluster_name: the name of the cluster as in self.names
        """
        path = self.session + self._clusters_base + cluster_name + self._clusters_sigproj_matrix_normal
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return Cluster(DataType.SIGNATURES, cluster_name, parse)

    def metadata(self, cluster_name: str) -> Cluster:
        """
        Returns the zscores and pvals of the cluster_name clustering over the signature data
            :cluster_name: the name of the cluster as in self.names
        """
        path = self.session + self._clusters_base + cluster_name + self._clusters_sigproj_matrix_meta
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return Cluster(DataType.METADATA, cluster_name, parse)

    def __repr__(self) -> str:
        return f"Clusters({self.session})"

class _LCA():
    """
    Main interface with the VISION Clusters API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._lca_normal = "/PearsonCorr/Normal"
        self._lca_proteins = "/PearsonCorr/Proteins"
        self._lca_meta = "/PearsonCorr/Meta"

        # data caches
        self._proteins: LatentClass = None
        self._signatures: LatentClass = None
        self._metadata: LatentClass = None

    @property
    def proteins(self) -> LatentClass:
        """
        Returns the zscores and pvalues of the latent class analysis over the protein data
        """
        if not self._proteins:
            path = self.session + self._lca_proteins
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            self._proteins = LatentClass(DataType.PROTEINS, parse)
        return self._proteins

    @property
    def signatures(self) -> LatentClass:
        """
        Returns the zscores and pvalues of the latent class analysis over the signatures data
        """
        if not self._signatures:
            path = self.session + self._lca_normal
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            self._signatures = LatentClass(DataType.SIGNATURES, parse)
        return self._signatures

    @property
    def metadata(self) -> LatentClass:
        """
        Returns the zscores and pvalues of the latent class analysis over the metadata data
        """
        if not self._metadata:
            path = self.session + self._lca_meta
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            self._metadata = LatentClass(DataType.METADATA, parse)
        return self._metadata

    def __repr__(self) -> str:
        return f"LCA({self.session})"

class _Projections():
    """
    Main interface with the VISION Projections API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._projections_coordinates_a = "/Projections/" # + projection_name
        self._projections_coordinates_b = "/coordinates/" # + projeciton_column
        self._projections_list = "/Projections/list"

        # cached data
        self._levels: Dict[str, List[str]] = dict()

    @property
    def names(self) -> List[str]:
        """
        Getter for the cluster_names. First time runs a request, afterwards caches the result
        This is equivalent to the metadata_names, for the labels see Clusters.levels attribute
        """
        if not self._levels:
            _ = self.levels
    
        return list(self._levels.keys())

    @property
    def levels(self) -> Dict[str, List[str]]:
        """
        Getter for the cluster levels of all cluster_names, sends a request for the specific info.
        This is equivalent to all labels from all metadata_levels. 
            :returns: a dictionary of all metadata names with a list of all levels of the specified entree
        """
        if not self._levels:
            path = self.session + self._projections_list
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)
            try: 
                parse["error"]
            except KeyError:
                pass
            else:
                raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

            temp = dict()
            for meta in parse:
                temp[meta] = parse[meta]
            self._levels = temp
        
        return self._levels

    def projection(self, projection_name) -> pd.DataFrame:
        """
        Returns x-coordinates of all dimensions of the projection_name (all defined in self.names)
        :returns: A pd.DataFrame of the projection value of each dimension vs cell_id
        """
        path_a = self.session + self._projections_coordinates_a + projection_name + self._projections_coordinates_b
        columns = self.levels[projection_name]
        
        output = dict()
        for column in columns:
            path = path_a + column
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)

            output[column] = pd.Series([x[0] for x in parse], index=[x[1] for x in parse])

        return pd.DataFrame(output)

    def __repr__(self) -> str:
        return f"Projections({self.session})"

class _Expression():
    """
    Main interface with the VISION Expression API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._expression_gene = "/Expression/Gene/" # + gene_name
        self._expression_gene_list = "/Expression/Genes/List"

        # data cache
        self._names: List[str] = []

    @property
    def names(self) -> List[str]:
        """
        Returns a list of all gene names
        """
        if not self._names:
            path = self.session + self._expression_gene_list
            try:
                request = requests.get(path)
                request.raise_for_status()
            except Exception as error:
                raise SessionError(f"Unable to send request to {path}") from error

            parse = json.loads(request.content)

            self._names = parse
        return self._names
    
    def expression(self, gene_name) -> pd.Series:
        """
        Getter for the gene expression data events/cells, sends a request for the specific info.
            :gene_name: the name of the gene as in self.names
            :returns: a dictionary of each cell/event with the gene expression
        """
        path = self.session + self._expression_gene + gene_name
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return pd.Series(parse["values"], index=parse["cells"])

    def __repr__(self) -> str:
        return f"Expression({self.session})"

class _Cell():
    """
    Main interface with the VISION Cell API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._cell_meta_a = "/Cell/" # + cell_id
        self._cell_meta_b = "/Meta"

    def cell(self, cell_id) -> Dict[str, Any]:
        """
        Requests the cell meta information
        """
        path = self.session + self._cell_meta_a + cell_id + self._cell_meta_b
        try:
            request = requests.get(path)
            request.raise_for_status()
        except Exception as error:
            raise SessionError(f"Unable to send request to {path}") from error

        parse = json.loads(request.content)
        try: 
            parse["error"]
        except KeyError:
            pass
        else:
            raise SessionError(f"Unable to parse reply from {path}. Is the path correct?")

        return parse

    def __repr__(self) -> str:
        return f"Cell({self.session})"

# Unimplemented
class _Tree():
    """
    Main interface with the VISION Tree API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._tree_projections: str = "/Tree/Projections/coordinates"
        self._tree_projections_list: str = "/Tree/Projections/list"
        self._tree_sigproj_matrix_meta: str = "/Tree/SigProjMatrix/Meta"
        self._tree_sigproj_matrix_normal: str = "/Tree/SigProjMatrix/Normal"
        self._tree_protein_matrix: str = "/Tree/ProteinMatrix"

    def __repr__(self) -> str:
        return f"Tree({self.session})"

class _Cells():
    """
    Main interface with the VISION Cells API
        :param session: the link to the vision session
    """
    def __init__(self, session: str) -> None:
        self.session: str = session

        # links to the internal Vision API
        self._cells_meta = "/Cells/Meta"
        self._cells_save_selection = "/Cells/Selections/" # + selection_name
        self._cells_get_selection = "/Cells/Selections" # + selection_name
        self._cells_list_selection = "/Cells/Selections"

    def __repr__(self) -> str:
        return f"Cells({self.session})"
