##############################################################################     ##    ######
#    A.J. Zwijnenburg                   2020-09-20           v1.1                 #  #      ##
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
Wraps the api in a convenience class for data management. Automatically requests
and caches the necessary data

:class: Data
Main storage class of a VISION dataset, provides:
.api            - a VISION API
.cells          - a list of all cell ids
.expression     - hook into the expression data of the dataset
.protein        - hook into the protein data of the dataset
.signature      - hook into the siganture data of the dataset
.meta           - hook into the metadata of the dataset
.projection     - hook into the projection data of the dataset

:class: _Expression
Storage class of expression data, provides:
.names      -   a list of all gene names
.data       -   a dataframe of all expression data indexed in cells vs genes (very slow)
.[]         -   provides access to the expression data of a (list of) gene(s)

:class: _Protein
Storage class of protein data, provides:
.names      -   a list of all gene names
.data       -   a dataframe of all protein data indexed in cells vs proteins (very slow)
.[]         -   provides access to the protein data of a (list of) protein(s)

:class: _Signature
Storage class of signature data, provides:
.names      -   a list of all signature names
.data       -   a dataframe of all signature data indexed in cells vs signature (very slow)
.[]         -   provides access to the signature scores of a (list of) signature(s)

:class: _Metadata
Storage class of metadata data, provides:
.names      -   a list of all gene names
.levels     -   a dictionary of the levels off all catagorical metadata parameters
.data       -   a dataframe of all metadata indexed in cells vs metadata (very slow)
.[]         -   provides access to the metadata of a (list of) metadata(s)

:class: _Projection
Storage class of projection data, provides:
.names      -   a list of all projections names
.dimension  -   provides a dictionary off projections to all dimension labels of the projection
.[]         -   a dataframe of a projection. The data indexed in cells vs projection dimension

"""

from __future__ import annotations
from .api import API

import pandas as pd

class Data:
    """
    Wraps the VISION api for convenient plotting
        :param session: the VISION session link
    """
    def __init__(self, session: str):
        self._api: API = API(session)

        self._cells: List[str] = None

        self._meta: _Metadata = _Metadata(self)
        self._expression: _Expression = _Expression(self)
        self._protein: _Protein = _Protein(self)
        self._signature: _Signature = _Signature(self)
        self._projection: _Projection = _Projection(self)

    @property
    def api(self) -> API:
        """
        Getter for the internal API
        """
        return self._api

    @property
    def cells(self) -> List[str]:
        """
        Getter for gene names
        """
        if not self._cells:
            temp_gene = self.expression.names[0]
            temp: pd.DataFrame = self._api.expression.expression(temp_gene)
            self._cells = list(temp.index)

        return self._cells

    @property
    def meta(self) -> _Metadata:
        """
        Getter for metadata
        """
        return self._meta

    @property
    def expression(self) -> _Expression:
        """
        Getter for expression
        """
        return self._expression

    @property
    def protein(self) -> _Protein:
        """
        Getter for protein
        """
        return self._protein
    
    @property
    def signature(self) -> _Signature:
        """
        Getter for protein
        """
        return self._signature

    @property
    def projection(self) -> _Projection:
        """
        Getter for projection
        """
        return self._projection
    
    def __getitem__(self, columns: Union[str, List[str]]) -> pd.DataFrame:
        """
        Looks through all data to find all columns of interest
        """
        if isinstance(columns, str):
            column = columns

            if column in self.expression.names:
                return self.expression[column]

            elif column in self.meta.names:
                return self.meta[column]

            elif column in self.protein.names:
                return self.protein[column]

            elif column in self.signature.names:
                return self.signature[column]

            else:
                raise ValueError(f"column parameter '{column}' couldnt be found in the dataset")

        elif isinstance(columns, list):
            warnings: List[str] = []
            output: pd.DataFrame = pd.DataFrame()
            for column in columns:
                try:
                    output[column] = self.__getitem__(column).iloc[:,0]
                except ValueError:
                    warnings.append(column)

            if warnings:
                raise ValueError(f"column parameters [{', '.join(warnings)}] couldnt be found in the dataset")

            return output
        else:
            raise TypeError(f"cannot extract expression with proteins defined as type {type(columns)}")

    def __contains__(self) -> None:
        """
        Accidental containment testing gives weird errors so raise error
        """
        raise AttributeError("for membership checking please use the names attribute of the relevant class")

    def __repr__(self) -> str:
        return f"Data({self.api.session})"

class _Expression:
    """
    Representation of expression data. Tries to limit io by caching data
        :param api: connects with API
    """
    def __init__(self, parent: Data):
        self._parent: Data = parent

        self._names: List[str] = None
        self._cache: pd.DataFrame = pd.DataFrame()

    @property
    def names(self) -> List[str]:
        """
        Getter for gene names
        """
        if not self._names:
            self._names = self._parent.api.expression.names()

        return self._names

    @property
    def data(self) -> pd.DataFrame:
        """
        Returns ALL expression data. This causes a lot of requests so use only when necessary
        Very slow -> pre-allocating enough memory for cache will likely improve performance a lot.
        """
        for gene in self.names:
            if not gene in self._cache:
                counts = self._parent.api.expression.expression(gene)
                self._cache[gene] = counts
        return self._cache

    def __len__(self):
        return len(self._names)

    def __getitem__(self, genes: Union[str, List[str]]) -> pd.DataFrame:
        """
        Gets the specified gene(s). Genes must be named (for proper requesting).
        If you want to see the full data set use the .data attribute which forces
        loading of all expression data
        """
        if isinstance(genes, str):
            gene = genes
            if gene not in self.names:
                raise ValueError(f"gene '{gene}' is not available in the data")

            # Check if gene is cached
            if not gene in self._cache:
                counts = self._parent.api.expression.expression(gene)
                self._cache[gene] = counts

            return self._cache.loc[:, [gene]]

        elif isinstance(genes, list):
            if not genes:
                raise ValueError("please define the genes you want to extract")

            for gene in genes:
                if not isinstance(gene, str):
                    raise ValueError("can only extract gene names (str)")

                if gene not in self.names:
                    raise ValueError(f"gene '{gene}' is not available in the data")

                # Check if gene is cached
                if not gene in self._cache:
                    counts = self._parent.api.expression.expression(gene)
                    self._cache[gene] = counts

            return self._cache[genes]
        else:
            raise TypeError(f"cannot extract expression with genes defined as type {type(genes)}")

    def __contains__(self, item) -> None:
        """
        Accidental containment testing gives weird errors so raise error
        """
        raise AttributeError("for membership checking please use the names class attribute")

class _Protein:
    """
    Representation of protein data. Tries to limit io by caching data
        :param api: connects with API
    """
    def __init__(self, parent: Data):
        self._parent: Data = parent

        self._names: List[str] = None
        self._cache: pd.DataFrame = pd.DataFrame()

    @property
    def names(self) -> List[str]:
        """
        Getter for protein names
        """
        if not self._names:
            self._names = self._parent.api.proteins.names()

        return self._names

    @property
    def data(self) -> pd.DataFrame:
        """
        Returns ALL protein data. This causes a lot of requests so use only when necessary.
        Very slow -> pre-allocating enough memory for cache will likely improve performance a lot.
        """
        for protein in self.names:
            if not protein in self._cache:
                counts = self._parent.api.proteins.value(protein)
                self._cache[protein] = counts
        return self._cache

    def __len__(self):
        return len(self._names)

    def __getitem__(self, proteins: Union[str, List[str]]) -> pd.DataFrame:
        """
        Gets the specified protein(s). Protein must be named (for proper requesting).
        If you want to see the full data set use the .data attribute which forces
        loading of all protein data
        """
        if isinstance(proteins, str):
            protein = proteins
            if protein not in self.names:
                raise ValueError(f"protein '{protein}' is not available in the data")

            # Check if protein is cached
            if not protein in self._cache:
                counts = self._parent.api.proteins.value(protein)
                self._cache[protein] = counts

            return self._cache.loc[:, [protein]]

        elif isinstance(proteins, list):
            if not proteins:
                raise ValueError("please define the proteins you want to extract")

            for protein in proteins:
                if not isinstance(protein, str):
                    raise ValueError("can only extract protein names (str)")

                if protein not in self.names:
                    raise ValueError(f"protein '{protein}' is not available in the data")

                # Check if protein is cached
                if not protein in self._cache:
                    counts = self._parent.api.proteins.value(protein)
                    self._cache[protein] = counts

            return self._cache[proteins]
        else:
            raise TypeError(f"cannot extract expression with proteins defined as type {type(proteins)}")

    def __contains__(self, item) -> None:
        """
        Accidental containment testing gives weird errors so raise error
        """
        raise AttributeError("for membership checking please use the names class attribute")

class _Signature:
    """
    Representation of signature data. Tries to limit io by caching data
        :param api: connects with API
    """
    def __init__(self, parent: Data):
        self._parent: Data = parent

        self._names: List[str] = None
        self._cache: pd.DataFrame = pd.DataFrame()

    @property
    def names(self) -> List[str]:
        """
        Getter for signature names
        """
        if not self._names:
            self._names = self._parent.api.signatures.names()

        return self._names

    @property
    def data(self) -> pd.DataFrame:
        """
        Returns ALL signature data. This causes a lot of requests so use only when necessary.
        Very slow -> pre-allocating enough memory for cache will likely improve performance a lot.
        """
        for signature in self.names:
            if not signature in self._cache:
                counts = self._parent.api.signatures.score(signature)
                self._cache[signature] = counts
        return self._cache

    def __len__(self):
        return len(self._names)

    def __getitem__(self, signatures: Union[str, List[str]]) -> pd.DataFrame:
        """
        Gets the specified signature(s). Signatures must be named (for proper requesting).
        If you want to see the full data set use the .data attribute which forces
        loading of all signature data
        """
        if isinstance(signatures, str):
            signature = signatures
            if signature not in self.names:
                raise ValueError(f"signature '{signature}' is not available in the data")

            # Check if signature is cached
            if not signature in self._cache:
                counts = self._parent.api.signatures.score(signature)
                self._cache[signature] = counts

            return self._cache.loc[:, [signature]]

        elif isinstance(signatures, list):
            if not signatures:
                raise ValueError("please define the signatures you want to extract")

            for signature in signatures:
                if not isinstance(signature, str):
                    raise ValueError("can only extract signature names (str)")

                if signature not in self.names:
                    raise ValueError(f"signature '{signature}' is not available in the data")

                # Check if signature is cached
                if not signature in self._cache:
                    counts = self._parent.api.signatures.score(signature)
                    self._cache[signature] = counts

            return self._cache[signatures]

        else:
            raise TypeError(f"cannot extract expression with signatures defined as type {type(signatures)}")

    def __contains__(self, item) -> None:
        """
        Accidental containment testing gives weird errors so raise error
        """
        raise AttributeError("for membership checking please use the names class attribute")

class _Metadata:
    """
    Representation of metadata data. Tries to limit io by caching data
        :param api: connects with API
    """
    def __init__(self, parent: Data):
        self._parent: Data = parent

        self._names: List[str] = None
        self._levels: Dict[str, List[str]] = None
        self._cache: pd.DataFrame = pd.DataFrame()

    @property
    def names(self) -> List[str]:
        """
        Getter for metadata names
        """
        if not self._names:
            self._names = list(self._parent.api.cell.cell(self._parent.api.expression.names()[0]).keys())

        return self._names

    @property
    def levels(self) -> Dict[str, List[str]]:
        """
        Getter for the levels of catagorized meta names
        """
        if not self._levels:
            self._levels = self._parent.api.clusters.levels()
        return self._levels

    @property
    def data(self) -> pd.DataFrame:
        """
        Returns meta data of all cells. This causes a lot of requests so use only when necessary
        Very slow -> pre-allocating enough memory for cache will likely improve performance a lot.
        """
        for meta in self.names:
            if not meta in self._cache:
                self._cache[meta] = self._parent.api.signatures._cells_metadata(meta)

        return self._cache

    def __len__(self):
        return len(self._names)

    def __getitem__(self, meta: Union[str, List[str]]) -> pd.DataFrame:
        """
        Gets the specified metadata. metadata must be named (for proper requesting).
        If you want to see the full data set use the .data attribute which forces
        loading of all metadata data
        """

        if isinstance(meta, str):
            if meta not in self.names:
                raise ValueError(f"meta '{meta}' is not available in the data")

            # Check if meta is cached
            if not meta in self._cache:
                counts = self._parent.api.signatures.metadata(meta)
                self._cache[meta] = counts

            return self._cache.loc[:, [meta]]

        elif isinstance(meta, list):
            if not meta:
                raise ValueError("please define the meta you want to extract")

            for label in meta:
                if not isinstance(label, str):
                    raise ValueError("can only extract metadata names of type str")

                if label not in self.names:
                    raise ValueError(f"metadata label '{label}' is not available in the data")

                # Check if label is cached
                if not label in self._cache:
                    counts = self._parent.api.signatures.metadata(label)
                    self._cache[label] = counts

            return self._cache[meta]
        else:
            raise TypeError(f"cannot extract expression with meta defined as type {type(meta)}")

    def __contains__(self, item) -> None:
        """
        Accidental containment testing gives weird errors so raise error
        """
        raise AttributeError("for membership checking please use the names/levels class attribute")

class _Projection:
    """
    Representation of projection data. Tries to limit io by caching data
        :param api: connects with API
    """
    def __init__(self, parent: Data):
        self._parent: Data = parent

        self._names: List[str] = None
        self._dimensions: Dict[str, List[str]] = None
        self._cache: Dict[str, pd.DataFrame] = dict()

    @property
    def names(self) -> List[str]:
        """
        Getter for gene names
        """
        if not self._names:
            self._names = self._parent.api.projections.names()

        return self._names

    @property
    def dimensions(self) -> Dict[str, List[str]]:
        """
        Getter for the dimensions of projections
        """
        if not self._dimensions:
            self._dimensions = self._parent.api.projections.levels()
        return self._dimensions

    def __len__(self):
        return len(self._names)

    def __getitem__(self, projection: str) -> pd.DataFrame:
        """
        Gets the specified projection. Projection must be named (for proper requesting).
        """

        if isinstance(projection, str):
            if projection not in self.names:
                raise ValueError(f"projection '{projection}' is not available in the data")

            # Check if meta is cached
            if not projection in self._cache:
                # Get all levels of the projection
                levels = self.dimensions[projection]
                data: List[pd.Series] = []
                for level in levels:
                    data.append(self._parent.api.projections.projection(projection, level))

                data = pd.concat(data, axis=1)
                data.columns = levels
                self._cache[projection] = data

            return self._cache[projection]

        else:
            raise TypeError(f"cannot extract projections with projection defined as type {type(meta)}, must be str")

    def __contains__(self, item) -> None:
        """
        Accidental containment testing gives weird errors so raise error
        """
        raise AttributeError("for membership checking please use the names/dimensions class attribute")
