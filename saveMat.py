from __future__ import absolute_import
import re
from collections import OrderedDict
from uuid import uuid4
from warnings import warn

import cobra
import scipy
from cobra.io.mat import _cell, create_mat_metabolite_id
from numpy import array, inf, isinf
from scipy.sparse import coo_matrix
from six import string_types

from cobra.core import Metabolite, Model, Reaction
from cobra.util import create_stoichiometric_matrix
from cobra.util.solver import set_objective


try:
    from scipy import io as scipy_io
    from scipy import sparse as scipy_sparse
except ImportError:
    scipy_sparse = None
    scipy_io = None


def new_load_matlab_model(infile_path, variable_name=None):
    data = scipy.io.loadmat(infile_path)
    possible_names = []
    meta_vars = {"__globals__", "__header__", "__version__"}
    possible_names = sorted(i for i in data if i not in meta_vars)
    if len(possible_names) == 1:
        variable_name = possible_names[0]
    if variable_name is not None:
        return new_from_mat_struct(data[variable_name], model_id=variable_name)

def new_save_matlab_model(model, file_name, varname=None):
    if not scipy_io:
        raise ImportError("load_matlab_model requires scipy")

    if varname is None:
        varname = (
            str(model.id)
            if model.id is not None and len(model.id) > 0
            else "exported_model"
        )
    mat = new_create_mat_dict(model)
    scipy_io.savemat(file_name, {varname: mat}, appendmat=True, oned_as="column")

def new_create_mat_dict(model):
    """Create a dict mapping model attributes to arrays."""
    rxns = model.reactions
    mets = model.metabolites
    mat = OrderedDict()
    mat["mets"] = _cell([met_id for met_id in create_mat_metabolite_id(model)])
    mat["metNames"] = _cell(mets.list_attr("name"))
    mat["metFormulas"] = _cell([str(m.formula) for m in mets])
    try:
        mat["metCharge"] = array(mets.list_attr("charge")) * 1.0
    except TypeError:
        # can't have any None entries for charge, or this will fail
        pass
    mat["genes"] = _cell(model.genes.list_attr("id"))
    # make a matrix for rxnGeneMat
    # reactions are rows, genes are columns
    rxn_gene = scipy_sparse.dok_matrix((len(model.reactions), len(model.genes)))
    if min(rxn_gene.shape) > 0:
        for i, reaction in enumerate(model.reactions):
            for gene in reaction.genes:
                rxn_gene[i, model.genes.index(gene)] = 1
        mat["rxnGeneMat"] = rxn_gene
    mat["grRules"] = _cell(rxns.list_attr("gene_reaction_rule"))
    mat["rxns"] = _cell(rxns.list_attr("id"))
    mat["rxnNames"] = _cell(rxns.list_attr("name"))
    mat["subSystems"] = _cell(rxns.list_attr("subsystem"))
    mat["csense"] = str("E" * len(model.metabolites))
    stoich_mat = create_stoichiometric_matrix(model)
    coo_stoich_mat = coo_matrix(stoich_mat)
    mat["S"] = coo_stoich_mat if coo_stoich_mat is not None else [[]]
    # multiply by 1 to convert to float, working around scipy bug
    # https://github.com/scipy/scipy/issues/4537
    mat["lb"] = array(rxns.list_attr("lower_bound")) * 1.0
    mat["ub"] = array(rxns.list_attr("upper_bound")) * 1.0
    mat["b"] = array(mets.list_attr("_bound")) * 1.0
    mat["c"] = array(rxns.list_attr("objective_coefficient")) * 1.0
    mat["rev"] = array(rxns.list_attr("reversibility")) * 1
    mat["description"] = str(model.id)
    return mat

def new_from_mat_struct(mat_struct, model_id=None):
    m = mat_struct
    if not {"rxns", "mets", "S", "lb", "ub"} <= set(m.dtype.names):
        raise ValueError("not a valid mat struct")
    if "c" in m.dtype.names:
        c_vec = m["c"][0, 0]
    model = Model()
    if model_id is not None:
        model.id = model_id
    for i, name in enumerate(m["mets"][0, 0]):
        print(i, str(name[0][0]))

if __name__ == '__main__':
    model = cobra.io.load_matlab_model('Model_iJB785_noSpace.mat')
    new_save_matlab_model(model, "temp.mat")