from sqlite3 import Connection
from typing import Dict, List, Literal, Tuple

from rdkit.Chem import MolFromSmiles, MolToInchiKey, MolToSmiles, rdmolops

from syn_retro.syn_types import Mol
from syn_retro.utils import decompose_mol


def fragment_compound(
    mol: Mol, rxn_smarts: Dict
) -> Dict[str, List[Tuple[Mol, Mol | None]]]:
    """fragment a compound into building blocks
    Iteratively match reaction smarts to the molecule
    Then, return the unique fragments per reaction

    Args:
        mol (Mol): rdkit mol object
        rxn_smarts (Dict): dictionary of reaction smarts
            keys: reaction names
            values: reaction smarts

    Returns:
        Dict: dictionary of fragments per suceed reaction
            keys: reaction names
            values: list of tuples of fragments as Mol objects
                each tuple contains one or two fragments, depends on the reaction
                Also, the fragments are unique
    """
    fragment_dict: Dict = {}
    for rxn_nm, smarts in rxn_smarts.items():
        bbs_all = decompose_mol(smarts, mol)
        if len(bbs_all) > 0:
            # remove duplicates
            joined_bbs = []
            for bbs in bbs_all:
                try:
                    bb_pairs_smiles = [
                        MolToSmiles(rdmolops.RemoveHs(MolFromSmiles(MolToSmiles(bb))))
                        for bb in bbs
                    ]
                    joined_bbs.append(".".join(bb_pairs_smiles))
                except:
                    print(f"skipped {rxn_nm}'s building block")
            joined_bbs = list(set(joined_bbs))
            fragment_dict[rxn_nm] = []
            for joined_smiles in joined_bbs:
                fragment_dict[rxn_nm].append(
                    tuple([MolFromSmiles(bb) for bb in joined_smiles.split(".")])
                )
    return fragment_dict


def search_fragment(
    fragments: List[Tuple[Mol, Mol | None]],
    type: Literal["exact", "similarity"] = "exact",
    connection: Connection | None = None,
    reactant_names: List[Tuple] | Tuple | None = None,
) -> List[Dict]:
    """Search fragments wrapper to search in the building block database
    exact search: search fragments with the same inchikey-first-14-characters
    similarity search: search fragments with the similar smiles tanioto coefficient > 0.8

    Args:
        fragments (List[Tuple[Mol, Mol | None]]): list of tuples of fragments
            each tuple contains one or two fragments (Mol object), depends on the reaction
        type (Literal["exact", "similarity"], optional): type of search. Defaults to "exact".
            exact: search_exact_fragment_in_db()
            similarity: search_similar_fragment_in_db() (not implemented yet)
        connection (Connection | None, optional): sqlite connection. Defaults to None.
        reactant_names (List[Tuple] | Tuple | None, optional): List of reactant names or None
            when None, all reactants are searched
            when Tuple, the list of fragments are searched within the same reactant
            when List[Tuple], the list of fragments are searched within the corresponding reactants
                len(reactant_names) == len(fragments). Defaults to None.

    Returns:
        List[Dict]: list of dictionary of fragments
            keys: parent_id, smiles, inchikey
            values: List of building block information
                for example, when two fragments are generated from one reaction,
                    the list contains two elements.
                    e.g., parent_id: [bb1_id, bb2_id], smiles: [bb1_smiles, bb2_smiles]
                when one fragment is generated from one reaction,
                    the list contains one element.
                    e.g., parent_id: [bb1_id], smiles: [bb1_smiles]
                when no building block is found from db, None is returned.
                    e.g., parent_id: [None, bb2_id], smiles: [None, bb2_smiles]
            when wrong type is given, empty list is returned
    """
    all_frag_dict: List = []
    if type == "exact" and connection is not None:
        all_frag_dict = search_exact_fragment_in_db(
            connection=connection, fragments=fragments, reactant_names=reactant_names
        )
    elif type == "similarity":
        search_similar_fragment_in_db()
    else:
        raise (ValueError("Type should be either 'exact' or 'similarity'"))
    return all_frag_dict


def search_similar_fragment_in_db():
    # placeholder
    raise (NotImplementedError("Only exact search is implemented at the moment"))


def search_exact_fragment_in_db(
    connection: Connection,
    fragments: List[Tuple[Mol, Mol | None]],
    reactant_names: List[Tuple] | Tuple | None = None,
) -> List[Dict]:
    """search exact fragments in the building block database

    Args:
        connection (Connection): sqlite connection
        fragments (List[Tuple[Mol, Mol | None]]): list of tuples of fragments
            each tuple contains one or two fragments (Mol object), depends on the reaction
        reactant_names (List[Tuple] | Tuple | None): List of reactant names or None
            when None, all reactants are searched
            when Tuple, the list of fragments are searched within the same reactant
            when List[Tuple], the list of fragments are searched within the corresponding reactants
                len(reactant_names) == len(fragments)


    Returns:
        List[Dict]: list of dictionary of fragments
            keys: parent_id, smiles, inchikey
            values: List of building block information
                for example, when two fragments are generated from one reaction,
                    the list contains two elements.
                    e.g., parent_id: [bb1_id, bb2_id], smiles: [bb1_smiles, bb2_smiles]
                when one fragment is generated from one reaction,
                    the list contains one element.
                    e.g., parent_id: [bb1_id], smiles: [bb1_smiles]
                when no building block is found from db, None is returned.
                    e.g., parent_id: [None, bb2_id], smiles: [None, bb2_smiles]
    """
    # if isinstance(reactant_names, list):
    #     assert (
    #         len(reactant_names) == len(fragments),
    #         "When reactant_names is a list, the length of the list should be the same as the length of fragments",
    #     )
    if isinstance(reactant_names, list):
        table_names = reactant_names
    elif isinstance(reactant_names, tuple):
        table_names = [reactant_names] * len(fragments)
    else:
        table_names = [("building_blocks", "building_blocks")] * len(fragments)
    cursor = connection.cursor()

    # search bb per fragment
    all_frag_dict: List = []
    for frag, tnm in zip(fragments, table_names):
        frag_dict: Dict = {"parent_id": [], "smiles": [], "inchikey": []}
        for bb, t in zip(frag, tnm):
            bb_inchikey = MolToInchiKey(bb)
            bb_sch = bb_inchikey.split("-")[0]
            sch_res = cursor.execute(
                f"SELECT id, cano_smiles, inchikey FROM {t} WHERE parent_inchikey LIKE '%{bb_sch}%'"
            ).fetchone()
            if sch_res is not None:
                frag_dict["parent_id"].append(sch_res[0])
                frag_dict["smiles"].append(sch_res[1])
                frag_dict["inchikey"].append(sch_res[2])
            else:
                frag_dict["parent_id"].append(None)
                frag_dict["smiles"].append(None)
                frag_dict["inchikey"].append(None)
        all_frag_dict.append(frag_dict)

    return all_frag_dict


def return_1_step_retro_plan(
    mol: Mol,
    rxn_smarts: Dict,
    connection: Connection,
    reactant_dict: Dict | None = None,
) -> List[Dict]:
    """Return 1-step retro plan of a molecule for given reactions

    Args:
        mol (Mol): rdkit mol object
        rxn_smarts (Dict): dictionary of reaction smarts
            keys: reaction names
            values: reaction smarts
        connection (Connection): sqlite connection to bb information
        reactant_dict (Dict): dictionary of reactant names, default None

    Returns:
        List[Dict]: list of dictionary of found building blocks
        when all building blocks are found "complete" is True
            each dict: {"rxn_nm": rxn_nm, "bb_info": bb_info, frag_info: [frag_mol, frag_mol|None], "complete": bool}
            bb_info: {
                "parent_id": [000, None] | [int or None],
                "smiles": [ZZZ, None] | [str or None],
                "inchikey": [III, None] | [str or None]
            }
    """
    fragment_dict = fragment_compound(mol=mol, rxn_smarts=rxn_smarts)
    retro_plans: List = []
    for rxn_nm, fragments in fragment_dict.items():
        if len(fragments) > 0:
            # for each reaction, returned fragments can be more than one
            # thus, returned search_info can include more than one elements
            reactants = tuple(reactant_dict[rxn_nm]) if reactant_dict else None
            search_info = search_fragment(
                connection=connection, fragments=fragments, reactant_names=reactants
            )
            for ii, res in enumerate(search_info):
                found_bb = [bb is not None for bb in res["parent_id"]]
                if any(found_bb):
                    plan_dict = {
                        "rxn_nm": rxn_nm,
                        "bb_info": res,
                        "frag_info": fragments[ii],
                        "complete": False,
                    }
                    if all(found_bb):
                        plan_dict["complete"] = True
                    retro_plans.append(plan_dict)
    return retro_plans


def combine_retro_plans(plan_1: Dict, plan_3: Dict) -> Dict:
    """Combine two retro plans into one
    Sequence is always 1 -> 2 (re_boc) -> 3
    Input plans can be incomplete or complete
    plan = {"rxn_nm": rxn_nm, "bb_info": res, "frag_info": fragment, "complete": bool}

    Args:
        plan_1 (Dict): Dictionary of first retro step plan, should be incomplete
            keys: rxn_nm, bb_info, frag_info, complete
        plan_3 (Dict): Dictionary of third retro step plan, should be complete or incomplete
            keys: rxn_nm, bb_info, frag_info, complete

    Returns:
        Dict: full retro plan
            keys: rxn_nm, bb_info, frag_info, complete
    """
    assert not plan_1["complete"], "plan_1 should be incomplete"
    found_lid = [
        lid for lid, bid in enumerate(plan_1["bb_info"]["parent_id"]) if bid is not None
    ]

    full_plan = {
        "rxn_nm": [plan_1["rxn_nm"], "re_boc", plan_3["rxn_nm"]],
        "bb_info": {
            k: [
                *[plan_1["bb_info"][k][lid] for lid in found_lid],
                *plan_3["bb_info"][k],
            ]
            for k in plan_1["bb_info"].keys()
        },
        "frag_info": [
            [plan_1["frag_info"][lid] for lid in found_lid],
            plan_3["frag_info"],
        ],
        "complete": True if plan_3["complete"] else False,
    }
    return full_plan
