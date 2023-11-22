import sqlite3
from pathlib import Path
from sqlite3 import Connection
from typing import Dict, List, Tuple

import pandas as pd
import yaml  # type: ignore
from rdkit import Chem  # type: ignore
from rdkit.Chem import Draw, rdChemReactions, rdmolops

from syn_retro.path import ASSET_PATH, DATA_PATH
from syn_retro.syn_types import Mol


def sanitize_mol(mol: Mol) -> Mol:
    """sanitize mol object
    Mol object should be sanitized after the reaction to avoid error in the next reaction

    Args:
        mol (Mol): rdkit mol object

    Returns:
        Mol: sanitized mol object
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(
        mol,
        Chem.SanitizeFlags.SANITIZE_FINDRADICALS
        | Chem.SanitizeFlags.SANITIZE_KEKULIZE
        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
        | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
        catchErrors=True,
    )
    return mol


def decompose_mol(rxn_smarts: str, mol: Mol) -> List[Tuple]:
    """decompose mol with reaction smarts

    Args:
        rxn_smarts (str): reaction smarts
        mol (Mol): rdkit mol object

    Returns:
        List[Tuple]: List of tuples of rdkit mol objects
            each tuple is a set of fragments
    """
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
    bbs = rxn.RunReactants((sanitize_mol(rdmolops.AddHs(mol)),))
    return bbs


def get_retro_rxn_smarts(
    asset_path: Path = ASSET_PATH / "rxn_smarts.yaml",
) -> Dict[str, str]:
    """get reaction smarts for retrosynthesis saved in asset file

    Args:
        asset_path (Path, optional): path to asset file. This yaml file contains the 'retro-rxn' section.
            It includes retro-reaction smarts. Defaults to ASSET_PATH/"rxn_smarts.yaml".

    Returns:
        Dict: dictionary of reaction name as keys and reaction smarts as values
            e.g., {"suzuki": "smarts"}
    """
    assert asset_path.exists(), f"{asset_path} does not exist"
    with open(asset_path, "r") as f:
        rxn_smarts = yaml.safe_load(f)
    return rxn_smarts["retro-rxn"]


def get_rxn_reactants(
    asset_path: Path = ASSET_PATH / "rxn_smarts.yaml",
) -> Dict[str, List[str]]:
    """get reactant names for each reaction

    Args:
        asset_path (Path, optional): path to asset file. This yaml file contains the 'reactants' section.
            It includes reactant names for each reaction. Defaults to ASSET_PATH/"rxn_smarts.yaml".

    Returns:
        Dict: dictionary of reaction name as keys and list of reactant names as values
            e.g., {"suzuki": ["boronic acid", "halide"]}
    """
    assert asset_path.exists(), f"{asset_path} does not exist"
    with open(asset_path, "r") as f:
        rxn_smarts = yaml.safe_load(f)
    reactant_dict: Dict = {}
    for k, v in rxn_smarts["reactants"].items():
        reactant_dict[k] = v.split(",")
    return reactant_dict


def get_substrate_smarts(
    asset_path: Path = ASSET_PATH / "sub_smarts.yaml",
) -> Dict[str, str]:
    """get substrate smarts for the outcome of the retrosynthesis

    Args:
        asset_path (Path, optional): path to asset file. Defaults to ASSET_PATH/"sub_smarts.yaml".
            this file contains the smarts for the substrates of each reaction followed by the substrate name.

    Returns:
        Dict[str, str]: dictionary of substrate name and smarts
    """
    assert asset_path.exists(), f"{asset_path} does not exist"
    with open(asset_path, "r") as f:
        smarts = yaml.safe_load(f)
    return smarts


def connect_to_sqlite(
    cache_bb_db_path: Path,
) -> Connection:
    """get sqlite connection and load chemicalite
    https://github.com/rvianello/chemicalite/blob/master/docs/tutorial_1st.rst

    Args:
        cache_bb_db_path (Path, optional): path to cached sqlite db. Defaults to DATA_PATH/"vlib_bb_cache_unwanted.db".

    Returns:
        Connection: sqlite connection
    """
    assert cache_bb_db_path.exists(), f"{cache_bb_db_path} does not exist"
    connection = sqlite3.connect(str(cache_bb_db_path))
    connection.enable_load_extension(True)
    connection.load_extension("chemicalite")
    connection.enable_load_extension(False)
    return connection


def plot_cmpds(
    mols: List[Mol], legends: List[str] | None = None, fignm: str = "diversity_pick.png"
):
    """Plot mols as a grid image of compounds

    Args:
        mols (List[Mol]): List of rdkit mol objects
        legends (List[str] | None, optional): legend of each molecule. Defaults to None.
        fignm (str, optional): file name od the saved fiture. Defaults to "diversity_pick.png".
    """
    assert legends is None or len(legends) == len(
        mols
    ), "legends should be None or the same length as mols"
    if legends is None:
        legends = [str(ii) for ii in range(len(mols))]
    # save img of mols
    img = Draw.MolsToGridImage(mols, molsPerRow=10, legends=legends)
    img.save(DATA_PATH / fignm)


def dict_to_df(full_retro_dict: Dict) -> pd.DataFrame:
    """convert full retro dict to dataframe
    dictionary should have the following format:
    {smi_ii: {
        "id": id,
        "smiles": smiles,
        "std_smiles": std_smiles,
        "plan": [
            {
                "rxn_nm": rxn_nm,
                "bb_info": bb_info,
                "complete": bool
            }
        ]
    }}

    Args:
        full_retro_dict (Dict): Dictionary of full retrosynthesis results

    Returns:
        pd.DataFrame: df with the following columns:
            id, smiles, std_smiles, complete, rxn_1, rxn_2, rxn_3,
            bb1_smiles, bb2_smiles, bb3_smiles, bb1_parent_id, bb2_parent_id, bb3_parent_id,
            bb1_inchikey, bb2_inchikey, bb3_inchikey
    """
    df = pd.DataFrame(
        columns=[
            "id",
            "smiles",
            "std_smiles",
            "complete",
            "rxn_1",
            "rxn_2",
            "rxn_3",
            "bb1_smiles",
            "bb2_smiles",
            "bb3_smiles",
            "bb1_parent_id",
            "bb2_parent_id",
            "bb3_parent_id",
            "bb1_inchikey",
            "bb2_inchikey",
            "bb3_inchikey",
        ]
    )
    row_num = 0
    for smi_ii in full_retro_dict.keys():
        compound = full_retro_dict[smi_ii]
        for plan in compound["plan"]:
            df.loc[row_num, "id"] = compound["id"]
            df.loc[row_num, "smiles"] = compound["smiles"]
            df.loc[row_num, "std_smiles"] = compound["std_smiles"]
            df.loc[row_num, "complete"] = plan["complete"]
            if isinstance(plan["rxn_nm"], list):
                # when three step retrosynthesis
                df.loc[row_num, "rxn_1"] = plan["rxn_nm"][-1]
                df.loc[row_num, "rxn_2"] = plan["rxn_nm"][1]
                df.loc[row_num, "rxn_3"] = plan["rxn_nm"][0]
                for bb_info_type in ["smiles", "parent_id", "inchikey"]:
                    df.loc[row_num, f"bb3_{bb_info_type}"] = plan["bb_info"][
                        bb_info_type
                    ][0]
                    df.loc[row_num, f"bb1_{bb_info_type}"] = plan["bb_info"][
                        bb_info_type
                    ][1]
                    df.loc[row_num, f"bb2_{bb_info_type}"] = plan["bb_info"][
                        bb_info_type
                    ][2]
            else:
                df.loc[row_num, "rxn_1"] = plan["rxn_nm"]
                for bb_info_type in ["smiles", "parent_id", "inchikey"]:
                    df.loc[row_num, f"bb1_{bb_info_type}"] = plan["bb_info"][
                        bb_info_type
                    ][0]
                    df.loc[row_num, f"bb2_{bb_info_type}"] = plan["bb_info"][
                        bb_info_type
                    ][1]

            row_num += 1

    return df
