from pathlib import Path
from typing import Dict, List

import yaml  # type: ignore
from pytest import fixture


# define fixtures here
@fixture
def prod_dict() -> Dict[str, List]:
    return {
        "bb1_parent_id": [507610, 507610, 527216],
        "bb2_parent_id": [439909, 479992, 530941],
        "bb3_parent_id": [29855, 475850, 533932],
        "smiles": [
            "CCC(C(C)C)N1CCC[C@@H]1CN(C)CC(O)c1ccc(O)cc1",
            "CCC(C(C)C)N1[C@H](CN2CC(C)OC(C)C2)COC1(C)C",
            "CC(CCCO)N1CCN(c2ccc(CN3CCN(c4ccc(C(F)(F)F)cc4)CC3)cc2)CC1",
        ],
    }


@fixture
def test_resources_path() -> Path:
    return Path(__file__).parent / "resources"


@fixture
def bb_csv_path(test_resources_path: Path):
    return test_resources_path / "test_bb.csv"


@fixture
def sub_smarts_path(test_resources_path: Path):
    return test_resources_path / "test_sub_smarts.yaml"


@fixture
def rxn_smarts_path(test_resources_path: Path):
    return test_resources_path / "test_rxn_smarts.yaml"


@fixture
def bb_bad_csv_path(test_resources_path: Path):
    return test_resources_path / "test_bad_bb.csv"


@fixture
def reactant_smarts(sub_smarts_path: Path) -> Dict:
    with open(sub_smarts_path, "r") as f:
        smarts = yaml.safe_load(f)
    return smarts


@fixture
def rxn_smarts(rxn_smarts_path: Path) -> Dict:
    with open(rxn_smarts_path, "r") as f:
        smarts = yaml.safe_load(f)
    return smarts["retro-rxn"]


@fixture
def reactant_dict(rxn_smarts_path: Path) -> Dict:
    with open(rxn_smarts_path, "r") as f:
        smarts = yaml.safe_load(f)
    reactant_dict: Dict = {}
    for k, v in smarts["reactants"].items():
        reactant_dict[k] = v.split(",")
    return reactant_dict
