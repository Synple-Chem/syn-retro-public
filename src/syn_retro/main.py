import argparse
import json
from pathlib import Path
from typing import List

import pandas as pd
from chembl_structure_pipeline import standardize_mol
from rdkit.Chem import MolFromSmiles, MolToSmiles

from syn_retro.path import DATA_PATH
from syn_retro.prod import apply_reboc, combine_retro_plans, return_1_step_retro_plan
from syn_retro.utils import RXN_SMARTS, connect_to_sqlite, dict_to_df


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-path",
        type=Path,
        required=True,
        help="path to the data csv file, should have columns: id, smiles",
    )
    parser.add_argument(
        "--save-path",
        type=Path,
        required=False,
        default=DATA_PATH / "test.json",
        help="path to save the results",
    )
    parser.add_argument("--db-path", type=Path, default=DATA_PATH / "test_bb.sqlite")

    return parser.parse_args()


def main():
    args = get_args()

    connection = connect_to_sqlite(cache_bb_db_path=args.db_path)

    rxn_defuilt_smarts = {k: v for k, v in RXN_SMARTS.items() if k not in ["re_boc"]}

    retro_plan_key_to_be_includes = ["rxn_nm", "bb_info", "complete"]

    df = pd.read_csv(args.data_path)  # df should have columns: id, smiles
    full_bb_dict = {}
    for smi_ii, smi in enumerate(df["smiles"].values):
        mol = MolFromSmiles(smi)
        std_mol = standardize_mol(mol)
        full_retro_plans: List = []
        # first step: fragment
        first_retro_plans = return_1_step_retro_plan(
            mol=std_mol, rxn_smarts=rxn_defuilt_smarts, connection=connection
        )
        for plan_1 in first_retro_plans:
            # plan = {"rxn_nm": rxn_nm, "bb_info": res, "frag_info": fragment, "complete": bool}
            if plan_1["complete"]:
                full_retro_plans.append(plan_1)
            else:
                # second step: re_boc
                not_found_lid = [
                    lid
                    for lid, bid in enumerate(plan_1["bb_info"]["parent_id"])
                    if bid is None
                ][0]
                mol_to_reboc = plan_1["frag_info"][not_found_lid]
                reboc_mols = apply_reboc(mol_to_reboc)
                for reboc_mol in reboc_mols:
                    # the third step: search again
                    third_retro_plans = return_1_step_retro_plan(
                        mol=reboc_mol,
                        rxn_smarts=rxn_defuilt_smarts,
                        connection=connection,
                    )

                    for plan_3 in third_retro_plans:
                        plan_combined = combine_retro_plans(
                            plan_1=plan_1, plan_3=plan_3
                        )
                        # append the full plan even if it is incomplete
                        full_retro_plans.append(plan_combined)

        full_bb_dict[smi_ii] = {
            "id": df["id"].values[smi_ii],
            "smiles": smi,
            "std_smiles": MolToSmiles(std_mol),
            "plan": [
                {k: v for k, v in p.items() if k in retro_plan_key_to_be_includes}
                for p in full_retro_plans
            ],
        }

    with open(args.save_path, "w") as fout:
        json_dumps_str = json.dumps(full_bb_dict, indent=4)
        print(json_dumps_str, file=fout)
    df = dict_to_df(full_bb_dict)
    df.to_csv(args.save_path.with_suffix(".csv"), index=False)


if __name__ == "__main__":
    main()
