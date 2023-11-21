import argparse
import logging
from pathlib import Path
from sqlite3 import Connection
from typing import Dict

import pandas as pd
from chembl_structure_pipeline import get_parent_mol, standardize_mol
from rdkit.Chem import MolFromSmiles, MolToInchiKey, MolToSmiles

from syn_retro.path import ASSET_PATH, DATA_PATH
from syn_retro.utils import connect_to_sqlite, get_substrate_smarts

REQUIRED_KEYS = [
    "id",
    "smiles",
    "inchikey",
    "cano_smiles",
    "parent_smiles",
    "parent_inchikey",
]


def query_create_table(table_name: str):
    return f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id INTEGER PRIMARY KEY,
        smiles TEXT,
        inchikey TEXT,
        cano_smiles TEXT,
        parent_smiles TEXT,
        parent_inchikey TEXT
    )
    """


def query_insert(table_name: str):
    return f"""
    INSERT INTO {table_name}
        (id, smiles, inchikey, cano_smiles, parent_smiles, parent_inchikey)
    """


def create_bb_db(db_path: Path) -> Connection:
    """create building block sqlite database
    When reagent_smarts is provided, the reagent specific tables will be created.
    SQLite database schema:
    building_blocks
        id INTEGER PRIMARY KEY,
        smiles TEXT,
        inchikey TEXT,
        cano_smiles TEXT,
        parent_smiles TEXT,
        parent_inchikey TEXT

    Args:
        db_path (Path): path to save sqlite database
    """

    connection = connect_to_sqlite(db_path)
    cursor = connection.cursor()
    cursor.execute(query_create_table("building_blocks"))
    connection.commit()
    return connection


def create_reagent_tables(connection: Connection, reagent_smarts: Dict):
    """Create reagent specific tables

    Args:
        connection (Connection): sqlite connection
        reagent_smarts (Dict): reagent smarts, keys are reagent class names
    """
    cursor = connection.cursor()
    # add molecule column
    cursor.execute("ALTER TABLE building_blocks ADD molecule MOL")
    cursor.execute("UPDATE building_blocks SET molecule=mol_from_smiles(smiles)")
    for k, v in reagent_smarts.items():
        cursor.execute(query_create_table(k))
        cursor.execute(
            query_insert(k)
            + """SELECT building_blocks.id, building_blocks.smiles, building_blocks.inchikey,
            building_blocks.cano_smiles, building_blocks.parent_smiles, building_blocks.parent_inchikey
            FROM building_blocks
            WHERE mol_is_substruct(building_blocks.molecule, mol_from_smarts(?1))""",
            (v,),
        )
        connection.commit()
    cursor.close()


def insert_df_to_db(df: pd.DataFrame, connection: Connection):
    """Insert building block information to sqlite database

    Args:
        df (pd.DataFrame): df including building block information, must include
            id, smiles, inchikey, cano_smiles, parent_smiles, parent_inchikey (REQUIRED_KEYS)
        connection (Connection): sqlite connection
    """
    assert all(
        [k in df.keys() for k in REQUIRED_KEYS]
    ), f"Building {', '.join(REQUIRED_KEYS)} must be provided"
    cursor = connection.cursor()
    cursor.executemany(
        query_insert("building_blocks") + "VALUES (?, ?, ?, ?, ?, ?)",
        df[REQUIRED_KEYS].values.tolist(),
    )
    connection.commit()
    cursor.close()


def add_mol_info(df: pd.DataFrame) -> pd.DataFrame:
    """Add mol, inchikey, cano_smiles, std_mol, parent_mol, parent_smiles, parent_inchikey columns to df

    Args:
        df (pd.DataFrame): dataframe including smiles and id columns

    Returns:
        pd.DataFrame: dataframe including mol, inchikey, cano_smiles, std_mol, parent_mol,
        parent_smiles, parent_inchikey columns
    """
    assert (
        "smiles" in df.keys() and "id" in df.keys()
    ), "Building block id and smiles must be provided"

    df["mol"] = df["smiles"].apply(MolFromSmiles)
    # remove invalid smiles
    df.dropna(subset=["mol"], inplace=True)
    df["inchikey"] = df["mol"].apply(MolToInchiKey)
    df["cano_smiles"] = df["mol"].apply(MolToSmiles)
    df["std_mol"] = df["mol"].apply(standardize_mol)
    df["parent_mol"] = None
    for ii, m in enumerate(df["std_mol"].values):
        if m is None:
            continue
        res, exc = get_parent_mol(m)
        if not exc:
            df.loc[ii, "parent_mol"] = res

    df["parent_smiles"] = df["parent_mol"].apply(
        lambda x: MolToSmiles(x) if x is not None else None
    )
    df["parent_inchikey"] = df["parent_mol"].apply(
        lambda x: MolToInchiKey(x) if x is not None else None
    )
    # remove invalid parent mols
    df.dropna(subset=["parent_inchikey"], inplace=True)
    return df


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create building block database using rdkit and chembl_structure_pipeline",
    )
    parser.add_argument("--csv-path", type=Path, required=True, help="Path to csv file")
    parser.add_argument(
        "--db-dir", type=Path, default=DATA_PATH, help="Dir to sqlite3 db file"
    )
    parser.add_argument(
        "--reagent-class-assets",
        type=Path,
        default=ASSET_PATH / "sub_smarts.yaml",
        help="Path to assets saving reagent smarts",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=10000,
        help="Number of rows to read in at a time",
    )
    return parser.parse_args()


def main():
    args = get_args()
    db_path = args.db_dir / (str(args.csv_path.name).split(".")[0] + ".sqlite")
    reagent_smarts = (
        get_substrate_smarts(args.reagent_class_assets)
        if args.reagent_class_assets.exists()
        else None
    )

    for ii, df in enumerate(pd.read_csv(args.csv_path, chunksize=args.chunksize)):
        if ii == 0:
            # create building block database
            connection = create_bb_db(db_path)

        logging.info(f"{df.shape[0]} building blocks are provided")
        df = add_mol_info(df)
        logging.info(f"{df.shape[0]} building blocks are valid")

        with connection:
            insert_df_to_db(df=df, connection=connection)
            if reagent_smarts is not None:
                create_reagent_tables(
                    connection=connection, reagent_smarts=reagent_smarts
                )

    logging.info(
        f"Building block database created at {args.db_dir / str(args.csv_path.name).split('.')[0]}.sqlite"
    )


if __name__ == "__main__":
    main()
