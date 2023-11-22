from pathlib import Path
from typing import Dict

import pandas as pd
import pytest

from syn_retro.prepare_bb_db import (
    REQUIRED_KEYS,
    add_mol_info,
    create_bb_db,
    create_reactant_tables,
    insert_df_to_db,
)


def test_create_bb_db(bb_csv_path: Path):
    db_path = bb_csv_path.parent / f"{str(bb_csv_path.name).split('.')[0]}.sqlite"
    connection = create_bb_db(db_path)
    cursor = connection.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    assert tables == [("building_blocks",)]
    # query column names of building_blocks table
    cursor.execute("PRAGMA table_info(building_blocks)")
    columns = cursor.fetchall()
    colnms = [c[1] for c in columns]
    assert all([k in colnms for k in REQUIRED_KEYS])
    cursor.close()
    connection.close()
    db_path.unlink()


def test_bad_df_add_mol_info(bb_bad_csv_path: Path):
    # test that add_mol_info raises an error when the df does not have the required keys
    df = pd.read_csv(bb_bad_csv_path)
    with pytest.raises(AssertionError):
        add_mol_info(df)


def test_add_mol_info(bb_csv_path: Path):
    df = pd.read_csv(bb_csv_path)
    df = add_mol_info(df)
    assert all([k in df.columns for k in REQUIRED_KEYS])


def test_insert_df_to_db(bb_csv_path: Path):
    db_path = bb_csv_path.parent / f"{str(bb_csv_path.name).split('.')[0]}.sqlite"
    connection = create_bb_db(db_path)
    df = pd.read_csv(bb_csv_path)
    df = add_mol_info(df)
    insert_df_to_db(df, connection)
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM building_blocks")
    rows = cursor.fetchall()
    assert len(rows) == len(df)
    for k in REQUIRED_KEYS:
        cursor.execute(f"SELECT {k} FROM building_blocks")
        rows = cursor.fetchall()
        assert all([r[0] in df[k].values for r in rows])
    cursor.close()
    connection.close()
    db_path.unlink()


def test_create_reactant_tables(bb_csv_path: Path, reactant_smarts: Dict):
    db_path = bb_csv_path.parent / f"{str(bb_csv_path.name).split('.')[0]}.sqlite"
    connection = create_bb_db(db_path)
    df = pd.read_csv(bb_csv_path)
    df = add_mol_info(df)
    insert_df_to_db(df, connection)
    create_reactant_tables(connection, reactant_smarts)
    cursor = connection.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    assert tables == [("building_blocks",)] + [(k,) for k in reactant_smarts.keys()]
    cursor.close()
    connection.close()
    db_path.unlink()
