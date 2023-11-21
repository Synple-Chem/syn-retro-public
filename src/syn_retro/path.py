from pathlib import Path

# project/package path
PROJ_PATH = Path(__file__).parent.parent.parent
TEST_PATH = PROJ_PATH / "tests"
# root path of where the project is exceuted
ROOT_PATH = Path(".").absolute()
ASSET_PATH = ROOT_PATH / "assets"
DATA_PATH = ROOT_PATH / "data"
