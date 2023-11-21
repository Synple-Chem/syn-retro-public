.PHONY: test # these are not real files
python=./env/bin/python

all: env

env: synutils.git chemblpipe.git
	conda env create -f ./environment.yaml -p ./env
	cd synutils.git && ../env/bin/python -m pip install .
	cd chemblpipe.git && ../env/bin/python -m pip install .
	${python} -m pip install -e .

env-dev: precommit synutils.git chemblpipe.git
	conda env create -f ./environment.yaml -p ./env
	cd synutils.git && ../env/bin/python -m pip install -e .
	cd chemblpipe.git && ../env/bin/python -m pip install -e .
	${python} -m pip install -e .

env-local: precommit chemblpipe.git
	conda env create -f ./environment.yaml -p ./env
	cd ../synple-utils && ../syn-retro-public/env/bin/python -m pip install -e .
	cd chemblpipe.git && ../env/bin/python -m pip install -e .
	${python} -m pip install -e .

chemblpipe.git:
	git clone https://github.com/Ohyeon5/ChEMBL_Structure_Pipeline.git chemblpipe.git --depth 1

synutils.git:
	git clone https://github.com/Synple-Chem/synple-utils.git synutils.git --depth 1

precommit:
	bash ./scripts/install_precommit.sh

prepare-bb-db:
	prepare_bb_db --csv-path ${CSV_PATH}

run-retro:
	run_retro --data-path ${CMPD_PATH}

test:
	${python} -m pytest ./tests
