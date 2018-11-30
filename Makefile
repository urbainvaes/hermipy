install:
	python setup.py install --prefix ~/.local

doc:
	make -C docs

push:
	rsync -rvut --exclude='/.git' --filter="dir-merge,- .gitignore" . urbain@fenec:phd/code/hermite

test:
	coverage run --source=hermipy -m unittest discover -v -f tests
	rm -f coverage.txt coverage.svg
	coverage report > coverage.txt
	coverage-badge -o coverage.svg
