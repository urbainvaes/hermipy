install:
	python setup.py install --prefix ~/.local

doc:
	make -C docs

push:
	rsync -rvut --exclude='/.git' --filter="dir-merge,- .gitignore" . urbain@fenec:phd/code/hermite

test:
	coverage run --source=hermipy -m unittest discover -v -f tests
	rm -f tests/.coverage tests/.coverage.svg
	coverage report > tests/.coverage
	coverage-badge -o tests/.coverage.svg
