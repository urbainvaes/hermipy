install:
	python setup.py install --prefix ~/.local

documentation:
	make -C doc

push:
	rsync -rvut --exclude='/.git' --filter="dir-merge,- .gitignore" . urbain@fenec:phd/code/hermite
