push:
	rsync -r --delete --exclude='/.git' --filter="dir-merge,- .gitignore" . urbain@155.198.193.89:phd/code/hermite
