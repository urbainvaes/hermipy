# Testing

To run the tests, run the following in the root directory.

```bash
$ python -m unittest discover -v -s tests
```

# Todo

- Tie position, directions, and function variables

- Duplicate split_operator

- Inconsistency with direction in project() and self.dirs

- Add function to lower degree

- Two passes in tensorize in 1d

- Unify interface for project(vector or int?)

- Is it possible to solve the projected eigenvalue problem?

- Rexpress project in terms of inner?

- In plot, add bounds of Hermite series

- Project: are directions indices, or actual directions (as in position.dirs)?

- Remove dim from attributes (only in Position)
