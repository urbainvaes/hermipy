# Testing

Run the following command in the root directory.

```bash
$ python -m unittest discover -v tests
```

# Todo

- Auto conversion from triangle -> cross when appropriate?

- Do varfd before tensorization

- Graph of index set with weights

- Add support for Hermite Fourier

- Tie position, directions, and function variables

- Inconsistency with direction in project() and self.dirs

- Project: are directions indices, or actual directions (as in position.dirs)?

- Two passes in tensorize in 1d

- In plot, add bounds of Hermite series

- Order of dim and degree not always consistent

# License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
