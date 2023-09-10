### CIS 3110 Assignment 1: Computing Human Evolution

Before you begin, read the [A1 writeup](https://canvas.cornell.edu/courses/54963/pages/a1-handout) in Canvas. It contains essential background material and specific instructions about what you need to do.

Some code is given to you, some you must supply. All the source code you need to write will go in the file `src/dna.ml`. You must also write OUnit2 tests, which will go in `test/main.ml`.

There are a few important Makefile targets you should be aware of. Type `make <target>` to execute them. They are:

- `make build` - compiles the project.
- `make utop` - starts the OCaml interactive interpreter with your code preloaded.
- `make test` - runs the OUnit2 test suite in `test/main.ml`.
- `make check` - checks your environment to ensure that all necessary OCaml components are installed and configured correctly.
- `make finalcheck` - does some extra checks to prepare for submission.
- `make zip` - creates a file `dna.zip` suitable for submission to CMS 
- `make clean` - removes compiled files and other compilation artifacts.

When you are finished, run `make zip` to create a file `dna.zip` and upload it to [CMS](https://cmsx.cs.cornell.edu/).
