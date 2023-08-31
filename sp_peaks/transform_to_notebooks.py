#!/usr/bin/env python

import os
import jupytext

def transform_to_notebooks():
    """Transform to Notebooks.

    Transform python scripts to jupyter notbooks.

    """
    notebook_dir = "./notebooks"

    for filename in os.listdir(notebook_dir):
        if filename.endswith(".py"):
            nb_file = os.path.join(
                notebook_dir,
                f"{os.path.splitext(filename)[0]}.ipynb",
            )
            if not os.path.exists(nb_file):
                py_file = os.path.join(notebook_dir, filename)
                jupytext.write(jupytext.read(py_file), nb_file)
                print(f"{py_file} -> {nb_file}")
            else:
                print(f"{nb_file} exists")


def main():
    transform_to_notebooks()

    return 0


if __name__ == "__main__":
    main()
