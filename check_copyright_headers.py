# This code is part of KQCircuits
# Copyright (C) 2021 IQM Finland Oy
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# https://www.gnu.org/licenses/gpl-3.0.html.
#
# The software distribution should follow IQM trademark policy for open-source software
# (meetiqm.com/developers/osstmpolicy). IQM welcomes contributions to the code. Please see our contribution agreements
# for individuals (meetiqm.com/developers/clas/individual) and organizations (meetiqm.com/developers/clas/organization).


"""Checks all files for existence of IQM copyright header.

The required copyright header is taken from a file named COPYRIGHT in the working directory.

Paths to exclude from the copyright check can be given using --exclude-paths. These paths should be relative to the
working directory.

Usage:
    python check_copyright_headers.py --exclude-paths path/to/exclude1 path/to/exclude2.py

Exits with code -1 if files without copyright were found, exits with code 0 otherwise. Prints the names of any files
without copyright.
"""

from argparse import ArgumentParser
from pathlib import Path
from sys import exit


if __name__ == "__main__":
    cwd = Path.cwd()

    parser = ArgumentParser()
    parser.add_argument("--exclude-paths", nargs="*", type=str, default=[])
    args = parser.parse_args()
    exclude_paths = [cwd/p for p in args.exclude_paths]

    print("Checking the existence of copyright header in all .py and .lym files...")

    file_paths_1 = list(cwd.glob("**/*.py")) + list(cwd.glob("**/*.lym"))
    file_paths_2 = []
    for path in file_paths_1:
        exclude = False
        for exclude_path in exclude_paths:
            if path == exclude_path or exclude_path in path.parents:
                exclude = True
                break
        if not exclude:
            file_paths_2.append(path)

    files_without_copyright = []
    copyright_string = open("COPYRIGHT").read()
    for file in file_paths_2:
        if copyright_string not in open(file, encoding="utf-8").read():
            files_without_copyright.append(file)
    if len(files_without_copyright) > 0:
        print("Files without copyright header:")
        [print(file) for file in files_without_copyright]
        exit(-1)
    else:
        print("No files without copyright header found.")
        exit(0)
