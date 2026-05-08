#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use this script to (mass) update either the name of a parameter or to update the value of a parameter.

The main purpose is to change multiple parameter files at once.
As such, the script recursively searches for all parameter files from
the calling directory and updates them. Can also be used on a single parameter
file or for parameter files with a given root.
"""


import argparse as ap
from pathlib import Path


def get_root_from_filepath(
    where, returnwhere=True
):
    """Get the root name of a Python simulation, extracting it from a file path.

    Parameters
    ----------
    where: str
        The directory path to a Python .pf file
    returnwhere: str [optional]
        Returns the directory containing the .pf file.

    Returns
    -------
    root: str
        The root name of the Python simulation
    where: str
        The directory path containing the provided Python .pf file"""

    if type(where) is not str:
        raise TypeError("expected a string as input for the file path, not whatever you put")

    dot = where.rfind(".")
    slash = where.rfind("/")

    root = where[slash:dot]
    where = where[:slash]
    if where == "":
        where = "."

    if returnwhere:
        return root, where
    else:
        return root


def find_parameter_files(
    updateroot=None, where="."
):
    """Search recursively for Python .pf files. This function will ignore
    py_wind.pf parameter files, as well as any root.out.pf files.

    Parameters
    ----------
    updateroot: str [optional]
        If given, only .pf files with the given root will be returned.
    where: str [optional]
        The directory to search for Python .pf files from

    Returns
    -------
    parameterfiles: list[str]
        The file paths for any Python .pf files founds"""

    parameterfiles = []

    for filepath in Path(where).glob("**/*.pf"):
        str_fp = str(filepath)
        if str_fp[0] == "/":  # Add dot to start of file path if it isn't there
            str_fp = "." + str_fp

        # Do not add anything with .out.pf or any py_wind.pf files

        if str_fp.find(".out.pf") != -1:
            continue
        elif str_fp.find("py_wind.pf") != -1:
            continue

        # If parameter files of given root name updateroot are only to be
        # updated, then check if the root name is correct

        if updateroot is not None:
            thisroot, thiswhere = get_root_from_filepath(str_fp)
            if thisroot != updateroot:
                continue

        parameterfiles.append(str_fp)

    # Sort the parameter files so the output is consistent between runs

    parameterfiles = sorted(parameterfiles, key=str.lower)

    return parameterfiles


def update_parameter(
    mode, where, parameter, newvalue
):
    """Change the name of a parameter in a Python parameter file. If the old and
    new parameter name are the same, the script will still update the parameter
    file anyway.

    Parameters
    ----------
    mode: str
        The run mode of the script, either "value" to update a parameter value
        or "name" to update a parameter name.
    where: str
        The path to the parameter file to update.
    parameter: str
        The name of the parameter to update.
    newvalue: str
        The new value of the parameter."""

    if where.find(".pf") == -1:
        raise IOError("provided parameter file path {} does not include a .pf".format(where))
    if mode not in ["value", "name"]:
        print("{} is an unknown mode for function".format(mode))

    new = None  # new is used to track to see if a change was made or not

    with open(where, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.find(parameter) != -1:
            if mode == "value":
                new = "{}{:20s}{}\n".format(parameter, " ", newvalue)
            else:
                parametervalue = line.split()[-1]
                new = "{}{:20s}{}\n".format(newvalue, " ", parametervalue)
            lines[i] = new
            break

    if not new:
        print("Haven't been able to update {} to {} for {}".format(parameter, newvalue, where))
        return

    with open(where, "w") as f:
        f.writelines(lines)

    return


def main():
    """The main function of the script: parses input from the command line using
    argparse and, most importantly, controls the flow for updating parameter
    files."""

    p = ap.ArgumentParser(description=__doc__)

    p.add_argument(
        "mode", choices=["value", "name"],
        help="Run mode of the script: either update a -value- or the -name- of a parameter."
    )
    p.add_argument(
        "update", help="The name of the parameter to update."
    )
    p.add_argument(
        "new", help="The new name/value of the parameter to update."
    )
    p.add_argument(
        "-r", "--root", default=None, help="Update parameter files with a given root name."
    )

    args = p.parse_args()

    for parameterfile in find_parameter_files(args.root):
        update_parameter(args.mode, parameterfile, args.update, args.new)

    return


if __name__ == "__main__":
    main()
