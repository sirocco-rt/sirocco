#!/usr/bin/env python
"""
Converts `.yml` format parameter descriptions to `.rst` files.

This script goes through the created yaml documentation and converts them to
.rst files, then writes them to an output directory

RST files to the git repository. It then calls
Sphinx to build the docs, in the documentation/html folder, and opens them in
the browser to test.

Usage: autogenerate_rtd_pages.py [-h] [-outdir whatever]

Arguments:
    None  Then the rst files are written to direcotry $Python/docs/rst/parameters
    -h    Print this documenation and quite
    -outdir whatever then write the rst files to that directoryk

Notes:
    The user is prevented (at least from the command line) from writing into
    the main sphinx directory, inorder to avoid overwriting the existing files

"""
from typing import TextIO
import os
import sys
from subprocess import call
import yaml


def write_header_by_level(output_file: TextIO, string: str, level: int = 0):
    """
    Writes a passed string as a header, with an appropriate type of underscore depending on
    nesting level (e.g. heading -> subheading)

    Arguments:
        output_file: The file to write to
        string: The text to write
        level: The level of heading, 0=heading, 1=subheading, 2=subsubheading, 3=paragraph heading
    """
    if level == 0:
        output_file.write("{}\n{}\n".format(string, ''.join('=' for i in range(len(string)))))
    elif level == 1:
        output_file.write("{}\n{}\n".format(string, ''.join('-' for i in range(len(string)))))
    elif level == 2:
        output_file.write("{}\n{}\n".format(string, ''.join('^' for i in range(len(string)))))
    else:
        output_file.write('**{}**\n{}\n'.format(string, ''.join('"' for i in range(len(string)+4))))


def write_str_indent(output_file: TextIO, string: str, indent: str = "  ", all: bool = False):
    """
    Writes a passed string to the output file. If it's split over multiple lines,
    indent the subsequent lines by the specified amount.

    Arguments:
        output_file: The file to write to
        string: The text to write
        indent: The indent to apply if the text splits over multiple lines
        all: Whether or not all lines should be indented
    """
    lines = string.splitlines()
    if all:
        output_file.write("{}{}\n".format(indent, lines[0]))
    else:
        output_file.write("{}\n".format(lines[0]))

    for line in lines[1:]:
        output_file.write("{}{}\n".format(indent, line))
    output_file.write("\n")
    return


def output_parameter(parameter: dict, output_file: TextIO):
    """
    Output a parameter to file

    Arguments:
        parameter: The parameter to write to file
        output_file: The file to write to
    """
    write_header_by_level(output_file, parameter['name'], 0)
    output_file.write("{}\n".format(parameter['description']))

    if parameter.get('type'):
        output_file.write("Type\n  {}\n\n".format(parameter['type']))
    if parameter.get('unit'):
        output_file.write("Unit\n  {}\n\n".format(parameter['unit']))

    if parameter.get('values'):
        if isinstance(parameter['values'], dict):
            output_file.write("Values\n")
            for key, value in parameter['values'].items():
                output_file.write("  {}\n".format(key))
                if isinstance(value, str):
                    write_str_indent(output_file, value.strip(), indent="    ", all=True)
                elif isinstance(value, list):
                    write_str_indent(
                        output_file,
                        ', '.join([str(x) for x in value]), indent="    ", all=True
                    )
                else:
                    output_file.write("    {}\n".format(value))

        elif isinstance(parameter['values'], list):
            # If this is a list of values, write each out as a bullet-point
            output_file.write("Values\n")
            for value in parameter['values'].items():
                write_str_indent(output_file, "  * {}".format(value), indent="  ")

        else:
            output_file.write("Values\n  {}\n".format(parameter['values']))
        output_file.write("\n")

    output_file.write(
        "File\n  `{} <https://github.com/agnwinds/python/blob/master/source/{}>`_\n\n\n".format(
            parameter['file'], parameter['file']
        )
    )

    if parameter.get('parent'):
        if isinstance(parameter['parent'], dict):
            output_file.write("Parent(s)\n")
            for key, value in parameter['parent'].items():
                if isinstance(value, str):
                    write_str_indent(output_file, "  * :ref:`{}`: {}".format(key, value), indent="    ")
                elif isinstance(value, list):
                    list_text = ', '.join(['``'+str(x)+'``' for x in value])
                    write_str_indent(output_file, "  * :ref:`{}`: {}".format(key, list_text, indent="    "))
                else:
                    output_file.write("  * :ref:`{}`: ``{}``\n\n".format(key, value))
            output_file.write("\n")

        elif isinstance(parameter['parent'], list):
            # If this is a list of parents, write each out as a bullet-point
            output_file.write("Parent(s)\n")
            for value in parameter['parent'].items():
                write_str_indent(output_file, "  * {}".format(value), indent="    ")
        else:
            output_file.write("Parent(s)\n  {}\n\n".format(parameter['parent']))

    # Go through all the children and output them
    if parameter.get('children'):
        output_file.write("Child(ren)\n")
        for key in parameter['children'].keys():
            write_str_indent(output_file, "  * :ref:`{}`".format(key), indent="    ")
    return


def read_yaml(input_folder: str) -> dict:
    """
    Read all the .yaml files in to dicts. Categorise by type, e.g. 'reverb.mode' is type 'reverb'

    Arguments:
        input_folder: Path to the folder of input yaml files.

    Returns:
        dict: Dictionary of parameter dicts
    """
    dox_all = {}

    # Read in all parameters from the input files
    for input_file in [filename for filename in os.listdir(input_folder) if filename.endswith(".yaml")]:
        with open(os.path.join(input_folder, input_file), "r") as file_object:
            # Open each yaml file and load from it.
            try:
                parameter = yaml.full_load(file_object)
                parameter_name = parameter['name']

                # Get the type from the name (e.g. reverb.mode and set no children)
                parameter_type = parameter['name'].split('.')[0]
                parameter['root_type'] = parameter_type.strip()
                parameter['children'] = {}

                # Register this on the full list of all parameters
                dox_all[parameter_name] = parameter

            except yaml.YAMLError as error:
                print(error)

    return dox_all


def write_rst(output_folder: str, dox_all: dict, dox_structured: dict):
    """
    Write the documentation out to rst files

    Arguments:
        output_folder: The folder to write to
        dox_structured: The tree structured documentation
    """

    output_file = open(os.path.join(output_folder, 'top_level_parameters.rst'), 'w')
    write_header_by_level(output_file, 'Top level parameters', 0)
    output_file.write("\n\n.. todo::\n    Fill in\n\n")
    for parameter in dox_all.values():
        if not 'parent' in parameter.keys():
            output_file.write(' * :ref:`{}`\n'.format(parameter['name']))
    output_file.close()

    for root, root_dict in dox_structured.items():
        if len(root_dict['all']) > 1:
            try:
                os.makedirs(os.path.join(output_folder, root))
            except OSError as e:
                pass

            output_file = open(os.path.join(output_folder, root+'.rst'), 'w')
            write_header_by_level(output_file, root, 0)
            output_file.write("\n\n.. todo::\n   Fill in\n")
            output_file.write("\n\n.. toctree::\n   :glob:\n\n   {}/*".format(root))
            output_file.close()

    # Go through the parameter tree and write it out
    for parameter in dox_all.values():
        if len(dox_structured[parameter['root_type']]['all']) == 1:
            output_file = open(
                os.path.join(output_folder, parameter['name']+'.rst'), 'w'
            )
        else:
            output_file = open(
                os.path.join(output_folder, parameter['root_type'], parameter['name']+'.rst'), 'w'
            )
        output_parameter(parameter, output_file)
        output_file.close()

def autogenerate_rtd_pages(output_folder):
    """
    Write the RTD files to disk and add them to git, then run sphinx-build to generate the docs.
    """
    # html_folder = os.path.join(os.environ["SIROCCO"], "docs", "html")
    docs_folder = os.path.join(os.environ["SIROCCO"], "docs")
    dox_all = {}

    par_folder=  os.path.join(os.environ["SIROCCO"], "docs", "parameters")
    print('Par folder',par_folder)

    dox_all = read_yaml(
        os.path.join(par_folder)
    )

    # This is a tree-structured dictionary e.g.
    # reverb: {
    #    type: {
    #        children: {
    #           mode: {},
    #           visualisation: {
    #               children: {
    #                   angular_bins: {}
    #               }
    #            }
    #        }
    #    },
    #    all: [type, mode, visualisation, angular_bins]
    # }
    dox_structured = {}

    # For each possible root type, set up a root dictionary for it
    for parameter_name, parameter in dox_all.items():
        dox_structured[parameter['root_type']] = {'all': []}

    # We want to build a tree with children etc. to respect parentage
    # Now, for the unstructured list of all parameters
    for parameter_name, parameter in dox_all.items():
        # For each parameter in the list of all parameters, check its parents
        parameter_type = parameter['root_type']

        if parameter.get('parent'):
            # If it has parents
            found_parent = False
            if isinstance(parameter['parent'], dict):
                for parent_name in parameter['parent'].keys():
                    # Check each parent.
                    if parent_name in dox_all.keys():
                        # If it's parent is another parameter (i.e. in .keys())
                        # We assign this parameter to the children's list of its parent
                        dox_all[parent_name]["children"][parameter_name] = parameter
                        found_parent = True

            if not found_parent:
                # If we didn't find the parameter's parent elsewhere,
                # then add it to the root.
                dox_structured[parameter_type][parameter_name] = parameter

        else:
            # If it has no parents, add it to the top layer
            dox_structured[parameter_type][parameter_name] = parameter

        dox_structured[parameter_type]['all'].append(parameter_name)


    write_rst(output_folder, dox_all, dox_structured)

def steer(argv):
    '''
    This is just a steering routine to enable better control of the program
    '''

    outdir= os.path.join(os.environ["SIROCCO"], "docs", "rst", "parameters")
    xdir= os.path.join(os.environ["SIROCCO"], "docs","sphinx")

    i=1
    while i< len(argv):
        if argv[0:2]=='-h':
            print(__doc__)
        elif argv[i]=='-outdir':
            i+=1
            outdir=argv[i]
            outdir=oa.path.abspath(outdir)

        i+=1

    if outdir.count(xdir)>0:
        print('Error: It is too dangeroust to write in a subdirector of %s, choose a different location' % xdir)
        return 

    if os.path.isdir(outdir)==False:
        print('Creating %s ' % outdir)
        os.makedirs(outdir,exist_ok=True)
    else:
        print('The output directory %s already exists; files will be replaced or added' % outdir)


    autogenerate_rtd_pages(outdir)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    steer(sys.argv)
