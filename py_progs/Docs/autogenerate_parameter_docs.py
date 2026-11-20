#!/usr/bin/env python
# -*- coding: <utf-8> -*-
"""
Generates `.yml` format descriptions of the parameters in the C files.

This script goes through the (non-blacklisted) input .c files, and identifies
the input variables, their types and any notes (e.g. units, choices).

It then compares this list to the existing documentation files, and lists
which existing documentation refers to parameters that are now deprecated,
and which parameters that are new and yet to be documented.

By default, the script will simply print the new and deprecated files to
the screen.

Arguments:
    -p / --print
        Print the full text of all .yaml files that would be created to the screen.

    -w / --write
        Move any deprecated parameters to `$SIROCCO/parameters/old/`,
        then write any new parameters to in `$SIROCCO/parameters/`

        Note, this will not over-write any parameters that have changed types
        but not names e.g. `rdflo('thing')` to `rdint('thing')`.

    -h / --help
        Prints this documentation.

After the program has been run $SIROCCO/parameters should contain a yaml file
for every possible input, and any input that has changed significantly should
be in $PYHON/parameters.old  One should normally be sure to add the new yaml
files to the repository.

Note:

    This routine does not autormatial add new yaml files that have been 
    created to the git repository, though it is possible that should 
    be the default.

    Therefore if one uses the -w option, it will create new files, 
    and these will remaing in the local repositiory, even if one changes 
    branches.  This can be confusing!!!

    The program checks which files in the directories it writes to
    are not committed, but it is up to the user to sort out what s/he
    wants to do.

    The command to list files in a directory that are not tracked is::

        git ls-files --other

    if one is in the directory in question.

    The command to remove files in a directory (from within that directory) is::

        git clean -f

    The recommendation is to:

    * to clean both $SIROCCO/parameters/, and $SIROCCO/parameters/old/ from your
      local directories before using writing files using this routine, and then
    * to add and comit  all of the files that are produced before going on to other
      stages of activities associated with documentation.

"""
import os
import sys
from typing import List, Union
from collections import OrderedDict, defaultdict
import yaml


# ==============================================================================
# CONSTANTS
# ==============================================================================

# Do not generate documentation for inputs in files containing these substrings
BLACKLIST = [
    "py_wind", "plot_roche", "py_grid", "synonyms",
    "t_bilinear", "rdpar", "parse", "para_update"
]

# Types of parameter read functions and the type they correspond to.
# Must include 'wrappers' e.g. get_spectype.
PARAM_TYPES = {
    "get_spectype": "Enumerator",
    "rdint": "Integer",
    "rdflo": "Float",
    "rddoub": "Double",
    "rdstr": "String",
    "rdchoice": "Enumerator",
    "rdline": "Line"
}

# Types of unit in the description, and what to write in the output.
# The first matching entry will be applied (e.g. msol/yr before msol).
UNIT_TYPES = OrderedDict([
    ("cm", "cm"),
    ("ergs/s", "ergs/s"),
    ("ev", "eV"),
    ("deg", 'Degrees'),
    ("msol/yr", u"M☉/year"),
    ("msol", u"M☉"),
    ("hr", "Hours"),
    ("r_star", "co.radius"),
    ("rstar", "co.radius"),
    ("wd_rad", "co.radius"),
    ("r_g", "co.gravitational_radius"),
    ("cgs", "cgs"),
    ("vescap", "Escape velocity"),
    ("k", "Kelvin"),
    ("K", "Kelvin"),
])

# Ways boolean variables are shown in the description and their types.
BOOL_TYPES = {"1=yes": "Boolean (1/0)",
              "y=1": "Boolean (1/0)",
              'anything_else': "Boolean (1/0)",
              'y/n': 'Boolean (Y/N)',
              'yes,no': 'Boolean (yes/no)'}


# ==============================================================================
# FUNCTIONS
# ==============================================================================

def parse_param_to_dict(found_parameters: List[dict], output_dict: dict,
                        existing_documentation: List[str],
                        text: str, input_file: str, param_type: str) -> bool:
    """
    Given a string that's the text of a c function call for a parameter,
    decide if this is a parameter and return so.

    Arguments:
        found_parameters: Array of found parameter entries so far.
        output_dict: The dict of parameter values.
        existing_documentation: List of already-found documentation.
        text: Parameter text e.g. "foo(bar, baz)".
        input_file: File the parameter is in e.g. "setup.c".
        param_type: Parameter type e.g. Integer, Enumerator.

    Returns:
        bool False: Probably unnecessary. Doesn't need to return a value.
    """
    param_dict = OrderedDict([
        ('name', 'Name'),
        ('description', 'Multi-line description, must keep indentation.\n'),
        ('type', param_type),
        ('unit', 'None'),
        ('values', 'Condition e.g. greater than 0 or list e.g. [1, 2, 5]'),
        ('parent', {"parameter": "Condition e.g. greater than 0 or list e.g. [1, 2, 5]"}),
        ('file', input_file)
    ])

    if "String" in param_type:
        # If this is a string parameter, it has no units or enumerator values
        param_dict.pop('unit', None)
        param_dict.pop('values', None)
    elif "Enumerator" in param_type:
        # If it's an enumerator, it has no unit
        param_dict.pop('unit', None)

    # If this is an 'advanced' parameter, record so
    if '@' in text[0]:
        param_dict["advanced"] = True
        text = text[1:]

    text_finished = False

    # If this parameter doesn't have a description, we have full name!
    if '(' not in text:
        param_dict["name"] = text.strip()
        text_finished = True
    # If this parameter *does* have a description. record the name
    else:
        if ')' in text:
            # If this is a properly formatted description with a )
            # We want to cut out the whole bracketed section e.g. "foo(bar,baz)_blah" -> "foo_blah"
            description = text[text.find('('):text.rfind(')')+1].strip()
        else:
            # Otherwise, this is a mangled name e.g. "foo(bar,baz" -> "foo"
            description = text[text.find('('):].strip()
        param_dict["name"] = text.replace(description, '').strip()
        text = description.replace('(', '').replace(')', '').strip()

    # Record we've found it
    found_parameters.append(param_dict["name"])

    # Exponent/power law parameters have no values
    if "exponent" in param_dict["name"] or "power_law" in param_dict["name"]:
        param_dict.pop("unit", None)

    # If this parameter isn't already documented, record!
    if param_dict["name"] not in existing_documentation:
        output_dict[param_dict["name"]] = param_dict
    # If this parameter is in the existing documentation, stop and find another
    else:
        return False

    # If this parameter had a name only, no description, stop here
    if text_finished:
        return False

    # Otherwise, we need to parse the remaining descriptive text e.g
    # the (bar,baz) from foo(bar,baz)
    for key, value in BOOL_TYPES.items():
        # Check to see if this is a de-facto boolean
        if key in text.lower():
            param_dict["type"] = value.strip()
            param_dict.pop("values", None)
            param_dict.pop("unit", None)
            text_finished = True
            break
    if text_finished:
        return False

    # So: Is this string (blahblahblah) or is it (foo,bar,baz)?
    if ',' not in text:
        # If this is not an enumerated list, it's a unit or description.
        for key, value in UNIT_TYPES.items():
            # Is the description a unit?
            if key in text.lower():
                param_dict["unit"] = value.strip()
                return False

        # If it's not a unit, and not an enumerator, it's just a description
        param_dict["description"] = text.strip()+'\n'

    else:
        # If it is enumerated, it has no units and a list of options
        param_dict.pop('unit', None)
        param_dict['values'] = {
            option: 'Multi-line description, must keep indentation.\n' for option in text.split(',')
        }

    return False


def list_existing_documentation(directory: str) -> List[str]:
    """
    Look through the output folder, and record every .yaml file in it.

    Arguments:
        directory: The location to look in.
    Returns:
        List[str]: Existing .yaml documentation files.
    """
    existing_documentation = []
    for filename in os.listdir(directory):
        parameter = os.path.splitext(filename)
        if 'yaml' in parameter[1]:
            existing_documentation.append(parameter[0])
    return existing_documentation


def list_input_files(directory: str) -> List[str]:
    """
    Look through the source folder, and track each file that isn't in the blacklist.

    Arguments:
        directory: The location to look in.
    Returns:
        List[str]: Input .c files.
    """
    input_files = []
    for filename in os.listdir(directory):
        # For each file in the input folder
        if filename.endswith(".c") and not any(black in filename for black in BLACKLIST):
            # This is a valid input file
            input_files.append(filename)
        else:
            continue
    return input_files


def read_parameters(input_folder: str, input_files: List[str],
                    existing_documentation: List[str]) -> Union[OrderedDict, List[str]]:
    """
    For each input file, process and add to found parameters list.

    Arguments:
        input_files: The .c files that need parsing for parameters.
    Returns:
        OrderedDict(): The parameter dict
        List[str]: The list of parameter names.
    """
    found_parameters = []
    output_dict = OrderedDict()

    for input_file in input_files:
        # For each input file with parameters in
        file_object = open(os.path.join(input_folder, input_file), "r")
        file_text = file_object.read()

        # Remove all the block comments
        while "/*" in file_text:
            start = file_text.find('/*')
            end = file_text.find('*/', start)
            file_text = file_text.replace(file_text[start:end+2], ' ')

        # Remove all the in-line comments
        while '//' in file_text:
            start = file_text.find('//')
            end = file_text.find('\n', start)
            file_text = file_text.replace(file_text[start:end], ' ')

        # Remove the linebreaks, then split based on command terminators.
        file_text = file_text.replace('\n', ' ')
        file_lines = [line.strip() for line in file_text.split(';')]

        for line in file_lines:
            # First, we check the dict of parameter read functions to see if this line
            # contains one of them.
            for rd_function, param_type in PARAM_TYPES.items():
                if rd_function in line.lower():
                    # We've found a rdpar line with a readable entry,
                    # now let's build the text to parse

                    if '(' not in line or '"' not in line:
                        # This isn't anything useful
                        continue
                    elif 'Error' in line:
                        continue

                    # Now we've got the full line, extract the parameter text
                    # e.g. rdchoice("var(foo,bar,baz)"", value, choices)
                    # we want to pass "var(foo,bar,baz)" for analysis.
                    start = line.find('"')
                    end = line.find('"', start+1)

                    parse_param_to_dict(found_parameters=found_parameters,
                                        output_dict=output_dict,
                                        existing_documentation=existing_documentation,
                                        text=line[start+1:end], input_file=input_file,
                                        param_type=param_type)

    return found_parameters, output_dict


def intersect_documentation(existing_documentation: List[str],
                            found_parameters: List[str]) -> Union[List[str], List[str]]:
    """
    Intersections of documentation

    Args:
        existing_documentation:
        found_parameters:

    Returns:
        List[str]:
        List[str]:
    """
    deprecated_documentation = sorted(list(set(existing_documentation) - set(found_parameters)),
                                      key=lambda s: s.lower())
    new_documentation = sorted(list(set(found_parameters) - set(existing_documentation)),
                               key=lambda s: s.lower())
    return deprecated_documentation, new_documentation


# https://stackoverflow.com/questions/8640959/how-can-i-control-what-scalar-form-pyyaml-uses-for-my-data
# Use literal outputs for multiline
def should_use_block(value: str) -> bool:
    """
    Whether or not this string is multi-line and should use the '|' literal string format

    Args:
        value: String to test
    Returns:
        bool: Whether or not this is a multi-line string
    """
    for character in u"\u000a\u000d\u001c\u001d\u001e\u0085\u2028\u2029":
        if character in value:
            return True
    return False


def my_represent_scalar(self, tag, value, style=None):
    """
    Replacement 'representer' for strings that intelligently switches between single- and
    multi-line output

    Args:
        tag (string): YAML tag
        value (string): YAML value to assess
    Returns:
        function: The representation function
    """
    if style is None:
        if should_use_block(value):
            style = '|'
        else:
            style = self.default_style

    node = yaml.representer.ScalarNode(tag, value, style=style)
    if self.alias_key is not None:
        self.represented_objects[self.alias_key] = node
    return node


def deprecate_documentation(deprecated_documentation: List[str], output_folder: str, output_old_folder: str):
    """
    Moves any deprecated documentation to the archive folder.

    Arguments:
        deprecated_documentation: The docs that need moving to old storage
        output_folder: The folder to write new docs to
        output_old_folder: The folder to move old docs to
    """
    for filename in deprecated_documentation:
        original_path = os.path.join(output_folder, '{}.yaml'.format(filename))
        moved_path = os.path.join(output_old_folder, '{}.yaml'.format(filename))
        os.rename(original_path, moved_path)


def yaml_output(output_dict: dict, new_documentation: List[str], output_folder: str,
                print_docs: bool = False, write_docs: bool = False):
    """
    Prints out new documentation to screen or writes it to file.

    Arguments:
        output_dict: The dictionary representation of the parameters found.
        new_documentation: The docs that need writing (from intersect_documentation)
        output_folder: The folder to write new docs to
        print_docs: Whether to print new documentation to the screen
        write_docs: Whether to write new documentation to disk/move old doc
    """

    # https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
    # Dump ordered dicts properly
    def dict_representer(dumper, data):
        """
        Function for representing ordered dictionaries like regular ones

        Args:
            dumper (class): Type of object to be dumped
            data(iterator): Output representation of class
        Returns:
            function: Way a dictionary should be represented
        """
        return dumper.represent_dict(data.items())

    yaml.representer.BaseRepresenter.represent_scalar = my_represent_scalar
    yaml.add_representer(OrderedDict, dict_representer)

    # Output to file and print to screen
    if print_docs:
        print("Printing new documentation to screen...\n")
    if write_docs:
        print("Writing out...")

    for key, value in output_dict.items():
        print("File: {}/{}.yaml".format(output_folder, key))
        if print_docs and key in new_documentation:
            print(yaml.dump(value, default_flow_style=False, allow_unicode=True))

        if write_docs and key in new_documentation:
            with open(os.path.join(output_folder, "{}.yaml".format(key)), "w") as file_object:
                # For each key, open a yaml file and write to it.
                yaml.dump(value, file_object, default_flow_style=False, allow_unicode=True)


def autogenerate_parameter_docs():
    """
    Function to autogenerate parameter documentation.
    """
    input_folder = os.path.join(os.environ["SIROCCO"], "source")
    output_folder = os.path.join(os.environ["SIROCCO"], "docs", "parameters")
    output_old_folder = os.path.join(os.environ["SIROCCO"], "docs", "parameters", "old")

    input_files = list_input_files(input_folder)
    existing_documentation = list_existing_documentation(output_folder)
    found_parameters, output_dict = read_parameters(input_folder, input_files, existing_documentation)
    deprecated_documentation, new_documentation = intersect_documentation(existing_documentation, found_parameters)

    duplicate_roots: List[str] = []
    parameter_roots: defaultdict = defaultdict(list)
    for parameter in found_parameters:
        parameter_roots[parameter.split('.')[0]].append(parameter)


    for index, root in enumerate(list(parameter_roots.keys())[:-1]):
        for root2 in list(parameter_roots.keys())[index+1:]:
            if root != root2 and root.casefold() == root2.casefold():
                if len(parameter_roots[root]) < len(parameter_roots[root2]):
                    duplicate_roots.append('- {} and {}: {}'.format(root, root2, ', '.join(parameter_roots[root])))
                else:
                    duplicate_roots.append('- {} and {}: {}'.format(root, root2, ', '.join(parameter_roots[root2])))


    if duplicate_roots:
        print("ERROR: There are parameter roots that only differ by case!")
        print('\n'.join(duplicate_roots))
        print("Cannot write documentation as it will not work in case-insensitive OSes (e.g. Mac, Windows)")
        return

    if len(sys.argv) == 1:
        # If we're not running in write mode
        print("Documentation for parameters that no longer exist:")
        for param in deprecated_documentation:
            print(" - {}".format(param))
        print("New, previously undocumented parameters:")
        for param in new_documentation:
            print(" - {}".format(param))
        print("Rerun with:\n" +
              "'-h/--help' to print the documentation for this script.\n" +
              "'-w/--write' to move old parameters to 'parameters/old' and write new files.\n" +
              "'-p/--print' to print new files to screen.")
    else:
        if '-w' in sys.argv or '--write' in sys.argv:
            deprecate_documentation(deprecated_documentation, output_folder, output_old_folder)

        yaml_output(
            output_dict, new_documentation, output_folder,
            print_docs=('-p' in sys.argv or '--print' in sys.argv),
            write_docs=('-w' in sys.argv or '--write' in sys.argv)
        )

    print('\n*** These are uncommitted files in the parameters directory ***\n')

    os.system('git ls-files --other %s' % output_folder)

    print('\n*** These are uncommitted files in the deprecated parameters directory ***\n')

    os.system('git ls-files --other %s' % output_old_folder)

    print('\n*** If either of the above showed uncommitted files, one should either commit them')
    print('or rm them (see the help for more info.)  ***\n')


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if '-h' in sys.argv or '--help' in sys.argv:
        print(__doc__)
    else:
        autogenerate_parameter_docs()
