import os
from lblcrn.surface_crn.surface_crns.readers.statements import *
import lblcrn.surface_crn.surface_crns.readers.transition_readers as transition_readers
import lblcrn.surface_crn.surface_crns.readers.colormap_readers as colormap_readers
import lblcrn.surface_crn.surface_crns.readers.grid_state_readers as grid_state_readers
import lblcrn.surface_crn.surface_crns.readers.totalistic_readers as totalistic_readers
'''
For external use
'''

def read_manifest(filename):
    '''
    Process a manifest file containing options for simulations, and returns a
    dictionary of those options. Options can be:

    1) Single-valued options of the form "FLAG = VALUE". These will be added
       directly to the options dictionary.

    or

    2) Special sections. Special sections include transition rule sections,
       initial state sections, and colormap sections. Sections start with a line
       starting with "!START_<section_name>" and ending with a line starting
       with "!END_<section_name>, where <section_name> is one of
       "TRANSITION_RULES", "INIT_STATE", or "COLORMAP", depending on the
       section. Each section has its own rules for parsing, as described in the
       functions below, and each returns a different data structure, which
       will be added to the options dictionary.

    or

    3) Include statements of the form "!INCLUDE FILENAME". Files with names
       declared in an include statement will be read as though they were in the
       manifest file at the position of the include statement.

    Lines beginning with "#" or "%" are considered comment lines.

    Return value: A dictionary of options.
    '''
    # TODO: if filename is only a stream, (not a string), then
    # think of it as if the comments or the include statements aren't present

    if isinstance(filename, str):

        file_location = filename[0:filename.rfind(os.sep)+1]
        with open(filename, 'rU') as manifest_file:
            manifest_stream = remove_comments(splice_includes(manifest_file,
                                                              file_location))
    else:
        manifest_stream = filename
    return parse_manifest_stream(manifest_stream)


def parse_manifest_stream(manifest_stream):
    '''
    Process a stream of text containing manifest information with include
    statements already flattened into the text.
    See documentation of read_manifest for information on the manifest file
    format.
    '''
    options = dict()
    for line in manifest_stream:
        if line.startswith('#') or line.startswith('%') or line.strip() == "":
            continue
        elif line.startswith(section_starts['totalistic_rule']):
            options['totalistic_rule'] = \
                totalistic_readers.parse_totalistic_rule_stream(manifest_stream)
        elif line.startswith(section_starts['transition_rules']):
            options['transition_rules'] = \
                transition_readers.parse_transition_rule_stream(manifest_stream)
        elif line.startswith(section_starts['init_state']):
            options['init_state'] = \
                     grid_state_readers.parse_grid_state_stream(manifest_stream)
        elif line.startswith(section_starts['colormap']):
            options['colormap'] = \
                         colormap_readers.parse_colormap_stream(manifest_stream)
        elif line.startswith('!'):
            raise Exception('Unrecognized section declaration ' + line +
                            ' in manifest file')
        else:
            parse_option(line, options)
    return options

def splice_includes(manifest_stream, file_location):
    '''
    Convert a stream of manifest file information that might have include
    statements into a stream with included files spliced in.

    See documentation of read_manifest for information on the manifest file
    format.
    '''
    for line in manifest_stream:
        if line.startswith(include_mark):
            filename = line[len(include_mark):].strip()
            if not os.path.isfile(filename):
                filename = file_location + filename
            with open(filename, 'rU') as include_file:
                for line in splice_includes(include_file, file_location):
                    yield line
        else:
            yield line

def remove_comments(manifest_stream):
    '''
    Removes all '#'s and everything after a '#' on the same line. 
    '''
    for line in manifest_stream:
        yield line.split('#')[0]

'''
Internal use only
'''

def parse_option(line, options):
    line_parts = list(map(lambda s: s.strip(), line.split("=")))
    if len(line_parts) != 2:
        raise Exception('Improperly formatted option "' + line + '"')
    options[line_parts[0]] = line_parts[1]
    return