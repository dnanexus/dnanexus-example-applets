#!/usr/bin/env python
"""
Script that generates site pages in /docs/pages

TODO: Multithreading logger the right way.
"""
from __future__ import print_function
import argparse
import os
import json
import re
import logging
import fnmatch
from multiprocessing import Pool
from itertools import izip
from datetime import date
from FrontMatterClass import FrontMatter
from SectionParserClass import SectionParser
from global_helper_vars import TUTORIAL_TYPES_SEARCH, SUPPORTED_INTERPRETERS, NUM_CORES, AppObj


def _get_section_parser(page_dict, logger):
    tutorial_parser = None
    if SUPPORTED_INTERPRETERS.get(page_dict["interpreter"]):
        tutorial_parser = SectionParser(
            code_file_path=page_dict["src_code"],
            ignore_comments=page_dict["ignore_comments"],
            interpreter=page_dict["interpreter"], logger=logger)
    return tutorial_parser


def _rehydrate_site(user_args):
    def update_dict(d1, d2):
        d1.update(d2)
        return d1

    def archived_applets():
        archiv_applets = os.path.join(os.path.abspath(os.path.join(user_args.tutorials_dir, '..')), 'Example', 'archived')
        return find_all_matches(archiv_applets, "dxapp.json")

    user_input_dict = {
        "overwrite": user_args.overwrite,
        "ignore_comments": user_args.skip_comments,
        "site_pages_dir": user_args.site_pages_dir
    }
    dxapp_files = find_all_matches(user_args.tutorials_dir, "dxapp.json")
    # Add non tutorial pages to page dict here
    dxapp_files.extend(archived_applets())
    # Added them
    page_dicts = [update_dict(d, user_input_dict) for d in _resolve_applets(dxapp_files)]
    md_maker = Pool(processes=NUM_CORES)
    status = md_maker.map(create_jekyll_markdown_tutorial, page_dicts)
    md_maker.close()
    md_maker.join()
    for i in xrange(len(status)):
        msg = "SUCCESS: {title}" if status[i][0] else "FAIL: {title}"
        print(msg.format(title=page_dicts[i]["title"]))
        if not status[i][0]:
            print(status[i][1])


def _resolve_applet(dxapp_path):
    """Moved to top level due to Pool.map() inability to use nested func or output dict"""
    with open(dxapp_path, "r") as dxf:
        dxapp_obj = json.load(dxf)
    src_code = os.path.join(dxapp_path[:-11], dxapp_obj["runSpec"]["file"])  # pretty sure this is a required field
    readme_md = os.path.join(dxapp_path[:-11], "Readme.md")
    if not os.path.exists(src_code) or not os.path.exists(readme_md):
        return None
    return AppObj(
        readme_md=os.path.join(dxapp_path[:-11], "Readme.md"),
        src_code=os.path.join(dxapp_path[:-11], dxapp_obj["runSpec"]["file"]),
        app_name=dxapp_obj["name"],
        title=dxapp_obj["title"],
        interpreter=dxapp_obj["runSpec"]["interpreter"])


def _resolve_applets(dxapp_paths):
    """Create app tutorial representations to be used downstream
    TODO: just using multiprocessing.processing since mapping pickling issue makes this complex
    TODO: Use namedtuple throughout code instead of converting to dict
    dictionary keys returned: "readme_md", "src_code", "name", "title", "interpreter"
    """
    workers = Pool(NUM_CORES)
    tutorial_tuples = [res for res in workers.map(_resolve_applet, dxapp_paths) if res is not None]
    workers.close()
    workers.join()
    mapping_list = [{k: v for k, v in izip(("readme_md", "src_code", "name", "title", "interpreter"), tupe)} for tupe in tutorial_tuples]
    return mapping_list


def _write_markdown(fh_md, readme_md_path, logger, section_parser):
    """Creates Kramdown webpage

    Notes: Add {match: func(line)} to dictionary in order to handle special matches
    """
    def _section_match(match):
        section = match.group(1).strip()
        logger.debug("Kramdown section in {location} {region_name}".format(
            location=readme_md_path, region_name=section))
        return section_dict[section]

    def _force_line_match(match):
        match = match.group(1).strip()
        logger.debug("Force line match: {}".format(match))
        return "\n{insert}\n".format(insert=match)

    def _func_match(match):
        func_name = match.group(1).strip()
        logger.debug("Adding Func: {}".format(func_name))
        return section_parser.get_func_code(func_name)

    def _opt_header2_match(match):
        fh_md.write('<hr>')

    section_dict = section_parser.parse_file().get_parse_dict()

    special_match = {  # Only one of these matches can be used. TODO flag multiple matches
        re.compile(r'.*<!--\s*SECTION:\s*\b(.*)-->'): _section_match,
        re.compile(r'.*<!--\s*INCLUDE:\s*(\S.*)-->'): _force_line_match,
        re.compile(r'.*<!--\s*FUNCTION:\s*(\S.*)-->'): _func_match
    }

    optional_match = {  # As many of these matches can be present in one line
        # re.compile(r'^\#\#\s.*|.*<!--\s*INCLUDE:\s*\#\#\s.*-->'): _opt_header2_match
    }

    with open(readme_md_path) as readme_md:
        for line in readme_md:
            for matcher, opt_func in optional_match.iteritems():
                match = matcher.match(line)
                if match:
                    opt_func(match)
            for matcher, func in special_match.iteritems():
                match = matcher.match(line)
                if match:
                    line = func(match)
                    break
            fh_md.write(line)


def find_all_matches(tutdir, filename, exclude_dirs=[]):
    """Return list of """
    matches = []
    for root, dirs, files in os.walk(tutdir):
        if filename in files:
            matches.append(os.path.join(root, filename))
        for dir_ex in exclude_dirs:
            if dir_ex in dirs:
                dirs.remove(dir_ex)
    return matches


def get_parser():
    """Return an argument parser for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--tutorials", help="Directory where tutorials are located", default=None, metavar="Tutorials directory", dest="tutorials_dir")
    parser.add_argument("--site-pages", help="Directory where site pages are located", default=None, metavar="Site Pages directory", dest="site_pages_dir")
    parser.add_argument("--keep-comments", help="Ignore other comments in code when parsing sections.", action="store_false", dest="skip_comments")
    parser.add_argument("--overwrite-files", help="Overwrites old files with generated files", action="store_true", dest="overwrite")

    return parser


def create_jekyll_markdown_tutorial(page_dict):
    """Creates a Jekyll site page

    Note:
        Generate a logger for each multithreaded process created.
        This can be improved with thread aware logging, TODO.

    Arguments:
        page_dict: Mapping with the following structure
            {
                readme_md: string path to Readme.md file location,
                src_code: string path to source code location
                name: string Page Name
                title: string Page Title
                interpreter: string Applet interpreter language
                overwrite: bool is okay to overwrite dest file if it exist?
                ignore_comments: bool is okay to ignore comments in files
                site_pages_dir: destination site page
            }

    TODO switch back to mapping once multiprocessing is implementended using Processing
    TODO support pages that don't need section parsers
    """
    page_basename = page_dict["name"].strip()
    setup_logger(
        logger_name=page_basename, log_dir='log_temp_dir',
        level=logging.DEBUG)
    proc_logger = logging.getLogger(page_basename)

    curr_date = date.today().isoformat()
    page_dict['date'] = curr_date
    target_file = os.path.join(
        page_dict["site_pages_dir"], curr_date + '-' + page_basename + ".md")

    for f in os.listdir('.'):  # TODO clean this loop up
        if fnmatch.fnmatch(f, '*{}*'.format(page_basename)):
            if page_dict["overwrite"]:
                proc_logger.debug("Removing: {0}".format(f))
                os.remove(f)
            else:
                proc_logger.info("Exist: {fn}".format(fn=f))
                return True, ""

    page_dict['ARCHIVE'] = 'archived' in page_dict['readme_md']
    section_parser = _get_section_parser(page_dict, logger=proc_logger)
    page_dict['isdocument'] = True  # Until Video support added leave this here

    try:
        with open(target_file, "w") as target_md:
            _write_front_matter(page_dict=page_dict, logger=proc_logger, file_handle=target_md)
            _write_markdown(
                fh_md=target_md,
                readme_md_path=page_dict["readme_md"],
                section_parser=section_parser,
                logger=proc_logger)
    except Exception as e:
        proc_logger.info("Exception: {msg}".format(msg=e))
        return False, "Failed with Error:\n{err}\n Review logs in log_temp_dir directory: {logname}.".format(
            err=e.message, logname=page_basename)

    return True, ""


def setup_logger(logger_name, log_dir, level=logging.INFO):
    """Create logger instance
    Inspiration: https://stackoverflow.com/questions/17035077/python-logging-to-multiple-log-files-from-different-classes
    """
    try:
        os.mkdir(log_dir)
    except OSError:
        pass
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    log_file = os.path.join(log_dir, logger_name)
    log_instance = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(asctime)s : %(message)s')
    fileHandler = logging.FileHandler(log_file, mode='w')
    fileHandler.setFormatter(formatter)
    # streamHandler = logging.StreamHandler()
    # streamHandler.setFormatter(formatter)

    log_instance.setLevel(level)
    log_instance.addHandler(fileHandler)
    # log_instance.addHandler(streamHandler)


def _write_front_matter(page_dict, logger, file_handle=None):
    """Generate and write liquid front matter"""
    proc_logger = logger
    tut_type = None
    for candidate_tut_type, re_compiled_obj in TUTORIAL_TYPES_SEARCH.iteritems():
        if re_compiled_obj.match(page_dict["name"]):
            tut_type = candidate_tut_type
            break

    frontmatter = FrontMatter(
        logger=logger,
        isdocument=page_dict['isdocument'])
        #source=page_dict["name"])

    frontmatter.add_field(field="date", value=page_dict["date"])
    frontmatter.add_field(field="title", value=page_dict["title"])
    language = SUPPORTED_INTERPRETERS.get(page_dict["interpreter"])
    if not page_dict['ARCHIVE'] and tut_type:
        frontmatter.add_field(field="categories", value=tut_type)
    if not page_dict['ARCHIVE'] and language:
        frontmatter.add_field(field="categories", value=language)
    if page_dict['ARCHIVE']:
        frontmatter.add_field('categories', value="Example Applet")

    if file_handle is None:
        return frontmatter.__str__()
    file_handle.write(frontmatter.__str__())
    proc_logger.info("front matter written")


def main():
    """Script Entry point."""
    parser = get_parser()
    args = parser.parse_args()

    if args.tutorials_dir is None:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        assumed_tutorial_path = os.path.join(script_dir, "..", "..", "Tutorials")
        args.tutorials_dir = os.path.normpath(assumed_tutorial_path)
    if args.site_pages_dir is None:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        assumed_site_page_path = os.path.join(script_dir, "..", "_posts")
        args.site_pages_dir = os.path.normpath(assumed_site_page_path)

    _rehydrate_site(args)


if __name__ == '__main__':
    main()
