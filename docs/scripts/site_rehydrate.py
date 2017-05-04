import dxpy
import argparse
import os
import json
import re
from __future__ import print_function
from multiprocessing import Pool
from itertools import izip
from helpers import FrontMatter, SectionParser
from global_helper_vars import TUTORIAL_TYPES_SEARCH, SUPPORTED_INTERPRETERS, NUM_CORES


def _rehydrate_site(user_args):
    user_input_dict = {
        "overwrite": user_args.overwrite,
        "ignore_comments": user_args.skip_comments,
        "site_pages_dir": user_args.site_pages_dir
    }
    dxapp_files = find_all_matches(user_args.tutorials_dir, "dxapp.json")
    tutorial_dicts = [d.update(user_input_dict) for d in _resolve_applets(dxapp_files)]
    md_maker = Pool(processes=NUM_CORES)
    status = md_maker.map(create_jekyll_markdown_tutorial, tutorial_dicts)
    md_maker.close()
    md_maker.join()

def _resolve_applet(dxapp_path):
    """Moved to top level due to Pool.map() inability to use nested func"""
    with open(dxapp_path, "r") as dxf:
        dxapp_obj = json.load(dxf)
    src_code = os.path.join(dxapp_path[:-11], dxapp_obj["runSpec"]["file"])  # pretty sure this is a required field
    readme_md = os.path.join(dxapp_path[:-11], "Readme.md")
    if os.path.exists(src_code) and os.path.exists(readme_md):
        app_name = dxapp_obj["name"]
        title = dxapp_obj["title"]
        interpreter = dxapp_obj["runSpec"]["interpreter"]
        return (readme_md, src_code, app_name, title, interpreter)
    return None


def _resolve_applets(dxapp_paths):
    """
    TODO: just using multiprocessing.processing since mapping pickling issue makes this complex
    Create app tutorial representations to be used downstream
    dictionary keys returned: "readme_md", "src_code", "name", "title", "interpreter"
    """
    workers = Pool(NUM_CORES)
    tutorial_tuples = [res for res in workers.map(_resolve_applet, dxapp_paths) if res is not None]
    workers.close()
    workers.join()
    mapping_list = [{k: v for k, v in izip(("readme_md", "src_code", "name", "title", "interpreter"), tupe)} for tupe in tutorial_tuples]
    return mapping_list


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
    parser.add_argument("tutorials_dir", help="Directory where tutorials are located", default=os.path.abspath('../../Tutorials'))
    parser.add_argument("site_pages_dir", help="Directory where site pages are located", default=os.path.abspath('../_tutorials'))
    parser.add_argument("--ignore-comments", help="Ignore other comments in code when parsing sections.", action="store_true", dest="skip_comments")
    parser.add_argument("--overwrite-files", help="Overwrites old files with generated files", action="store_true", dest="overwrite")
    return parser


def create_jekyll_markdown_tutorial(tutorial_dict):
    """
    :param dxapp_obj: tuple describing an applet.
        dict keys: "readme_md", "src_code", "name", "title"
                    "interpreter" "overwrite" "ignore_comments" "site_pages_dir"
    :type dxapp_obj: dict
    :overwrite: if target markdown already exist determins if it is overriden.
    :type overwrite: boolean
    TODO switch back to mapping once multiprocessing is implementended using Processing
    """
    class InvalidSection(Exception):
        pass

    def _write_front_matter(dxapp_obj, file_handle):
        """Generate and write liquid front matter"""
        frontmatter = FrontMatter(
            {
                "title": dxapp_obj["title"],
                "source": dxapp_obj["name"]
            })
        for k, v in TUTORIAL_TYPES_SEARCH:
            if v.match(dxapp_obj["name"]):
                frontmatter.add_field(k, v)
                break
        lang = SUPPORTED_INTERPRETERS.get(dxapp_obj["interpreter"], "none")
        frontmatter.add_field("language", lang)
        file_handle.write(frontmatter.__str__())

    def _write_kmarkdown(fh_md, readme_md_path, section_parser={}):
        """Creates Kramdown webpage"""
        section_match = re.compile(r'')
        with open(readme_md_path):
            for line in readme_md_path:

        pass

    target_file = os.path.join(
        tutorial_dict["site_pages_dir"], tutorial_dict["name"].strip() + ".md")
    if os.path.exists(target_file):
        if tutorial_dict["overwrite"]:
            os.remove(target_file)
        else:
            return True

    # Handle JSON for assets at some point, should be easy.
    code_block_dict = {}
    if SUPPORTED_INTERPRETERS.get(tutorial_dict["interpreter"]) is not None:
        code_block_dict = SectionParser(
            code_file_path=tutorial_dict["src_code"],
            ignore_comments=tutorial_dict["ignore_comments"],
            interpreter=tutorial_dict["interpreter"]).parse_file().get_section_dict()
    try:
        with open(target_file, "w") as tutorial_md:
            _write_front_matter(tutorial_dict, tutorial_md)
            _write_kmarkdown(
                fh_md=tutorial_md,
                readme_md_path=tutorial_dict["readme_md"],
                section_dict=code_block_dict)
    except InvalidSection as IE:
        # Change to logging
        print(IE.message)
        return False

    return True


def main():
    """Script Entry point."""
    parser = get_parser()
    args = parser.parse_args()
    _rehydrate_site(args)


if __name__ == '__main__':
    main()
