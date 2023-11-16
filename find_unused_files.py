#!/usr/bin/env python3

from pathlib import Path
import os
import argparse
from fnmatch import fnmatch
import re
import sys


# ex = re.compile(r"(\b(?<!std::)size_t)\b")
#
# github = "GITHUB_ACTIONS" in os.environ


def main():
    # p = argparse.ArgumentParser()
    # p.add_argument("input")
    # p.add_argument(
    #     "--fix", action="store_true", help="Attempt to fix any license issues found."
    # )
    # p.add_argument("--exclude", "-e", action="append", default=[])
    #
    # args = p.parse_args()

    # walk over all files
    exit = 0

    # Collect all suffix for investigation
    suffix_dict = {}
    for root, _, files in os.walk("."):
        root = Path(root)

        # Skip "Scripts"-folder
        if str(root).find("Scripts") != -1:
            continue

        # Skip "thirdparty"-folder
        if str(root).find("thirdparty") != -1:
            continue

        # Skip "CI"-folder
        if str(root).find("CI") != -1:
            continue

        # Skip "git"-folders
        if str(root).find("git") != -1:
            continue

        # Skip "cmake"-folders
        if str(root).find("cmake") != -1:
            continue

        # Skip base-directory
        if str(root) == ".":
            continue

        # Skip ".idea"-folders
        # TODO remove after testing
        if str(root).find(".idea") != -1:
            continue

        for filename in files:

            if str(filename) == ".gitignore":
                continue

            # get the full path of the file
            filepath = root / filename
            suffix_dict[filepath.suffix] = ""


            # Check header files and remove
            if filepath.suffix in (".hpp", ".h", ".cuh"):
                # continue
                cmd = 'grep -IR "' + filename + '" Alignment Core Examples Fatras Plugins Tests > unused_files_tmp'
                os.system(cmd)
                output = os.popen('cat unused_files_tmp').read()
                if output.count('\n') == 0:
                    print(f"Remove file\n{filename}\n{filepath}\n")
                    remove_cmd = 'rm ' + str(filepath)
                    os.system(remove_cmd)


            # Check source files and remove
            if filepath.suffix in (".cpp", ".c", ".C", ".cu", ".ipp"):
                # continue
                cmd = 'grep -IR "' + filename + '" Alignment Core Examples Fatras Plugins Tests > unused_files_tmp'
                os.system(cmd)
                output = os.popen('cat unused_files_tmp').read()
                if output.count('\n') == 0:
                    print(f"Remove file\n{filename}\n{filepath}\n")
                    remove_cmd = 'rm ' + str(filepath)
                    os.system(remove_cmd)


            # Check images and remove
            # TODO make jpeg unnecessary
            if filepath.suffix in (".png", ".jpg", ".svg", ".jpeg", ".gif"):
                # continue
                # Skip "acts_logo"
                if str(filename).find("acts_logo") != -1:
                    continue

                # Skip "white-paper-figures"
                if str(root).find("white_papers/figures") != -1:
                    continue

                cmd = 'grep -IR "' + filename + '" Alignment Core Examples Fatras Plugins Tests docs > unused_files_tmp'
                os.system(cmd)
                output = os.popen('cat unused_files_tmp').read()
                if output.count('\n') == 0:
                    print(f"Remove file\n{filename}\n{filepath}\n")
                    remove_cmd = 'rm ' + str(filepath)
                    os.system(remove_cmd)


            # Check and print other files
            if filepath.suffix in (
                    '.yml',
                    '.txt',
                    '.ini',
                    '',
                    '.toml',
                    '.cff',
                    '.json',
                    '.in',
                    '.patch',
                    '.sh',
                    '.xsl',
                    '.lock',
                    '.imp',
                    '.yaml',
                    '.root',
                    '.ipynb',
                    '.csv',
                    '.j2',
                    '.css',
                    '.gdml',
                    '.hepmc3',
                    '.onnx',
                    '.idx',
                    '.pack',
                    '.sample',
                    '.xml',
                    '.iml',
                    ):
                # Skip "vertexing_event_mu20" because they don't appear in full-text
                if str(filename).find("vertexing_event_mu20_") != -1:
                    continue

                cmd = 'grep -IR "' + filename + '" Alignment Core Examples Fatras Plugins Tests docs > unused_files_tmp'
                os.system(cmd)
                output = os.popen('cat unused_files_tmp').read()
                if output.count('\n') == 0:
                    print(f"Remove file\n{filename}\n{filepath}\n")
                    # remove_cmd = 'rm ' + str(filepath)
                    # os.system(remove_cmd)


            # Not implemented tests
            if filepath.suffix in (".py"):
                continue

            # Check documentation files (weak tests)
            if filepath.suffix in (".md", ".rst"):
                cmd = 'grep -IR "' + filepath.stem + '" Alignment Core Examples Fatras Plugins Tests docs > unused_files_tmp'
                os.system(cmd)
                output = os.popen('cat unused_files_tmp').read()
                if output.count('\n') == 0:
                    print(f"Remove file\n{filename}\n{filepath}\n")
                    # remove_cmd = 'rm ' + str(filepath)
                    # os.system(remove_cmd)

    print(suffix_dict)

    return exit



if "__main__" == __name__:
    sys.exit(main())
