#! /usr/bin/env python
import os
import sys
import getopt
import frontmatter


def main():
    opts, _ = getopt.getopt(sys.argv[1:], 'b:')
    with open('headers/' + opts[0][1]) as f:
        m, _ = frontmatter.parse(f.read())

        code = """
        pandoc headers/{} -i {} --bibliography=./mscl_refs.bib --filter=pandoc-eqnos --columns 6 --filter=pandoc-crossref -o {}
        """.format(m['header'], m['include'], m['name'])
        os.system(code)


if __name__ == '__main__':
    main()
    print('document successfully compiled')
