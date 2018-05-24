#! /usr/bin/env python
import os
import sys
import getopt 
import frontmatter
import glob

def main():
    if os.path.exists('_html') == False:
        os.mkdir('_html')
    opts, _ =  getopt.getopt(sys.argv[1:], 'b:')
    with open('headers/' + opts[0][1]) as f:
        m, _ = frontmatter.parse(f.read())
        num_sects = len([x for x in m['include'].split(' ')])
        sect_ids = [x for x in m['include'].split(' ')]
        code = """
        pandoc headers/{} -i {} --filter pandoc-crossref  --filter pandoc-eqnos --filter pandoc-citeproc  --bibliography=./mscl_refs.bib  --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML -H headers/header.txt  -o _html/_{} 
        """.format(m['header'], m['include'], m['name'])
        os.system(code)

        with open("_html/_{}".format(m['name']), 'rt') as _ht:
            # Set up the empty list of sections.
            sections = [[] for i in range(num_sects + 2)]
            sect_id = 0
            for line in _ht:
                if '<h2' in line:
                    sect_id += 1
                else:
                    sections[sect_id].append(line)
 
        for i in range(1, num_sects + 1):
            new_name = '2018-05-{:02d}-{}.html'.format(i, sect_ids[i-1].split('.md')[0].split('_')[-1])
            with open("_html/{}".format(new_name), 'wt') as ht:
                for line in sections[i]:
                    if line != None:
                        if ('.pdf' in line) or ('../figs/' in line):
                            out = line.replace('.pdf', '.png')
                            out = out.replace('../figs/', 'assets/')
                        else:
                            out = line
                        ht.write('{}{}{}'.format('{% raw %}',out,'{% endraw %}'))
                    
        os.remove('_html/_{}'.format(m['name']))
    
if __name__ == '__main__':
    main()
    print('html document compiled.')



