#!/usr/bin/env python3

from . import tree
import re

# TODO: Add tokens for internal nodes, branch-lengths etc.
class NewickIO(object):
    """Parse and write Newick trees"""

    def __init__(self):
        self._grammar = {
                'addbranch': self._add_branch,
                'closebranch': self._close_branch,
                'addleaf': self._add_leaf
                }
        self._tree = tree.Tree() 
        self._queue = []
        self._internalcount = 0

    ############################################################################
    ## Public reader methods
    ############################################################################

    def read(self, filename):
        with open(filename) as infile:
            nwkstring = infile.readline().strip()
        return self.parse(nwkstring)

    def parse(self, tokens=None):
        if type(tokens) == str:
            self._tree = tree.Tree()
            tokens = self._tokenize(tokens)
        
        currtoken = tokens[0]
        if currtoken == '(':
            treeop = self._grammar['addbranch']
        elif currtoken == ')':
            treeop = self._grammar['closebranch']
        elif re.match('\)\:\d+\.\d+', currtoken):
            treeop = self._grammar['closebranch']
        elif re.match('\)\w+\:\d+\.\d+', currtoken):
            treeop = self._grammar['closebranch']
        elif re.match('\)\w+', currtoken):
            treeop = self._grammar['closebranch']
        elif re.match(r'\w+', currtoken):
            treeop = self._grammar['addleaf']
        elif currtoken == ';':
            return
        
        treeop(currtoken)
        self.parse(tokens[1:])
        
        return self._tree
    
    ############################################################################
    ## Public writer methods
    ############################################################################
    
    def write(self, tree, filename):
        with open(filename, 'w') as outfile:
            nwkstring = self.to_string(tree)
            outfile.write(f'{nwkstring}\n')

    def to_string(self, tree):
        nwkstring = self._to_string(tree) + ';'
        return nwkstring

    ############################################################################
    ## Private parser subroutines
    ############################################################################
    
    def _to_string(self, tree, node=None, nwkstring=''):
        if not node:
            node = tree.get_root()
        
        if len(node.children) == 0:
            if node.brlen:
                token = f'{node.name}:{node.brlen}'
            else:
                token = node.name
            nwkstring += token
        else:
            nwkstring += '('
            n = len(node.children)
            for i, child in enumerate(node.children):
                nwkstring = self._to_string(tree, child, nwkstring)
                if i < n - 1:
                    nwkstring += ','
            if not node.name.startswith('_') and node.brlen:
                token = f'){node.name}:{node.brlen}'
            elif node.brlen:
                token = f'):{node.brlen}'
            elif not node.name.startswith('_'):
                token = f'){node.name}'
            else:
                token = ')'
            nwkstring += token
        
        return nwkstring
        
    def _tokenize(self, newick_string):
        """Splits newick string into set of tokens for parsing."""
        if newick_string.count('(') != newick_string.count(')'):
            raise ValueError('parentheses are unbalanced')
        elif not newick_string.endswith(';'):
            raise ValueError('missing tree terminator: ;;')

        # Token list
        # TODO: This list of tokens could be improved. Maybe make it a class
        # atttribute?
        leaf = r'\w+'
        leaf_br = r'\w+\:\d+\.\d+'
        intl = r'\('
        intl_c = r'\)'
        intl_c_n = r'\)\w+'
        intl_c_n_br = r'\)\w+\:\d+\.\d+'
        intl_c_br = r'\)\:\d+\.\d+'
        
        tokpat = f'{intl}|{leaf_br}|{leaf}|{intl_c_n_br}|{intl_c_n}|{intl_c_br}|{intl_c}|;'
        return re.findall(tokpat, newick_string)

    def _add_branch(self, token):
        """Adds internal node from which a new branch will grow.

        '(' -> open branch
        """
        meta = re.findall(r'\w+|\(', token)
        if len(meta) == 2:
            name = meta[0]
        else:
            name = f'_{self._internalcount}'
        if len(self._queue) == 0:
            newinternal = tree.Node(name=name)
        else:
            newinternal = tree.Node(name=name, parent=self._queue[-1])
            self._queue[-1].add_child(newinternal)
        self._queue.append(newinternal)
        self._tree.add_node(newinternal)
        self._internalcount += 1

    def _close_branch(self, token):
        """Returns to internal node from which branch originates.

        ')' -> close branch
        """
        meta = re.findall(r'\)|\w+|\:\d+\.\d+', token)
        if len(meta) == 2:
            if meta[1].startswith(':'):
                self._queue[-1].set_brlen(float(meta[1].strip(':')))
            else:
                self._queue[-1].name = meta[1]  # TODO: Write setter for node name
        elif len(meta) == 3:
            self._queue[-1].name = meta[1]  # TODO: Write setter for node name
            self._queue[-1].set_brlen(float(meta[2].strip(':')))
        self._queue.pop()

    def _add_leaf(self, token):
        """Adds terminal leaf without descendents.

        <alphanumeric leaf name> -> Add leaf
        """
        meta = token.split(':')
        if len(meta) == 1:
            name, brlen = meta[0], None
        elif len(meta) == 2:
            name, brlen = meta
        newleaf = tree.Node(name=name, parent=self._queue[-1], brlen=brlen)
        self._queue[-1].add_child(newleaf)
        self._tree.add_node(newleaf)


if __name__ == '__main__':
    NwkParser = NewickIO()
    nwkstrs = ['((A,B),(C,(D,E)));',
               '((A:0.2,B:0.1),(C:0.3,(D:0.1,E:0.1)));',
               '((A:0.2,B:0.1):0.6,(C:0.3,(D:0.1,E:0.1):0.5):0.7);',
               '((A,B)i1,(C,(D,E)i3)i2)i0;',
               '((A:0.2,B:0.1)i1:0.6,(C:0.3,(D:0.1,E:0.1)i3:0.5)i2:0.7)i0;']
    for nwkstr in nwkstrs:
        thing = NwkParser.parse(nwkstr)
        newnwkstr = NwkParser.to_string(thing)
        assert newnwkstr == nwkstr
        print(nwkstr, newnwkstr, sep='\n', end='\n\n')
