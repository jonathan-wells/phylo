#!/usr/bin/env python3

import re


class Node(object):

    def __init__(self, name=None, parent=None):
        self.name = name
        self.parent = parent
        self.children = []

    def add_child(self, child):
        self.children.append(child)
    
    def add_parent(self, parent):
        self.parent = parent

    def get_parent(self):
        return self.parent

    def get_children(self):
        return self.children

    def __repr__(self):
        name = self.name
        if self.parent == None:
            parent = 'None'
        else:
            parent = self.parent.name
        children = ', '.join([child.name for child in self.children]) 

        return f'{name}\n\tParent: {parent}\n\tChildren: [{children}]'


class NewickParser(object):
    """Parse Newick trees"""

    def __init__(self, newick_string=None):
        self.grammar = {
                '(': self._add_branch,
                ')': self._close_branch
                }
        
        self.tokens = None
        self.currnode = Node('_0')
        self.inode = 0
        self.tree = {'_0': self.currnode}
        if newick_string:
            self._tokenize(newick_string)

    def parse(self, tokens=None):
        if type(tokens) == str:
            self._tokenize(tokens)

        token = self.tokens[0]
        if token == ';':
            return self.tree

        action = self.grammar.get(token, self._add_leaf)
        action(token)

        self.tokens = self.tokens[1:]
        return self.parse(self.tokens)

    def _tokenize(self, newick_string):
        if newick_string.count('(') != newick_string.count(')'):
            raise ValueError('parentheses are unbalanced')
        self.tokens = re.findall(r'\w+|\(|\)|\;', newick_string)

    def _add_branch(self, *argv):
        self.inode += 1
        newinternal = Node(f'_{self.inode}', self.currnode)
        self.currnode.add_child(newinternal)
        self.currnode = newinternal
        self.tree[f'_{self.inode}'] = self.currnode

    def _close_branch(self, *argv):
        self.inode -= 1
        self.currnode = self.tree[f'_{self.inode}']

    def _add_leaf(self, *argv):
        newleaf = Node(argv[-1], self.currnode)
        self.currnode.add_child(newleaf)
        self.tree[argv[-1]] = newleaf


class BFS(object):
    def __init__(self, tree):
        pass

    def search(node, condition):
        pass
        


if __name__ == '__main__':
    tree = '(B,(A,C,E),D);'
    NwkParser = NewickParser(tree)
    tree = NwkParser.parse()
    for key, val in tree.items():
        print(val)
