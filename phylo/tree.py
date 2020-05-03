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

        return f'Name: {name}\n\tParent: {parent}\n\tChildren: [{children}]'


class Tree(object):
    def __init__(self, nodes):
        self.root = None
        if nodes:
            self.nodes = nodes
        else:
            self.nodes = {}

    def add_node(self, name, node):
        self.nodes[name] = node

    def get_node(self, name):
        return self.nodes[name]

    def list_nodes(self):
        return list(self.nodes.items())


class NewickIO(object):
    """Parse and write Newick trees"""

    def __init__(self):
        self.grammar = {
                '(': self._add_branch,
                ')': self._close_branch
                }
        
        self.tokens = None
        self.currnode = Node('_0')
        self.inode = 0
        self.tree = Tree({'_0': self.currnode})

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
        self.tree.add_node(f'_{self.inode}', self.currnode)

    def _close_branch(self, *argv):
        self.inode -= 1
        self.currnode = self.tree.get_node(f'_{self.inode}')

    def _add_leaf(self, *argv):
        newleaf = Node(argv[-1], self.currnode)
        self.currnode.add_child(newleaf)
        self.tree.add_node(argv[-1], newleaf)
        


if __name__ == '__main__':
    NwkParser = NewickIO()
    tree = NwkParser.parse('(B,(A,C,E),D);')
    for name, node in tree.list_nodes():
        print(node)
