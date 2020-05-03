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
    """Simple tree structure - represents tree as linked list/dict."""
    def __init__(self, nodes=None):
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
        return list(self.nodes.values())

    def get_root(self):
        for node in self.nodes.values():
            if node.parent == None:
                return node

    def post_order_traversal(self):
        root = self.get_root()
        return self._passdown_child(root, [])


    def _passdown_child(self, node, pot):
        for child in node.children:
            self._passdown_child(child, pot)
        pot.append(node)
        return pot
    
    def __repr__(self):
        nodes = []
        for node in self.list_nodes():
            nodes.append(str(node))
        return '\n'.join(nodes)


# TODO: Add tokens for internal nodes, branch-lengths etc.
class NewickIO(object):
    """Parse and write Newick trees"""

    def __init__(self):
        self.grammar = {
                '(': self._add_branch,
                ')': self._close_branch
                }
        
        self.tokens = None
        self.currnode = Node('_0')
        self.currnode = None
        # self.inode = 0
        # self.tree = Tree({'_0': self.currnode})
        self.tree = Tree()

    def parse(self, tokens=None):
        """Parses simplest format newick string and returns Tree structure.
        
        By default, internal nodes are numbered in order of encounter and
            prefixed with an underscore.
        """
        if type(tokens) == str:
            self._tokenize(tokens)

        # Base case - close if terminal ';' is reached.
        token = self.tokens[0]
        if token == ';':
            return self.tree

        # Parse tokens recursively
        action = self.grammar.get(token, self._add_leaf)
        action(token)
        self.tokens = self.tokens[1:]

        return self.parse(self.tokens)

    def _tokenize(self, newick_string):
        """Splits newick string into set of tokens for parsing."""
        if newick_string.count('(') != newick_string.count(')'):
            raise ValueError('parentheses are unbalanced')
        self.tokens = re.findall(r'\w+|\(|\)|\;', newick_string)

    def _add_branch(self, *tokens):
        """Adds internal node from which a new branch will grow.

        '(' -> open branch
        """
        self.inode += 1
        newinternal = Node(f'_{self.inode}', self.currnode)
        self.currnode.add_child(newinternal)
        self.currnode = newinternal
        self.tree.add_node(f'_{self.inode}', self.currnode)

    def _close_branch(self, *tokens):
        """Returns to internal node from which branch originates.

        ')' -> close branch
        """
        self.inode -= 1
        self.currnode = self.tree.get_node(f'_{self.inode}')

    def _add_leaf(self, *tokens):
        """Adds terminal leaf without descendents.

        <alphanumeric leaf name> -> Add leaf
        """
        newleaf = Node(tokens[-1], self.currnode)
        self.currnode.add_child(newleaf)
        self.tree.add_node(tokens[-1], newleaf)
        


if __name__ == '__main__':
    NwkParser = NewickIO()
    tree = NwkParser.parse('((A,B),(H,(F,G)));')
    print(tree)
    # test = tree.post_order_traversal()
    # for i in test:
    #     print(i.name)
