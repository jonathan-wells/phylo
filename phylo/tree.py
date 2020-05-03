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
<<<<<<< HEAD
    """Simple tree structure - represents tree as linked list/dict."""
=======
    """Simple self._tree structure - represents self._tree as linked list/dict."""
>>>>>>> fixparser
    def __init__(self, nodes=None):
        self.root = None
        if nodes:
            self.nodes = nodes
        else:
            self.nodes = []

    def add_node(self, node):
        self.nodes.append(node)

    def get_root(self):
        for node in self.nodes:
            if node.parent == None:
                return node

<<<<<<< HEAD
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

=======
    def post_order_traversal(self):
        root = self.get_root()
        return self._passdown_child(root, [])

    def _passdown_child(self, node, pot):
        pot.append(node)
        for child in node.children:
            self._passdown_child(child, pot)
        return pot
    
    def __repr__(self):
        nodes = []
        for node in self.nodes:
            nodes.append(str(node))
        return '\n'.join(nodes)

>>>>>>> fixparser

# TODO: Add tokens for internal nodes, branch-lengths etc.
class NewickIO(object):
    """Parse and write Newick self._trees"""

    def __init__(self):
        self._grammar = {
                'addbranch': self._add_branch,
                'closebranch': self._close_branch,
                'addleaf': self._add_leaf
                }
<<<<<<< HEAD
        
        self.tokens = None
        self.currnode = Node('_0')
        self.currnode = None
        # self.inode = 0
        # self.tree = Tree({'_0': self.currnode})
        self.tree = Tree()
=======
        self._tree = Tree() 
        self._queue = []
        self._internalcount = 0
>>>>>>> fixparser

    def parse(self, tokens=None):
        if type(tokens) == str:
            tokens = self._tokenize(tokens)
        
        currtoken = tokens[0]
        if currtoken == '(':
            self._grammar['addbranch']()
        elif currtoken == ')':
            self._grammar['closebranch']()
        elif re.match(r'\w+', currtoken):
            self._grammar['addleaf'](currtoken)
        elif currtoken == ';':
            return
        
        self.parse(tokens[1:])
        return self._tree

    def _tokenize(self, newick_string):
        """Splits newick string into set of tokens for parsing."""
        if newick_string.count('(') != newick_string.count(')'):
            raise ValueError('parentheses are unbalanced')
        return re.findall(r'\w+|\(|\)|\;', newick_string)

    def _add_branch(self):
        """Adds internal node from which a new branch will grow.

        '(' -> open branch
        """
        if len(self._queue) == 0:
            newinternal = Node(f'_{self._internalcount}')
        else:
            newinternal = Node(f'_{self._internalcount}', self._queue[-1])
            self._queue[-1].add_child(newinternal)
        self._queue.append(newinternal)
        self._tree.add_node(newinternal)
        self._internalcount += 1

    def _close_branch(self):
        """Returns to internal node from which branch originates.

        ')' -> close branch
        """
        self._queue.pop()

    def _add_leaf(self, name):
        """Adds terminal leaf without descendents.

        <alphanumeric leaf name> -> Add leaf
        """
        newleaf = Node(name, self._queue[-1])
        self._queue[-1].add_child(newleaf)
        self._tree.add_node(newleaf)


if __name__ == '__main__':
    NwkParser = NewickIO()
<<<<<<< HEAD
    tree = NwkParser.parse('((A,B),(H,(F,G)));')
    print(tree)
    # test = tree.post_order_traversal()
    # for i in test:
    #     print(i.name)
=======
    tree = NwkParser.parse('((A,B),(H,(G,F)));')
    # print(tree)
    # print(tree.get_root())

    test = tree.post_order_traversal()
    for i in test:
        print(i.name)
>>>>>>> fixparser
