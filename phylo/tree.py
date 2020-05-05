#!/usr/bin/env python3


class Node(object):

    def __init__(self, name=None, parent=None, brlen=None):
        self.name = name
        self.parent = parent
        self.children = []
        self.brlen = brlen

    def add_child(self, child):
        self.children.append(child)
    
    def add_parent(self, parent):
        self.parent = parent

    def set_brlen(self, brlen):
        self.brlen = brlen

    def get_parent(self):
        return self.parent

    def get_children(self):
        return self.children
    
    @property
    def isleaf(self):
        if self.degree == 1:
            return True
        else:
            return False

    @property 
    def degree(self):
        degree = 0
        if self.parent:
            degree += 1
        degree += len(self.children)
        return degree

    def __repr__(self):
        name = f'Name: {self.name}'
        if self.parent == None:
            parent = 'Parent: None'
        else:
            parent = f'Parent: {self.parent.name}'
        children = 'Children: ['+', '.join([c.name for c in self.children])+']'
        brlen = f'Branch length: {self.brlen}'

        return f'{name}\n\t{parent}\n\t{children}\n\t{brlen}'


class Tree(object):
    """Simple self._tree structure - represents self._tree as linked list/dict."""
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

    def post_order_traversal(self):
        root = self.get_root()
        return self._passdown_child(root, [])

    def _passdown_child(self, node, pot):
        pot.append(node)
        for child in node.children:
            self._passdown_child(child, pot)
        return pot
    
    @property
    def isbinary(self):
        for node in self.nodes:
            if node.degree > 3:
                return False
        return True

    def __repr__(self):
        nodes = []
        for node in self.nodes:
            nodes.append(str(node))
        return '\n'.join(nodes)


if __name__ == '__main__':
    root = Node(name='_0')
    node = Node(name='A', parent=root, brlen=0.1)
    root.add_child(node)
    tree = Tree()
    tree.add_node(root)
    tree.add_node(node)
    print(tree)
    print(root.degree)
    print(node.degree)
    print(root.isleaf)
