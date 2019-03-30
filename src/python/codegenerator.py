#!/usr/bin/python3

import sympy as sp
from sympy.printing.codeprinter import Assignment
from sympy.printing.cxxcode import CXX17CodePrinter

class CxxMatrixPrinter(CXX17CodePrinter):
    def _print_list(self, list_of_exprs):
        if all(isinstance(x, sp.Eq) for x in list_of_exprs):
            list_of_lhs = [x.lhs for x in list_of_exprs]
            list_of_rhs = [x.rhs for x in list_of_exprs]
            common_variables, list_of_rhs = sp.cse(list_of_rhs, symbols=sp.numbered_symbols(prefix='com'))
            lines = []
            for variable, expression in common_variables:
                ass = Assignment(variable, expression)
                lines.append('const double ' + self._print(ass))
            for i in range(len(list_of_rhs)):
                ass = Assignment(list_of_lhs[i], list_of_rhs[i])
                lines.append(self._print(ass))
        return '\n'.join(lines)

from sympy.printing.ccode import _as_macro_if_defined
from sympy.printing.ccode import precedence
class EigenMatrixPrinter(CxxMatrixPrinter):
    def _print_MatrixElement(self, expr):
        return "{0}({1},{2})".format(
            self.parenthesize(expr.parent, sp.printing.precedence.PRECEDENCE["Atom"], strict=True),
            expr.i, expr.j)
    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for j in range(cols) for i in range(rows))

    # @_as_macro_if_defined
    # def _print_Pow(self, expr):
        
    #     # if "Pow" in self.known_functions:
    #     #     return self._print_Function(expr)
    #     PREC = precedence(expr)
    #     # suffix = self._get_func_suffix(real)
    #     return '%s.pow(%s)' % (self.parenthesize(expr.base, PREC), self._print(expr.exp))

    #     # if expr.exp == -1:
    #     #     return '1.0%s/%s' % (suffix.upper(), self.parenthesize(expr.base, PREC))
    #     # elif expr.exp == 0.5:
    #     #     return '%ssqrt%s(%s)' % (self._ns, suffix, self._print(expr.base))
    #     # elif expr.exp == S.One/3 and self.standard != 'C89':
    #     #     return '%scbrt%s(%s)' % (self._ns, suffix, self._print(expr.base))
    #     # else:
    #     #     return '%spow%s(%s, %s)' % (self._ns, suffix, self._print(expr.base),
    #     #                            self._print(expr.exp))

def codegen(expressions, template_path, output_path):
    p = EigenMatrixPrinter()
    p._print_Pow
    with open(template_path, 'r') as f:
        template = f.read()
    
    code = template.format(code=p.doprint(expressions))     
    with open(output_path, 'w') as f:
        f.write(code)